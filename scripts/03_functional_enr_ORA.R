# ==== Loading libraries ====
# Suppress logs
pkgs <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
          "enrichplot", "pheatmap", "dplyr", "ggplot2", 
          "stringr", "gridExtra", "patchwork", "msigdbr", 
          "ggplotify", "pathview", "KEGGREST")
invisible(lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
}))

# ==== Loading previous steps ====
dds <- readRDS(snakemake@input[["dds_results"]])
res <- results(dds, contrast = c("condition", "treated", "control"))
resOrdered <- res[order(res$padj), ]

# ==== Gene Ontology (GO) Enrichment
# clean Ensembl IDs
ensembl_ids <- gsub("\\..*", "", rownames(resOrdered))

# Filter significant genes (padj < 0.05 and |LFC| > 1)
sig_mask <- !is.na(resOrdered$padj) & resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) > 1
sig_genes <- ensembl_ids[sig_mask]

# Convert Ensembl ID -> Entrez ID (for clusterProfiler)
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# enrichGO Biological Processes
go_enrich <- enrichGO(gene          = entrez_ids,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

# GO graphs
if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
    # Reduction of GO terms redundancy
    go_simple <- simplify(go_enrich, cutoff = 0.7, by = "p.adjust", select_fun = min)
    
    # dotplot
    p_go <- dotplot(go_simple, showCategory = 10) +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
            theme_minimal(base_size = 15) +
            labs(title = 'GO Enrichment Analysis', subtitle = 'Biological Processes')
    
    ggsave(snakemake@output[["go_plot"]], plot = p_go, width = 14, height = 8, dpi = 150)

    # barplot
    #p_bar_go <- barplot(go_simple, showCategory = 10) + 
    #            scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    #            theme_minimal(base_size = 10)
    #ggsave(snakemake@output[["go_barplot"]], plot = p_bar_go, width = 8, height = 7, dpi = 150)
    
    # combined plot
    #combined_plot_go <- grid.arrange(p_go, p_bar_go, ncol = 2)
    #combined_plot_go <- (p_go | p_bar_go) + 
    #  plot_annotation(title = 'GO Enrichment Analysis - Biological Processes',
    #                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) +
    #  plot_layout(guides = "collect")
    #ggsave(snakemake@output[["go_combined_plot"]], plot = combined_plot_go, width = 14, height = 8, dpi = 150)
    
} else {
    # Creates empty files if no enriched pathways to avoid Snakemake failure
    #file.create(snakemake@output[["go_plot"]], snakemake@output[["go_barplot"]])
    file.create(snakemake@output[["go_plot"]])
}

# ==== KEGG Enrichment ====
kegg_ora <- enrichKEGG(gene         = entrez_ids,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)

p_kegg <- dotplot(kegg_ora, showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_minimal(base_size = 15) +
  labs(title = 'KEGG Enrichment Analysis', subtitle = 'Biological Processes')
p_kegg

ggsave(snakemake@output[["kegg_plot"]], plot = p_go, width = 14, height = 8, dpi = 150)


# hsa04110 è il Cell Cycle. Usa l'ID del tuo top pathway.
pv.out <- pathview(gene.data  = entrez_ids, 
         pathway.id = kegg_ora$ID[1], 
         species    = "hsa", 
         limit      = list(gene=2, cpd=1)) # color limit to LFC 2
# Visualize in plot panel in interactive mode
#img <- readPNG(paste0(kegg_ora$ID[1],".pathview.png"))
#grid.raster(img)
# For Snakemake:
file.rename(from = paste0(kegg_ora$ID[1], ".pathview.png"), 
            to = snakemake@output[["kegg_map"]])

# Cnetplot (network pathway genes)
# Usiamo i Symbol per rendere il grafico leggibile
kegg_symbol <- setReadable(kegg_ora, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
cnetplot(kegg_symbol, categorySize="pvalue", foldChange=entrez_ids)

