# ==== Loading libraries ====
# Suppress logs
pkgs <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
          "enrichplot", "pheatmap", "dplyr", "ggplot2", 
          "stringr", "gridExtra", "patchwork", "msigdbr", 
          "ggplotify")
invisible(lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
}))

# ==== Loading previous steps ====
dds <- readRDS(snakemake@input[["dds_results"]])
res <- results(dds, contrast = c("condition", "treated", "control"))
resOrdered <- res[order(res$padj), ]

# ==== Gene Set Enrichment Analysis (GSEA) ====
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Prepare gene list (removing NAs)
res_gsea <- res[!is.na(res$padj) & !is.na(res$entrez) & !is.na(res$pvalue), ]
# Ranked list of genes for GSEA
# Statistics for GSEA: sign(log2FC) * -log10(pvalue) -> nothing enriched
res_gsea$stat_gsea <- sign(res_gsea$log2FoldChange) * -log10(res_gsea$pvalue)
gene_list <- res_gsea$stat_gsea
names(gene_list) <- res_gsea$entrez
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

# New statistics: log2FC
#gene_list <- res_gsea$log2FoldChange
#gene_list
#names(gene_list) <- res_gsea$entrez
#gene_list <- sort(gene_list, decreasing = TRUE)
# Remove entrez duplicates if present
#gene_list <- gene_list[!duplicated(names(gene_list))]

# Check
message("Total number of genes for GSEA: ", length(gene_list))

# gseGO execution
gsea_res <- gseGO(geneList      = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  minGSSize     = 10,       
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.25,    
                  pAdjustMethod = "BH")

message("Total number of rows from GSEA results: ", nrow(gsea_res))

p_titles <- c("Four most significant pathways", "First 10 activated and suppressed pathways", 
              "Ridgeplot", "Genes guiding enrichment in top pathways")

if(nrow(gsea_res) > 0) {
  # Transforms IDs in readable gene names
  gsea_res <- setReadable(gsea_res, OrgDb = org.Hs.eg.db)
  
  # Gsea plot (4 most enriched pathways)
  as_grob <- function(plot_obj) {
    grid::grid.grabExpr(print(plot_obj))
  }
  get_gsea_title <- function(res, i) {
    paste0(str_wrap(res$Description[i], 40), 
           "\n(p-adj: ", formatC(res$p.adjust[i], format = "e", digits = 2), ")")
  }
  p1 <- gseaplot2(gsea_res, geneSetID = 1, title = get_gsea_title(gsea_res, 1), 
                  ES_geom = "line", pvalue_table = FALSE)
  p2 <- gseaplot2(gsea_res, geneSetID = 2, title = get_gsea_title(gsea_res, 2), 
                  ES_geom = "line", pvalue_table = FALSE)
  p3 <- gseaplot2(gsea_res, geneSetID = 3, title = get_gsea_title(gsea_res, 3), 
                  ES_geom = "line", pvalue_table = FALSE)
  p4 <- gseaplot2(gsea_res, geneSetID = 4, title = get_gsea_title(gsea_res, 4), 
                  ES_geom = "line", pvalue_table = FALSE)
  # Convert plots to grob for patchwork
  g1 <- as_grob(p1)
  g2 <- as_grob(p2)
  g3 <- as_grob(p3)
  g4 <- as_grob(p4)
  p_combined <- wrap_elements(g1) + wrap_elements(g2) + 
    wrap_elements(g3) + wrap_elements(g4) + 
    plot_layout(ncol = 2, nrow = 2) +
    theme_minimal(base_size = 15) +
    plot_annotation(title = "Top 4 GSEA Enriched Pathways",
                    theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)))
  ggsave(snakemake@output[["gsea_plot"]], 
         plot = p_combined, 
         width = 14, height = 8, dpi = 150)
  
  # Dotplot split by sign (first 10 activated and suppressed pathways)
  p_dot_gsea <- dotplot(gsea_res, showCategory = 10, split = ".sign") + 
    facet_grid(.~.sign) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
    theme_minimal(base_size = 10) +
    labs(title = "GSEA Analysis Results", subtitle= p_titles[2])
  ggsave(snakemake@output[["gsea_dotplot"]], plot = p_dot_gsea, 
         width = 14, height = 8, dpi = 150)
  
  # Ridgeplot
  gsea_plot_data <- gsea_res
  # filter for pathways with at least 10 genes
  valid_rows <- sapply(strsplit(gsea_plot_data@result$core_enrichment, "/"), length) >= 10
  gsea_plot_data@result <- gsea_plot_data@result[valid_rows, ]
  p_ridge <- ridgeplot(gsea_plot_data, showCategory = 20) + labs(x = "Enrichment distribution") + 
    theme_minimal(base_size = 10) +
    labs(title = "GSEA Analysis Results", subtitle = p_titles[4]) + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
    theme(plot.title = element_text(face = "bold")) +
    geom_vline(xintercept = 0, 
               color = "grey20", 
               linewidth = 0.8)
  ggsave(snakemake@output[["gsea_ridge"]], plot = p_ridge, 
         width = 14, height = 8, dpi = 150)
  
  # Heatmap of core genes (genes guiding the enrichment in first pathway)
  pathway_id <- gsea_res$ID[1]
  core_genes <- strsplit(gsea_res@result[pathway_id, "core_enrichment"], "/")[[1]]
  
  # Transform vst for heatmap
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  # Mapping to gene names
  rownames(vsd_mat) <- mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", rownames(vsd_mat)), 
                              column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  # Filter only core genes for pathway
  counts_core <- vsd_mat[intersect(core_genes, rownames(vsd_mat)), ]
  
  p_heat_core <- pheatmap(counts_core, 
                          scale = "row", 
                          annotation_col = as.data.frame(colData(dds)[, "condition", drop=FALSE]),
                          main = gsea_res$Description[1],
                          silent = TRUE)
  
  p_heat_gg <- as.ggplot(p_heat_core)
  
  # Heatmap of core genes (genes guiding the enrichment in second pathway)
  pathway_id_2 <- gsea_res$ID[2]
  core_genes_2 <- strsplit(gsea_res@result[pathway_id_2, "core_enrichment"], "/")[[1]]
  # Filter only core genes for pathway
  counts_core_2 <- vsd_mat[intersect(core_genes_2, rownames(vsd_mat)), ]
  
  p_heat_core_2 <- pheatmap(counts_core_2, 
                            scale = "row", 
                            annotation_col = as.data.frame(colData(dds)[, "condition", drop=FALSE]),
                            main = gsea_res$Description[2],
                            silent = TRUE)
  
  p_heat_gg_2 <- as.ggplot(p_heat_core_2)
  
  p_combined_heatmaps <- p_heat_gg | p_heat_gg_2
  p_combined_heatmaps <- p_combined_heatmaps + 
    plot_annotation(
      title = "Top 2 GSEA Pathways - Core Enrichment Heatmaps",
      subtitle = "Z-score of normalized counts for leading edge genes",
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
                    plot.subtitle = element_text(size = 15, hjust = 0.5))
    )
  
  ggsave(snakemake@output[["heatmap_core"]], 
         plot = p_combined_heatmaps, 
         width = 18,
         height = 10, 
         dpi = 150)
  
  #ggsave(snakemake@output[["heatmap_core"]], 
  #       plot = p_heat_gg, 
  #       width = 14, height = 8, dpi = 150)
} else {
  # Fallback if no results for GSEA
  message("No GSEA terms enriched. Creating empty placeholder files.")
  
  invisible(lapply(seq_along(output_files), function(i) {
    f <- output_files[i]
    title_text <- p_titles[i]
    
    p_empty <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No significant GSEA enrichment found\n(p-adjust > 0.05)", 
               size = 6, fontface = "italic") + 
      theme_void() + 
      labs(title = paste("GSEA Analysis Results -", title_text)) + 
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", 
                                  margin = margin(t = 60, b = 20))
      )
    
    if (!is.null(snakemake@output[[f]])) {
      ggsave(snakemake@output[[f]], plot = p_empty, width = 8, height = 6)
    }
  }))
}

message("Success: Functional analysis successfully completed")
