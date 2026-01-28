# ==== Loading libraries ====
# Suppress logs
pkgs <- c("DESeq2", "org.Hs.eg.db", "ggplot2", "matrixStats", 
          "dplyr", "factoextra", "GGally", "gridExtra")
invisible(lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
}))

# ==== Loading output from previous script ====
input_data <- readRDS(snakemake@input[["rds_import"]])
txi <- input_data$txi
samples <- input_data$metadata

# ==== DESeq2 setup and analysis ====
samples$condition <- factor(ifelse(samples$Treatment == "Untreated", "control", "treated"))
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = samples, 
                                design = ~condition) 

# Pre-filtering: keep genes with at least 10 reads in total
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq analysis
dds <- DESeq(dds)

# Annotation of gene names
rowData(dds)$gene <- mapIds(org.Hs.eg.db, 
                            keys = rownames(dds),
                            column = "SYMBOL", 
                            keytype = "ENSEMBL", 
                            multiVals = "first")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# ==== Exploratory Data Analysis (EDA) ====

# PCA Plot
p_pca <- plotPCA(vsd, intgroup = "condition") + 
  geom_point(size = 1.5) + 
  # Remove 'a' from legend
  guides(
    color = guide_legend(override.aes = list(label = "")),
    shape = guide_legend(override.aes = list(label = ""))
  ) +
  theme_bw(base_size = 5) +
  ggtitle(paste("PCA:", vsd$condition[1], "vs", vsd$condition[2])) +
  geom_text(aes(label = name), vjust = 2) +
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20))
ggsave(snakemake@output[["pca_plot"]], p_pca, width = 4, height = 4, units = "in", dpi = 300)

# Variance explained (Scree plot)
rv <- rowVars(assay(vsd))
# select first 500 components (500 most variables genes)
select <- order(rv, decreasing = T)[seq_len(min(500, length(rv)))]
pca_res <- prcomp(t(assay(vsd)[select, ]))
# explained variance for each principal component
var_expl <- pca_res$sdev^2 / sum(pca_res$sdev^2)
var_perc <- round(var_expl * 100, 1)

var_df <- data.frame(PC = paste0("PC", 1:length(var_perc)), var = var_perc)
var_df$PC <- factor(var_df$PC, levels = var_df$PC)

p_scree <- ggplot(var_df, aes(x = PC, y = var)) + 
                geom_bar(stat = "identity", fill = "steelblue") + 
                geom_text(aes(label = paste0(var, "%")), vjust = -0.5) +
                theme_minimal(base_size = 7) +
                labs(title = "Scree plot", x = "", y = "variance %") +
                scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
                theme(plot.title = element_text(margin = margin(b = 20))) +
                theme_minimal()
ggsave(snakemake@output[["scree_plot"]], p_scree, width = 4, height = 4, units = "in", dpi = 200)

# PCA component pairs plots (first 4 components)
pca_df <- as.data.frame(pca_res$x[, 1:4])
pca_df$condition <- vsd$condition
pca_df$cell.line <- vsd$Cell.Line
#pca_df$sample_name <- colnames(vsd)

p_pca_pairs <- ggpairs(pca_df, columns = 1:4,  
                       aes(color = condition, shape = cell.line), 
                       upper = list(continuous = wrap("cor", size = 2.5)), # correlation text
                       lower = list(continuous = wrap("points", size = 1.5))) + # points 
  theme_minimal(base_size = 7) +
  theme(
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "grey80"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold", colour = "black"), 
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Add samples labels to points
#for(i in 2:4) {
#  for(j in 1:(i-1)) {
#    p_pca_pairs[i,j] <- p_pca_pairs[i,j] + 
#      geom_text(aes(label = sample_name), 
#                vjust = 1.5, size = 2, show.legend = FALSE)
#  }
#}

ggsave(snakemake@output[["pca_pairs"]], p_pca_pairs, width = 6, height = 6, units = "in", dpi = 150)

# PCA loadings ("guide" genes)
loadings <- as.data.frame(pca_res$rotation)
# Assign corresponding gene name to genes selected for PCA
loadings$gene <- rowData(dds)[rownames(loadings), "gene"]
loadings$geneID <- rownames(loadings)
# Handle missing gene names
loadings$gene <- ifelse(is.na(loadings$gene), loadings$geneID, loadings$gene)

# "Positive" guide genes (genes having high expression in samples on the right in PCA plot)
# "Negative" guide genes (genes having high expression in samples on the left in the PCA plot)
top_loadings <- rbind(loadings %>% arrange(desc(PC1)) %>% head(15) , 
                       loadings %>% arrange(PC1) %>% head(15))
                       
p_load <- ggplot(top_loadings, aes(y = reorder(gene, PC1), x = PC1, 
                   fill = PC1 > 0)) + 
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values = c("blue", "red"), guide = "none") +
                  labs(y = "loadings PC1") + theme_bw() + 
                  labs(title = "PCA loadings - Top 15 genes", color = "Condition")
  
#ggsave(snakemake@output[["pca_load"]], p_load)

# Biplot (top 15 most influent genes)
gene_names_pca <- rowData(dds)$gene[select]
gene_names_pca <- rowData(dds)$gene[select]
gene_names_pca[is.na(gene_names_pca)] <- names(gene_names_pca)[is.na(gene_names_pca)]
pca_res_labeled <- pca_res
rownames(pca_res_labeled$rotation) <- gene_names_pca

p_biplot <- fviz_pca_biplot(pca_res_labeled, label = "var",              
                repel = TRUE, select.var = list(contrib = 15), 
                habillage = vsd$condition, palette = "jco",
                geom.ind = "point", pointsize = 3, pointshape = 16,
                arrowsize = 0.5, labelsize = 4) + theme_minimal() +
                labs(title = "PCA Biplot - Top 20 genes", color = "Condition") +
  guides(colour = guide_legend(title = "Condition", override.aes = list(shape = 16)),
         fill = "none", shape = "none")
#ggsave(snakemake@output[["biplot"]], p_biplot)

combined_plot_load <- grid.arrange(p_load, p_biplot, ncol = 2)
ggsave(snakemake@output[["combined_plot_load"]], combined_plot_load, width = 12, height = 6)

# ==== Analysis of DESeq results ====
res <- results(dds, contrast = c("condition", "treated", "control"))
# order results by padj
resOrdered <- res[order(res$padj), ]

# write complete results
write.csv(as.data.frame(resOrdered), snakemake@output[["csv_all"]])

# Filter and save significant genes from Differential Expression Analysis
# at different log2FC thresholds
thresholds <- snakemake@config[["fc_thresholds"]]
output_csvs <- snakemake@output[["csv_sign"]]

for (i in seq_along(thresholds)) {
    fc_val <- thresholds[i]
    # Filtro: padj < 0.05 e Log2FoldChange > threshold
    res_sig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > fc_val)
    write.csv(as.data.frame(res_sig), file = output_csvs[i])
}

# Saving final dds object
saveRDS(dds, snakemake@output[["dds_final"]])

message("Success: Results of Differential Expression Analysis saved in ", snakemake@output[["dds_final"]])
