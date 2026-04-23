setwd("/Users/okendojo/Downloads/kplots/sexdetermination/variantcalling/heterozygosity/sar4_dosage_out/rna_nadia_coverage/volcanoPlots")

## =========================================================
## Volcano plot with EnhancedVolcano
## Input:
##   - sarblockGene.txt
##   - meta_txt (sample / sex)
## Output:
##   - sex_determination_volcano.png
##   - sex_determination_volcano.pdf
## =========================================================

## ---- packages ----
pkgs <- c("DESeq2", "EnhancedVolcano", "dplyr", "readr", "tibble", "ggplot2")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

## ---- metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

## ---- read count matrix ----
counts_file <- "sarblockGene.txt"   # change this if needed
dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

count_mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Make sure samples are in the same order as the metadata
count_mat <- count_mat[, rownames(meta), drop = FALSE]

## DESeq2 expects integer counts.
## This file contains decimals, so round them before analysis.
storage.mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

## Optional filtering: keep genes with some expression
keep <- rowSums(count_mat) >= 0
count_mat <- count_mat[keep, , drop = FALSE]

## ---- DESeq2 differential expression ----
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ sex
)

dds <- DESeq(dds)

## Compare Male vs Female
res <- results(dds, contrast = c("sex", "Male", "Female"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id")

## Remove genes with missing statistics
res_df <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue))

## ---- choose labels ----
## Label known sex-determination genes if they are present,
## plus the most significant genes in this dataset.
key_genes <- c(
  "LOC103911019", "zgc:113119", "zgc:153116", "psmb10", "ptchd3l", "LOC110437810")

present_key_genes <- res_df$gene_id[toupper(res_df$gene_id) %in% toupper(key_genes)]

top_sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 15) %>%
  pull(gene_id)

label_genes <- unique(c(present_key_genes, top_sig_genes))

## If nothing passes the filter, label the top 15 by adjusted P value
if (length(label_genes) == 0) {
  label_genes <- res_df %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 15) %>%
    pull(gene_id)
}

## ---- EnhancedVolcano plot ----
p <- EnhancedVolcano(
  res_df,
  lab = res_df$gene_id,
  x = "log2FoldChange",
  y = "pvalue",
  selectLab = label_genes,
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Male vs Female",
  subtitle = " ",
  caption = "Positive log2FC = higher in Male; negative log2FC = higher in Female",
  pointSize = 2.5,
  labSize = 5,
  col = c("grey75", "#4E79A7", "#E15759", "#D62728"),
  colAlpha = 0.85,
  drawConnectors = T,
  widthConnectors = 0.5,
  legendPosition = "right",
  legendLabSize = 11,
  legendIconSize = 3
) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(p)

## ---- save high-resolution outputs ----
ggsave(
  filename = "sex_determination_volcano.pdf",
  plot = p,
  width = 9,
  height = 7
)

## ---- optional: write DE results ----
#write.csv(res_df, "DESeq2_Male_vs_Female_results.csv", row.names = FALSE)


## =========================================================
## ComplexHeatmap for sex-labeled samples
## - Shows all genes
## - Groups columns by Female / Male
## - Uses large output size and small row labels for readability
## =========================================================

## ---- packages ----
pkgs <- c("ComplexHeatmap", "circlize", "grid", "tibble", "dplyr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
  library(dplyr)
})

## ---- input metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

## ---- read expression matrix ----
counts_file <- "sarblockGene.txt"   # adjust path if needed

dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Reorder columns to match metadata exactly
mat <- mat[, rownames(meta), drop = FALSE]

## The file contains decimal values, so keep numeric and transform for display
mode(mat) <- "numeric"

## Log-transform for a cleaner display
mat_log <- log2(mat + 1)

## Optional: row z-score scaling so each gene pattern is easier to compare
mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0
mat_z[is.infinite(mat_z)] <- 0

## Order columns by sex for visual grouping
col_order <- order(meta$sex)
mat_z <- mat_z[, col_order, drop = FALSE]
meta <- meta[col_order, , drop = FALSE]

## ---- colors ----
sex_cols <- c(
  Female = "#E76F9A",
  Male   = "#4C78A8"
)

## Continuous heatmap colors
heat_cols <- colorRamp2(
  c(-2, 0, 2),
  c("#2C7BB6", "white", "#D7191C")
)

## ---- annotations ----
ha_top <- HeatmapAnnotation(
  Sex = meta$sex,
  col = list(Sex = sex_cols),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(5, "mm"),
  show_annotation_name = TRUE
)

## ---- make the heatmap ----
ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = heat_cols,
  
  ## sample grouping
  top_annotation = ha_top,
  column_split = meta$sex,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  
  ## gene labels
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "left",
  row_names_max_width = unit(8, "cm"),
  
  ## sample labels
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_side = "bottom",
  
  ## clustering
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  
  ## readability / performance
  use_raster = TRUE,
  raster_quality = 4,
  border = TRUE,
  
  ## titles
  column_title = "Samples grouped by sex",
  row_title = "All genes",
  heatmap_legend_param = list(title = "Row\nZ-score")
)

## ---- draw to screen ----
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

## ---- save high-resolution output ----
## Make the figure tall enough for all genes
out_h <- max(10, 0.18 * nrow(mat_z) + 3)   # dynamic height
out_w <- 11

pdf("sex_grouped_complexheatmap2.pdf", width = out_w, height = out_h)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


#Consider genez with z-score >= 3
## =========================================================
## ComplexHeatmap for sex-labeled samples
## - Keeps only genes with |z-score| >= 2 in at least one sample
## - Groups columns by Female / Male
## - Uses large output size and small row labels for readability
## =========================================================

## ---- packages ----
pkgs <- c("ComplexHeatmap", "circlize", "grid", "tibble", "dplyr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
  library(dplyr)
})

## ---- input metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

## ---- read expression matrix ----
counts_file <- "sarblockGene.txt"   # adjust path if needed

dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Reorder columns to match metadata exactly
mat <- mat[, rownames(meta), drop = FALSE]

## The file contains decimal values, so keep numeric and transform for display
mode(mat) <- "numeric"

## Log-transform for a cleaner display
mat_log <- log2(mat + 1)

## Row z-score scaling so each gene pattern is easier to compare
mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0
mat_z[is.infinite(mat_z)] <- 0

## ---- filter genes by z-score threshold ----
## Keep genes with at least one sample having z-score >= 1 or <= -1
keep_genes <- apply(mat_z, 1, function(x) any(x >= 1.5 | x <= -1.5))

mat_z <- mat_z[keep_genes, , drop = FALSE]

cat("Genes before filtering:", nrow(mat), "\n")
cat("Genes after filtering (|z-score| >= 2):", nrow(mat_z), "\n")

if (nrow(mat_z) == 0) {
  stop("No genes passed the |z-score| >= 2 filter.")
}

## Order columns by sex for visual grouping
col_order <- order(meta$sex)
mat_z <- mat_z[, col_order, drop = FALSE]
meta <- meta[col_order, , drop = FALSE]

## ---- colors ----
sex_cols <- c(
  Female = "#E76F9A",
  Male   = "#4C78A8"
)

heat_cols <- colorRamp2(
  c(-2, 0, 2),
  c("#2C7BB6", "white", "#D7191C")
)

## ---- annotations ----
ha_top <- HeatmapAnnotation(
  Sex = meta$sex,
  col = list(Sex = sex_cols),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(5, "mm"),
  show_annotation_name = TRUE
)

## ---- make the heatmap ----
ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = heat_cols,
  
  ## sample grouping
  top_annotation = ha_top,
  column_split = meta$sex,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  
  ## gene labels
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "left",
  row_names_max_width = unit(8, "cm"),
  
  ## sample labels
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_side = "bottom",
  
  ## clustering
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  
  ## readability / performance
  use_raster = TRUE,
  raster_quality = 4,
  border = TRUE,
  
  ## titles
  column_title = "",
  row_title = "Genes with sigLFC",
  heatmap_legend_param = list(title = "Row\nZ-score")
)

## ---- draw to screen ----
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

## ---- save high-resolution output ----
out_h <- max(10, 0.18 * nrow(mat_z) + 3)
out_w <- 11

pdf("sex_grouped_complexheatmap_filtered2.pdf", width = out_w, height = out_h)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

png("sex_grouped_complexheatmap_filtered.png",
    width = 3200,
    height = as.integer(out_h * 300),
    res = 300)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


#sar4 block stratified by sex
## =========================================================
## ComplexHeatmap with DESeq2-based male/female gene blocks
## - Keeps only genes with |z-score| >= 2 in at least one sample
## - Splits rows into Female-biased / Male-biased / House keeping
## - Adds a row annotation for direction from DESeq2
## - Groups columns by Female / Male
## =========================================================

## ---- packages ----
pkgs <- c(
  "DESeq2", "ComplexHeatmap", "circlize", "grid",
  "tibble", "dplyr"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
  library(dplyr)
})

## ---- input metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

## ---- read expression matrix ----
counts_file <- "sarblockGene.txt"   # adjust path if needed

dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Reorder columns to match metadata exactly
mat <- mat[, rownames(meta), drop = FALSE]

## The file contains decimal values, so keep numeric and transform for display
mode(mat) <- "numeric"

## ---- DESeq2 for direction labels ----
## Use the same samples as the heatmap
dds <- DESeqDataSetFromMatrix(
  countData = round(mat),
  colData   = meta,
  design    = ~ sex
)

dds <- DESeq(dds)

## Male vs Female:
## positive log2FC = higher in Male
## negative log2FC = higher in Female
res <- results(dds, contrast = c("sex", "Male", "Female"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id")

## Keep only genes present in both objects
common_genes <- intersect(rownames(mat), res_df$gene_id)
mat <- mat[common_genes, , drop = FALSE]
res_df <- res_df[match(common_genes, res_df$gene_id), , drop = FALSE]

## ---- transform for heatmap ----
mat_log <- log2(mat + 1)

## Row z-score scaling
mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0
mat_z[is.infinite(mat_z)] <- 0

## ---- filter genes by z-score threshold ----
## Keep genes with at least one sample having z-score >= 2 or <= -2
keep_genes <- apply(mat_z, 1, function(x) any(x >= 1 | x <= -1))

mat_z <- mat_z[keep_genes, , drop = FALSE]
res_df <- res_df[keep_genes, , drop = FALSE]

cat("Genes before filtering:", nrow(mat), "\n")
cat("Genes after filtering (|z-score| >= 2):", nrow(mat_z), "\n")

if (nrow(mat_z) == 0) {
  stop("No genes passed the |z-score| >= 2 filter.")
}

## ---- define gene direction blocks from DESeq2 ----
## You can tighten/loosen thresholds here if needed
lfc_threshold <- 2
padj_threshold <- 0.05

gene_class <- rep("House keeping", nrow(res_df))
gene_class[!is.na(res_df$padj) & res_df$padj < padj_threshold & res_df$log2FoldChange >=  lfc_threshold] <- "Male-biased"
gene_class[!is.na(res_df$padj) & res_df$padj < padj_threshold & res_df$log2FoldChange <= -lfc_threshold] <- "Female-biased"

gene_class <- factor(
  gene_class,
  levels = c("Female-biased", "Male-biased", "House keeping")
)

rownames(gene_class) <- rownames(mat_z)

## Order rows so the blocks appear in a readable order
row_order <- order(gene_class, -abs(res_df$log2FoldChange), na.last = TRUE)
mat_z <- mat_z[row_order, , drop = FALSE]
res_df <- res_df[row_order, , drop = FALSE]
gene_class <- gene_class[row_order]

## ---- order columns by sex for visual grouping ----
col_order <- order(meta$sex)
mat_z <- mat_z[, col_order, drop = FALSE]
meta <- meta[col_order, , drop = FALSE]

## ---- colors ----
sex_cols <- c(
  Female = "#E76F9A",
  Male   = "#4C78A8"
)

class_cols <- c(
  "Female-biased" = "#E76F9A",
  "Male-biased" = "#4C78A8",
  "House keeping" = "grey70"
)

heat_cols <- colorRamp2(
  c(-3, 0, 3),
  c("#2C7BB6", "white", "#D7191C")
)

## ---- annotations ----
ha_top <- HeatmapAnnotation(
  Sex = meta$sex,
  col = list(Sex = sex_cols),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(5, "mm"),
  show_annotation_name = TRUE
)

## Row annotation showing DESeq2 direction
ra_left <- rowAnnotation(
  Class = gene_class,
  col = list(Class = class_cols),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(5, "mm")
)

## ---- make the heatmap ----
ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = heat_cols,
  
  ## group columns by sex
  top_annotation = ha_top,
  column_split = meta$sex,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  
  ## group rows by DESeq2 direction
  row_split = gene_class,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  show_parent_dend_line = TRUE,
  
  ## row annotation
  left_annotation = ra_left,
  
  ## gene labels
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "left",
  row_names_max_width = unit(8, "cm"),
  
  ## sample labels
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_side = "bottom",
  
  ## readability / performance
  use_raster = TRUE,
  raster_quality = 4,
  border = TRUE,
  row_gap = unit(2, "mm"),
  
  ## titles
  column_title = "",
  row_title = "%s",
  heatmap_legend_param = list(title = "Row\nZ-score")
)

## ---- draw to screen ----
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

## ---- save high-resolution output ----
out_h <- max(10, 0.18 * nrow(mat_z) + 4)
out_w <- 12

pdf("sex_grouped_complexheatmap_DESeq2_blocks.pdf", width = out_w, height = out_h)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


## ---- optional export of filtered results ----
write.csv(
  data.frame(
    gene_id = rownames(res_df),
    gene_class = as.character(gene_class),
    log2FoldChange = res_df$log2FoldChange,
    padj = res_df$padj
  ),
  "gene_classification_from_DESeq2.csv",
  row.names = FALSE
)

#Print zero genes
## ---- read expression matrix ----
dat <- read.delim(
  "sarblockGene.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## Convert to matrix with gene_id as rownames
mat <- dat %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

mode(mat) <- "numeric"

## ---- extract genes with all zeros ----
zero_genes <- rownames(mat)[apply(mat, 1, function(x) all(x == 0))]

## Print results
cat("Number of genes with all zero expression:", length(zero_genes), "\n")

## View gene names
print(zero_genes)

## ---- save to file ----
write.table(
  zero_genes,
  file = "all_zero_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#Genes exclusively zero in female
## =========================================================
## Extract genes with:
## - Expression = 0 in ALL males
## - Expression > 0 in at least one female
## =========================================================

## Identify sample groups
male_samples <- rownames(meta)[meta$sex == "Male"]
female_samples <- rownames(meta)[meta$sex == "Female"]

## Subset matrices
mat_male <- mat[, male_samples, drop = FALSE]
mat_female <- mat[, female_samples, drop = FALSE]

## Condition:
## 1. All male values == 0
## 2. At least one female value > 0
female_specific_idx <- apply(mat_male, 1, function(x) all(x > 0)) &
  apply(mat_female, 1, function(x) any(x == 0))

female_specific_genes <- rownames(mat)[female_specific_idx]

## View results
length(female_specific_genes)
head(female_specific_genes)

## Save to file
write.csv(
  data.frame(gene_id = female_specific_genes),
  "female_specific_genes_male_zero.csv",
  row.names = FALSE
)



#Complete R Code: Pearson + Spearman intersection analysis
## =========================================================
## Robust correlation analysis (Pearson + Spearman)
## - Select genes significant in BOTH methods
## - Visualize intersection candidates
## =========================================================

## ---- packages ----
pkgs <- c("dplyr", "tibble", "ggplot2", "ggrepel", "ComplexHeatmap", "circlize", "grid")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## ---- metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(text = meta_txt, sep = "\t")
meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

sex_num <- ifelse(meta$sex == "Male", 1, 0)

## ---- read data ----
dat <- read.delim("sarblockGene.txt", check.names = FALSE)
mat <- dat %>% column_to_rownames("gene_id") %>% as.matrix()
mat <- mat[, rownames(meta)]
mode(mat) <- "numeric"

## ---- transform ----
mat_log <- log2(mat + 1)

## ---- correlation (Pearson + Spearman) ----
cor_res <- lapply(rownames(mat_log), function(g) {
  x <- mat_log[g, ]
  
  p <- suppressWarnings(cor.test(x, sex_num, method = "pearson"))
  s <- suppressWarnings(cor.test(x, sex_num, method = "spearman"))
  
  data.frame(
    gene_id = g,
    pearson_cor = p$estimate,
    pearson_p = p$p.value,
    spearman_cor = s$estimate,
    spearman_p = s$p.value
  )
}) %>% bind_rows()

## Adjust p-values
cor_df <- cor_res %>%
  mutate(
    pearson_padj = p.adjust(pearson_p, method = "BH"),
    spearman_padj = p.adjust(spearman_p, method = "BH")
  )

## ---- define strict candidates (intersection) ----
cor_thresh <- 0.7
padj_thresh <- 0.05

cor_df <- cor_df %>%
  mutate(
    pearson_sig = abs(pearson_cor) >= cor_thresh & pearson_padj < padj_thresh,
    spearman_sig = abs(spearman_cor) >= cor_thresh & spearman_padj < padj_thresh,
    consensus = pearson_sig & spearman_sig,
    direction = case_when(
      pearson_cor > 0 & spearman_cor > 0 ~ "Male-associated",
      pearson_cor < 0 & spearman_cor < 0 ~ "Female-associated",
      TRUE ~ "Mixed"
    )
  )

## Extract strict candidates
candidates <- cor_df %>% filter(consensus)

#write.csv(cor_df, "correlation_all_results.csv", row.names = FALSE)
#write.csv(candidates, "correlation_strict_candidates.csv", row.names = FALSE)

cat("Number of strict candidates:", nrow(candidates), "\n")

## ---- plot 1: Pearson vs Spearman agreement ----
p1 <- ggplot(cor_df, aes(x = pearson_cor, y = spearman_cor)) +
  geom_point(aes(color = consensus), alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = " ")

ggsave("pearson_vs_spearman.pdf", p1, width = 6, height = 5)

## ---- plot 2: volcano-style correlation plot ----
p2 <- ggplot(cor_df, aes(x = pearson_cor, y = -log10(pearson_padj))) +
  geom_point(aes(color = consensus), alpha = 0.7) +
  geom_vline(xintercept = c(-cor_thresh, cor_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
  theme_classic() +
  labs(title = "Correlation with sex (Pearson)")

ggsave("correlation_volcano.pdf", p2, width = 7, height = 5)

## ---- heatmap of strict candidates ----
if (nrow(candidates) > 1) {
  
  genes <- candidates$gene_id
  mat_hm <- mat_log[genes, ]
  
  ## z-score
  mat_z <- t(scale(t(mat_hm)))
  mat_z[is.na(mat_z)] <- 0
  
  ## order samples
  ord <- order(meta$sex)
  mat_z <- mat_z[, ord]
  meta_ord <- meta[ord, ]
  
  ## colors
  sex_cols <- c(Female = "#E76F9A", Male = "#4C78A8")
  
  ha <- HeatmapAnnotation(
    Sex = meta_ord$sex,
    col = list(Sex = sex_cols)
  )
  
  ht <- Heatmap(
    mat_z,
    name = "Z",
    top_annotation = ha,
    column_split = meta_ord$sex,
    col = colorRamp2(c(-2,0,2), c("blue","white","red")),
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_rot = 45
  )
  
  pdf("strict_candidates_heatmap.pdf", width = 10, height = 10)
  draw(ht)
  dev.off()
}

## ---- label top genes ----
top_genes <- candidates %>%
  arrange(desc(abs(pearson_cor))) %>%
  slice_head(n = 30)

p3 <- ggplot(cor_df, aes(x = pearson_cor, y = -log10(pearson_padj))) +
  geom_point(color = "grey70") +
  geom_point(data = candidates, color = "red") +
  ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = gene_id),
    size = 3
  ) +
  theme_classic() +
  labs(title = "Top consensus sex-associated genes")

ggsave("top_candidates_labeled.pdf", p3, width = 7, height = 5)


#build a regulatory network (co-expression graph)
## =========================================================
## Sex-associated co-expression network
## - Uses Pearson + Spearman correlation with sex
## - Keeps intersection as stricter candidate genes
## - Builds a co-expression graph among candidates
## - Saves node/edge tables and a network plot
## =========================================================

## ---- packages ----
pkgs <- c(
  "dplyr", "tibble", "ggplot2", "ggrepel",
  "igraph", "ggraph", "tidygraph", "scales"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(scales)
})

## ---- metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

sex_num <- ifelse(meta$sex == "Male", 1, 0)
names(sex_num) <- rownames(meta)

## ---- read counts matrix ----
counts_file <- "sarblockGene.txt"   # adjust path if needed

dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Keep only the known samples, in the same order
mat <- mat[, rownames(meta), drop = FALSE]
mode(mat) <- "numeric"

## ---- transform ----
mat_log <- log2(mat + 1)

## =========================================================
## 1) Gene-sex correlation analysis
## =========================================================
cor_df <- lapply(rownames(mat_log), function(g) {
  x <- as.numeric(mat_log[g, ])
  
  p <- suppressWarnings(cor.test(x, sex_num, method = "pearson"))
  s <- suppressWarnings(cor.test(x, sex_num, method = "spearman"))
  
  data.frame(
    gene_id = g,
    pearson_cor = unname(p$estimate),
    pearson_p = p$p.value,
    spearman_cor = unname(s$estimate),
    spearman_p = s$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

cor_df <- cor_df %>%
  mutate(
    pearson_padj = p.adjust(pearson_p, method = "BH"),
    spearman_padj = p.adjust(spearman_p, method = "BH"),
    pearson_sig = abs(pearson_cor) >= 0.70 & pearson_padj < 0.05,
    spearman_sig = abs(spearman_cor) >= 0.70 & spearman_padj < 0.05,
    consensus = pearson_sig & spearman_sig,
    direction = case_when(
      pearson_cor > 0 & spearman_cor > 0 ~ "Male-associated",
      pearson_cor < 0 & spearman_cor < 0 ~ "Female-associated",
      TRUE ~ "Mixed"
    )
  )

## Stricter candidate genes = intersection of Pearson + Spearman hits
candidates <- cor_df %>%
  filter(consensus) %>%
  arrange(desc(abs(pearson_cor) + abs(spearman_cor)))

## Fallback if too few genes pass strict cutoff
if (nrow(candidates) < 5) {
  candidates <- cor_df %>%
    arrange(desc(abs(pearson_cor) + abs(spearman_cor))) %>%
    slice_head(n = 30)
}

#write.csv(cor_df, "sex_correlation_all_genes.csv", row.names = FALSE)
#write.csv(candidates, "sex_correlation_consensus_candidates.csv", row.names = FALSE)

cat("Consensus candidate genes:", nrow(candidates), "\n")

## =========================================================
## 2) Build co-expression network among candidate genes
## =========================================================
cand_genes <- candidates$gene_id
expr_sub <- mat_log[cand_genes, , drop = FALSE]

## Gene-gene correlation across samples
pear_mat <- cor(t(expr_sub), method = "pearson", use = "pairwise.complete.obs")
spear_mat <- cor(t(expr_sub), method = "spearman", use = "pairwise.complete.obs")

## Build edge list from upper triangle
idx <- which(upper.tri(pear_mat), arr.ind = TRUE)

edges <- data.frame(
  from = rownames(pear_mat)[idx[, 1]],
  to   = colnames(pear_mat)[idx[, 2]],
  pearson = pear_mat[idx],
  spearman = spear_mat[idx],
  stringsAsFactors = FALSE
) %>%
  mutate(
    same_direction = sign(pearson) == sign(spearman),
    abs_mean_cor = (abs(pearson) + abs(spearman)) / 2,
    edge_sign = case_when(
      pearson > 0 & spearman > 0 ~ "Positive",
      pearson < 0 & spearman < 0 ~ "Negative",
      TRUE ~ "Mixed"
    )
  ) %>%
  filter(
    same_direction,
    abs(pearson) >= 0.85,
    abs(spearman) >= 0.85
  ) %>%
  arrange(desc(abs_mean_cor))

## If the network is empty, relax thresholds
if (nrow(edges) == 0) {
  message("No edges passed the strict threshold; relaxing to |cor| >= 0.75")
  edges <- data.frame(
    from = rownames(pear_mat)[idx[, 1]],
    to   = colnames(pear_mat)[idx[, 2]],
    pearson = pear_mat[idx],
    spearman = spear_mat[idx],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      same_direction = sign(pearson) == sign(spearman),
      abs_mean_cor = (abs(pearson) + abs(spearman)) / 2,
      edge_sign = case_when(
        pearson > 0 & spearman > 0 ~ "Positive",
        pearson < 0 & spearman < 0 ~ "Negative",
        TRUE ~ "Mixed"
      )
    ) %>%
    filter(
      same_direction,
      abs(pearson) >= 0.75,
      abs(spearman) >= 0.75
    ) %>%
    arrange(desc(abs_mean_cor))
}

write.csv(edges, "sex_coexpression_edges.csv", row.names = FALSE)

## =========================================================
## 3) Node table with sex association and degree
## =========================================================
nodes <- candidates %>%
  select(
    gene_id,
    pearson_cor, pearson_padj,
    spearman_cor, spearman_padj,
    consensus, direction
  ) %>%
  rename(name = gene_id)

if (nrow(edges) > 0) {
  g_tmp <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  nodes$degree <- degree(g_tmp)
} else {
  nodes$degree <- 0
}

#write.csv(nodes, "sex_coexpression_nodes.csv", row.names = FALSE)

## =========================================================
## 4) Plot network
## =========================================================
if (nrow(edges) > 0) {
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  ## Keep only the largest connected component for a cleaner plot
  comps <- components(g)
  largest_comp <- which.max(comps$csize)
  g <- induced_subgraph(g, vids = which(comps$membership == largest_comp))
  
  ## Recompute degree after subsetting
  V(g)$degree <- degree(g)
  
  ## Node aesthetics
  V(g)$node_group <- ifelse(V(g)$direction == "Male-associated", "Male-associated",
                            ifelse(V(g)$direction == "Female-associated", "Female-associated", "Neutral"))
  
  V(g)$node_size <- rescale(V(g)$degree, to = c(3, 10))
  
  ## Label the strongest hubs
  hub_cutoff <- if (vcount(g) >= 10) {
    sort(V(g)$degree, decreasing = TRUE)[min(10, vcount(g))]
  } else {
    1
  }
  
  V(g)$label <- ifelse(V(g)$degree >= hub_cutoff, V(g)$name, "")
  
  node_cols <- c(
    "Female-associated" = "#E76F9A",
    "Male-associated"   = "#4C78A8",
    "Neutral"           = "grey70"
  )
  
  edge_cols <- c(
    "Positive" = "#B2182B",
    "Negative" = "#2166AC"
  )
  
  p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(width = abs_mean_cor, color = edge_sign),
      alpha = 0.6,
      lineend = "round"
    ) +
    geom_node_point(
      aes(size = node_size, color = node_group),
      alpha = 0.95
    ) +
    geom_node_text(
      aes(label = label),
      repel = TRUE,
      size = 3,
      max.overlaps = Inf
    ) +
    scale_edge_color_manual(values = edge_cols) +
    scale_color_manual(values = node_cols) +
    scale_edge_width(range = c(0.3, 2.5)) +
    scale_size_identity() +
    theme_void(base_size = 12) +
    labs(
      title = "Sex-associated co-expression network",
      subtitle = "Consensus genes supported by Pearson + Spearman correlation with sex"
    )
  
  ggsave("sex_coexpression_network.pdf", p_net, width = 11, height = 8.5)
  #ggsave("sex_coexpression_network.png", p_net, width = 11, height = 8.5, dpi = 300)
  
  print(p_net)
} else {
  message("No network edges passed the filtering threshold.")
}

## =========================================================
## 5) Optional: simple candidate summary
## =========================================================
summary_tbl <- candidates %>%
  mutate(
    abs_consensus = (abs(pearson_cor) + abs(spearman_cor)) / 2
  ) %>%
  arrange(desc(abs_consensus))

write.csv(summary_tbl, "sex_coexpression_candidate_summary.csv", row.names = FALSE)

cat("\nTop candidate genes:\n")
print(summary_tbl %>% select(gene_id, pearson_cor, spearman_cor, direction) %>% head(20))


#Weighted Gene Co-expression Network Analysis (WGCNA)
## =========================================================
## WGCNA for sex-determination candidate genes
## Input:
##   - sarblockGene.txt (gene x sample count matrix)
##   - meta_txt (sample / sex)
##
## Output:
##   - WGCNA diagnostic plots
##   - module-trait heatmap
##   - module membership vs gene significance plots
##   - hub gene table
## =========================================================

## ---- packages ----
pkgs <- c("WGCNA", "dplyr", "tibble", "ggplot2", "ggrepel")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(impute)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## ---- metadata ----
meta_txt <- "
sample\tsex
SRR8176748_ZW\tFemale
SRR8176749_ZZ\tMale
SRR8176750_ZZ\tMale
SRR8176751_ZW\tFemale
SRR8176752_ZZ\tMale
SRR8176755_ZW\tFemale
SRR8176757_ZZ\tMale
SRR8176758_ZZ\tMale
SRR8176759_ZW\tFemale
SRR8176760_ZW\tFemale
"

meta <- read.delim(
  text = meta_txt,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
rownames(meta) <- meta$sample

## Encode sex as 0/1 for module-trait correlation
sex_num <- ifelse(meta$sex == "Male", 1, 0)
names(sex_num) <- rownames(meta)

## ---- read count matrix ----
counts_file <- "sarblockGene.txt"   # adjust path if needed

dat <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("gene_id" %in% colnames(dat))

mat <- dat %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## Keep only samples in metadata, in the same order
mat <- mat[, rownames(meta), drop = FALSE]
mode(mat) <- "numeric"

## ---- transform ----
## WGCNA works well on transformed expression values.
## With count-like data, log2(x + 1) is a reasonable starting point.
expr_gene_by_sample <- log2(mat + 1)

## WGCNA expects samples in rows and genes in columns
datExpr0 <- t(expr_gene_by_sample)

## Optional: check for missing values / bad genes
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## Reorder metadata to match WGCNA matrix rows
meta2 <- meta[rownames(datExpr0), , drop = FALSE]
sex_num2 <- sex_num[rownames(datExpr0)]

cat("Samples:", nrow(datExpr0), "\n")
cat("Genes:", ncol(datExpr0), "\n")

## ---- sample clustering ----
sampleTree <- hclust(dist(datExpr0), method = "average")
pdf("wgcna_sample_clustering.pdf", width = 8, height = 5)
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex = 0.9)
abline(h = 50, col = "red", lty = 2)
dev.off()

## ---- choose soft-threshold power ----
powers <- c(1:10, seq(12, 20, 2))

sft <- pickSoftThreshold(
  datExpr0,
  powerVector = powers,
  networkType = "signed",
  corFnc = "cor",
  corOptions = list(use = "p"),
  verbose = 5
)

pdf("wgcna_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "blue"
)
abline(h = 0.8, col = "red", lty = 2)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "blue"
)
dev.off()

## Pick a power:
## first power with signed R^2 >= 0.8; if none, fall back to 6
fitTbl <- sft$fitIndices
fitTbl$R2 <- -sign(fitTbl[, 3]) * fitTbl[, 2]

softPower <- fitTbl$Power[which(fitTbl$R2 >= 0.8)[1]]
if (is.na(softPower)) softPower <- 6

cat("Chosen soft power:", softPower, "\n")

## ---- network construction and module detection ----
## With few samples, keep modules small and avoid over-merging.
net <- blockwiseModules(
  datExpr0,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  maxBlockSize = ncol(datExpr0) + 1,
  minModuleSize = 5,
  deepSplit = 2,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "WGCNA_TOM",
  verbose = 3
)

moduleColors <- labels2colors(net$colors)

## ---- dendrogram with module colors ----
pdf("wgcna_gene_dendrogram_modules.pdf", width = 12, height = 7)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

## ---- module eigengenes ----
MEs0 <- moduleEigengenes(datExpr0, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

## ---- module-trait relationship ----
moduleTraitCor <- cor(MEs, sex_num2, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr0))

## Save module-trait table
moduleTraitTbl <- data.frame(
  module = rownames(moduleTraitCor),
  cor = as.numeric(moduleTraitCor[, 1]),
  pvalue = as.numeric(moduleTraitPvalue[, 1]),
  padj = p.adjust(as.numeric(moduleTraitPvalue[, 1]), method = "BH")
)
write.csv(moduleTraitTbl, "wgcna_module_trait_table.csv", row.names = FALSE)

## Heatmap plot
pdf("wgcna_module_trait_heatmap.pdf", width = 7, height = 8)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Male (1) vs Female (0)",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = paste0(
    signif(moduleTraitCor, 2), "\n(",
    signif(moduleTraitPvalue, 1), ")"
  ),
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

## ---- gene significance and module membership ----
## Gene significance for sex
GS.sex <- as.numeric(cor(datExpr0, sex_num2, use = "p"))
names(GS.sex) <- colnames(datExpr0)

## Module membership
MEs_num <- as.data.frame(MEs)
moduleMembership <- as.data.frame(cor(datExpr0, MEs_num, use = "p"))
moduleMembershipP <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples = nrow(datExpr0)))

## Build a gene-level summary table
geneSummary <- data.frame(
  gene_id = colnames(datExpr0),
  module = moduleColors,
  GS_sex = GS.sex,
  stringsAsFactors = FALSE
)

## Add top module membership value for each gene
mm_max <- apply(moduleMembership, 1, function(x) max(abs(x), na.rm = TRUE))
geneSummary$MM_max <- mm_max

## Identify hub genes per module
hubGenes <- chooseTopHubInEachModule(datExpr0, moduleColors, power = softPower, type = "signed")
hubTbl <- data.frame(
  module = names(hubGenes),
  hub_gene = as.vector(hubGenes),
  stringsAsFactors = FALSE
)

#write.csv(geneSummary, "wgcna_gene_summary.csv", row.names = FALSE)
#write.csv(hubTbl, "wgcna_hub_genes.csv", row.names = FALSE)

## ---- plot GS vs MM for the most sex-associated module ----
if (nrow(moduleTraitTbl) > 0) {
  bestModule <- moduleTraitTbl %>%
    arrange(desc(abs(cor))) %>%
    slice(1) %>%
    pull(module)
  
  bestModuleColor <- sub("^ME", "", bestModule)
  
  inModule <- moduleColors == bestModuleColor
  mm_use <- as.numeric(cor(datExpr0[, inModule, drop = FALSE], MEs[, bestModule, drop = FALSE], use = "p"))
  gs_use <- as.numeric(cor(datExpr0[, inModule, drop = FALSE], sex_num2, use = "p"))
  
  plotDF <- data.frame(
    gene = colnames(datExpr0)[inModule],
    MM = mm_use,
    GS = gs_use
  )
  
  p_mmgs <- ggplot(plotDF, aes(x = MM, y = GS)) +
    geom_point(alpha = 0.8) +
    ggrepel::geom_text_repel(
      data = subset(plotDF, abs(MM) > 0.8 | abs(GS) > 0.8),
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Module membership vs gene significance: ", bestModule),
      x = "Module membership (kME)",
      y = "Gene significance for sex"
    )
  
  ggsave("wgcna_MM_vs_GS_best_module.pdf", p_mmgs, width = 7, height = 6)
  ggsave("wgcna_MM_vs_GS_best_module.png", p_mmgs, width = 7, height = 6, dpi = 300)
  print(p_mmgs)
}

## ---- optionally export module colors for all genes ----
moduleColorTbl <- data.frame(
  gene_id = colnames(datExpr0),
  module_color = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(moduleColorTbl, "wgcna_module_colors.csv", row.names = FALSE)

cat("\nWGCNA completed.\n")
cat("Soft power used:", softPower, "\n")
cat("Modules found:", paste(unique(moduleColors), collapse = ", "), "\n")
