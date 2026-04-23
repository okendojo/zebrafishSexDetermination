#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(vcfR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Inputs
# ============================================================
bed_file <- "genes_on_NC_133179.1.bed"
vcf_file <- "merged.samples_snps.vcf.gz"
outdir <- "sex_linked_snp_density_results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "vcf_filtered"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Fixed metadata
# ============================================================
meta_txt <- "
sample\tWT_line\tsex\thaps
10WW\tNadia\tFemale\tWW
1ZZ\tNadia\tMale\tZZ
2ZZ\tNadia\tMale\tZZ
3ZZ\tNadia\tMale\tZZ
4ZZ\tNadia\tMale\tZZ
6WW\tNadia\tFemale\tWW
8WW\tNadia\tFemale\tWW
9WW\tNadia\tFemale\tWW
CBW1_4\tCB\tMale\tWW
CBZW_1\tCB\tFemale\tZW
CBZW_11\tCB\tFemale\tZW
CBZW_12\tCB\tFemale\tZW
CBZW_9\tCB\tFemale\tZW
W2_13\tCB\tFemale\tZW
W2_14\tCB\tFemale\tZW
W2_15\tCB\tFemale\tZW
W2_16\tCB\tFemale\tZW
ww_cbw1_1\tCB\tMale\tWW
ww_cbw1_3\tCB\tMale\tWW
"

meta <- fread(text = meta_txt)
stopifnot(nrow(meta) == 19)

# ============================================================
# Distinctive color palettes
# ============================================================
sex_cols <- c(
  Female = "#D81B60",
  Male   = "#1E88E5"
)

comp_cols <- c(
  "Nadia: Male - Female" = "#E69F00",
  "CB: Male - Female"    = "#0072B2"
)

line_cols <- c(
  Nadia = "#6A3D9A",
  CB    = "#33A02C"
)

hap_cols <- c(
  WW = "#4DAF4A",
  ZW = "#984EA3",
  ZZ = "#FF7F00"
)

# ============================================================
# Read BED and standardize intervals
# BED is assumed to be 0-based, half-open.
# ============================================================
bed <- fread(bed_file, header = FALSE)
if (ncol(bed) < 3) stop("BED file must have at least 3 columns: chrom, start, end.")

setnames(bed, 1:3, c("chrom", "start0", "end"))

if (ncol(bed) >= 4) {
  setnames(bed, 4, "gene_id")
} else {
  bed[, gene_id := paste0(chrom, ":", start0, "-", end)]
}

bed[, `:=`(
  start = as.integer(start0) + 1L,
  end = as.integer(end),
  len_bp = as.integer(end) - as.integer(start0),
  mid_kb = ((as.integer(start0) + as.integer(end)) / 2) / 1000
)]

bed[, interval_key := paste0(gene_id, " | ", chrom, ":", start, "-", end)]
bed <- bed[, .(chrom, start, end, gene_id, interval_key, len_bp, mid_kb)]
bed[, interval_idx := .I]

# ============================================================
# Filter VCF to BED intervals using vcftools
# ============================================================
filtered_prefix <- file.path(outdir, "vcf_filtered", "genes_only")
filtered_vcf <- paste0(filtered_prefix, ".recode.vcf")

if (!file.exists(filtered_vcf)) {
  message("Running vcftools BED filtering...")
  args <- c(
    "--gzvcf", vcf_file,
    "--bed", bed_file,
    "--recode",
    "--recode-INFO-all",
    "--out", filtered_prefix
  )
  system2("vcftools", args = args)
  if (!file.exists(filtered_vcf)) {
    stop("vcftools did not produce the expected filtered VCF: ", filtered_vcf)
  }
}

# ============================================================
# Read filtered VCF
# ============================================================
vcf <- read.vcfR(filtered_vcf, verbose = FALSE)
fixed <- as.data.table(vcf@fix)

variants <- data.table(
  var_idx = seq_len(nrow(fixed)),
  chrom = fixed$CHROM,
  pos = as.integer(fixed$POS)
)

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
vcf_samples <- colnames(gt)
keep_samples <- intersect(meta$sample, vcf_samples)

if (length(keep_samples) == 0) stop("No metadata samples were found in the VCF.")

meta <- meta[match(keep_samples, meta$sample), ]
if (anyNA(meta$sample)) stop("Some keep_samples are missing in metadata.")
if (!all(meta$sample == keep_samples)) stop("Metadata ordering mismatch after reordering.")

# ============================================================
# Map variants to BED intervals using non-equi join
# ============================================================
hits <- bed[variants,
  on = .(chrom, start <= pos, end >= pos),
  nomatch = 0,
  allow.cartesian = TRUE
]

variant_to_interval <- unique(hits[, .(var_idx, interval_idx)])

map_interval <- rep(NA_integer_, nrow(fixed))
map_interval[variant_to_interval$var_idx] <- variant_to_interval$interval_idx

# ============================================================
# Extract genotypes and count SNPs per sample per interval
# ============================================================
gt <- gt[, keep_samples, drop = FALSE]

is_variant <- !is.na(gt) & gt != "./." & gt != ".|." & !grepl("^0[\\/|]0$", gt)

count_mat <- matrix(
  0L,
  nrow = nrow(bed),
  ncol = length(keep_samples),
  dimnames = list(bed$interval_key, keep_samples)
)

for (j in seq_along(keep_samples)) {
  idx <- map_interval[which(is_variant[, j])]
  idx <- idx[!is.na(idx)]
  if (length(idx)) {
    count_mat[, j] <- tabulate(idx, nbins = nrow(bed))
  }
}

# ============================================================
# Sample-level SNP density table
# ============================================================
sample_density <- as.data.table(as.table(count_mat))
setnames(sample_density, c("interval_key", "sample", "snp_count"))
sample_density[, snp_count := as.integer(snp_count)]
sample_density <- merge(sample_density, bed, by = "interval_key", all.x = TRUE)
sample_density <- merge(sample_density, meta, by = "sample", all.x = TRUE)
sample_density[, snp_density := snp_count / (len_bp / 1000)]

fwrite(sample_density, file.path(outdir, "sample_interval_snp_density.tsv"), sep = "\t")

# ============================================================
# Helper functions
# ============================================================
safe_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
}

safe_kruskal_p <- function(df, group_col, value_col) {
  g <- df[[group_col]]
  x <- df[[value_col]]
  keep <- is.finite(x) & !is.na(g)
  g <- droplevels(factor(g[keep]))
  x <- x[keep]
  if (nlevels(g) < 2) return(NA_real_)
  tryCatch(kruskal.test(x ~ g)$p.value, error = function(e) NA_real_)
}

# ============================================================
# Interval summaries
# ============================================================
line_sex_summary <- sample_density[, .(
  n_samples = .N,
  mean_density = mean(snp_density),
  median_density = median(snp_density),
  sd_density = sd(snp_density),
  total_snps = sum(snp_count)
), by = .(interval_idx, interval_key, gene_id, chrom, start, end, len_bp, mid_kb, WT_line, sex)]

haps_summary <- sample_density[, .(
  n_samples = .N,
  mean_density = mean(snp_density),
  median_density = median(snp_density),
  sd_density = sd(snp_density),
  total_snps = sum(snp_count)
), by = .(interval_idx, interval_key, gene_id, chrom, start, end, len_bp, mid_kb, haps)]

# ============================================================
# Per-interval testing
# ============================================================
stats_list <- lapply(split(sample_density, by = "interval_idx"), function(d) {
  d <- as.data.table(d)
  nadia <- d[WT_line == "Nadia"]
  cb <- d[WT_line == "CB"]

  male_nadia <- nadia[sex == "Male", snp_density]
  female_nadia <- nadia[sex == "Female", snp_density]
  male_cb <- cb[sex == "Male", snp_density]
  female_cb <- cb[sex == "Female", snp_density]

  haps_means <- d[, .(mean_density = mean(snp_density)), by = haps]
  haps_w <- haps_means[haps == "WW", mean_density]
  haps_zw <- haps_means[haps == "ZW", mean_density]
  haps_zz <- haps_means[haps == "ZZ", mean_density]

  d[1, .(
    interval_idx = interval_idx,
    interval_key = interval_key,
    gene_id = gene_id,
    chrom = chrom,
    start = start,
    end = end,
    len_bp = len_bp,
    mid_kb = mid_kb,

    mean_density_all = mean(snp_density),
    median_density_all = median(snp_density),

    mean_male_nadia = if (length(male_nadia)) mean(male_nadia) else NA_real_,
    mean_female_nadia = if (length(female_nadia)) mean(female_nadia) else NA_real_,
    delta_male_female_nadia = if (length(male_nadia) && length(female_nadia)) mean(male_nadia) - mean(female_nadia) else NA_real_,
    p_male_female_nadia = safe_wilcox_p(male_nadia, female_nadia),

    mean_male_cb = if (length(male_cb)) mean(male_cb) else NA_real_,
    mean_female_cb = if (length(female_cb)) mean(female_cb) else NA_real_,
    delta_male_female_cb = if (length(male_cb) && length(female_cb)) mean(male_cb) - mean(female_cb) else NA_real_,
    p_male_female_cb = safe_wilcox_p(male_cb, female_cb),

    mean_WW = if (length(haps_w)) mean(haps_w) else NA_real_,
    mean_ZW = if (length(haps_zw)) mean(haps_zw) else NA_real_,
    mean_ZZ = if (length(haps_zz)) mean(haps_zz) else NA_real_,

    p_haps = safe_kruskal_p(d, "haps", "snp_density"),
    p_sex = safe_kruskal_p(d, "sex", "snp_density"),
    p_line = safe_kruskal_p(d, "WT_line", "snp_density")
  )]
})

stats <- rbindlist(stats_list, fill = TRUE)
stats[, q_male_female_nadia := p.adjust(p_male_female_nadia, method = "BH")]
stats[, q_male_female_cb := p.adjust(p_male_female_cb, method = "BH")]
stats[, q_haps := p.adjust(p_haps, method = "BH")]

stats[, consistent_direction := !is.na(delta_male_female_nadia) &
        !is.na(delta_male_female_cb) &
        sign(delta_male_female_nadia) == sign(delta_male_female_cb)]

stats[, sex_linked_flag := (
  !is.na(q_male_female_nadia) & q_male_female_nadia < 0.05 &
  !is.na(q_male_female_cb) & q_male_female_cb < 0.05 &
  consistent_direction
)]

stats[, score := (
  -log10(pmax(q_male_female_nadia, 1e-300)) +
  -log10(pmax(q_male_female_cb, 1e-300)) +
  abs(delta_male_female_nadia) +
  abs(delta_male_female_cb)
)]

setorder(stats, -score)

fwrite(stats, file.path(outdir, "interval_stats.tsv"), sep = "\t")
fwrite(haps_summary, file.path(outdir, "haps_interval_summary.tsv"), sep = "\t")
fwrite(line_sex_summary, file.path(outdir, "line_sex_interval_summary.tsv"), sep = "\t")

# ============================================================
# Plot data
# ============================================================
scatter_df <- rbindlist(list(
  stats[, .(
    WT_line = "Nadia",
    female_density = mean_female_nadia,
    male_density = mean_male_nadia,
    gene_id,
    interval_key
  )],
  stats[, .(
    WT_line = "CB",
    female_density = mean_female_cb,
    male_density = mean_male_cb,
    gene_id,
    interval_key
  )]
), fill = TRUE)

# Haplotype summary for line plot with CI
hap_line_summary <- sample_density[, .(
  mean_density = mean(snp_density),
  sd_density = if (.N > 1) sd(snp_density) else 0,
  n_samples = .N
), by = .(WT_line, haps, mid_kb, interval_idx, interval_key, gene_id, chrom)]

hap_line_summary[, se_density := ifelse(n_samples > 1, sd_density / sqrt(n_samples), 0)]
hap_line_summary[, ci95 := 1.96 * se_density]
hap_line_summary[, ymin := mean_density - ci95]
hap_line_summary[, ymax := mean_density + ci95]

# ============================================================
# Plots
# ============================================================

# 1) Mean SNP density by sex across intervals
p1 <- ggplot(
  line_sex_summary,
  aes(x = mid_kb, y = mean_density, color = sex, group = sex)
) +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~ WT_line, ncol = 1, scales = "free_x") +
  scale_color_manual(values = sex_cols) +
  labs(
    x = "Genomic position (kb)",
    y = "Mean SNP density (SNPs / kb)",
    color = "Sex",
    title = "SNP density across BED-restricted gene intervals"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

# 2) Sex-associated SNP density differences
delta_df <- stats %>%
  select(interval_key, gene_id, chrom, start, end, mid_kb,
         delta_male_female_nadia, delta_male_female_cb) %>%
  pivot_longer(
    cols = c(delta_male_female_nadia, delta_male_female_cb),
    names_to = "comparison",
    values_to = "delta_density"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      delta_male_female_nadia = "Nadia: Male - Female",
      delta_male_female_cb = "CB: Male - Female"
    )
  )

annot_df <- stats[1:min(.N, 10), .(
  gene_id,
  mid_kb,
  delta_male_female_nadia,
  delta_male_female_cb
)]

annot_long <- melt(
  annot_df,
  id.vars = c("gene_id", "mid_kb"),
  variable.name = "comparison",
  value.name = "delta_density"
)

annot_long[, comparison := recode(
  comparison,
  delta_male_female_nadia = "Nadia: Male - Female",
  delta_male_female_cb = "CB: Male - Female"
)]

p2 <- ggplot(delta_df, aes(x = mid_kb, y = delta_density, color = comparison)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey40") +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~ comparison, ncol = 1, scales = "free_x") +
  scale_color_manual(values = comp_cols, guide = "none") +
  geom_text_repel(
    data = annot_long,
    aes(x = mid_kb, y = delta_density, label = gene_id),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = 20,
    box.padding = 0.2,
    point.padding = 0.2
  ) +
  labs(
    x = "Genomic position (kb)",
    y = "Male - Female mean SNP density",
    title = "Sex-associated SNP density differences"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

# 3) Female vs male scatter by WT line
p3 <- ggplot(scatter_df, aes(x = female_density, y = male_density, color = WT_line)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5, color = "grey45") +
  geom_point(size = 2.0, alpha = 0.88) +
  facet_wrap(~ WT_line) +
  scale_color_manual(values = line_cols) +
  labs(
    x = "Female mean SNP density (SNPs / kb)",
    y = "Male mean SNP density (SNPs / kb)",
    color = "WT line",
    title = "Female vs male SNP density scatter"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# 4) Haplotype line plot with 95% confidence interval
p3_haps <- ggplot(
  hap_line_summary,
  aes(x = mid_kb, y = mean_density, color = haps, fill = haps, group = haps)
) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~ WT_line, ncol = 1, scales = "free_x") +
  scale_color_manual(values = hap_cols, drop = FALSE) +
  scale_fill_manual(values = hap_cols, drop = FALSE) +
  labs(
    x = "Genomic position (kb)",
    y = "Mean SNP density (SNPs / kb)",
    color = "Haplotype",
    fill = "Haplotype",
    title = "Haplotype-specific SNP density with 95% confidence interval"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

combined <- p1 / p2 / p3 / p3_haps + plot_layout(heights = c(1.2, 1.2, 1.2, 1.5))

ggsave(
  filename = file.path(outdir, "plots", "sex_linked_snp_density_panels.pdf"),
  plot = combined, width = 14, height = 20
)
ggsave(
  filename = file.path(outdir, "plots", "sex_linked_snp_density_panels.png"),
  plot = combined, width = 14, height = 20, dpi = 300
)

ggsave(
  filename = file.path(outdir, "plots", "male_female_density_track.pdf"),
  plot = p1, width = 14, height = 8
)
ggsave(
  filename = file.path(outdir, "plots", "male_female_density_track.png"),
  plot = p1, width = 14, height = 8, dpi = 300
)

ggsave(
  filename = file.path(outdir, "plots", "male_female_delta_track.pdf"),
  plot = p2, width = 14, height = 8
)
ggsave(
  filename = file.path(outdir, "plots", "male_female_delta_track.png"),
  plot = p2, width = 14, height = 8, dpi = 300
)

ggsave(
  filename = file.path(outdir, "plots", "female_vs_male_scatter.pdf"),
  plot = p3, width = 12.5, height = 6
)
ggsave(
  filename = file.path(outdir, "plots", "female_vs_male_scatter.png"),
  plot = p3, width = 12.5, height = 6, dpi = 300
)

ggsave(
  filename = file.path(outdir, "plots", "haplotype_density_line_ci.pdf"),
  plot = p3_haps, width = 12.5, height = 8
)
ggsave(
  filename = file.path(outdir, "plots", "haplotype_density_line_ci.png"),
  plot = p3_haps, width = 12.5, height = 8, dpi = 300
)

# ============================================================
# Ranked candidate intervals
# ============================================================
candidate_cols <- c(
  "interval_key", "gene_id", "chrom", "start", "end", "len_bp", "mid_kb",
  "mean_male_nadia", "mean_female_nadia", "delta_male_female_nadia", "q_male_female_nadia",
  "mean_male_cb", "mean_female_cb", "delta_male_female_cb", "q_male_female_cb",
  "mean_WW", "mean_ZW", "mean_ZZ", "p_haps", "sex_linked_flag", "score"
)

candidates <- stats[, ..candidate_cols]
fwrite(candidates, file.path(outdir, "ranked_candidate_intervals.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("Done. Results written to: ", normalizePath(outdir))
