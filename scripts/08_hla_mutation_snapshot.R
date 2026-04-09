# ==========================================================
# Script: 08_hla_mutation_snapshot.R
# Purpose:
#   Download TCGA-LIHC somatic mutation data, build a MAF
#   object, summarize non-synonymous mutations in HLA class I
#   and full antigen processing machinery (APM) genes, and
#   generate oncoplots, barplots, lollipop plots, and
#   supplementary mutation tables.
#
# Input:
#   Downloaded from GDC using TCGAbiolinks:
#     - TCGA-LIHC Masked Somatic Mutation
#
# Output:
#   results/mutation/tables/
#     - TCGA_LIHC_MaskedSomaticMutation_raw.csv
#     - HLA_fullAPM_mutation_frequency_TCGA_LIHC.csv
#     - HLA_fullAPM_nonsyn_variants_TCGA_LIHC.csv
#
#   results/mutation/figures/
#     - Oncoplot_HLA_fullAPM_TCGA_LIHC.pdf
#     - Oncoplot_HLA_fullAPM_TCGA_LIHC_600dpi.tiff
#     - MutationFreq_Barplot_HLA_fullAPM_TCGA_LIHC.pdf
#     - MutationFreq_Barplot_HLA_fullAPM_TCGA_LIHC_600dpi.tiff
#     - Lollipop_HLA-A_TCGA_LIHC.pdf
#     - Lollipop_HLA-A_TCGA_LIHC_600dpi.tiff
#     - Lollipop_HLA-B_TCGA_LIHC.pdf
#     - Lollipop_HLA-B_TCGA_LIHC_600dpi.tiff
#     - Lollipop_HLA-C_TCGA_LIHC.pdf
#     - Lollipop_HLA-C_TCGA_LIHC_600dpi.tiff
#     - Combined_Oncoplot_Barplot_HLA_fullAPM_TCGA_LIHC.png
#
#   results/mutation/logs/
#     - maf_summary.txt
#     - mutation_sessionInfo.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "TCGAbiolinks", "maftools", "ggplot2",
  "data.table", "cowplot", "png", "grid"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running the script."
  )
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(maftools)
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(grid)
  library(png)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."

data_dir <- file.path(project_root, "data", "mutation")
results_root <- file.path(project_root, "results", "mutation")

figures_dir <- file.path(results_root, "figures")
tables_dir <- file.path(results_root, "tables")
logs_dir <- file.path(results_root, "logs")

output_dirs <- c(
  data_dir,
  results_root,
  figures_dir,
  tables_dir,
  logs_dir
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---------------------------
# 2) Helpers
# ---------------------------
save_pdf <- function(file, expr, width = 8, height = 6) {
  pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

save_tiff <- function(file, expr, width = 8, height = 6, dpi = 600) {
  tiff(
    filename = file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

save_png <- function(file, expr, width = 8, height = 6, dpi = 600) {
  png(
    filename = file,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# ---------------------------
# 3) Target genes
# ---------------------------
hla_genes <- c("HLA-A", "HLA-B", "HLA-C")

apm_genes <- c(
  "PSMB8", "PSMB9", "PSMB10",
  "TAP1", "TAP2", "TAPBP",
  "CALR", "CANX", "PDIA3",
  "ERAP1", "ERAP2",
  "NLRC5",
  "B2M"
)

genes_of_interest <- c(hla_genes, apm_genes)

# ---------------------------
# 4) Download TCGA-LIHC somatic mutations
# ---------------------------
query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query)
maf_tbl <- GDCprepare(query)

write.csv(
  maf_tbl,
  file.path(tables_dir, "TCGA_LIHC_MaskedSomaticMutation_raw.csv"),
  row.names = FALSE
)

# ---------------------------
# 5) Convert to MAF object
# ---------------------------
maf <- read.maf(maf = maf_tbl, isTCGA = TRUE)

sink(file.path(logs_dir, "maf_summary.txt"))
print(maf)
sink()

# ---------------------------
# 6) Filter non-synonymous mutations
# ---------------------------
nonsyn_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Translation_Start_Site"
)

m <- maf@data

m_sub <- subset(
  m,
  Hugo_Symbol %in% genes_of_interest &
    Variant_Classification %in% nonsyn_classes
)

# ---------------------------
# 7) Mutation frequency table
# ---------------------------
total_samples <- length(getSampleSummary(maf)$Tumor_Sample_Barcode)

freq_df <- aggregate(
  Tumor_Sample_Barcode ~ Hugo_Symbol,
  data = unique(m_sub[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]),
  FUN = function(x) length(unique(x))
)

colnames(freq_df) <- c("Gene", "Mutated_Samples")
freq_df$Total_Samples <- total_samples
freq_df$Freq <- round(100 * freq_df$Mutated_Samples / freq_df$Total_Samples, 2)
freq_df <- freq_df[order(match(freq_df$Gene, genes_of_interest)), ]

write.csv(
  freq_df,
  file.path(tables_dir, "HLA_fullAPM_mutation_frequency_TCGA_LIHC.csv"),
  row.names = FALSE
)

# ---------------------------
# 8) Oncoplot
# ---------------------------
save_pdf(
  file.path(figures_dir, "Oncoplot_HLA_fullAPM_TCGA_LIHC.pdf"),
  {
    oncoplot(
      maf = maf,
      genes = genes_of_interest,
      drawRowBar = TRUE,
      drawColBar = TRUE,
      showTumorSampleBarcodes = FALSE,
      sortByAnnotation = FALSE
    )
  },
  width = 10,
  height = 7
)

save_tiff(
  file.path(figures_dir, "Oncoplot_HLA_fullAPM_TCGA_LIHC_600dpi.tiff"),
  {
    oncoplot(
      maf = maf,
      genes = genes_of_interest,
      drawRowBar = TRUE,
      drawColBar = TRUE,
      showTumorSampleBarcodes = FALSE,
      sortByAnnotation = FALSE
    )
  },
  width = 10,
  height = 7,
  dpi = 600
)

# ---------------------------
# 9) Mutation frequency barplot
# ---------------------------
freq_df$Gene <- factor(freq_df$Gene, levels = genes_of_interest)

p_bar <- ggplot(freq_df, aes(x = Gene, y = Freq)) +
  geom_col(width = 0.7, fill = "#377EB8") +
  geom_text(
    aes(label = paste0(Freq, "% (", Mutated_Samples, "/", Total_Samples, ")")),
    vjust = -0.4,
    size = 3
  ) +
  labs(
    x = NULL,
    y = "Mutation frequency (%)",
    title = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_pdf(
  file.path(figures_dir, "MutationFreq_Barplot_HLA_fullAPM_TCGA_LIHC.pdf"),
  {
    print(p_bar)
  },
  width = 11,
  height = 6
)

save_tiff(
  file.path(figures_dir, "MutationFreq_Barplot_HLA_fullAPM_TCGA_LIHC_600dpi.tiff"),
  {
    print(p_bar)
  },
  width = 11,
  height = 6,
  dpi = 600
)

# ---------------------------
# 10) Supplementary table with detailed variants
# ---------------------------
supp_cols <- c(
  "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome",
  "Start_Position", "End_Position", "Reference_Allele",
  "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type"
)

keep <- intersect(supp_cols, colnames(m_sub))
supp_tbl <- m_sub[, ..keep]
data.table::setcolorder(supp_tbl, keep)

data.table::fwrite(
  supp_tbl,
  file.path(tables_dir, "HLA_fullAPM_nonsyn_variants_TCGA_LIHC.csv")
)

# ---------------------------
# 11) Lollipop plots for HLA-I genes
# ---------------------------
gene_summary <- getGeneSummary(maf)

for (g in hla_genes) {
  if (g %in% gene_summary$Hugo_Symbol) {
    save_pdf(
      file.path(figures_dir, paste0("Lollipop_", g, "_TCGA_LIHC.pdf")),
      {
        lollipopPlot(
          maf = maf,
          gene = g,
          AACol = "HGVSp_Short",
          showMutationRate = TRUE
        )
      },
      width = 10,
      height = 5
    )
    
    save_tiff(
      file.path(figures_dir, paste0("Lollipop_", g, "_TCGA_LIHC_600dpi.tiff")),
      {
        lollipopPlot(
          maf = maf,
          gene = g,
          AACol = "HGVSp_Short",
          showMutationRate = TRUE
        )
      },
      width = 10,
      height = 5,
      dpi = 600
    )
  }
}

# ---------------------------
# 12) Combined figure: Oncoplot (A) + Barplot (B)
# ---------------------------
tmp_oncoplot_png <- file.path(figures_dir, "tmp_oncoplot_panel.png")

save_png(
  tmp_oncoplot_png,
  {
    oncoplot(
      maf = maf,
      genes = genes_of_interest,
      drawRowBar = TRUE,
      drawColBar = TRUE,
      showTumorSampleBarcodes = FALSE,
      sortByAnnotation = FALSE
    )
  },
  width = 12,
  height = 7,
  dpi = 300
)

oncoplot_img <- png::readPNG(tmp_oncoplot_png)

p_onco <- ggplot() +
  annotation_raster(
    oncoplot_img,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf
  ) +
  theme_void()

combined_plot <- cowplot::plot_grid(
  p_onco,
  p_bar + theme(plot.margin = margin(5.5, 20, 5.5, 5.5)),
  labels = c("A", "B"),
  ncol = 1,
  rel_heights = c(1.7, 1)
)

ggsave(
  filename = file.path(figures_dir, "Combined_Oncoplot_Barplot_HLA_fullAPM_TCGA_LIHC.png"),
  plot = combined_plot,
  width = 11,
  height = 10,
  dpi = 1200,
  bg = "white"
)

unlink(tmp_oncoplot_png)

# ---------------------------
# 13) Session info log
# ---------------------------
sink(file.path(logs_dir, "mutation_sessionInfo.txt"))
sessionInfo()
sink()

cat("Done. All mutation outputs were written to:\n", results_root, "\n")