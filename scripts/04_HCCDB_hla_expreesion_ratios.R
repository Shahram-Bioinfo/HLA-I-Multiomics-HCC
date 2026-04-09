# ==========================================================
# Script: 04_HCCDB_hla_expreesion_ratios.R
# Purpose:
#   Generate boxplots for HLA-A, HLA-B, HLA-C, and their
#   expression ratios across HCCDB datasets, then export
#   single-panel and multi-panel high-resolution figures.
#
# Input:
#   data/bulk_RNA_seq/boxplot/Other-Data-set/3 Group/
#     - Data_Subset_HCCDB*.xlsx
#
#   data/bulk_RNA_seq/
#     - Dataset.ID.csv
#
# Output:
#   results/hccdb_boxplots_hla_ratios/figures/
#     - HLA_A.pdf
#     - HLA_A_600dpi.tiff
#     - HLA_A_600dpi.png
#     - HLA_B.pdf
#     - HLA_B_600dpi.tiff
#     - HLA_B_600dpi.png
#     - HLA_C.pdf
#     - HLA_C_600dpi.tiff
#     - HLA_C_600dpi.png
#     - AB.pdf
#     - AB_600dpi.tiff
#     - AB_600dpi.png
#     - AC.pdf
#     - AC_600dpi.tiff
#     - AC_600dpi.png
#     - BC.pdf
#     - BC_600dpi.tiff
#     - BC_600dpi.png
#     - Plots.pdf
#     - Plots_vector.pdf
#     - Plots_600dpi.tiff
#
#   results/hccdb_boxplots_hla_ratios/logs/
#     - 09_hccdb_boxplots_hla_ratios_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "readxl", "ggplot2", "dplyr", "ggsignif",
  "ggpubr", "gridExtra"
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
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(ggsignif)
  library(ggpubr)
  library(gridExtra)
})

has_cairo <- requireNamespace("Cairo", quietly = TRUE)

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "G:/HLA-I-Multiomics-HCC"

input_dir <- file.path(
  project_root, "data", "bulk_RNA_seq", "boxplot", "Other-Data-set", "3 Group"
)

mapping_file <- file.path(
  project_root, "data", "bulk_RNA_seq", "Dataset.ID.csv"
)

results_root <- file.path(project_root, "results", "hccdb_boxplots_hla_ratios")
fig_dir <- file.path(results_root, "figures")
log_dir <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop("Input folder does not exist: ", input_dir)
}

if (!file.exists(mapping_file)) {
  stop("The mapping file does not exist: ", mapping_file)
}

file_paths <- list.files(
  path = input_dir,
  pattern = "Data_Subset_HCCDB.*\\.xlsx",
  full.names = TRUE
)

# remove temp + hidden + corrupted
file_paths <- file_paths[
  !grepl("^~\\$", basename(file_paths)) &
    file.exists(file_paths)
]

if (length(file_paths) == 0) {
  stop("No matching Excel files found in: ", input_dir)
}

# ---------------------------
# 2) Load input data
# ---------------------------
data_list <- lapply(file_paths, read_excel)
mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)

name_replacements <- c(
  "HCCDB1" = "GSE22058",
  "HCCDB4" = "GSE36376",
  "HCCDB8" = "GSE9843",
  "HCCDB11" = "GSE46444",
  "HCCDB12" = "GSE54236",
  "HCCDB13" = "GSE63898",
  "HCCDB14" = "GSE43619",
  "HCCDB15" = "TCGA-LIHC",
  "HCCDB17" = "GSE76427",
  "HCCDB18" = "ICGC-LIRI-JP",
  "HCCDB19" = "GSE124751",
  "HCCDB20" = "GSE164760",
  "HCCDB21" = "GSE134568",
  "HCCDB22" = "GSE87630",
  "HCCDB24" = "GSE76427",
  "HCCDB25" = "OEP000321",
  "HCCDB26" = "GSE112790",
  "HCCDB27" = "GSE121248",
  "HCCDB28" = "GSE109211",
  "HCCDB29" = "GSE89377",
  "HCCDB30" = "GSE148355"
)

dataset_names <- sapply(file_paths, function(x) {
  dataset_id <- sub(".*HCCDB\\.(\\d+).*", "HCCDB\\1", x)
  name_replacements[dataset_id]
})

names(data_list) <- dataset_names

# ---------------------------
# 3) Helpers
# ---------------------------
plot_pal <- c("Adjacent" = "#18cc50", "HCC" = "#f91e0d")

save_plot_set <- function(plot_obj, base_name, width = 7, height = 5) {
  ggsave(
    filename = file.path(fig_dir, paste0(base_name, ".pdf")),
    plot = plot_obj,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(base_name, "_600dpi.tiff")),
    plot = plot_obj,
    device = "tiff",
    dpi = 600,
    compression = "lzw",
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(base_name, "_600dpi.png")),
    plot = plot_obj,
    dpi = 600,
    width = width,
    height = height,
    units = "in"
  )
}

make_combined_plot_data <- function(data_list_obj, value_col) {
  extractor <- function(df) {
    df %>% dplyr::select(Group = 1, Type = 3, Value = all_of(value_col))
  }
  
  data_list_tmp <- lapply(data_list_obj, extractor)
  bind_rows(data_list_tmp, .id = "Dataset")
}

make_boxplot <- function(df, panel_title, y_label) {
  ggplot(df, aes(x = Dataset, y = Value, fill = Type)) +
    geom_boxplot() +
    stat_compare_means(
      aes(group = Type),
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 5,
      color = "#4518a6",
      vjust = 1
    ) +
    theme_bw() +
    labs(title = panel_title, x = "", y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = plot_pal)
}

# ---------------------------
# 4) HLA-A
# ---------------------------
combined_data_A <- make_combined_plot_data(data_list, 4)

p_A <- make_boxplot(
  combined_data_A,
  panel_title = "A",
  y_label = "HLA-A Expression"
)

p1 <- p_A
print(p1)

save_plot_set(p1, "HLA_A")

# ---------------------------
# 5) HLA-B
# ---------------------------
combined_data_B <- make_combined_plot_data(data_list, 5)

p_B <- make_boxplot(
  combined_data_B,
  panel_title = "B",
  y_label = "HLA-B Expression"
)

p2 <- p_B
print(p2)

save_plot_set(p2, "HLA_B")

# ---------------------------
# 6) HLA-C
# ---------------------------
combined_data_C <- make_combined_plot_data(data_list, 6)

p_C <- make_boxplot(
  combined_data_C,
  panel_title = "C",
  y_label = "HLA-C Expression"
)

p3 <- p_C
print(p3)

save_plot_set(p3, "HLA_C")

# ---------------------------
# 7) AB ratio
# ---------------------------
combined_data_AB <- make_combined_plot_data(data_list, 7)

p_AB <- make_boxplot(
  combined_data_AB,
  panel_title = "D",
  y_label = "HLA-A to HLA-B Expression ratio"
)

p4 <- p_AB
print(p4)

save_plot_set(p4, "AB")

# ---------------------------
# 8) AC ratio
# ---------------------------
combined_data_AC <- make_combined_plot_data(data_list, 8)

p_AC <- make_boxplot(
  combined_data_AC,
  panel_title = "E",
  y_label = "HLA-A to HLA-C Expression ratio"
)

p5 <- p_AC
print(p5)

save_plot_set(p5, "AC")

# ---------------------------
# 9) BC ratio
# ---------------------------
combined_data_BC <- make_combined_plot_data(data_list, 9)

p_BC <- make_boxplot(
  combined_data_BC,
  panel_title = "F",
  y_label = "HLA-B to HLA-C Expression ratio"
)

p6 <- p_BC
print(p6)

save_plot_set(p6, "BC")

# ---------------------------
# 10) Multi-panel grid
# ---------------------------
pdf(file.path(fig_dir, "Plots.pdf"), width = 14, height = 10)
grid.arrange(p1, p4, p2, p5, p3, p6, ncol = 2, nrow = 3)
dev.off()

if (has_cairo) {
  Cairo::CairoPDF(file.path(fig_dir, "Plots_vector.pdf"), width = 14, height = 10)
  grid.arrange(p1, p4, p2, p5, p3, p6, ncol = 2, nrow = 3)
  dev.off()
}

tiff(
  filename = file.path(fig_dir, "Plots_600dpi.tiff"),
  width = 14,
  height = 10,
  units = "in",
  res = 600,
  compression = "lzw"
)
grid.arrange(p1, p4, p2, p5, p3, p6, ncol = 2, nrow = 3)
dev.off()

# ---------------------------
# 11) Log
# ---------------------------
sink(file.path(log_dir, "09_hccdb_boxplots_hla_ratios_log.txt"))
cat("09_hccdb_boxplots_hla_ratios.R completed\n")
cat("Input directory: ", input_dir, "\n", sep = "")
cat("Mapping file: ", mapping_file, "\n", sep = "")
cat("Number of input files: ", length(file_paths), "\n", sep = "")
cat("Figures directory: ", fig_dir, "\n", sep = "")
cat("Cairo available: ", has_cairo, "\n", sep = "")
sink()

cat("Done: 09_hccdb_boxplots_hla_ratios.R\n")

cat("Done. Updated Supplementary Tables 2, 3, and 4 were written to:\n", results_root, "\n")