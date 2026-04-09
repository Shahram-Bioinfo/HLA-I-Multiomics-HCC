# ==========================================================
# Script: 01_allele_frequency.R
# Purpose:
#   Calculate HLA-I allele frequencies from the
#   genotype-expression file and identify prevalent alleles
#   with frequency > 5% for downstream analyses.
#
# Input:
#   data/HLA_genotype_expression.csv
#
# Output:
#   results/allele_specific_neoantigen/tables/
#     - allele_frequency_table.csv
#     - prevalent_alleles_by_locus.csv
#   results/allele_specific_neoantigen/figures/
#     - Figure_11A_Allele_Frequency.png
#     - Figure_11A_Allele_Frequency.pdf
#   results/allele_specific_neoantigen/logs/
#     - 01_allele_frequency_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c("dplyr", "stringr", "ggplot2", "cowplot", "readr")

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(readr)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."   #set your directory

ase_file <- file.path(
  project_root, "data", "genotype_allele_specific_expression",
  "HLA_genotype_expression.csv"
)

results_root <- file.path(project_root, "results", "allele_specific_neoantigen")
tab_dir      <- file.path(results_root, "tables")
fig_dir      <- file.path(results_root, "figures")
log_dir      <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ase_file)) stop("ASE/HLA file not found: ", ase_file)

# ---------------------------
# 2) Helpers
# ---------------------------
sample_to_case <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) >= 12, substr(x, 1, 12), x)
}

normalize_display_allele <- function(x, locus = NULL) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  
  # HLA-A02:01 -> A*02:01
  x <- stringr::str_replace(x, "^HLA-([ABC])", "\\1")
  x <- ifelse(
    !is.na(x) & stringr::str_detect(x, "^[ABC][0-9]"),
    stringr::str_replace(x, "^([ABC])", "\\1*"),
    x
  )
  
  if (!is.null(locus)) {
    x <- ifelse(!is.na(x) & !stringr::str_detect(x, paste0("^", locus, "\\*")), NA, x)
  }
  
  x
}

make_frequency_plot <- function(df, locus_label) {
  d <- df %>% dplyr::filter(locus == locus_label, is_prevalent)
  
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No prevalent HLA-", locus_label, " alleles", sep = "")))
  }
  
  ggplot(d, aes(x = allele, y = allele_frequency * 100, fill = allele)) +
    geom_col() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = paste0("HLA-", locus_label, " alleles"),
      y = "Frequency (%)"
    )
}

# ---------------------------
# 3) Read genotype-expression file
# ---------------------------
ase_raw <- read.csv(ase_file, check.names = FALSE, stringsAsFactors = FALSE)

required_ase_cols <- c(
  "sample_id",
  "A1_allele", "A2-allele",
  "B1_allele", "B2-allele",
  "C1_allele", "C2-allele"
)

missing_ase_cols <- setdiff(required_ase_cols, colnames(ase_raw))
if (length(missing_ase_cols) > 0) {
  stop("Missing ASE columns: ", paste(missing_ase_cols, collapse = ", "))
}

# ---------------------------
# 4) Build allele long table
# ---------------------------
allele_long <- dplyr::bind_rows(
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "A",
      allele = normalize_display_allele(A1_allele, "A")
    ),
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "A",
      allele = normalize_display_allele(`A2-allele`, "A")
    ),
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "B",
      allele = normalize_display_allele(B1_allele, "B")
    ),
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "B",
      allele = normalize_display_allele(`B2-allele`, "B")
    ),
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "C",
      allele = normalize_display_allele(C1_allele, "C")
    ),
  ase_raw %>%
    dplyr::transmute(
      sample_id = as.character(sample_id),
      patientBarcode = sample_to_case(sample_id),
      locus = "C",
      allele = normalize_display_allele(`C2-allele`, "C")
    )
) %>%
  dplyr::filter(!is.na(patientBarcode), patientBarcode != "") %>%
  dplyr::filter(!is.na(allele), allele != "")

# ---------------------------
# 5) Calculate allele frequencies
# ---------------------------
allele_freq <- allele_long %>%
  dplyr::group_by(locus, allele) %>%
  dplyr::summarise(
    allele_calls = dplyr::n(),
    n_patients = dplyr::n_distinct(patientBarcode),
    .groups = "drop"
  ) %>%
  dplyr::group_by(locus) %>%
  dplyr::mutate(
    total_allele_calls = sum(allele_calls),
    allele_frequency = allele_calls / total_allele_calls,
    is_prevalent = allele_frequency > 0.05
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(locus, dplyr::desc(allele_frequency))

prevalent_tbl <- allele_freq %>%
  dplyr::filter(is_prevalent) %>%
  dplyr::select(locus, allele, allele_calls, n_patients, total_allele_calls, allele_frequency)

write.csv(allele_freq, file.path(tab_dir, "allele_frequency_table.csv"), row.names = FALSE)
write.csv(prevalent_tbl, file.path(tab_dir, "prevalent_alleles_by_locus.csv"), row.names = FALSE)

# ---------------------------
# 6) Plot Figure 11A
# ---------------------------
pA <- make_frequency_plot(allele_freq, "A")
pB <- make_frequency_plot(allele_freq, "B")
pC <- make_frequency_plot(allele_freq, "C")

figure11A <- cowplot::plot_grid(
  pA, pB, pC,
  labels = c("A1", "A2", "A3"),
  ncol = 3,
  label_size = 14
)

ggsave(
  filename = file.path(fig_dir, "Figure_11A_Allele_Frequency.png"),
  plot = figure11A,
  width = 14,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure_11A_Allele_Frequency.pdf"),
  plot = figure11A,
  width = 14,
  height = 5,
  bg = "white"
)

# ---------------------------
# 7) Log
# ---------------------------
sink(file.path(log_dir, "01_allele_frequency_log.txt"))
cat("01_allele_frequency.R completed\n")
cat("Input ASE file: ", ase_file, "\n", sep = "")
cat("Total allele rows: ", nrow(allele_long), "\n", sep = "")
cat("Total prevalent alleles: ", nrow(prevalent_tbl), "\n", sep = "")
sink()

cat("Done: 01_allele_frequency.R\n")