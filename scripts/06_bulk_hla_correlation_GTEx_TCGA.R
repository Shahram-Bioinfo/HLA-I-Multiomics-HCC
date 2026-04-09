# ==========================================================
# Script: 06_bulk_hla_correlation_GTEx_TCGA.R
# Purpose:
#   Perform differential correlation analysis for HLA-A,
#   HLA-B, and HLA-C across GTEx liver, TCGA-LIHC tumor,
#   and TCGA-LIHC normal samples, then export correlation
#   tables, differential Z-score comparisons, heatmap text
#   files, and common genes across HLA-A/B/C.
#
# Input:
#   data/differential_correlation/
#     - TcgaTargetGtex_rsem_gene_tpm
#     - SampleID-LIHC-GTEx.csv
#
# Output:
#   results/differential_hla_correlation/tables/
#     - Cor-A-GTEx.xlsx
#     - Cor-A-TCGA-T.xlsx
#     - Cor-A-TCGA-N.xlsx
#     - Cor-B-GTEx.xlsx
#     - Cor-B-TCGA-T.xlsx
#     - Cor-B-TCGA-N.xlsx
#     - Cor-C-GTEx.xlsx
#     - Cor-C-TCGA-T.xlsx
#     - Cor-C-TCGA-N.xlsx
#     - A_results.xlsx
#     - B_results.xlsx
#     - C_results.xlsx
#     - Merged.Data.A.xlsx
#     - Merged.Data.B.xlsx
#     - Merged.Data.C.xlsx
#     - Common_genes_HLA_ABC.xlsx
#     - merged_zscore_results_A.csv
#     - merged_zscore_results_B.csv
#     - merged_zscore_results_C.csv
#     - Heatmap-A.txt
#     - Heatmap-B.txt
#     - Heatmap-C.txt
#
#   results/differential_hla_correlation/logs/
#     - 11_differential_hla_correlation_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c("dplyr", "openxlsx")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running the script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."

data_dir <- file.path(project_root, "data", "differential_correlation")

expr_file <- file.path(data_dir, "TcgaTargetGtex_rsem_gene_tpm")
sample_file <- file.path(data_dir, "SampleID-LIHC-GTEx.csv")

results_root <- file.path(project_root, "results", "differential_hla_correlation")
table_dir <- file.path(results_root, "tables")
log_dir <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_file)) {
  stop("Expression file not found: ", expr_file)
}

if (!file.exists(sample_file)) {
  stop("Sample annotation file not found: ", sample_file)
}

# ---------------------------
# 2) Read expression matrix
# ---------------------------
expr_raw <- read.table(
  expr_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)

cat("Expression matrix dimensions:", dim(expr_raw), "\n")

# Remove Ensembl version suffix
rownames(expr_raw) <- gsub("\\.[0-9]+$", "", rownames(expr_raw))

# ---------------------------
# 3) Read sample annotation
# ---------------------------
annot <- read.csv(sample_file, stringsAsFactors = FALSE)

# Standardize sample IDs
annot$Sample.ID <- gsub("-", ".", annot$Sample.ID)
colnames(expr_raw) <- gsub("-", ".", colnames(expr_raw))

# ---------------------------
# 4) Extract sample groups
# ---------------------------
tumor_samples  <- annot$Sample.ID[annot$Group == "Primary Tumor"]
normal_samples <- annot$Sample.ID[annot$Group == "Solid Tissue Normal"]
gtex_samples   <- annot$Sample.ID[annot$Group == "GTEx(Liver)"]

tumor_common  <- intersect(tumor_samples, colnames(expr_raw))
normal_common <- intersect(normal_samples, colnames(expr_raw))
gtex_common   <- intersect(gtex_samples, colnames(expr_raw))

cat("GTEx liver samples found:", length(gtex_common), "\n")
cat("TCGA tumor samples found:", length(tumor_common), "\n")
cat("TCGA normal samples found:", length(normal_common), "\n")

# Transpose so that samples are rows and genes are columns
mat_gtex   <- t(expr_raw[, gtex_common, drop = FALSE])
mat_tumor  <- t(expr_raw[, tumor_common, drop = FALSE])
mat_normal <- t(expr_raw[, normal_common, drop = FALSE])

# Force numeric storage
storage.mode(mat_gtex) <- "numeric"
storage.mode(mat_tumor) <- "numeric"
storage.mode(mat_normal) <- "numeric"

# ---------------------------
# 5) Define HLA target genes
# ---------------------------
hla_targets <- list(
  A = "ENSG00000206503",
  B = "ENSG00000234745",
  C = "ENSG00000204525"
)

# ---------------------------
# 6) Helper functions
# ---------------------------
clip_r <- function(r, eps = 1e-6) {
  r[r >= 1] <- 1 - eps
  r[r <= -1] <- -1 + eps
  r
}

compute_cor_table <- function(mat, target_id) {
  if (!target_id %in% colnames(mat)) {
    stop("Target gene not found in matrix: ", target_id)
  }
  
  target_vec <- mat[, target_id]
  genes <- colnames(mat)
  
  r <- suppressWarnings(cor(mat, target_vec, use = "pairwise.complete.obs", method = "pearson"))
  r <- as.numeric(r)
  
  n_complete <- colSums(is.finite(mat) & is.finite(target_vec))
  t_stat <- r * sqrt((n_complete - 2) / pmax(1 - r^2, .Machine$double.eps))
  p_val <- 2 * pt(-abs(t_stat), df = pmax(n_complete - 2, 1))
  
  p_val[is.na(r)] <- NA
  p_val[n_complete < 3] <- NA
  
  p_adj <- p.adjust(p_val, method = "BH")
  
  r_for_z <- clip_r(r)
  fisher_z <- 0.5 * log((1 + r_for_z) / (1 - r_for_z))
  fisher_z[is.na(r)] <- NA
  
  data.frame(
    Gene.ID = genes,
    Correlation = r,
    P_Value = p_val,
    P_Value_Adjusted = p_adj,
    Fisher_Z = fisher_z,
    N = n_complete,
    stringsAsFactors = FALSE
  )
}

compare_groups <- function(gtex_res, tumor_res, target_id, label_prefix, out_dir, n_gtex, n_tumor) {
  merged <- merge(
    gtex_res,
    tumor_res,
    by = "Gene.ID",
    suffixes = c("_GTEx", "_Tumor")
  )
  
  merged <- merged[merged$Gene.ID != target_id, , drop = FALSE]
  
  merged$Correlation_GTEx_Sig <- ifelse(
    !is.na(merged$P_Value_Adjusted_GTEx) & merged$P_Value_Adjusted_GTEx < 0.05,
    merged$Correlation_GTEx,
    0
  )
  
  merged$Correlation_Tumor_Sig <- ifelse(
    !is.na(merged$P_Value_Adjusted_Tumor) & merged$P_Value_Adjusted_Tumor < 0.05,
    merged$Correlation_Tumor,
    0
  )
  
  merged$Correlation_GTEx_Sig <- clip_r(merged$Correlation_GTEx_Sig)
  merged$Correlation_Tumor_Sig <- clip_r(merged$Correlation_Tumor_Sig)
  
  merged$Z_GTEx <- 0.5 * log((1 + merged$Correlation_GTEx_Sig) / (1 - merged$Correlation_GTEx_Sig))
  merged$Z_Tumor <- 0.5 * log((1 + merged$Correlation_Tumor_Sig) / (1 - merged$Correlation_Tumor_Sig))
  
  merged$Z_diff <- merged$Z_GTEx - merged$Z_Tumor
  merged$SE_diff <- sqrt(1 / (n_gtex - 3) + 1 / (n_tumor - 3))
  merged$Z_stat <- merged$Z_diff / merged$SE_diff
  merged$P_diff <- 2 * pnorm(-abs(merged$Z_stat))
  merged$P_diff_adj <- p.adjust(merged$P_diff, method = "BH")
  
  merged$AbsCorMax <- pmax(
    abs(merged$Correlation_GTEx_Sig),
    abs(merged$Correlation_Tumor_Sig),
    na.rm = TRUE
  )
  
  merged_filtered <- merged %>%
    filter(!is.na(P_diff_adj)) %>%
    filter(P_diff_adj < 0.05) %>%
    filter(abs(Correlation_GTEx_Sig) > 0.7 | abs(Correlation_Tumor_Sig) > 0.7) %>%
    arrange(desc(abs(Z_diff)))
  
  heatmap_df <- merged_filtered %>%
    transmute(
      Gene.ID = Gene.ID,
      TCGA = round(Correlation_Tumor_Sig, 2),
      GTEx = round(Correlation_GTEx_Sig, 2)
    )
  
  merged_export <- merged %>%
    select(
      Gene.ID,
      Correlation_GTEx,
      P_Value_GTEx,
      P_Value_Adjusted_GTEx,
      Fisher_Z_GTEx,
      N_GTEx,
      Correlation_Tumor,
      P_Value_Tumor,
      P_Value_Adjusted_Tumor,
      Fisher_Z_Tumor,
      N_Tumor,
      Correlation_GTEx_Sig,
      Correlation_Tumor_Sig,
      Z_GTEx,
      Z_Tumor,
      Z_diff,
      SE_diff,
      Z_stat,
      P_diff,
      P_diff_adj,
      AbsCorMax
    )
  
  summary_export <- merged_filtered %>%
    transmute(
      `Gene Id` = Gene.ID,
      !!paste0("HLA-", label_prefix, " (GTEx-Liver)") := Correlation_GTEx_Sig,
      `Adjusted P GTEx` = P_Value_Adjusted_GTEx,
      !!paste0("HLA-", label_prefix, " (TCGA-LIHC-Normal)") := NA_real_,
      `Adjusted P Normal` = NA_real_,
      !!paste0("HLA-", label_prefix, " (TCGA-LIHC-Tumor)") := Correlation_Tumor_Sig,
      `Adjusted P Tumor` = P_Value_Adjusted_Tumor,
      `Z GTEx` = Z_GTEx,
      `Z Tumor` = Z_Tumor,
      `Z diff` = Z_diff,
      `P diff` = P_diff,
      `P diff adjusted` = P_diff_adj
    )
  
  write.csv(
    merged_export,
    file.path(out_dir, paste0("merged_zscore_results_", label_prefix, ".csv")),
    row.names = FALSE
  )
  
  write.table(
    heatmap_df,
    file.path(out_dir, paste0("Heatmap-", label_prefix, ".txt")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  write.xlsx(
    summary_export,
    file.path(out_dir, paste0(label_prefix, "_results.xlsx")),
    overwrite = TRUE
  )
  
  write.xlsx(
    merged_export,
    file.path(out_dir, paste0("Merged.Data.", label_prefix, ".xlsx")),
    overwrite = TRUE
  )
  
  list(
    all = merged_export,
    filtered = merged_filtered,
    heatmap = heatmap_df
  )
}

# ---------------------------
# 7) Compute correlation tables
# ---------------------------
results_GTEx_A <- compute_cor_table(mat_gtex, hla_targets$A)
results_TCGA_T_A <- compute_cor_table(mat_tumor, hla_targets$A)
results_TCGA_N_A <- compute_cor_table(mat_normal, hla_targets$A)

results_GTEx_B <- compute_cor_table(mat_gtex, hla_targets$B)
results_TCGA_T_B <- compute_cor_table(mat_tumor, hla_targets$B)
results_TCGA_N_B <- compute_cor_table(mat_normal, hla_targets$B)

results_GTEx_C <- compute_cor_table(mat_gtex, hla_targets$C)
results_TCGA_T_C <- compute_cor_table(mat_tumor, hla_targets$C)
results_TCGA_N_C <- compute_cor_table(mat_normal, hla_targets$C)

# ---------------------------
# 8) Save correlation tables
# ---------------------------
write.xlsx(results_GTEx_A, file.path(table_dir, "Cor-A-GTEx.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_T_A, file.path(table_dir, "Cor-A-TCGA-T.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_N_A, file.path(table_dir, "Cor-A-TCGA-N.xlsx"), overwrite = TRUE)

write.xlsx(results_GTEx_B, file.path(table_dir, "Cor-B-GTEx.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_T_B, file.path(table_dir, "Cor-B-TCGA-T.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_N_B, file.path(table_dir, "Cor-B-TCGA-N.xlsx"), overwrite = TRUE)

write.xlsx(results_GTEx_C, file.path(table_dir, "Cor-C-GTEx.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_T_C, file.path(table_dir, "Cor-C-TCGA-T.xlsx"), overwrite = TRUE)
write.xlsx(results_TCGA_N_C, file.path(table_dir, "Cor-C-TCGA-N.xlsx"), overwrite = TRUE)

# ---------------------------
# 9) Differential comparison between GTEx and TCGA tumor
# ---------------------------
n_gtex <- nrow(mat_gtex)
n_tumor <- nrow(mat_tumor)

cat("GTEx sample size used for differential comparison:", n_gtex, "\n")
cat("TCGA tumor sample size used for differential comparison:", n_tumor, "\n")

res_A <- compare_groups(
  gtex_res = results_GTEx_A,
  tumor_res = results_TCGA_T_A,
  target_id = hla_targets$A,
  label_prefix = "A",
  out_dir = table_dir,
  n_gtex = n_gtex,
  n_tumor = n_tumor
)

res_B <- compare_groups(
  gtex_res = results_GTEx_B,
  tumor_res = results_TCGA_T_B,
  target_id = hla_targets$B,
  label_prefix = "B",
  out_dir = table_dir,
  n_gtex = n_gtex,
  n_tumor = n_tumor
)

res_C <- compare_groups(
  gtex_res = results_GTEx_C,
  tumor_res = results_TCGA_T_C,
  target_id = hla_targets$C,
  label_prefix = "C",
  out_dir = table_dir,
  n_gtex = n_gtex,
  n_tumor = n_tumor
)

# ---------------------------
# 10) Common genes across HLA-A, HLA-B, and HLA-C
# ---------------------------
common_genes <- Reduce(
  intersect,
  list(
    res_A$filtered$Gene.ID,
    res_B$filtered$Gene.ID,
    res_C$filtered$Gene.ID
  )
)

common_df <- data.frame(Gene.ID = common_genes, stringsAsFactors = FALSE)

if (length(common_genes) > 0) {
  tmp_A <- res_A$filtered %>%
    filter(Gene.ID %in% common_genes) %>%
    select(Gene.ID, Correlation_GTEx_Sig, Correlation_Tumor_Sig, Z_diff, P_diff_adj)
  
  tmp_B <- res_B$filtered %>%
    filter(Gene.ID %in% common_genes) %>%
    select(Gene.ID, Correlation_GTEx_Sig, Correlation_Tumor_Sig, Z_diff, P_diff_adj)
  
  tmp_C <- res_C$filtered %>%
    filter(Gene.ID %in% common_genes) %>%
    select(Gene.ID, Correlation_GTEx_Sig, Correlation_Tumor_Sig, Z_diff, P_diff_adj)
  
  colnames(tmp_A)[-1] <- paste0("A_", colnames(tmp_A)[-1])
  colnames(tmp_B)[-1] <- paste0("B_", colnames(tmp_B)[-1])
  colnames(tmp_C)[-1] <- paste0("C_", colnames(tmp_C)[-1])
  
  common_df <- common_df %>%
    left_join(tmp_A, by = "Gene.ID") %>%
    left_join(tmp_B, by = "Gene.ID") %>%
    left_join(tmp_C, by = "Gene.ID")
}

write.xlsx(
  common_df,
  file.path(table_dir, "Common_genes_HLA_ABC.xlsx"),
  overwrite = TRUE
)

# ---------------------------
# 11) Log
# ---------------------------
sink(file.path(log_dir, "11_differential_hla_correlation_log.txt"))
cat("11_differential_hla_correlation.R completed\n")
cat("Expression file: ", expr_file, "\n", sep = "")
cat("Sample annotation file: ", sample_file, "\n", sep = "")
cat("Tables directory: ", table_dir, "\n", sep = "")
cat("GTEx samples: ", n_gtex, "\n", sep = "")
cat("TCGA tumor samples: ", n_tumor, "\n", sep = "")
cat("TCGA normal samples: ", nrow(mat_normal), "\n", sep = "")
cat("Filtered genes HLA-A: ", nrow(res_A$filtered), "\n", sep = "")
cat("Filtered genes HLA-B: ", nrow(res_B$filtered), "\n", sep = "")
cat("Filtered genes HLA-C: ", nrow(res_C$filtered), "\n", sep = "")
cat("Common genes across HLA-A/B/C: ", length(common_genes), "\n", sep = "")
sink()

# ---------------------------
# 12) Console summary
# ---------------------------
cat("\nSummary of filtered genes\n")
cat("HLA-A:", nrow(res_A$filtered), "\n")
cat("HLA-B:", nrow(res_B$filtered), "\n")
cat("HLA-C:", nrow(res_C$filtered), "\n")
cat("Common genes across A, B, and C:", length(common_genes), "\n")
cat("\nAll analyses completed successfully.\n")
cat("Results saved in:", results_root, "\n")