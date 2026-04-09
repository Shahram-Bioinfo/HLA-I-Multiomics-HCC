# ==========================================================
# Script: 02_allele_expression.R
# Purpose:
#   Analyze allele-specific expression of prevalent HLA-I
#   alleles (>5%) using Shapiro-Wilk, IQR-based outlier
#   removal, ANOVA, Tukey post hoc, and one-vs-rest t-tests.
#
# Input:
#   data/HLA_genotype_expression.csv
#   results/allele_specific_neoantigen/tables/
#     - prevalent_alleles_by_locus.csv
#
# Output:
#   results/allele_specific_neoantigen/tables/
#     - ASE_long_raw.csv
#     - expression_prevalent_only.csv
#     - expression_clean.csv
#     - expression_shapiro_raw.csv
#     - expression_shapiro_clean.csv
#     - expression_anova_by_locus.csv
#     - expression_tukey_posthoc.csv
#     - expression_vs_rest_tests.csv
#   results/allele_specific_neoantigen/figures/
#     - Figure_11B_D_Allele_Expression.png
#     - Figure_11B_D_Allele_Expression.pdf
#   results/allele_specific_neoantigen/logs/
#     - 02_allele_expression_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "dplyr", "stringr", "ggplot2", "ggpubr", "cowplot", "readr"
)

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(readr)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."     #set your directory

ase_file <- file.path(
  project_root, "data", "genotype_allele_specific_expression",
  "HLA_genotype_expression.csv"
)
prev_file <- file.path(project_root, "results", "allele_specific_neoantigen", "tables", "prevalent_alleles_by_locus.csv")

results_root <- file.path(project_root, "results", "allele_specific_neoantigen")
tab_dir      <- file.path(results_root, "tables")
fig_dir      <- file.path(results_root, "figures")
log_dir      <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ase_file)) stop("ASE/HLA file not found: ", ase_file)
if (!file.exists(prev_file)) stop("Prevalent allele file not found: ", prev_file)

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

remove_iqr_outliers <- function(x) {
  q1 <- stats::quantile(x, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(x, 0.75, na.rm = TRUE)
  iqr_val <- stats::IQR(x, na.rm = TRUE)
  lower <- q1 - 1.5 * iqr_val
  upper <- q3 + 1.5 * iqr_val
  !(x < lower | x > upper)
}

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (length(x) > 5000) x <- sample(x, 5000)
  tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
}

safe_anova <- function(df, y_col = "value") {
  if (dplyr::n_distinct(df$group) < 2) return(NA_real_)
  form <- stats::as.formula(paste(y_col, "~ group"))
  tryCatch({
    fit <- stats::aov(form, data = df)
    summary(fit)[[1]][["Pr(>F)"]][1]
  }, error = function(e) NA_real_)
}

safe_tukey <- function(df, y_col = "value") {
  if (dplyr::n_distinct(df$group) < 2) return(NULL)
  form <- stats::as.formula(paste(y_col, "~ group"))
  tryCatch({
    fit <- stats::aov(form, data = df)
    out <- stats::TukeyHSD(fit)$group
    out_df <- as.data.frame(out)
    out_df$comparison <- rownames(out_df)
    rownames(out_df) <- NULL
    out_df
  }, error = function(e) NULL)
}

safe_ttest_vs_all <- function(df, y_col = "value") {
  out <- lapply(unique(as.character(df$group)), function(g) {
    x <- df[df$group == g, y_col, drop = TRUE]
    y <- df[df$group != g, y_col, drop = TRUE]
    p <- tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
    data.frame(
      group = g,
      mean_group = mean(x, na.rm = TRUE),
      mean_all = mean(df[[y_col]], na.rm = TRUE),
      p_value = p,
      signif = ifelse(
        is.na(p), "ns",
        ifelse(p < 0.001, "***",
               ifelse(p < 0.01, "**",
                      ifelse(p < 0.05, "*", "ns")))
      ),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

make_compare_plot <- function(data_long, xlab_text, ylab_text, label_y = NULL) {
  data_long <- stats::na.omit(data_long)
  
  if (nrow(data_long) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }
  
  if (is.null(label_y)) {
    label_y <- max(data_long$value, na.rm = TRUE) * 1.10
  }
  
  ggpubr::ggboxplot(
    data_long,
    x = "group",
    y = "value",
    color = "group",
    add = "jitter",
    legend = "none"
  ) +
    ggpubr::rotate_x_text(angle = 45) +
    ggplot2::geom_hline(yintercept = mean(data_long$value, na.rm = TRUE), linetype = 5) +
    ggpubr::stat_compare_means(method = "anova", label.y = label_y) +
    ggpubr::stat_compare_means(
      label = "p.signif",
      method = "t.test",
      ref.group = ".all.",
      size = 5
    ) +
    ggplot2::labs(x = xlab_text, y = ylab_text) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 13, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 10),
      axis.text.y  = ggplot2::element_text(size = 10),
      text = ggplot2::element_text(size = 12)
    )
}

# ---------------------------
# 3) Read input files
# ---------------------------
ase_raw  <- read.csv(ase_file, check.names = FALSE, stringsAsFactors = FALSE)
prev_tbl <- read.csv(prev_file, stringsAsFactors = FALSE)

required_ase_cols <- c(
  "sample_id",
  "A1_allele", "A1_allele_exp",
  "A2-allele", "A2_allele_exp",
  "B1_allele", "B1_allele_exp",
  "B2-allele", "B2_allele_exp",
  "C1_allele", "C1_allele_exp",
  "C2-allele", "C2_allele_exp"
)

missing_ase_cols <- setdiff(required_ase_cols, colnames(ase_raw))
if (length(missing_ase_cols) > 0) {
  stop("Missing ASE columns: ", paste(missing_ase_cols, collapse = ", "))
}

prev_A <- prev_tbl %>% dplyr::filter(locus == "A") %>% dplyr::pull(allele)
prev_B <- prev_tbl %>% dplyr::filter(locus == "B") %>% dplyr::pull(allele)
prev_C <- prev_tbl %>% dplyr::filter(locus == "C") %>% dplyr::pull(allele)

# ---------------------------
# 4) Build expression long table
# ---------------------------
expr_long_raw <- dplyr::bind_rows(
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "A", group = normalize_display_allele(A1_allele, "A"), value = as.numeric(A1_allele_exp)),
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "A", group = normalize_display_allele(`A2-allele`, "A"), value = as.numeric(A2_allele_exp)),
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "B", group = normalize_display_allele(B1_allele, "B"), value = as.numeric(B1_allele_exp)),
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "B", group = normalize_display_allele(`B2-allele`, "B"), value = as.numeric(B2_allele_exp)),
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "C", group = normalize_display_allele(C1_allele, "C"), value = as.numeric(C1_allele_exp)),
  ase_raw %>%
    dplyr::transmute(sample_id = as.character(sample_id), patientBarcode = sample_to_case(sample_id), locus = "C", group = normalize_display_allele(`C2-allele`, "C"), value = as.numeric(C2_allele_exp))
) %>%
  dplyr::filter(!is.na(patientBarcode), patientBarcode != "") %>%
  dplyr::filter(!is.na(group), group != "") %>%
  dplyr::filter(!is.na(value))

write.csv(expr_long_raw, file.path(tab_dir, "ASE_long_raw.csv"), row.names = FALSE)

# ---------------------------
# 5) Keep only prevalent alleles
# ---------------------------
expr_prev <- expr_long_raw %>%
  dplyr::filter(
    (locus == "A" & group %in% prev_A) |
      (locus == "B" & group %in% prev_B) |
      (locus == "C" & group %in% prev_C)
  )

write.csv(expr_prev, file.path(tab_dir, "expression_prevalent_only.csv"), row.names = FALSE)

# ---------------------------
# 6) Normality and outlier handling
# ---------------------------
expr_shapiro_raw <- expr_prev %>%
  dplyr::group_by(locus, group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    shapiro_p = safe_shapiro(value),
    .groups = "drop"
  )

expr_clean <- expr_prev %>%
  dplyr::group_by(locus, group) %>%
  dplyr::mutate(keep = remove_iqr_outliers(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(keep)

write.csv(expr_clean, file.path(tab_dir, "expression_clean.csv"), row.names = FALSE)

expr_shapiro_clean <- expr_clean %>%
  dplyr::group_by(locus, group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    shapiro_p = safe_shapiro(value),
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------
# 7) Statistical analysis
# ---------------------------
expr_anova_tbl <- expr_clean %>%
  dplyr::group_by(locus) %>%
  dplyr::group_modify(~data.frame(anova_p = safe_anova(.x, "value"))) %>%
  dplyr::ungroup()

expr_tukey_tbl <- dplyr::bind_rows(lapply(unique(expr_clean$locus), function(loc) {
  d <- expr_clean %>% dplyr::filter(locus == loc)
  x <- safe_tukey(d, "value")
  if (is.null(x)) return(NULL)
  x$locus <- loc
  x
}))

expr_vs_all_tbl <- dplyr::bind_rows(lapply(unique(expr_clean$locus), function(loc) {
  d <- expr_clean %>% dplyr::filter(locus == loc)
  x <- safe_ttest_vs_all(d, "value")
  x$locus <- loc
  x
}))

write.csv(expr_shapiro_raw,   file.path(tab_dir, "expression_shapiro_raw.csv"), row.names = FALSE)
write.csv(expr_shapiro_clean, file.path(tab_dir, "expression_shapiro_clean.csv"), row.names = FALSE)
write.csv(expr_anova_tbl,     file.path(tab_dir, "expression_anova_by_locus.csv"), row.names = FALSE)
write.csv(expr_tukey_tbl,     file.path(tab_dir, "expression_tukey_posthoc.csv"), row.names = FALSE)
write.csv(expr_vs_all_tbl,    file.path(tab_dir, "expression_vs_rest_tests.csv"), row.names = FALSE)

# ---------------------------
# 8) Plot Figure 11B-D
# ---------------------------
pB <- make_compare_plot(
  expr_clean %>% dplyr::filter(locus == "A"),
  xlab_text = "HLA-A Alleles",
  ylab_text = "Allele Expression",
  label_y = max(expr_clean %>% dplyr::filter(locus == "A") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
)

pC <- make_compare_plot(
  expr_clean %>% dplyr::filter(locus == "B"),
  xlab_text = "HLA-B Alleles",
  ylab_text = "Allele Expression",
  label_y = max(expr_clean %>% dplyr::filter(locus == "B") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
)

pD <- make_compare_plot(
  expr_clean %>% dplyr::filter(locus == "C"),
  xlab_text = "HLA-C Alleles",
  ylab_text = "Allele Expression",
  label_y = max(expr_clean %>% dplyr::filter(locus == "C") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
)

figure11_BD <- cowplot::plot_grid(
  pB, pC, pD,
  labels = c("B", "C", "D"),
  ncol = 3,
  label_size = 16
)

ggsave(
  filename = file.path(fig_dir, "Figure_11B_D_Allele_Expression.png"),
  plot = figure11_BD,
  width = 14,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure_11B_D_Allele_Expression.pdf"),
  plot = figure11_BD,
  width = 14,
  height = 6,
  bg = "white"
)

# ---------------------------
# 9) Log
# ---------------------------
sink(file.path(log_dir, "02_allele_expression_log.txt"))
cat("02_allele_expression.R completed\n")
cat("Input ASE file: ", ase_file, "\n", sep = "")
cat("Input prevalent allele file: ", prev_file, "\n", sep = "")
cat("Rows in expr_long_raw: ", nrow(expr_long_raw), "\n", sep = "")
cat("Rows in expr_prev: ", nrow(expr_prev), "\n", sep = "")
cat("Rows in expr_clean: ", nrow(expr_clean), "\n", sep = "")
sink()

cat("Done: 02_allele_expression.R\n")