# ==========================================================
# Script: 03_neoantigen_affinity.R
# Purpose:
#   Prepare neoantigen peptide files, generate netMHCpan
#   Linux shell scripts, parse raw netMHCpan outputs, build
#   compare CSVs, and analyze peptide-HLA binding affinities
#   for prevalent HLA-I alleles (>5%).
#
# Input:
#   data/HLA_genotype_expression.csv
#   data/neoantigen/TCIA-NeoantigensData.tsv
#   results/allele_specific_neoantigen/tables/
#     - prevalent_alleles_by_locus.csv
#
# Raw netMHCpan outputs:
#   Priority 1:
#     results/neoantigen/A
#     results/neoantigen/B
#     results/neoantigen/C
#
#   Fallback packaged outputs:
#     data/neoantigen/A
#     data/neoantigen/B
#     data/neoantigen/C
#
# Output:
#   results/neoantigen/
#     - peptides/*.pep.text
#     - linux_scripts/Code-Linux-HLA-A.sh
#     - linux_scripts/Code-Linux-HLA-B.sh
#     - linux_scripts/Code-Linux-HLA-C.sh
#
#   results/allele_specific_neoantigen/tables/
#     - neoantigen_clean.csv
#     - HLA_Genotype_for_netMHCpan.csv
#     - peptide_file_qc.csv
#     - manifest.csv
#     - parsed_netMHCpan_affinity.csv
#     - HLA-A-compare.csv
#     - HLA-B-compare.csv
#     - HLA-C-compare.csv
#     - HLA-Total-compare.csv
#     - affinity_long_from_compare_csvs.csv
#     - affinity_anova_by_locus.csv
#     - affinity_dunn_posthoc.csv
#     - affinity_vs_rest_tests.csv
#     - affinity_summary_by_allele.csv
#
#   results/allele_specific_neoantigen/figures/
#     - Figure_11E_G_Affinity.png
#     - Figure_11E_G_Affinity.pdf
#
#   results/allele_specific_neoantigen/logs/
#     - 03_neoantigen_affinity_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "ggplot2",
  "ggpubr", "reshape2", "dunn.test", "cowplot"
)

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(dunn.test)
  library(cowplot)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."  

ase_file <- file.path(
  project_root, "data", "genotype_allele_specific_expression",
  "HLA_genotype_expression.csv"
)
neo_file  <- file.path(project_root, "data", "neoantigen", "TCIA-NeoantigensData.tsv")
prev_file <- file.path(project_root, "results", "allele_specific_neoantigen", "tables", "prevalent_alleles_by_locus.csv")

results_root <- file.path(project_root, "results", "allele_specific_neoantigen")
tab_dir      <- file.path(results_root, "tables")
fig_dir      <- file.path(results_root, "figures")
log_dir      <- file.path(results_root, "logs")

# netMHCpan workflow folders
nmhc_root      <- file.path(project_root, "results", "neoantigen")
peptide_dir    <- file.path(nmhc_root, "peptides")
script_dir     <- file.path(nmhc_root, "linux_scripts")
raw_A_dir      <- file.path(nmhc_root, "A")
raw_B_dir      <- file.path(nmhc_root, "B")
raw_C_dir      <- file.path(nmhc_root, "C")

# fallback packaged raw outputs
neo_data_root  <- file.path(project_root, "data", "neoantigen")
fallback_A_dir <- file.path(neo_data_root, "A")
fallback_B_dir <- file.path(neo_data_root, "B")
fallback_C_dir <- file.path(neo_data_root, "C")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(nmhc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(peptide_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_A_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_B_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_C_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ase_file)) stop("ASE/HLA file not found: ", ase_file)
if (!file.exists(neo_file)) stop("Neoantigen file not found: ", neo_file)
if (!file.exists(prev_file)) stop("Prevalent allele file not found: ", prev_file)

# ---------------------------
# 2) Helpers
# ---------------------------
sample_to_case <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) >= 12, substr(x, 1, 12), x)
}

sanitize_id <- function(x) {
  stringr::str_replace_all(as.character(x), "-", "_")
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

normalize_netmhc_allele <- function(x, locus = NULL) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  x <- stringr::str_replace(x, "^HLA-", "")
  x <- stringr::str_replace_all(x, "\\*", "")
  x <- stringr::str_replace_all(x, "\\s+", "")
  
  if (!is.null(locus)) {
    x <- ifelse(!is.na(x) & !stringr::str_detect(x, paste0("^", locus)), paste0(locus, x), x)
  }
  
  ifelse(is.na(x), NA, paste0("HLA-", x))
}

write_peptide_file <- function(peptides, out_file) {
  write.table(
    peptides,
    file = out_file,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

build_netmhcpan_cmd <- function(alleles, peptide_file_rel, output_file_rel) {
  paste(
    "../netMHCpan",
    "-a", paste(alleles, collapse = ","),
    "-p", peptide_file_rel,
    "-BA -w",
    paste0(">", output_file_rel)
  )
}

parse_netmhcpan_file <- function(file_path, locus_label) {
  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[stringr::str_detect(lines, "^\\s*[0-9]+\\s+")]
  
  if (length(lines) == 0) return(NULL)
  
  out <- vector("list", length(lines))
  
  for (i in seq_along(lines)) {
    tokens <- stringr::str_split(stringr::str_squish(lines[i]), "\\s+")[[1]]
    
    if (length(tokens) < 8) {
      out[[i]] <- NULL
      next
    }
    
    n <- length(tokens)
    bind_level <- if (tokens[n] %in% c("SB", "WB", "NB")) tokens[n] else NA_character_
    idx_last_numeric <- if (!is.na(bind_level)) n - 1 else n
    
    rank_ba  <- suppressWarnings(as.numeric(tokens[idx_last_numeric]))
    affinity <- suppressWarnings(as.numeric(tokens[idx_last_numeric - 1]))
    score_ba <- suppressWarnings(as.numeric(tokens[idx_last_numeric - 2]))
    
    mhc_raw <- tokens[2]
    peptide <- tokens[3]
    
    out[[i]] <- data.frame(
      source_file = basename(file_path),
      patientBarcode = sample_to_case(stringr::str_replace_all(basename(file_path), "_", "-")),
      locus = locus_label,
      allele_raw = mhc_raw,
      allele_display = normalize_display_allele(mhc_raw, locus_label),
      allele_netmhc = normalize_netmhc_allele(mhc_raw, locus_label),
      peptide = peptide,
      score_ba = score_ba,
      affinity_nM = affinity,
      rank_ba = rank_ba,
      bind_level = bind_level,
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(out)
}

pivot_compare_csv <- function(parsed_df, prevalent_display, locus_label) {
  d <- parsed_df %>%
    dplyr::filter(locus == locus_label, allele_display %in% prevalent_display) %>%
    dplyr::group_by(patientBarcode, allele_display) %>%
    dplyr::summarise(value = mean(affinity_nM, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(group = allele_display)
  
  if (nrow(d) == 0) return(NULL)
  
  wide <- tidyr::pivot_wider(d, names_from = group, values_from = value)
  as.data.frame(wide)
}

safe_anova <- function(df, y_col = "value") {
  if (dplyr::n_distinct(df$group) < 2) return(NA_real_)
  form <- stats::as.formula(paste(y_col, "~ group"))
  tryCatch({
    fit <- stats::aov(form, data = df)
    summary(fit)[[1]][["Pr(>F)"]][1]
  }, error = function(e) NA_real_)
}

safe_dunn <- function(df, y_col = "value") {
  if (dplyr::n_distinct(df$group) < 2) return(NULL)
  tryCatch({
    d <- dunn.test::dunn.test(df[[y_col]], df$group, method = "bh", kw = FALSE, list = TRUE)
    data.frame(
      comparison = d$comparisons,
      z = d$Z,
      p_raw = d$P,
      p_adj = d$P.adjusted,
      stringsAsFactors = FALSE
    )
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

prepare_affinity_long_from_compare <- function(compare_df, locus_label) {
  if (is.null(compare_df)) return(NULL)
  if (ncol(compare_df) == 0) return(NULL)
  
  if ("patientBarcode" %in% colnames(compare_df)) {
    compare_df <- compare_df[, setdiff(colnames(compare_df), "patientBarcode"), drop = FALSE]
  }
  
  long <- reshape2::melt(compare_df, variable.name = "group", value.name = "value")
  long <- stats::na.omit(long)
  long$locus <- locus_label
  long
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
# 3) Read inputs
# ---------------------------
ase_raw  <- read.csv(ase_file, check.names = FALSE, stringsAsFactors = FALSE)
neo_raw  <- read.delim(neo_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
prev_tbl <- read.csv(prev_file, stringsAsFactors = FALSE)

prev_A <- prev_tbl %>% dplyr::filter(locus == "A") %>% dplyr::pull(allele)
prev_B <- prev_tbl %>% dplyr::filter(locus == "B") %>% dplyr::pull(allele)
prev_C <- prev_tbl %>% dplyr::filter(locus == "C") %>% dplyr::pull(allele)

required_neo_cols <- c("patientBarcode", "peptide")
missing_neo_cols <- setdiff(required_neo_cols, colnames(neo_raw))
if (length(missing_neo_cols) > 0) {
  stop("Missing neoantigen columns: ", paste(missing_neo_cols, collapse = ", "))
}

# ---------------------------
# 4) Clean neoantigen table
# ---------------------------
neo_tbl <- neo_raw %>%
  dplyr::transmute(
    patientBarcode = as.character(patientBarcode),
    disease = if ("disease" %in% colnames(neo_raw)) as.character(disease) else NA_character_,
    gene    = if ("gene" %in% colnames(neo_raw)) as.character(gene) else NA_character_,
    peptide = as.character(peptide)
  ) %>%
  dplyr::mutate(
    patientBarcode = stringr::str_trim(patientBarcode),
    peptide = stringr::str_trim(peptide)
  ) %>%
  dplyr::filter(!is.na(patientBarcode), patientBarcode != "") %>%
  dplyr::filter(!is.na(peptide), peptide != "") %>%
  dplyr::distinct()

write.csv(neo_tbl, file.path(tab_dir, "neoantigen_clean.csv"), row.names = FALSE)

# ---------------------------
# 5) Build genotype table for netMHCpan
# ---------------------------
genotype_tbl <- ase_raw %>%
  dplyr::transmute(
    sample_id = as.character(sample_id),
    patientBarcode = sample_to_case(sample_id),
    HLA_A_1 = normalize_netmhc_allele(A1_allele, "A"),
    HLA_A_2 = normalize_netmhc_allele(`A2-allele`, "A"),
    HLA_B_1 = normalize_netmhc_allele(B1_allele, "B"),
    HLA_B_2 = normalize_netmhc_allele(`B2-allele`, "B"),
    HLA_C_1 = normalize_netmhc_allele(C1_allele, "C"),
    HLA_C_2 = normalize_netmhc_allele(`C2-allele`, "C")
  ) %>%
  dplyr::distinct()

write.csv(genotype_tbl, file.path(tab_dir, "HLA_Genotype_for_netMHCpan.csv"), row.names = FALSE)

# ---------------------------
# 6) Build peptide files
# ---------------------------
patients_all <- sort(unique(neo_tbl$patientBarcode))

peptide_qc <- lapply(patients_all, function(pid) {
  pid_safe <- sanitize_id(pid)
  
  peps <- neo_tbl %>%
    dplyr::filter(patientBarcode == pid) %>%
    dplyr::pull(peptide) %>%
    unique()
  
  pep_file <- file.path(peptide_dir, paste0(pid_safe, ".pep.text"))
  
  if (length(peps) > 0) {
    write_peptide_file(peps, pep_file)
  }
  
  data.frame(
    patientBarcode = pid,
    safe_id = pid_safe,
    n_peptides = length(peps),
    peptide_file = pep_file,
    peptide_file_exists = file.exists(pep_file),
    stringsAsFactors = FALSE
  )
})

peptide_qc_tbl <- dplyr::bind_rows(peptide_qc)
write.csv(peptide_qc_tbl, file.path(tab_dir, "peptide_file_qc.csv"), row.names = FALSE)

# ---------------------------
# 7) Build netMHCpan shell scripts
# ---------------------------
manifest <- peptide_qc_tbl %>%
  dplyr::left_join(genotype_tbl, by = "patientBarcode") %>%
  dplyr::mutate(
    has_A = !is.na(HLA_A_1) & !is.na(HLA_A_2),
    has_B = !is.na(HLA_B_1) & !is.na(HLA_B_2),
    has_C = !is.na(HLA_C_1) & !is.na(HLA_C_2)
  )

write.csv(manifest, file.path(tab_dir, "manifest.csv"), row.names = FALSE)

cmd_A <- c()
cmd_B <- c()
cmd_C <- c()

for (i in seq_len(nrow(manifest))) {
  x <- manifest[i, ]
  
  if (!isTRUE(x$peptide_file_exists) || x$n_peptides == 0) next
  
  pid_safe <- x$safe_id
  pep_rel  <- paste0("./peptides/", basename(x$peptide_file))
  
  if (isTRUE(x$has_A)) {
    cmd_A <- c(cmd_A, build_netmhcpan_cmd(c(x$HLA_A_1, x$HLA_A_2), pep_rel, paste0("./A/", pid_safe)))
  }
  if (isTRUE(x$has_B)) {
    cmd_B <- c(cmd_B, build_netmhcpan_cmd(c(x$HLA_B_1, x$HLA_B_2), pep_rel, paste0("./B/", pid_safe)))
  }
  if (isTRUE(x$has_C)) {
    cmd_C <- c(cmd_C, build_netmhcpan_cmd(c(x$HLA_C_1, x$HLA_C_2), pep_rel, paste0("./C/", pid_safe)))
  }
}

writeLines(cmd_A, file.path(script_dir, "Code-Linux-HLA-A.sh"))
writeLines(cmd_B, file.path(script_dir, "Code-Linux-HLA-B.sh"))
writeLines(cmd_C, file.path(script_dir, "Code-Linux-HLA-C.sh"))

# ---------------------------
# 8) Locate raw netMHCpan outputs
# Priority:
#   1) results/neoantigen/A-B-C
#   2) data/neoantigen/A-B-C
# ---------------------------
raw_A_files <- list.files(raw_A_dir, full.names = TRUE)
raw_B_files <- list.files(raw_B_dir, full.names = TRUE)
raw_C_files <- list.files(raw_C_dir, full.names = TRUE)

source_A_dir <- raw_A_dir
source_B_dir <- raw_B_dir
source_C_dir <- raw_C_dir

if (length(raw_A_files) == 0 && dir.exists(fallback_A_dir)) {
  raw_A_files <- list.files(fallback_A_dir, full.names = TRUE)
  source_A_dir <- fallback_A_dir
}

if (length(raw_B_files) == 0 && dir.exists(fallback_B_dir)) {
  raw_B_files <- list.files(fallback_B_dir, full.names = TRUE)
  source_B_dir <- fallback_B_dir
}

if (length(raw_C_files) == 0 && dir.exists(fallback_C_dir)) {
  raw_C_files <- list.files(fallback_C_dir, full.names = TRUE)
  source_C_dir <- fallback_C_dir
}

# ---------------------------
# 9) Parse raw netMHCpan outputs if available
# ---------------------------
parsed_affinity <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(raw_A_files, parse_netmhcpan_file, locus_label = "A")),
  dplyr::bind_rows(lapply(raw_B_files, parse_netmhcpan_file, locus_label = "B")),
  dplyr::bind_rows(lapply(raw_C_files, parse_netmhcpan_file, locus_label = "C"))
)

compare_A <- NULL
compare_B <- NULL
compare_C <- NULL
compare_T <- NULL

if (nrow(parsed_affinity) > 0) {
  write.csv(parsed_affinity, file.path(tab_dir, "parsed_netMHCpan_affinity.csv"), row.names = FALSE)
  
  compare_A <- pivot_compare_csv(parsed_affinity, prev_A, "A")
  compare_B <- pivot_compare_csv(parsed_affinity, prev_B, "B")
  compare_C <- pivot_compare_csv(parsed_affinity, prev_C, "C")
  
  if (!is.null(compare_A)) {
    write.csv(compare_A, file.path(tab_dir, "HLA-A-compare.csv"), row.names = FALSE)
  }
  if (!is.null(compare_B)) {
    write.csv(compare_B, file.path(tab_dir, "HLA-B-compare.csv"), row.names = FALSE)
  }
  if (!is.null(compare_C)) {
    write.csv(compare_C, file.path(tab_dir, "HLA-C-compare.csv"), row.names = FALSE)
  }
  
  compare_T <- parsed_affinity %>%
    dplyr::filter(
      (locus == "A" & allele_display %in% prev_A) |
        (locus == "B" & allele_display %in% prev_B) |
        (locus == "C" & allele_display %in% prev_C)
    ) %>%
    dplyr::group_by(patientBarcode, allele_display) %>%
    dplyr::summarise(value = mean(affinity_nM, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(group = allele_display) %>%
    tidyr::pivot_wider(names_from = group, values_from = value) %>%
    as.data.frame()
  
  write.csv(compare_T, file.path(tab_dir, "HLA-Total-compare.csv"), row.names = FALSE)
}

# ---------------------------
# 10) Affinity analysis
# ---------------------------
aff_A_long <- prepare_affinity_long_from_compare(compare_A, "A")
aff_B_long <- prepare_affinity_long_from_compare(compare_B, "B")
aff_C_long <- prepare_affinity_long_from_compare(compare_C, "C")

affinity_long <- dplyr::bind_rows(aff_A_long, aff_B_long, aff_C_long)

if (nrow(affinity_long) > 0) {
  write.csv(affinity_long, file.path(tab_dir, "affinity_long_from_compare_csvs.csv"), row.names = FALSE)
  
  aff_anova_tbl <- affinity_long %>%
    dplyr::group_by(locus) %>%
    dplyr::group_modify(~data.frame(anova_p = safe_anova(.x, "value"))) %>%
    dplyr::ungroup()
  
  aff_dunn_tbl <- dplyr::bind_rows(lapply(unique(affinity_long$locus), function(loc) {
    d <- affinity_long %>% dplyr::filter(locus == loc)
    x <- safe_dunn(d, "value")
    if (is.null(x)) return(NULL)
    x$locus <- loc
    x
  }))
  
  aff_vs_all_tbl <- dplyr::bind_rows(lapply(unique(affinity_long$locus), function(loc) {
    d <- affinity_long %>% dplyr::filter(locus == loc)
    x <- safe_ttest_vs_all(d, "value")
    x$locus <- loc
    x
  }))
  
  aff_summary_tbl <- affinity_long %>%
    dplyr::group_by(locus, group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(aff_anova_tbl,   file.path(tab_dir, "affinity_anova_by_locus.csv"), row.names = FALSE)
  write.csv(aff_dunn_tbl,    file.path(tab_dir, "affinity_dunn_posthoc.csv"), row.names = FALSE)
  write.csv(aff_vs_all_tbl,  file.path(tab_dir, "affinity_vs_rest_tests.csv"), row.names = FALSE)
  write.csv(aff_summary_tbl, file.path(tab_dir, "affinity_summary_by_allele.csv"), row.names = FALSE)
}

# ---------------------------
# 11) Plot Figure 11E-G
# ---------------------------
if (nrow(affinity_long) > 0) {
  pE <- make_compare_plot(
    affinity_long %>% dplyr::filter(locus == "A"),
    xlab_text = "HLA-A Alleles",
    ylab_text = "Affinity binding",
    label_y = max(affinity_long %>% dplyr::filter(locus == "A") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
  )
  
  pF <- make_compare_plot(
    affinity_long %>% dplyr::filter(locus == "B"),
    xlab_text = "HLA-B Alleles",
    ylab_text = "Affinity binding",
    label_y = max(affinity_long %>% dplyr::filter(locus == "B") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
  )
  
  pG <- make_compare_plot(
    affinity_long %>% dplyr::filter(locus == "C"),
    xlab_text = "HLA-C Alleles",
    ylab_text = "Affinity binding",
    label_y = max(affinity_long %>% dplyr::filter(locus == "C") %>% dplyr::pull(value), na.rm = TRUE) * 1.12
  )
} else {
  pE <- ggplot() + theme_void() + labs(title = "No affinity data yet")
  pF <- ggplot() + theme_void() + labs(title = "No affinity data yet")
  pG <- ggplot() + theme_void() + labs(title = "No affinity data yet")
}

figure11_EG <- cowplot::plot_grid(
  pE, pF, pG,
  labels = c("E", "F", "G"),
  ncol = 3,
  label_size = 16
)

ggsave(
  filename = file.path(fig_dir, "Figure_11E_G_Affinity.png"),
  plot = figure11_EG,
  width = 14,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure_11E_G_Affinity.pdf"),
  plot = figure11_EG,
  width = 14,
  height = 6,
  bg = "white"
)

# ---------------------------
# 12) Log
# ---------------------------
sink(file.path(log_dir, "03_neoantigen_affinity_log.txt"))
cat("03_neoantigen_affinity.R completed\n")
cat("ASE input: ", ase_file, "\n", sep = "")
cat("Neoantigen input: ", neo_file, "\n", sep = "")
cat("Prevalent allele file: ", prev_file, "\n", sep = "")
cat("Peptide files dir: ", peptide_dir, "\n", sep = "")
cat("Linux scripts dir: ", script_dir, "\n", sep = "")
cat("Fallback raw A dir: ", fallback_A_dir, "\n", sep = "")
cat("Fallback raw B dir: ", fallback_B_dir, "\n", sep = "")
cat("Fallback raw C dir: ", fallback_C_dir, "\n", sep = "")
cat("Raw A source dir used: ", source_A_dir, "\n", sep = "")
cat("Raw B source dir used: ", source_B_dir, "\n", sep = "")
cat("Raw C source dir used: ", source_C_dir, "\n", sep = "")
cat("Raw A outputs found: ", length(raw_A_files), "\n", sep = "")
cat("Raw B outputs found: ", length(raw_B_files), "\n", sep = "")
cat("Raw C outputs found: ", length(raw_C_files), "\n", sep = "")
cat("Affinity compare CSVs built: ", ifelse(nrow(parsed_affinity) > 0, "YES", "NO - no raw netMHCpan outputs found"), "\n", sep = "")
sink()

cat("Done: 03_neoantigen_affinity.R\n")