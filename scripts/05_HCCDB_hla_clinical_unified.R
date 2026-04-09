# ==========================================================
# Script: 05_HCCDB_hla_clinical_unified.R
# Purpose:
#   Unified clinical + survival + HCC-vs-Adjacent analysis
#   using the previous logic with updated style and paths.
#
# Output:
#   results/clinical_and_survival_unified/
#     - tables/
#         - HCC_vs_Adjacent_Summary.csv
#         - Clinical_Continuous_Associations.csv
#         - Clinical_Grouped_Associations.csv
#         - Survival_OS_RFS_Summary.csv
#         - Clinical_and_Survival_Unified.xlsx
#     - logs/
#         - 10_clinical_and_survival_unified_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "readxl", "openxlsx", "dplyr", "stringr", "survival"
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
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(survival)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "G:/HLA-I-Multiomics-HCC/"

clinical_dir <- file.path(project_root, "data", "bulk_RNA_seq", "clinical")
boxplot_dir  <- file.path(project_root, "data", "bulk_RNA_seq", "boxplot", "Other-Data-set", "3 Group")

file_survival_os  <- file.path(clinical_dir, "Survival_Analysis1.xlsx")
file_survival_rfs <- file.path(clinical_dir, "Survival_Analysis2.xlsx")

file_bclc    <- file.path(clinical_dir, "BCLC.xlsx")
file_gender  <- file.path(clinical_dir, "Gender.xlsx")
file_tnm     <- file.path(clinical_dir, "Data-TNM-AJCC.xlsx")
file_virus   <- file.path(clinical_dir, "Virus.xlsx")

file_age     <- file.path(clinical_dir, "Age.xlsx")
file_afp     <- file.path(clinical_dir, "AFP.xlsx")
file_maxsize <- file.path(clinical_dir, "Maxmumsize.xlsx")
file_ploidy  <- file.path(clinical_dir, "poloidy.xlsx")
file_purity  <- file.path(clinical_dir, "purity.xlsx")
file_tmb     <- file.path(clinical_dir, "TMB.xlsx")
file_tdt     <- file.path(clinical_dir, "Tumor_double_time.xlsx")

results_root <- file.path(project_root, "results", "clinical_and_survival_unified")
table_dir <- file.path(results_root, "tables")
log_dir <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 2) Common definitions
# ---------------------------
hla_vars_survival <- c(
  "HLA-A.Expression",
  "HLA-B.Expression",
  "HLA-C.Expression",
  "Log.Ratio.(HLA-A.to.HLA-B.Expression)",
  "Log.Ratio.(HLA-A.to.HLA-C.Expression)",
  "Log.Ratio.(HLA-B.to.HLA-C.Expression)"
)

hla_vars_clinical <- c("HLA-A", "HLA-B", "HLA-C", "AB", "AC", "BC")

hla_label_map <- c(
  "HLA-A.Expression" = "HLA-A",
  "HLA-B.Expression" = "HLA-B",
  "HLA-C.Expression" = "HLA-C",
  "Log.Ratio.(HLA-A.to.HLA-B.Expression)" = "AB",
  "Log.Ratio.(HLA-A.to.HLA-C.Expression)" = "AC",
  "Log.Ratio.(HLA-B.to.HLA-C.Expression)" = "BC",
  "HLA-A" = "HLA-A",
  "HLA-B" = "HLA-B",
  "HLA-C" = "HLA-C",
  "AB" = "AB",
  "AC" = "AC",
  "BC" = "BC"
)

dataset_sources <- c(
  "HCCDB.1"  = "GSE22058",
  "HCCDB.4"  = "GSE36376",
  "HCCDB.8"  = "GSE9843",
  "HCCDB.11" = "GSE46444",
  "HCCDB.12" = "GSE54236",
  "HCCDB.13" = "GSE63898",
  "HCCDB.14" = "GSE43619",
  "HCCDB.15" = "TCGA-LIHC",
  "HCCDB.17" = "GSE76427",
  "HCCDB.18" = "ICGC-LIRI-JP",
  "HCCDB.19" = "GSE124751",
  "HCCDB.20" = "GSE164760",
  "HCCDB.21" = "GSE134568",
  "HCCDB.22" = "GSE87630",
  "HCCDB.24" = "GSE76427",
  "HCCDB.25" = "OEP000321",
  "HCCDB.26" = "GSE112790",
  "HCCDB.27" = "GSE121248",
  "HCCDB.28" = "GSE109211",
  "HCCDB.29" = "GSE89377",
  "HCCDB.30" = "GSE148355",
  "HCCDB15"  = "TCGA-LIHC",
  "HCCDB17"  = "GSE76427",
  "HCCDB18"  = "ICGC-LIRI-JP",
  "HCCDB19"  = "GSE124751",
  "HCCDB24"  = "GSE76427",
  "HCCDB25"  = "OEP000321",
  "HCCDB30"  = "GSE148355"
)

# ---------------------------
# 3) Helpers
# ---------------------------
safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) {
    return(data.frame(Correlation = NA_real_, P.value = NA_real_))
  }
  test <- tryCatch(cor.test(x, y, method = method), error = function(e) NULL)
  if (is.null(test)) {
    return(data.frame(Correlation = NA_real_, P.value = NA_real_))
  }
  data.frame(
    Correlation = unname(test$estimate),
    P.value = test$p.value
  )
}

safe_source <- function(x) {
  ifelse(x %in% names(dataset_sources), dataset_sources[[x]], x)
}

save_csv_if_not_null <- function(df, path) {
  if (!is.null(df) && nrow(df) > 0) {
    write.csv(df, path, row.names = FALSE)
  }
}

# ---------------------------
# 4) HCC vs Adjacent
# ---------------------------
build_hcc_adjacent_table <- function(folder_path) {
  file_names <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  file_names <- file_names[!grepl("^~\\$", basename(file_names))]
  
  if (length(file_names) == 0) {
    return(NULL)
  }
  
  marker_map <- list(
    list(column_index = 4, marker_name = "HLA-A"),
    list(column_index = 5, marker_name = "HLA-B"),
    list(column_index = 6, marker_name = "HLA-C"),
    list(column_index = 7, marker_name = "AB"),
    list(column_index = 8, marker_name = "AC"),
    list(column_index = 9, marker_name = "BC")
  )
  
  process_file <- function(file_path, column_index, marker_name) {
    dat <- read.xlsx(file_path)
    ds <- str_extract(basename(file_path), "HCCDB\\.[0-9]{1,2}")
    
    hcc_vals <- dat[dat$Type == "HCC", column_index]
    adj_vals <- dat[dat$Type == "Adjacent", column_index]
    
    data.frame(
      Section = "HCC_vs_Adjacent",
      Metric = marker_name,
      Dataset = ds,
      Source = safe_source(ds),
      N_Adjacent = sum(dat$Type == "Adjacent", na.rm = TRUE),
      N_HCC = sum(dat$Type == "HCC", na.rm = TRUE),
      Mean_Adjacent = mean(adj_vals, na.rm = TRUE),
      Mean_HCC = mean(hcc_vals, na.rm = TRUE),
      P.value = safe_wilcox(adj_vals, hcc_vals),
      Status = ifelse(
        mean(hcc_vals, na.rm = TRUE) > mean(adj_vals, na.rm = TRUE),
        "Upregulated",
        "Downregulated"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(lapply(marker_map, function(mk) {
    bind_rows(lapply(file_names, process_file, mk$column_index, mk$marker_name))
  }))
}

# ---------------------------
# 5) Continuous clinical features
# ---------------------------
run_continuous_feature <- function(file_path, feature_name, sheet_feature_name = NULL) {
  if (!file.exists(file_path)) return(NULL)
  
  sheets <- excel_sheets(file_path)
  
  out <- lapply(sheets, function(sheet) {
    df <- read.xlsx(file_path, sheet = sheet)
    
    feature_col <- if (!is.null(sheet_feature_name) && sheet_feature_name %in% colnames(df)) {
      sheet_feature_name
    } else {
      tail(colnames(df), 1)
    }
    
    bind_rows(lapply(hla_vars_clinical, function(metric) {
      if (!(metric %in% colnames(df)) || !(feature_col %in% colnames(df))) return(NULL)
      res <- safe_cor(df[[feature_col]], df[[metric]], method = "pearson")
      data.frame(
        Section = "Clinical_Continuous",
        Feature = feature_name,
        Dataset = sheet,
        Source = safe_source(sheet),
        Metric = metric,
        Feature_Column = feature_col,
        Correlation = res$Correlation,
        P.value = res$P.value,
        stringsAsFactors = FALSE
      )
    }))
  })
  
  bind_rows(out)
}

# ---------------------------
# 6) Grouped clinical features
# ---------------------------
run_grouped_feature <- function(file_path, feature_name, groups) {
  if (!file.exists(file_path)) return(NULL)
  
  sheets <- excel_sheets(file_path)
  data_columns <- 2:7
  
  out <- lapply(sheets, function(sheet) {
    df <- read_excel(file_path, sheet = sheet)
    if (ncol(df) < 8) return(NULL)
    
    group_column <- df[[8]]
    
    bind_rows(lapply(data_columns, function(data_col_idx) {
      data_column <- df[[data_col_idx]]
      metric_name <- names(metric_columns <- c(
        "HLA-A" = 2, "HLA-B" = 3, "HLA-C" = 4,
        "AB" = 5, "AC" = 6, "BC" = 7
      ))[match(data_col_idx, unname(metric_columns))]
      
      if (is.na(metric_name)) return(NULL)
      
      counts <- sapply(groups, function(g) sum(group_column == g, na.rm = TRUE))
      means  <- sapply(groups, function(g) mean(data_column[group_column == g], na.rm = TRUE))
      
      comps <- combn(groups, 2, simplify = FALSE)
      
      bind_rows(lapply(comps, function(comp) {
        data.frame(
          Section = "Clinical_Grouped",
          Feature = feature_name,
          Dataset = sheet,
          Source = safe_source(sheet),
          Metric = metric_name,
          Group1 = comp[1],
          Group2 = comp[2],
          N_Group1 = counts[[comp[1]]],
          N_Group2 = counts[[comp[2]]],
          Mean_Group1 = means[[comp[1]]],
          Mean_Group2 = means[[comp[2]]],
          Difference = means[[comp[2]]] - means[[comp[1]]],
          P.value = safe_wilcox(
            data_column[group_column == comp[1]],
            data_column[group_column == comp[2]]
          ),
          stringsAsFactors = FALSE
        )
      }))
    }))
  })
  
  bind_rows(out)
}

# ---------------------------
# 7) Virus feature
# ---------------------------
run_virus_feature <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  
  sheets <- c("HCCDB15", "HCCDB18", "HCCDB19", "HCCDB30")
  
  out <- lapply(sheets, function(sheet) {
    df <- read.xlsx(file_path, sheet = sheet)
    
    bind_rows(lapply(names(c("HLA-A" = 3, "HLA-B" = 4, "HLA-C" = 5, "AB" = 6, "AC" = 7, "BC" = 8)), function(metric) {
      col_idx <- c("HLA-A" = 3, "HLA-B" = 4, "HLA-C" = 5, "AB" = 6, "AC" = 7, "BC" = 8)[[metric]]
      
      hbv  <- df[df$Virus == "HBV", col_idx]
      hcv  <- df[df$Virus == "HCV", col_idx]
      free <- df[df$Virus == "Free", col_idx]
      
      bind_rows(
        data.frame(
          Section = "Clinical_Grouped",
          Feature = "Virus",
          Dataset = sheet,
          Source = safe_source(sheet),
          Metric = metric,
          Group1 = "HBV",
          Group2 = "HCV",
          N_Group1 = sum(df$Virus == "HBV", na.rm = TRUE),
          N_Group2 = sum(df$Virus == "HCV", na.rm = TRUE),
          Mean_Group1 = mean(hbv, na.rm = TRUE),
          Mean_Group2 = mean(hcv, na.rm = TRUE),
          Difference = mean(hcv, na.rm = TRUE) - mean(hbv, na.rm = TRUE),
          P.value = safe_wilcox(hbv, hcv),
          stringsAsFactors = FALSE
        ),
        data.frame(
          Section = "Clinical_Grouped",
          Feature = "Virus",
          Dataset = sheet,
          Source = safe_source(sheet),
          Metric = metric,
          Group1 = "HBV",
          Group2 = "Free",
          N_Group1 = sum(df$Virus == "HBV", na.rm = TRUE),
          N_Group2 = sum(df$Virus == "Free", na.rm = TRUE),
          Mean_Group1 = mean(hbv, na.rm = TRUE),
          Mean_Group2 = mean(free, na.rm = TRUE),
          Difference = mean(free, na.rm = TRUE) - mean(hbv, na.rm = TRUE),
          P.value = safe_wilcox(hbv, free),
          stringsAsFactors = FALSE
        ),
        data.frame(
          Section = "Clinical_Grouped",
          Feature = "Virus",
          Dataset = sheet,
          Source = safe_source(sheet),
          Metric = metric,
          Group1 = "HCV",
          Group2 = "Free",
          N_Group1 = sum(df$Virus == "HCV", na.rm = TRUE),
          N_Group2 = sum(df$Virus == "Free", na.rm = TRUE),
          Mean_Group1 = mean(hcv, na.rm = TRUE),
          Mean_Group2 = mean(free, na.rm = TRUE),
          Difference = mean(free, na.rm = TRUE) - mean(hcv, na.rm = TRUE),
          P.value = safe_wilcox(hcv, free),
          stringsAsFactors = FALSE
        )
      )
    }))
  })
  
  bind_rows(out)
}

# ---------------------------
# 8) Survival
# ---------------------------
run_survival_summary <- function(data, dataset_name, variables, time_col, status_col, endpoint_name) {
  results <- list()
  
  for (variable in variables) {
    median_value <- median(data[[variable]], na.rm = TRUE)
    
    grouped_data <- data %>%
      mutate(Group = ifelse(.data[[variable]] > median_value, "High", "Low"))
    
    if (endpoint_name == "OS") {
      HCC <- grouped_data[, c(time_col, status_col, "Group")]
    } else {
      HCC <- grouped_data[, c(time_col, status_col, "Group")]
      HCC[[status_col]] <- ifelse(HCC[[status_col]] == "Yes", 1, HCC[[status_col]])
    }
    
    fit <- survfit(Surv(HCC[[time_col]], HCC[[status_col]]) ~ Group, data = HCC)
    cox_model <- coxph(Surv(HCC[[time_col]], HCC[[status_col]]) ~ Group, data = HCC)
    summary_cox <- summary(cox_model)
    
    results[[variable]] <- data.frame(
      Section = "Survival",
      Endpoint = endpoint_name,
      Metric = hla_label_map[[variable]],
      Dataset = dataset_name,
      Source = safe_source(dataset_name),
      LogHR = round(coef(summary_cox)[1], 3),
      P.value = round(summary_cox$coefficients[, "Pr(>|z|)"][1], 3),
      CI.Lower = round(confint(cox_model)[1, 1], 3),
      CI.Upper = round(confint(cox_model)[1, 2], 3),
      N.High = sum(grouped_data$Group == "High"),
      N.Low = sum(grouped_data$Group == "Low"),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(results)
}

build_survival_table <- function(file_os, file_rfs) {
  os_datasets  <- c("HCCDB15", "HCCDB17", "HCCDB18", "HCCDB19", "HCCDB24", "HCCDB25", "HCCDB30")
  rfs_datasets <- c("HCCDB19", "HCCDB24", "HCCDB25")
  
  os_res <- lapply(os_datasets, function(ds) {
    dat <- read.xlsx(file_os, sheet = ds)
    run_survival_summary(dat, ds, hla_vars_survival, "OS(Month)", "Status", "OS")
  }) %>% bind_rows()
  
  rfs_res <- lapply(rfs_datasets, function(ds) {
    dat <- read.xlsx(file_rfs, sheet = ds)
    dat$`RFS(Status)` <- ifelse(dat$`RFS(Status)` == "Yes", 1, 0)
    run_survival_summary(dat, ds, hla_vars_survival, "RFS", "RFS(Status)", "RFS")
  }) %>% bind_rows()
  
  bind_rows(os_res, rfs_res)
}

# ---------------------------
# 9) Run all sections
# ---------------------------
hcc_adjacent_results <- build_hcc_adjacent_table(boxplot_dir)

continuous_results <- bind_rows(
  run_continuous_feature(file_age, "Age", "Age"),
  run_continuous_feature(file_afp, "AFP", "AFP"),
  run_continuous_feature(file_purity, "Purity", "Purity"),
  run_continuous_feature(file_ploidy, "Ploidy", "Ploidy"),
  run_continuous_feature(file_tmb, "TMB", "TMB"),
  run_continuous_feature(file_tdt, "Tumor_Double_Time", "TUMOR DOUBLE TIME"),
  run_continuous_feature(file_maxsize, "Maximum_Tumor_Size", "Maxmumsize")
)

grouped_results <- bind_rows(
  run_grouped_feature(file_gender, "Gender", c("F", "M")),
  run_grouped_feature(file_bclc, "BCLC", c("0", "A", "B", "C")),
  run_grouped_feature(file_tnm, "TNM", c("I", "II", "III")),
  run_virus_feature(file_virus)
)

survival_results <- build_survival_table(file_survival_os, file_survival_rfs)

# ---------------------------
# 10) Save outputs
# ---------------------------
save_csv_if_not_null(hcc_adjacent_results, file.path(table_dir, "HCC_vs_Adjacent_Summary.csv"))
save_csv_if_not_null(continuous_results, file.path(table_dir, "Clinical_Continuous_Associations.csv"))
save_csv_if_not_null(grouped_results, file.path(table_dir, "Clinical_Grouped_Associations.csv"))
save_csv_if_not_null(survival_results, file.path(table_dir, "Survival_OS_RFS_Summary.csv"))

wb <- createWorkbook()

if (!is.null(hcc_adjacent_results) && nrow(hcc_adjacent_results) > 0) {
  addWorksheet(wb, "HCC_vs_Adjacent")
  writeData(wb, "HCC_vs_Adjacent", hcc_adjacent_results)
}

if (!is.null(continuous_results) && nrow(continuous_results) > 0) {
  addWorksheet(wb, "Clinical_Continuous")
  writeData(wb, "Clinical_Continuous", continuous_results)
}

if (!is.null(grouped_results) && nrow(grouped_results) > 0) {
  addWorksheet(wb, "Clinical_Grouped")
  writeData(wb, "Clinical_Grouped", grouped_results)
}

if (!is.null(survival_results) && nrow(survival_results) > 0) {
  addWorksheet(wb, "Survival")
  writeData(wb, "Survival", survival_results)
}

saveWorkbook(
  wb,
  file.path(table_dir, "Clinical_and_Survival_Unified.xlsx"),
  overwrite = TRUE
)

# ---------------------------
# 11) Log
# ---------------------------
sink(file.path(log_dir, "10_clinical_and_survival_unified_log.txt"))
cat("10_clinical_and_survival_unified.R completed\n")
cat("Clinical directory: ", clinical_dir, "\n", sep = "")
cat("Boxplot directory: ", boxplot_dir, "\n", sep = "")
cat("Results directory: ", results_root, "\n", sep = "")
cat("HCC vs Adjacent rows: ", ifelse(is.null(hcc_adjacent_results), 0, nrow(hcc_adjacent_results)), "\n", sep = "")
cat("Continuous rows: ", ifelse(is.null(continuous_results), 0, nrow(continuous_results)), "\n", sep = "")
cat("Grouped rows: ", ifelse(is.null(grouped_results), 0, nrow(grouped_results)), "\n", sep = "")
cat("Survival rows: ", ifelse(is.null(survival_results), 0, nrow(survival_results)), "\n", sep = "")
sink()

cat("Done: 10_clinical_and_survival_unified.R\n")