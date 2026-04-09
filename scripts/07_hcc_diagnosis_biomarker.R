# ==========================================================
# Script: 07_hcc_diagnosis_biomarker.R
# Purpose:
#   Build and evaluate logistic regression models for HCC
#   diagnosis using HLA-I expression ratios, including ROC,
#   PR-AUC, calibration, Brier score, decision curve analysis,
#   threshold-based operating points, Figure 5 single-ratio
#   models, and external dataset evaluation.
#
# Input:
#   data/model/normalized_batch_corrected_data.xlsx
#   data/model/external/*.xlsx
#
# Output:
#   results/hcc_diagnosis_logistic_regression/tables/
#     - model_fit/*
#     - thresholds/*
#     - calibration/*
#     - decision_curve/*
#     - external_validation/*
#     - figure5/*
#
#   results/hcc_diagnosis_logistic_regression/figures/
#     - roc_pr/*
#     - calibration/*
#     - decision_curve/*
#     - external_validation/*
#     - figure5/*
#
#   results/hcc_diagnosis_logistic_regression/logs/
#     - sessionInfo_model.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "caret", "readxl", "dplyr", "pROC", "PRROC",
  "ggplot2", "patchwork"
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
  library(caret)
  library(readxl)
  library(dplyr)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."

data_dir <- file.path(project_root, "data", "model")
external_dir <- file.path(data_dir, "external")
input_file <- file.path(data_dir, "model_data.xlsx")

results_root <- file.path(project_root, "results", "hcc_diagnosis_logistic_regression")
tables_root <- file.path(results_root, "tables")
figures_root <- file.path(results_root, "figures")
logs_root <- file.path(results_root, "logs")

tables_model_fit_dir <- file.path(tables_root, "model_fit")
tables_thresholds_dir <- file.path(tables_root, "thresholds")
tables_calibration_dir <- file.path(tables_root, "calibration")
tables_decision_curve_dir <- file.path(tables_root, "decision_curve")
tables_external_dir <- file.path(tables_root, "external_validation")
tables_figure5_dir <- file.path(tables_root, "figure5")

figures_roc_pr_dir <- file.path(figures_root, "roc_pr")
figures_calibration_dir <- file.path(figures_root, "calibration")
figures_decision_curve_dir <- file.path(figures_root, "decision_curve")
figures_external_dir <- file.path(figures_root, "external_validation")
figures_figure5_dir <- file.path(figures_root, "figure5")

output_dirs <- c(
  results_root,
  tables_root,
  figures_root,
  logs_root,
  tables_model_fit_dir,
  tables_thresholds_dir,
  tables_calibration_dir,
  tables_decision_curve_dir,
  tables_external_dir,
  tables_figure5_dir,
  figures_roc_pr_dir,
  figures_calibration_dir,
  figures_decision_curve_dir,
  figures_external_dir,
  figures_figure5_dir
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

if (!dir.exists(external_dir)) {
  warning("External directory not found: ", external_dir)
}

external_files <- sort(list.files(
  external_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE
))

# ---------------------------
# 2) Reproducibility
# ---------------------------
set.seed(42)

# ---------------------------
# 3) Validation helpers
# ---------------------------
check_required_columns <- function(df, required_cols, df_name = "data frame") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", df_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

get_positive_class <- function(y_factor) {
  lvl <- levels(y_factor)
  if ("Tumor" %in% lvl) return("Tumor")
  if ("HCC" %in% lvl) return("HCC")
  return(lvl[2])
}

safe_factor_match <- function(y, ref_levels) {
  factor(as.character(y), levels = ref_levels)
}

# ---------------------------
# 4) Modeling helpers
# ---------------------------
fit_logistic_model <- function(x, y) {
  x <- as.data.frame(x)
  
  train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  train(
    x = x,
    y = y,
    method = "glm",
    family = "binomial",
    trControl = train_control,
    preProc = c("center", "scale"),
    metric = "ROC"
  )
}

compute_binary_metrics <- function(model, x_data, y_data, positive_class) {
  x_data <- as.data.frame(x_data)
  
  predictions <- predict(model, x_data)
  
  conf_mat <- confusionMatrix(
    data = predictions,
    reference = y_data,
    positive = positive_class
  )
  
  probabilities <- predict(model, x_data, type = "prob")[, positive_class]
  
  roc_obj <- roc(
    response = y_data,
    predictor = probabilities,
    levels = levels(y_data),
    direction = "<"
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  labels <- ifelse(y_data == positive_class, 1, 0)
  
  pr_obj <- pr.curve(
    scores.class0 = probabilities[labels == 1],
    scores.class1 = probabilities[labels == 0],
    curve = TRUE
  )
  
  pr_auc <- as.numeric(pr_obj$auc.integral)
  brier <- mean((labels - probabilities)^2)
  
  sens <- as.numeric(conf_mat$byClass["Sensitivity"])
  spec <- as.numeric(conf_mat$byClass["Specificity"])
  precision <- as.numeric(conf_mat$byClass["Precision"])
  recall <- as.numeric(conf_mat$byClass["Recall"])
  f1 <- as.numeric(conf_mat$byClass["F1"])
  bal_acc <- as.numeric(conf_mat$byClass["Balanced Accuracy"])
  accuracy <- as.numeric(conf_mat$overall["Accuracy"])
  kappa <- as.numeric(conf_mat$overall["Kappa"])
  
  list(
    roc = roc_obj,
    auc = auc_val,
    pr = pr_obj,
    pr_auc = pr_auc,
    brier = brier,
    sensitivity = sens,
    specificity = spec,
    precision = precision,
    recall = recall,
    f1 = f1,
    balanced_accuracy = bal_acc,
    accuracy = accuracy,
    kappa = kappa,
    confusion_matrix = conf_mat,
    probabilities = probabilities,
    labels = labels
  )
}

find_best_threshold_youden <- function(probabilities, labels) {
  roc_obj <- pROC::roc(labels, probabilities, direction = "<")
  
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  data.frame(
    Threshold = as.numeric(best[, "threshold"]),
    Sensitivity = as.numeric(best[, "sensitivity"]),
    Specificity = as.numeric(best[, "specificity"])
  )
}

compute_metrics_at_fixed_threshold <- function(probabilities, labels, threshold) {
  labels <- as.numeric(labels)
  pred <- ifelse(probabilities >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & labels == 1, na.rm = TRUE)
  tn <- sum(pred == 0 & labels == 0, na.rm = TRUE)
  fp <- sum(pred == 1 & labels == 0, na.rm = TRUE)
  fn <- sum(pred == 0 & labels == 1, na.rm = TRUE)
  
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- ppv
  recall <- sensitivity
  f1 <- ifelse(
    !is.na(precision) && !is.na(recall) && (precision + recall) > 0,
    2 * precision * recall / (precision + recall),
    NA
  )
  
  data.frame(
    Threshold = threshold,
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    NPV = npv,
    Accuracy = accuracy,
    F1 = f1
  )
}

compute_metrics_multiple_fixed_thresholds <- function(probabilities, labels, thresholds) {
  bind_rows(
    lapply(thresholds, function(thr) {
      compute_metrics_at_fixed_threshold(probabilities, labels, thr)
    })
  )
}

# ---------------------------
# 5) Plot helpers
# ---------------------------
add_metrics_annotation <- function(auc, pr_auc, sens, spec) {
  usr <- par("usr")
  
  x_pos <- usr[1] + 0.58 * (usr[2] - usr[1])
  y_pos <- usr[3] + 0.16 * (usr[4] - usr[3])
  
  label_text <- paste0(
    "ROC AUC = ", sprintf("%.3f", auc), "\n",
    "PR AUC = ", sprintf("%.3f", pr_auc), "\n",
    "Sens = ", sprintf("%.2f", sens), "\n",
    "Spec = ", sprintf("%.2f", spec)
  )
  
  text(
    x = x_pos,
    y = y_pos,
    labels = label_text,
    adj = c(0, 0),
    cex = 0.92,
    col = "blue"
  )
}

save_single_roc_plot <- function(roc_obj, auc_val, pr_auc, sens, spec, file_name, main_title) {
  tiff(
    filename = file_name,
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  par(mar = c(5, 5, 3, 2))
  plot(
    roc_obj,
    main = main_title,
    col = "#C98500",
    lwd = 3,
    legacy.axes = TRUE,
    xlab = "Specificity",
    ylab = "Sensitivity"
  )
  
  add_metrics_annotation(
    auc = auc_val,
    pr_auc = pr_auc,
    sens = sens,
    spec = spec
  )
  
  dev.off()
}

save_single_pr_plot <- function(pr_obj, file_name, main_title) {
  tiff(
    filename = file_name,
    width = 7,
    height = 7,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  par(mar = c(5, 5, 3, 2))
  plot(pr_obj, main = main_title)
  dev.off()
}

build_calibration_summary <- function(probabilities, labels, n_bins = 10) {
  calibration_df <- data.frame(
    Observed = labels,
    Predicted = probabilities
  )
  
  quantile_breaks <- unique(
    quantile(
      calibration_df$Predicted,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    )
  )
  
  if (length(quantile_breaks) < 3) {
    quantile_breaks <- seq(0, 1, length.out = n_bins + 1)
  }
  
  calibration_df$Bin <- cut(
    calibration_df$Predicted,
    breaks = quantile_breaks,
    include.lowest = TRUE
  )
  
  calibration_df %>%
    group_by(Bin) %>%
    summarise(
      Mean_Predicted = mean(Predicted, na.rm = TRUE),
      Observed_Fraction = mean(Observed, na.rm = TRUE),
      Count = n(),
      .groups = "drop"
    )
}

save_calibration_plot <- function(calibration_summary, file_prefix, plot_title) {
  calibration_plot <- ggplot(calibration_summary, aes(x = Mean_Predicted, y = Observed_Fraction)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = plot_title,
      x = "Mean Predicted Probability",
      y = "Observed Fraction of Tumor"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = calibration_plot,
    width = 7,
    height = 7,
    dpi = 600
  )
  
  ggsave(
    filename = paste0(file_prefix, ".tiff"),
    plot = calibration_plot,
    width = 7,
    height = 7,
    dpi = 600,
    compression = "lzw"
  )
}

build_decision_curve_data <- function(
    probabilities,
    labels,
    thresholds = seq(0.01, 0.50, by = 0.01),
    n_boot = 200,
    seed = 42
) {
  set.seed(seed)
  
  n <- length(labels)
  prevalence <- mean(labels == 1, na.rm = TRUE)
  
  model_df <- lapply(thresholds, function(pt) {
    pred_positive <- ifelse(probabilities >= pt, 1, 0)
    
    tp <- sum(pred_positive == 1 & labels == 1, na.rm = TRUE)
    fp <- sum(pred_positive == 1 & labels == 0, na.rm = TRUE)
    
    net_benefit_model <- (tp / n) - (fp / n) * (pt / (1 - pt))
    net_benefit_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
    net_benefit_none <- 0
    
    data.frame(
      Threshold = pt,
      Model = net_benefit_model,
      Treat_All = net_benefit_all,
      Treat_None = net_benefit_none
    )
  }) %>%
    bind_rows()
  
  boot_mat <- sapply(thresholds, function(pt) {
    replicate(n_boot, {
      idx <- sample(seq_len(n), size = n, replace = TRUE)
      probs_b <- probabilities[idx]
      labels_b <- labels[idx]
      
      pred_positive <- ifelse(probs_b >= pt, 1, 0)
      tp <- sum(pred_positive == 1 & labels_b == 1, na.rm = TRUE)
      fp <- sum(pred_positive == 1 & labels_b == 0, na.rm = TRUE)
      
      (tp / n) - (fp / n) * (pt / (1 - pt))
    })
  })
  
  model_df$Model_Lower <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  model_df$Model_Upper <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  model_df
}

save_decision_curve_plot <- function(dca_df, file_prefix, plot_title) {
  dca_plot <- ggplot(dca_df, aes(x = Threshold)) +
    geom_ribbon(
      aes(ymin = Model_Lower, ymax = Model_Upper),
      fill = "#1b9e77",
      alpha = 0.25
    ) +
    geom_line(aes(y = Model, color = "Model"), linewidth = 1.6) +
    geom_line(aes(y = Treat_All, color = "Treat all"), linewidth = 1.4) +
    geom_line(aes(y = Treat_None, color = "Treat none"), linewidth = 1.4, linetype = "dashed") +
    scale_color_manual(
      values = c(
        "Model" = "#1b9e77",
        "Treat all" = "black",
        "Treat none" = "gray40"
      )
    ) +
    scale_x_continuous(
      limits = c(0, 0.50),
      breaks = seq(0, 0.50, by = 0.10),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    labs(
      title = plot_title,
      x = "Decision threshold",
      y = "Net Benefit",
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 12),
      panel.grid.minor = element_line(color = "gray90"),
      panel.grid.major = element_line(color = "gray85")
    )
  
  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = dca_plot,
    width = 10,
    height = 6,
    dpi = 600
  )
  
  ggsave(
    filename = paste0(file_prefix, ".tiff"),
    plot = dca_plot,
    width = 10,
    height = 6,
    dpi = 600,
    compression = "lzw"
  )
  
  invisible(dca_plot)
}

summarize_net_benefit <- function(dca_df, selected_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  dca_df %>%
    filter(round(Threshold, 2) %in% selected_thresholds)
}

save_multi_panel_roc <- function(results_list, output_file) {
  tiff(
    filename = output_file,
    width = 16,
    height = 9,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  par(
    mfrow = c(2, 4),
    mar = c(5.0, 5.0, 2.2, 1.5),
    oma = c(0.2, 0.2, 0.2, 0.2),
    cex.axis = 1.05,
    cex.lab = 1.10
  )
  
  for (nm in names(results_list)) {
    res <- results_list[[nm]]
    
    plot(
      res$roc,
      main = "",
      col = "#C98500",
      lwd = 3,
      legacy.axes = TRUE,
      xlab = "Specificity",
      ylab = "Sensitivity"
    )
    
    mtext(
      text = nm,
      side = 3,
      adj = 0.00,
      line = 1.00,
      cex = 1.35,
      font = 2
    )
    
    add_metrics_annotation(
      auc = res$auc,
      pr_auc = res$pr_auc,
      sens = res$sensitivity,
      spec = res$specificity
    )
  }
  
  dev.off()
}

# ---------------------------
# 6) External evaluation helper
# ---------------------------
evaluate_external_dataset <- function(file_path, model, selected_normals_df, train_levels, figures_external_dir, output_prefix = NULL) {
  if (!file.exists(file_path)) {
    warning("External dataset not found: ", file_path)
    return(NULL)
  }
  
  external_data <- readxl::read_xlsx(file_path) %>%
    as.data.frame()
  
  required_cols_external <- c("AB", "AC", "BC")
  missing_external <- setdiff(required_cols_external, colnames(external_data))
  if (length(missing_external) > 0) {
    stop(
      paste0(
        "Missing required columns in external dataset ",
        basename(file_path),
        ": ",
        paste(missing_external, collapse = ", ")
      )
    )
  }
  
  external_data$Group <- "Tumor"
  
  tumor_part <- external_data[, c("AB", "AC", "BC", "Group"), drop = FALSE]
  normal_part <- selected_normals_df[, c("AB", "AC", "BC", "Group"), drop = FALSE]
  
  combined_test <- rbind(tumor_part, normal_part)
  
  y_ext <- as.factor(combined_test$Group)
  y_ext <- factor(y_ext, levels = train_levels)
  positive_class_ext <- get_positive_class(y_ext)
  
  x_ext <- combined_test[, c("AB", "AC", "BC"), drop = FALSE]
  
  ext_results <- compute_binary_metrics(
    model = model,
    x_data = x_ext,
    y_data = y_ext,
    positive_class = positive_class_ext
  )
  
  result_row <- data.frame(
    Dataset = basename(file_path),
    Accuracy = ext_results$accuracy,
    Kappa = ext_results$kappa,
    Sensitivity = ext_results$sensitivity,
    Specificity = ext_results$specificity,
    Precision = ext_results$precision,
    Recall = ext_results$recall,
    F1 = ext_results$f1,
    Balanced_Accuracy = ext_results$balanced_accuracy,
    ROC_AUC = ext_results$auc,
    PR_AUC = ext_results$pr_auc,
    Brier_Score = ext_results$brier,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(output_prefix)) {
    save_single_roc_plot(
      roc_obj = ext_results$roc,
      auc_val = ext_results$auc,
      pr_auc = ext_results$pr_auc,
      sens = ext_results$sensitivity,
      spec = ext_results$specificity,
      file_name = file.path(figures_external_dir, paste0(output_prefix, "_ROC.tiff")),
      main_title = paste("ROC Curve -", basename(file_path))
    )
    
    save_single_pr_plot(
      pr_obj = ext_results$pr,
      file_name = file.path(figures_external_dir, paste0(output_prefix, "_PR.tiff")),
      main_title = paste("PR Curve -", basename(file_path))
    )
  }
  
  result_row
}

# ---------------------------
# 7) Read input data
# ---------------------------
TCGA_GTEx <- read_excel(input_file, sheet = "TCGA-GTEx") %>%
  as.data.frame()

LIRI_test <- read_excel(input_file, sheet = "LIRJ-test") %>%
  as.data.frame()

required_cols <- c("AB", "AC", "BC", "Group")
required_join_cols <- c("HLA-A", "HLA-B", "HLA-C", "AB", "AC", "BC", "Group")

check_required_columns(TCGA_GTEx, required_cols, "TCGA_GTEx")
check_required_columns(LIRI_test, required_cols, "LIRI_test")
check_required_columns(TCGA_GTEx, required_join_cols, "TCGA_GTEx")

# ---------------------------
# 8) Hold out 50 normal samples for independent validation
# ---------------------------
selected_normals <- TCGA_GTEx %>%
  filter(Group == "Normal") %>%
  sample_n(50) %>%
  as.data.frame()

remaining_data <- anti_join(
  TCGA_GTEx,
  selected_normals,
  by = c("HLA-A", "HLA-B", "HLA-C", "AB", "AC", "BC", "Group")
) %>%
  as.data.frame()

# ---------------------------
# 9) Prepare training data
# ---------------------------
x_train_all <- remaining_data[, c("AB", "AC", "BC"), drop = FALSE]
y_train <- as.factor(remaining_data$Group)

positive_class_train <- get_positive_class(y_train)

if (length(levels(y_train)) == 2 && levels(y_train)[2] != positive_class_train) {
  other_class <- setdiff(levels(y_train), positive_class_train)
  y_train <- factor(y_train, levels = c(other_class, positive_class_train))
}

# ---------------------------
# 10) Fit combined model
# ---------------------------
model_combined <- fit_logistic_model(x_train_all, y_train)
print(model_combined)

coef_vec <- coef(model_combined$finalModel)
coefficients_df <- data.frame(
  Term = names(coef_vec),
  Coefficient = as.numeric(coef_vec),
  row.names = NULL
)

performance_metrics_df <- as.data.frame(model_combined$results, stringsAsFactors = FALSE)

write.csv(coefficients_df, file.path(tables_model_fit_dir, "model_coefficients.csv"), row.names = FALSE)
write.csv(performance_metrics_df, file.path(tables_model_fit_dir, "model_cv_results.csv"), row.names = FALSE)

# ---------------------------
# 11) Training evaluation
# ---------------------------
train_results <- compute_binary_metrics(
  model = model_combined,
  x_data = x_train_all,
  y_data = y_train,
  positive_class = positive_class_train
)

print(paste("Training ROC AUC:", round(train_results$auc, 4)))
print(paste("Training PR AUC:", round(train_results$pr_auc, 4)))
print(paste("Training Brier Score:", round(train_results$brier, 4)))

save_single_roc_plot(
  roc_obj = train_results$roc,
  auc_val = train_results$auc,
  pr_auc = train_results$pr_auc,
  sens = train_results$sensitivity,
  spec = train_results$specificity,
  file_name = file.path(figures_roc_pr_dir, "ROC_Curve_Training.tiff"),
  main_title = "Training Performance: ROC Curve"
)

save_single_pr_plot(
  pr_obj = train_results$pr,
  file_name = file.path(figures_roc_pr_dir, "PR_Curve_Training.tiff"),
  main_title = "Training Performance: PR Curve"
)

# ---------------------------
# 12) Prepare independent test set
# ---------------------------
test_data_main <- LIRI_test[, c("AB", "AC", "BC", "Group"), drop = FALSE]
test_data_normals <- selected_normals[, c("AB", "AC", "BC", "Group"), drop = FALSE]
test_data_independent <- rbind(test_data_main, test_data_normals)

x_test_all <- test_data_independent[, c("AB", "AC", "BC"), drop = FALSE]
y_test <- as.factor(test_data_independent$Group)
y_test <- safe_factor_match(y_test, levels(y_train))
positive_class_test <- get_positive_class(y_test)

# ---------------------------
# 13) Independent test evaluation
# ---------------------------
test_results <- compute_binary_metrics(
  model = model_combined,
  x_data = x_test_all,
  y_data = y_test,
  positive_class = positive_class_test
)

print(test_results$confusion_matrix)
print(paste("Independent Test ROC AUC:", round(test_results$auc, 4)))
print(paste("Independent Test PR AUC:", round(test_results$pr_auc, 4)))
print(paste("Independent Test Brier Score:", round(test_results$brier, 4)))

save_single_roc_plot(
  roc_obj = test_results$roc,
  auc_val = test_results$auc,
  pr_auc = test_results$pr_auc,
  sens = test_results$sensitivity,
  spec = test_results$specificity,
  file_name = file.path(figures_roc_pr_dir, "ROC_Curve_Independent_Test.tiff"),
  main_title = "Independent Dataset Test Performance: ROC Curve"
)

save_single_pr_plot(
  pr_obj = test_results$pr,
  file_name = file.path(figures_roc_pr_dir, "PR_Curve_Independent_Test.tiff"),
  main_title = "Independent Dataset Test Performance: PR Curve"
)

# ---------------------------
# 14) Threshold selection on training set
# ---------------------------
best_threshold_train <- find_best_threshold_youden(
  probabilities = train_results$probabilities,
  labels = train_results$labels
)

training_optimal_threshold <- best_threshold_train$Threshold[1]

write.csv(
  best_threshold_train,
  file.path(tables_thresholds_dir, "best_threshold_train_youden.csv"),
  row.names = FALSE
)

print("Best threshold by Youden index - Training:")
print(best_threshold_train)

train_operating_point <- compute_metrics_at_fixed_threshold(
  probabilities = train_results$probabilities,
  labels = train_results$labels,
  threshold = training_optimal_threshold
)

test_operating_point <- compute_metrics_at_fixed_threshold(
  probabilities = test_results$probabilities,
  labels = test_results$labels,
  threshold = training_optimal_threshold
)

write.csv(
  train_operating_point,
  file.path(tables_thresholds_dir, "training_operating_point_at_youden_threshold.csv"),
  row.names = FALSE
)

write.csv(
  test_operating_point,
  file.path(tables_thresholds_dir, "independent_test_operating_point_at_training_threshold.csv"),
  row.names = FALSE
)

fixed_thresholds <- c(0.20, 0.30, 0.50, training_optimal_threshold)

train_fixed_threshold_metrics <- compute_metrics_multiple_fixed_thresholds(
  probabilities = train_results$probabilities,
  labels = train_results$labels,
  thresholds = fixed_thresholds
)

test_fixed_threshold_metrics <- compute_metrics_multiple_fixed_thresholds(
  probabilities = test_results$probabilities,
  labels = test_results$labels,
  thresholds = fixed_thresholds
)

write.csv(
  train_fixed_threshold_metrics,
  file.path(tables_thresholds_dir, "training_fixed_threshold_metrics.csv"),
  row.names = FALSE
)

write.csv(
  test_fixed_threshold_metrics,
  file.path(tables_thresholds_dir, "independent_test_fixed_threshold_metrics.csv"),
  row.names = FALSE
)

supp_threshold_table <- bind_rows(
  mutate(train_fixed_threshold_metrics, Dataset = "Training"),
  mutate(test_fixed_threshold_metrics, Dataset = "Independent")
) %>%
  select(
    Dataset, Threshold, Sensitivity, Specificity, PPV, NPV, Accuracy, F1, TP, TN, FP, FN
  )

write.csv(
  supp_threshold_table,
  file.path(tables_thresholds_dir, "Supplementary_Table_threshold_operating_points.csv"),
  row.names = FALSE
)

print(supp_threshold_table)

# ---------------------------
# 15) Calibration
# ---------------------------
calibration_summary_train <- build_calibration_summary(
  probabilities = train_results$probabilities,
  labels = train_results$labels,
  n_bins = 10
)

calibration_summary_test <- build_calibration_summary(
  probabilities = test_results$probabilities,
  labels = test_results$labels,
  n_bins = 10
)

write.csv(
  calibration_summary_train,
  file.path(tables_calibration_dir, "calibration_summary_training.csv"),
  row.names = FALSE
)

write.csv(
  calibration_summary_test,
  file.path(tables_calibration_dir, "calibration_summary_independent_test.csv"),
  row.names = FALSE
)

save_calibration_plot(
  calibration_summary = calibration_summary_train,
  file_prefix = file.path(figures_calibration_dir, "Calibration_Curve_Training"),
  plot_title = "Calibration Curve - Training"
)

save_calibration_plot(
  calibration_summary = calibration_summary_test,
  file_prefix = file.path(figures_calibration_dir, "Calibration_Curve_Independent_Test"),
  plot_title = "Calibration Curve - Independent Test"
)

# ---------------------------
# 16) Decision curve analysis
# ---------------------------
dca_train <- build_decision_curve_data(
  probabilities = train_results$probabilities,
  labels = train_results$labels,
  thresholds = seq(0.01, 0.50, by = 0.01),
  n_boot = 200,
  seed = 42
)

dca_test <- build_decision_curve_data(
  probabilities = test_results$probabilities,
  labels = test_results$labels,
  thresholds = seq(0.01, 0.50, by = 0.01),
  n_boot = 200,
  seed = 42
)

write.csv(dca_train, file.path(tables_decision_curve_dir, "decision_curve_training.csv"), row.names = FALSE)
write.csv(dca_test, file.path(tables_decision_curve_dir, "decision_curve_independent_test.csv"), row.names = FALSE)

dca_plot_train <- save_decision_curve_plot(
  dca_df = dca_train,
  file_prefix = file.path(figures_decision_curve_dir, "Decision_Curve_Training"),
  plot_title = "Decision Curve Analysis - Training"
)

dca_plot_test <- save_decision_curve_plot(
  dca_df = dca_test,
  file_prefix = file.path(figures_decision_curve_dir, "Decision_Curve_Independent_Test"),
  plot_title = "Decision Curve Analysis - Independent Test"
)

dca_train_summary <- summarize_net_benefit(dca_train)
dca_test_summary <- summarize_net_benefit(dca_test)

write.csv(
  dca_train_summary,
  file.path(tables_decision_curve_dir, "decision_curve_training_summary.csv"),
  row.names = FALSE
)

write.csv(
  dca_test_summary,
  file.path(tables_decision_curve_dir, "decision_curve_independent_test_summary.csv"),
  row.names = FALSE
)

dca_combined_plot <- dca_plot_train + dca_plot_test + plot_layout(ncol = 2)

ggsave(
  filename = file.path(figures_decision_curve_dir, "Decision_Curve_Combined_AB.png"),
  plot = dca_combined_plot,
  width = 12,
  height = 5.5,
  dpi = 600
)

ggsave(
  filename = file.path(figures_decision_curve_dir, "Decision_Curve_Combined_AB.tiff"),
  plot = dca_combined_plot,
  width = 12,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ---------------------------
# 17) Save independent test metrics
# ---------------------------
independent_metrics <- data.frame(
  Accuracy = test_results$accuracy,
  Kappa = test_results$kappa,
  Sensitivity = test_results$sensitivity,
  Specificity = test_results$specificity,
  Precision = test_results$precision,
  Recall = test_results$recall,
  F1 = test_results$f1,
  Balanced_Accuracy = test_results$balanced_accuracy,
  ROC_AUC = test_results$auc,
  PR_AUC = test_results$pr_auc,
  Brier_Score = test_results$brier
)

write.csv(
  independent_metrics,
  file.path(tables_model_fit_dir, "independent_test_metrics.csv"),
  row.names = FALSE
)

# ---------------------------
# 18) Fit single-ratio models for Figure 5
# ---------------------------
x_train_AB <- remaining_data[, "AB", drop = FALSE]
x_train_AC <- remaining_data[, "AC", drop = FALSE]
x_train_BC <- remaining_data[, "BC", drop = FALSE]

x_test_AB <- test_data_independent[, "AB", drop = FALSE]
x_test_AC <- test_data_independent[, "AC", drop = FALSE]
x_test_BC <- test_data_independent[, "BC", drop = FALSE]

model_AB <- fit_logistic_model(x_train_AB, y_train)
model_AC <- fit_logistic_model(x_train_AC, y_train)
model_BC <- fit_logistic_model(x_train_BC, y_train)

res_AB_train <- compute_binary_metrics(model_AB, x_train_AB, y_train, positive_class_train)
res_AB_test  <- compute_binary_metrics(model_AB, x_test_AB, y_test, positive_class_test)

res_AC_train <- compute_binary_metrics(model_AC, x_train_AC, y_train, positive_class_train)
res_AC_test  <- compute_binary_metrics(model_AC, x_test_AC, y_test, positive_class_test)

res_BC_train <- compute_binary_metrics(model_BC, x_train_BC, y_train, positive_class_train)
res_BC_test  <- compute_binary_metrics(model_BC, x_test_BC, y_test, positive_class_test)

res_ALL_train <- compute_binary_metrics(model_combined, x_train_all, y_train, positive_class_train)
res_ALL_test  <- compute_binary_metrics(model_combined, x_test_all, y_test, positive_class_test)

roc_panel_results <- list(
  "A" = res_AB_train,
  "B" = res_AB_test,
  "C" = res_AC_train,
  "D" = res_AC_test,
  "E" = res_BC_train,
  "F" = res_BC_test,
  "G" = res_ALL_train,
  "H" = res_ALL_test
)

save_multi_panel_roc(
  results_list = roc_panel_results,
  output_file = file.path(figures_figure5_dir, "Figure5_ROC_with_metrics.tiff")
)

figure5_summary <- data.frame(
  Panel = c("A", "B", "C", "D", "E", "F", "G", "H"),
  Model = c(
    "AB training", "AB independent",
    "AC training", "AC independent",
    "BC training", "BC independent",
    "Combined training", "Combined independent"
  ),
  ROC_AUC = c(
    res_AB_train$auc, res_AB_test$auc,
    res_AC_train$auc, res_AC_test$auc,
    res_BC_train$auc, res_BC_test$auc,
    res_ALL_train$auc, res_ALL_test$auc
  ),
  PR_AUC = c(
    res_AB_train$pr_auc, res_AB_test$pr_auc,
    res_AC_train$pr_auc, res_AC_test$pr_auc,
    res_BC_train$pr_auc, res_BC_test$pr_auc,
    res_ALL_train$pr_auc, res_ALL_test$pr_auc
  ),
  Sensitivity = c(
    res_AB_train$sensitivity, res_AB_test$sensitivity,
    res_AC_train$sensitivity, res_AC_test$sensitivity,
    res_BC_train$sensitivity, res_BC_test$sensitivity,
    res_ALL_train$sensitivity, res_ALL_test$sensitivity
  ),
  Specificity = c(
    res_AB_train$specificity, res_AB_test$specificity,
    res_AC_train$specificity, res_AC_test$specificity,
    res_BC_train$specificity, res_BC_test$specificity,
    res_ALL_train$specificity, res_ALL_test$specificity
  ),
  Brier_Score = c(
    res_AB_train$brier, res_AB_test$brier,
    res_AC_train$brier, res_AC_test$brier,
    res_BC_train$brier, res_BC_test$brier,
    res_ALL_train$brier, res_ALL_test$brier
  )
)

write.csv(
  figure5_summary,
  file.path(tables_figure5_dir, "Figure5_metrics_summary.csv"),
  row.names = FALSE
)

print(figure5_summary)

# ---------------------------
# 19) External dataset evaluation
# ---------------------------
external_results <- lapply(seq_along(external_files), function(i) {
  evaluate_external_dataset(
    file_path = external_files[i],
    model = model_combined,
    selected_normals_df = selected_normals,
    train_levels = levels(y_train),
    figures_external_dir = figures_external_dir,
    output_prefix = paste0("External_", sprintf("%02d", i))
  )
})

external_results <- external_results[!vapply(external_results, is.null, logical(1))]
external_results_df <- bind_rows(external_results)

write.csv(
  external_results_df,
  file.path(tables_external_dir, "external_dataset_evaluation_summary.csv"),
  row.names = FALSE
)

print(external_results_df)

# ---------------------------
# 20) Combined calibration + Brier figure
# ---------------------------
calibration_plot_train <- ggplot(
  calibration_summary_train,
  aes(x = Mean_Predicted, y = Observed_Fraction)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "",
    x = "Mean Predicted Probability",
    y = "Observed Fraction"
  ) +
  theme_classic(base_size = 12)

calibration_plot_test <- ggplot(
  calibration_summary_test,
  aes(x = Mean_Predicted, y = Observed_Fraction)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "",
    x = "Mean Predicted Probability",
    y = "Observed Fraction"
  ) +
  theme_classic(base_size = 12)

brier_df <- data.frame(
  Dataset = c("Training", "Independent Test"),
  Brier = c(train_results$brier, test_results$brier)
)

brier_plot <- ggplot(brier_df, aes(x = Dataset, y = Brier)) +
  geom_col(width = 0.6, fill = "#4C72B0") +
  geom_text(aes(label = sprintf("%.3f", Brier)), vjust = -0.5, size = 4) +
  ylim(0, max(brier_df$Brier) + 0.05) +
  labs(
    title = "",
    x = "",
    y = "Brier Score"
  ) +
  theme_classic(base_size = 12)

final_plot <- calibration_plot_train + calibration_plot_test + brier_plot +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(figures_calibration_dir, "Calibration_Brier_Combined.png"),
  final_plot,
  width = 12,
  height = 4.5,
  dpi = 600
)

ggsave(
  file.path(figures_calibration_dir, "Calibration_Brier_Combined.tiff"),
  final_plot,
  width = 12,
  height = 4.5,
  dpi = 600,
  compression = "lzw"
)

# ---------------------------
# 21) Save session info
# ---------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(logs_root, "sessionInfo_model.txt")
)

message("Workflow completed successfully.")