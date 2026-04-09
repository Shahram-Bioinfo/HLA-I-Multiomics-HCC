# ==========================================================
# Script: 10_scRNAseq_hla_workflow.R
# Purpose:
#   Perform a comprehensive single-cell workflow for HLA-I
#   analysis in HCC, including preprocessing, integration,
#   expression profiling, stage analysis, co-expression
#   analysis, and marker extraction.
#
# Input:
#   data/scrna/
#     - Ma2021.download.Rds
#     - Guilliams2022.download.Rds
#     - Stage_metadata.xlsx   [optional]
#
# Output:
#   results/scrna_hla_workflow/tables/
#   results/scrna_hla_workflow/figures/
#   results/scrna_hla_workflow/objects/
#   results/scrna_hla_workflow/logs/
# ==========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(harmony)
  library(readxl)
  library(writexl)
  library(openxlsx)
  library(cowplot)
  library(tibble)
  library(purrr)
  library(stringr)
  library(forcats)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."

data_dir <- file.path(project_root, "data", "scRNAseq")

tumor_rds_path  <- file.path(data_dir, "Ma2021.download.Rds")
normal_rds_path <- file.path(data_dir, "Guilliams2022.download.Rds")
stage_xlsx_path <- file.path(data_dir, "Stage_metadata.xlsx")

results_root <- file.path(project_root, "results", "scrna_hla_workflow")
table_dir <- file.path(results_root, "tables")
fig_dir <- file.path(results_root, "figures")
object_dir <- file.path(results_root, "objects")
log_dir <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 2) Parameters
# ---------------------------
min_features <- 200
max_features <- 2500
max_percent_mt <- 5

norm_method <- "LogNormalize"
scale_factor <- 10000
n_hvg <- 2000

global_pcs <- 20
subset_pcs <- 10

hla_genes <- c("HLA-A", "HLA-B", "HLA-C")

min_nonmalignant_fraction <- 0.05
min_cells_per_patient_celltype <- 10
min_patients_per_celltype <- 3

# ---------------------------
# 3) Helper functions
# ---------------------------
pick_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Could not find any of these columns: ", paste(candidates, collapse = ", "))
  }
  NULL
}

ensure_percent_mt <- function(seu) {
  if (!"percent.mt" %in% colnames(seu@meta.data)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  }
  seu
}

run_basic_qc_preprocess <- function(seu, dataset_label) {
  seu <- ensure_percent_mt(seu)
  
  qc_before <- seu@meta.data %>%
    transmute(
      cell = rownames(seu@meta.data),
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      percent.mt = percent.mt
    )
  
  write.csv(
    qc_before,
    file.path(table_dir, paste0(dataset_label, "_QC_before_filtering.csv")),
    row.names = FALSE
  )
  
  seu <- subset(
    seu,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
  
  qc_after <- seu@meta.data %>%
    transmute(
      cell = rownames(seu@meta.data),
      nFeature_RNA = nFeature_RNA,
      nCount_RNA = nCount_RNA,
      percent.mt = percent.mt
    )
  
  write.csv(
    qc_after,
    file.path(table_dir, paste0(dataset_label, "_QC_after_filtering.csv")),
    row.names = FALSE
  )
  
  seu <- NormalizeData(
    seu,
    normalization.method = norm_method,
    scale.factor = scale_factor,
    verbose = FALSE
  )
  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = n_hvg,
    verbose = FALSE
  )
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
  
  return(seu)
}

standardize_celltypes <- function(x) {
  x0 <- as.character(x)
  x0 <- ifelse(grepl("^b( |$)|\\bb[-_ ]?cells?\\b", x0, ignore.case = TRUE), "B cells", x0)
  x0 <- ifelse(grepl("plasma.*b|^plasma\\b", x0, ignore.case = TRUE), "Plasma B cells", x0)
  x0 <- ifelse(grepl("cd8", x0, ignore.case = TRUE) & grepl("t", x0, ignore.case = TRUE), "CD8+ T cells", x0)
  x0 <- ifelse(grepl("cd4", x0, ignore.case = TRUE) & grepl("t", x0, ignore.case = TRUE), "CD4+ T cells", x0)
  x0 <- ifelse(grepl("treg", x0, ignore.case = TRUE), "Tregs", x0)
  x0 <- ifelse(grepl("\\bc?dc1\\b", x0, ignore.case = TRUE), "cDC1s", x0)
  x0 <- ifelse(grepl("\\bc?dc2\\b", x0, ignore.case = TRUE), "cDC2s", x0)
  x0 <- ifelse(grepl("\\bpdc\\b", x0, ignore.case = TRUE), "pDCs", x0)
  x0 <- ifelse(grepl("circulating.*nk|nk/nkt", x0, ignore.case = TRUE), "Circulating NK/NKT", x0)
  x0 <- ifelse(grepl("resident.*nk", x0, ignore.case = TRUE), "Resident NK cell", x0)
  x0 <- ifelse(grepl("endothel|tec", x0, ignore.case = TRUE), "Endothelial cells", x0)
  x0 <- ifelse(grepl("macroph", x0, ignore.case = TRUE), "Macrophages", x0)
  x0 <- ifelse(grepl("mono", x0, ignore.case = TRUE), "Monocytes", x0)
  x0 <- ifelse(grepl("malign|tumou?r", x0, ignore.case = TRUE), "Malignant cells", x0)
  x0 <- ifelse(grepl("strom", x0, ignore.case = TRUE), "Stromal cells", x0)
  x0 <- ifelse(grepl("cholangi", x0, ignore.case = TRUE), "Cholangiocytes", x0)
  x0 <- ifelse(grepl("hepatocyte", x0, ignore.case = TRUE), "Hepatocytes", x0)
  x0
}

safe_sheet_names <- function(x) {
  s <- substr(x, 1, 31)
  dup <- names(table(s))[table(s) > 1]
  for (d in dup) {
    idx <- which(s == d)
    base <- substr(d, 1, 28)
    s[idx] <- paste0(base, sprintf("_%02d", seq_along(idx)))
  }
  s
}

fisher_z <- function(r) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  0.5 * log((1 + r) / (1 - r))
}

fast_p_from_r <- function(r, n) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tval <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  2 * pt(abs(tval), df = n - 2, lower.tail = FALSE)
}

compare_z_between_groups <- function(z1, n1, z2, n2) {
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  (z1 - z2) / se
}

# ---------------------------
# 4) Load raw data
# ---------------------------
tumor_raw  <- readRDS(tumor_rds_path)
normal_raw <- readRDS(normal_rds_path)

tumor_raw$dataset_label  <- "Tumor"
normal_raw$dataset_label <- "Normal"

tumor_ct_col  <- pick_col(tumor_raw@meta.data,  c("Celltype", "CellType", "celltype"), required = FALSE)
normal_ct_col <- pick_col(normal_raw@meta.data, c("Celltype", "CellType", "celltype"), required = FALSE)

if (!is.null(tumor_ct_col)) {
  tumor_raw$Celltype <- standardize_celltypes(tumor_raw@meta.data[[tumor_ct_col]])
}
if (!is.null(normal_ct_col)) {
  normal_raw$Celltype <- standardize_celltypes(normal_raw@meta.data[[normal_ct_col]])
}

# ---------------------------
# 5) QC and preprocess each dataset
# ---------------------------
tumor_pre  <- run_basic_qc_preprocess(tumor_raw,  "Tumor")
normal_pre <- run_basic_qc_preprocess(normal_raw, "Normal")

saveRDS(tumor_pre,  file.path(object_dir, "Tumor_preprocessed.rds"))
saveRDS(normal_pre, file.path(object_dir, "Normal_preprocessed.rds"))

# ---------------------------
# 6) Merge and Harmony batch correction
# ---------------------------
combined_data <- merge(
  tumor_pre,
  y = normal_pre,
  add.cell.ids = c("Batch1", "Batch2")
)

if (!"orig.ident" %in% colnames(combined_data@meta.data)) {
  combined_data$orig.ident <- ifelse(grepl("^Batch1_", colnames(combined_data)), "Tumor", "Normal")
}

combined_data <- RunUMAP(
  combined_data,
  dims = 1:global_pcs,
  reduction = "pca",
  reduction.name = "umap_before",
  verbose = FALSE
)

p_before <- DimPlot(
  combined_data,
  reduction = "umap_before",
  group.by = "orig.ident"
) + ggtitle("Before Batch Correction")

ggsave(
  file.path(fig_dir, "UMAP_before_batch_correction.png"),
  p_before, width = 7, height = 6, dpi = 300
)

combined_data <- RunHarmony(
  combined_data,
  group.by.vars = "orig.ident",
  verbose = FALSE
)

combined_data <- RunUMAP(
  combined_data,
  reduction = "harmony",
  dims = 1:global_pcs,
  reduction.name = "umap_after",
  verbose = FALSE
)

p_after <- DimPlot(
  combined_data,
  reduction = "umap_after",
  group.by = "orig.ident"
) + ggtitle("After Batch Correction")

ggsave(
  file.path(fig_dir, "UMAP_after_batch_correction.png"),
  p_after, width = 7, height = 6, dpi = 300
)

ggsave(
  file.path(fig_dir, "batch_correction_plots.png"),
  p_before + p_after, width = 12, height = 6, dpi = 300
)

saveRDS(combined_data, file.path(object_dir, "corrected_combined_data.rds"))

# ---------------------------
# 7) Split corrected object into tumor and normal
# ---------------------------
cell_names <- colnames(combined_data)
tumor_cells  <- grep("^Batch1_", cell_names, value = TRUE)
normal_cells <- grep("^Batch2_", cell_names, value = TRUE)

tumor_obj  <- subset(combined_data, cells = tumor_cells)
normal_obj <- subset(combined_data, cells = normal_cells)

tumor_obj  <- FindNeighbors(tumor_obj, dims = 1:subset_pcs, reduction = "harmony", verbose = FALSE)
tumor_obj  <- FindClusters(tumor_obj, resolution = 0.5, verbose = FALSE)
tumor_obj  <- RunUMAP(tumor_obj, dims = 1:subset_pcs, reduction = "harmony", verbose = FALSE)
tumor_obj  <- RunTSNE(tumor_obj, dims = 1:subset_pcs, reduction = "harmony", check_duplicates = FALSE, verbose = FALSE)

normal_obj <- FindNeighbors(normal_obj, dims = 1:subset_pcs, reduction = "harmony", verbose = FALSE)
normal_obj <- FindClusters(normal_obj, resolution = 0.5, verbose = FALSE)
normal_obj <- RunUMAP(normal_obj, dims = 1:subset_pcs, reduction = "harmony", verbose = FALSE)
normal_obj <- RunTSNE(normal_obj, dims = 1:subset_pcs, reduction = "harmony", check_duplicates = FALSE, verbose = FALSE)

saveRDS(tumor_obj,  file.path(object_dir, "Tumor_corrected_subset.rds"))
saveRDS(normal_obj, file.path(object_dir, "Normal_corrected_subset.rds"))

# ---------------------------
# 8) Extract HLA expression tables
# ---------------------------
extract_hla_table <- function(seu, group_label) {
  stopifnot(all(hla_genes %in% rownames(seu[["RNA"]]@data)))
  
  df <- data.frame(
    Cell = colnames(seu),
    Celltype = as.character(seu$Celltype),
    HLA_A_Expression = as.numeric(GetAssayData(seu, assay = "RNA", slot = "data")["HLA-A", ]),
    HLA_B_Expression = as.numeric(GetAssayData(seu, assay = "RNA", slot = "data")["HLA-B", ]),
    HLA_C_Expression = as.numeric(GetAssayData(seu, assay = "RNA", slot = "data")["HLA-C", ]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Group = group_label,
      Total = HLA_A_Expression + HLA_B_Expression + HLA_C_Expression
    )
  
  return(df)
}

hla_tumor_df  <- extract_hla_table(tumor_obj,  "HCC")
hla_normal_df <- extract_hla_table(normal_obj, "Normal")

write_xlsx(
  list(Exp = hla_tumor_df),
  file.path(table_dir, "HLA-Tumor.xlsx")
)

write_xlsx(
  list(Exp = hla_normal_df),
  file.path(table_dir, "HLA-Normal.xlsx")
)

write.csv(hla_tumor_df,  file.path(table_dir, "HLA-T.csv"), row.names = FALSE)
write.csv(hla_normal_df, file.path(table_dir, "HLA-N.csv"), row.names = FALSE)

combined_hla_df <- bind_rows(hla_normal_df, hla_tumor_df)

# ---------------------------
# 9) Global transcriptional separation analysis
# ---------------------------
fig8a_input <- combined_hla_df %>%
  mutate(Celltype = standardize_celltypes(Celltype)) %>%
  mutate(
    CellClass = case_when(
      Celltype == "Malignant cells" ~ "Malignant cells",
      Celltype %in% c("Hepatocytes", "Cholangiocytes") ~ Celltype,
      Celltype %in% c("Stromal cells", "Endothelial cells") ~ "Other stromal cells",
      TRUE ~ "Immune cells"
    )
  )

pca_mat <- fig8a_input %>%
  select(HLA_A_Expression, HLA_B_Expression, HLA_C_Expression) %>%
  as.matrix()

pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- fig8a_input %>%
  mutate(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )

p_fig8a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Celltype)) +
  geom_point(size = 0.6, alpha = 0.7) +
  stat_ellipse(aes(group = Celltype), linewidth = 0.4, alpha = 0.6) +
  theme_bw() +
  labs(
    title = "Global transcriptional separation",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")
  )

ggsave(
  file.path(fig_dir, "PCA_HLA_expression.png"),
  p_fig8a, width = 10, height = 7, dpi = 300
)

write.csv(pca_df, file.path(table_dir, "PCA_coordinates.csv"), row.names = FALSE)

# ---------------------------
# 10) Zero-expression frequency analysis
# ---------------------------
tumor_comp <- hla_tumor_df %>%
  count(Celltype, name = "n_cells") %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  arrange(desc(freq))

write.csv(tumor_comp, file.path(table_dir, "Tumor_celltype_composition.csv"), row.names = FALSE)

top_nonmalignant <- tumor_comp %>%
  filter(Celltype != "Malignant cells") %>%
  filter(freq > min_nonmalignant_fraction) %>%
  pull(Celltype)

top_nonmalignant <- head(top_nonmalignant, 4)

fig8b_celltypes <- c("Malignant cells", top_nonmalignant)

fig8b_df <- hla_tumor_df %>%
  filter(Celltype %in% fig8b_celltypes) %>%
  mutate(
    zero_HLA_A = HLA_A_Expression == 0,
    zero_HLA_B = HLA_B_Expression == 0,
    zero_HLA_C = HLA_C_Expression == 0,
    zero_all_three = HLA_A_Expression == 0 & HLA_B_Expression == 0 & HLA_C_Expression == 0
  ) %>%
  group_by(Celltype) %>%
  summarise(
    HLA_A = mean(zero_HLA_A) * 100,
    HLA_B = mean(zero_HLA_B) * 100,
    HLA_C = mean(zero_HLA_C) * 100,
    All_Three = mean(zero_all_three) * 100,
    .groups = "drop"
  )

write.csv(fig8b_df, file.path(table_dir, "Zero_expression_frequencies.csv"), row.names = FALSE)

fig8b_long <- fig8b_df %>%
  pivot_longer(cols = c(HLA_A, HLA_B, HLA_C, All_Three),
               names_to = "Measure",
               values_to = "Percent")

p_fig8b <- ggplot(fig8b_long, aes(x = Celltype, y = Percent, fill = Measure)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(
    title = "Zero-expression frequencies",
    x = "",
    y = "Frequency of zero-expression cells (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fig_dir, "Zero_expression_barplot.png"),
  p_fig8b, width = 10, height = 6, dpi = 300
)

# ---------------------------
# 11) Overall comparison between normal and tumor
# ---------------------------
plot_overall_violin <- function(df, gene_col, ylab, label_letter) {
  p <- ggplot(df, aes(x = Group, y = .data[[gene_col]], color = Group, fill = Group)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    stat_compare_means(
      comparisons = list(c("Normal", "HCC")),
      method = "wilcox.test",
      label = "p.signif",
      size = 6
    ) +
    labs(x = NULL, y = ylab) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p_lab <- ggdraw(p) + draw_plot_label(label_letter, x = 0, y = 1,
                                       hjust = -0.1, vjust = 1.1,
                                       size = 14, fontface = "bold")
  p_lab
}

p8c <- plot_overall_violin(combined_hla_df, "HLA_A_Expression", "HLA-A Expression", "C")
p8d <- plot_overall_violin(combined_hla_df, "HLA_B_Expression", "HLA-B Expression", "D")
p8e <- plot_overall_violin(combined_hla_df, "HLA_C_Expression", "HLA-C Expression", "E")

ggsave(
  file.path(fig_dir, "Overall_HLA_expression_comparison.png"),
  plot_grid(p8c, p8d, p8e, ncol = 1),
  width = 6, height = 14, dpi = 300
)

# ---------------------------
# 12) Cell-type specific comparison between normal and tumor
# ---------------------------
plot_celltype_group_violin <- function(df, gene_col, ylab, label_letter, p_adjust_method = "BH") {
  d <- df %>%
    mutate(Celltype = standardize_celltypes(Celltype)) %>%
    group_by(Celltype, Group) %>%
    mutate(n_in_group = n()) %>%
    ungroup()
  
  keep_ct <- d %>%
    count(Celltype, Group) %>%
    group_by(Celltype) %>%
    filter(n() == 2) %>%
    pull(Celltype)
  
  d <- d %>% filter(Celltype %in% keep_ct)
  
  pv <- compare_means(
    formula = as.formula(paste(gene_col, "~ Group")),
    group.by = "Celltype",
    data = d,
    method = "wilcox.test",
    p.adjust.method = p_adjust_method
  )
  
  ypos <- d %>%
    group_by(Celltype) %>%
    summarise(
      ymax = max(.data[[gene_col]], na.rm = TRUE),
      yrng = diff(range(.data[[gene_col]], na.rm = TRUE)),
      y.position = ymax + ifelse(yrng > 0, 0.1 * yrng, 0.3),
      .groups = "drop"
    )
  
  pv <- left_join(pv, ypos, by = "Celltype") %>%
    mutate(xmin = Celltype, xmax = Celltype)
  
  p <- ggplot(d, aes(x = Celltype, y = .data[[gene_col]], fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    stat_pvalue_manual(pv, label = "p.adj.signif", hide.ns = TRUE, tip.length = 0, size = 4) +
    theme_bw() +
    labs(x = "Celltype", y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_lab <- ggdraw(p) + draw_plot_label(label_letter, x = 0, y = 1,
                                       hjust = -0.1, vjust = 1.1,
                                       size = 14, fontface = "bold")
  p_lab
}

p8f <- plot_celltype_group_violin(combined_hla_df, "HLA_A_Expression", "HLA-A Expression", "F")
p8g <- plot_celltype_group_violin(combined_hla_df, "HLA_B_Expression", "HLA-B Expression", "G")
p8h <- plot_celltype_group_violin(combined_hla_df, "HLA_C_Expression", "HLA-C Expression", "H")

ggsave(
  file.path(fig_dir, "Celltype_specific_HLA_expression_comparison.png"),
  plot_grid(p8f, p8g, p8h, ncol = 1),
  width = 14, height = 18, dpi = 300
)

# ---------------------------
# 13) Save supporting tables
# ---------------------------
write_xlsx(
  list(
    Overall_HLA = combined_hla_df,
    Zero_Frequency = fig8b_df,
    PCA_Coordinates = pca_df,
    Tumor_Composition = tumor_comp
  ),
  file.path(table_dir, "HLA_supporting_tables.xlsx")
)

# ---------------------------
# 14) Stage analysis
# ---------------------------
if (file.exists(stage_xlsx_path)) {
  stage_map <- read_excel(stage_xlsx_path)
  
  if (!all(c("Cell", "PatientID", "StageGroup") %in% colnames(stage_map))) {
    stop("Stage file must contain at least: Cell, PatientID, StageGroup")
  }
  
  stage_df <- hla_tumor_df %>%
    left_join(
      stage_map %>%
        transmute(
          Cell = as.character(Cell),
          PatientID = as.character(PatientID),
          StageGroup = as.character(StageGroup),
          Stage = ifelse(grepl("IV", StageGroup, ignore.case = TRUE), "Stage IV", "Stage I-IIIB")
        ),
      by = c("Cell" = "Cell")
    ) %>%
    filter(!is.na(Stage), !is.na(PatientID)) %>%
    mutate(
      Stage = factor(Stage, levels = c("Stage I-IIIB", "Stage IV")),
      Celltype = standardize_celltypes(Celltype)
    )
  
  write.csv(stage_df, file.path(table_dir, "Stage_joined_cell_level_table.csv"), row.names = FALSE)
  
  stage_pb <- stage_df %>%
    group_by(PatientID, Celltype, Stage) %>%
    summarise(
      HLA_A = mean(HLA_A_Expression, na.rm = TRUE),
      HLA_B = mean(HLA_B_Expression, na.rm = TRUE),
      HLA_C = mean(HLA_C_Expression, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells_per_patient_celltype)
  
  write.csv(stage_pb, file.path(table_dir, "Stage_patient_celltype_pseudobulk.csv"), row.names = FALSE)
  
  valid_ct <- stage_pb %>%
    group_by(Celltype, Stage) %>%
    summarise(n_patients = n_distinct(PatientID), .groups = "drop") %>%
    group_by(Celltype) %>%
    summarise(
      min_patients = min(n_patients),
      n_groups = n(),
      .groups = "drop"
    ) %>%
    filter(n_groups == 2, min_patients >= min_patients_per_celltype) %>%
    pull(Celltype)
  
  stage_pb_valid <- stage_pb %>% filter(Celltype %in% valid_ct)
  
  stage_overall <- stage_pb_valid %>%
    group_by(PatientID, Stage) %>%
    summarise(
      HLA_A = mean(HLA_A, na.rm = TRUE),
      HLA_B = mean(HLA_B, na.rm = TRUE),
      HLA_C = mean(HLA_C, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_stage_overall <- function(df, gene_col, ylab, letter) {
    p <- ggplot(df, aes(x = Stage, y = .data[[gene_col]], color = Stage, fill = Stage)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
      stat_compare_means(
        comparisons = list(c("Stage I-IIIB", "Stage IV")),
        method = "wilcox.test",
        label = "p.signif",
        size = 6
      ) +
      labs(x = NULL, y = ylab) +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggdraw(p) + draw_plot_label(letter, x = 0, y = 1,
                                hjust = -0.1, vjust = 1.1,
                                size = 14, fontface = "bold")
  }
  
  p9a <- plot_stage_overall(stage_overall, "HLA_A", "HLA-A expression", "A")
  p9b <- plot_stage_overall(stage_overall, "HLA_B", "HLA-B expression", "B")
  p9c <- plot_stage_overall(stage_overall, "HLA_C", "HLA-C expression", "C")
  
  ggsave(
    file.path(fig_dir, "Stage_overall_expression_comparison.png"),
    plot_grid(p9a, p9b, p9c, ncol = 1),
    width = 6, height = 14, dpi = 300
  )
  
  plot_stage_by_celltype <- function(df, gene_col, ylab, letter) {
    d <- df %>%
      transmute(
        PatientID = PatientID,
        Celltype = Celltype,
        Stage = Stage,
        Expression = .data[[gene_col]]
      )
    
    pv <- compare_means(
      Expression ~ Stage,
      group.by = "Celltype",
      data = d,
      method = "wilcox.test",
      p.adjust.method = "BH"
    )
    
    ypos <- d %>%
      group_by(Celltype) %>%
      summarise(
        ymax = max(Expression, na.rm = TRUE),
        yrng = diff(range(Expression, na.rm = TRUE)),
        y.position = ymax + ifelse(yrng > 0, 0.1 * yrng, 0.3),
        .groups = "drop"
      )
    
    pv <- left_join(pv, ypos, by = "Celltype") %>%
      mutate(xmin = Celltype, xmax = Celltype)
    
    p <- ggplot(d, aes(x = Celltype, y = Expression, fill = Stage)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
      stat_pvalue_manual(pv, label = "p.adj.signif", hide.ns = FALSE, tip.length = 0, size = 4) +
      theme_bw() +
      labs(x = "Celltype in HCC", y = ylab) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggdraw(p) + draw_plot_label(letter, x = 0, y = 1,
                                hjust = -0.1, vjust = 1.1,
                                size = 14, fontface = "bold")
  }
  
  p9d <- plot_stage_by_celltype(stage_pb_valid, "HLA_A", "HLA-A expression", "D")
  p9e <- plot_stage_by_celltype(stage_pb_valid, "HLA_B", "HLA-B expression", "E")
  p9f <- plot_stage_by_celltype(stage_pb_valid, "HLA_C", "HLA-C expression", "F")
  
  ggsave(
    file.path(fig_dir, "Stage_celltype_specific_expression_comparison.png"),
    plot_grid(p9d, p9e, p9f, ncol = 1),
    width = 14, height = 18, dpi = 300
  )
  
  ggsave(
    file.path(fig_dir, "Stage_full_expression_comparison.png"),
    plot_grid(
      plot_grid(p9a, p9b, p9c, ncol = 1),
      plot_grid(p9d, p9e, p9f, ncol = 1),
      ncol = 2, rel_widths = c(1, 2.2)
    ),
    width = 18, height = 18, dpi = 300
  )
  
  write_xlsx(
    list(
      Cell_level = stage_df,
      Patient_Celltype_Pseudobulk = stage_pb_valid,
      Patient_Overall = stage_overall
    ),
    file.path(table_dir, "Stage_supporting_tables.xlsx")
  )
}

# ---------------------------
# 15) Single-cell co-expression analysis
# ---------------------------
tumor_md <- tumor_obj@meta.data

patient_col <- pick_col(
  tumor_md,
  c("PatientID", "Patient_ID", "patient", "patient_id", "SampleID", "sample", "Donor", "DonorID"),
  required = FALSE
)

celltype_col <- pick_col(
  tumor_md,
  c("Celltype", "CellType", "celltype"),
  required = TRUE
)

if (!is.null(patient_col)) {
  tumor_meta <- tumor_md %>%
    transmute(
      Cell = rownames(tumor_md),
      PatientID = as.character(.data[[patient_col]]),
      Celltype = standardize_celltypes(.data[[celltype_col]])
    ) %>%
    filter(!is.na(PatientID), !is.na(Celltype))
  
  expr <- GetAssayData(tumor_obj, assay = "RNA", slot = "data")
  
  tumor_meta$group_id <- paste(tumor_meta$PatientID, tumor_meta$Celltype, sep = "|")
  grp_counts <- table(tumor_meta$group_id)
  keep_groups <- names(grp_counts)[grp_counts >= min_cells_per_patient_celltype]
  keep_cells <- tumor_meta$Cell[tumor_meta$group_id %in% keep_groups]
  
  expr_valid <- expr[, keep_cells, drop = FALSE]
  meta_valid <- tumor_meta %>% filter(Cell %in% keep_cells)
  
  group_ids <- meta_valid$group_id
  design <- model.matrix(~ 0 + factor(group_ids))
  colnames(design) <- sub("^factor\\(group_ids\\)", "", colnames(design))
  
  group_sizes <- colSums(design)
  design <- sweep(design, 2, group_sizes, FUN = "/")
  pseudo <- as.matrix(expr_valid) %*% design
  
  group_info <- tibble(
    Group = colnames(pseudo),
    PatientID = sub("\\|.*$", "", Group),
    Celltype = sub("^.*\\|", "", Group)
  )
  
  write.csv(group_info, file.path(table_dir, "Tumor_pseudobulk_group_info.csv"), row.names = FALSE)
  
  common_nonmalignant <- tumor_comp %>%
    filter(Celltype != "Malignant cells", freq > min_nonmalignant_fraction) %>%
    pull(Celltype) %>%
    head(4)
  
  target_celltypes <- c("Malignant cells", common_nonmalignant)
  
  for (anchor in hla_genes) {
    if (!(anchor %in% rownames(pseudo))) next
    
    wb <- createWorkbook()
    ct_results <- list()
    
    for (ct in target_celltypes) {
      cols <- which(group_info$Celltype == ct)
      if (length(cols) < min_patients_per_celltype) next
      
      anchor_vec <- as.numeric(pseudo[anchor, cols, drop = TRUE])
      if (sd(anchor_vec) == 0) next
      
      X <- t(pseudo[, cols, drop = FALSE])
      r <- suppressWarnings(cor(X, anchor_vec, method = "pearson"))
      p <- fast_p_from_r(r, n = length(anchor_vec))
      q <- p.adjust(p, method = "BH")
      
      res <- tibble(
        Gene = rownames(pseudo),
        Correlation = as.numeric(r),
        P_value = as.numeric(p),
        FDR = as.numeric(q),
        Z = fisher_z(as.numeric(r)),
        N = length(anchor_vec),
        Celltype = ct,
        Anchor = anchor
      ) %>%
        filter(Gene != anchor)
      
      ct_results[[ct]] <- res
      
      sname <- safe_sheet_names(ct)
      addWorksheet(wb, sname)
      writeData(wb, sname, res)
    }
    
    if (length(ct_results) > 0) {
      saveWorkbook(
        wb,
        file.path(table_dir, paste0(anchor, "_celltype_correlations.xlsx")),
        overwrite = TRUE
      )
      
      merged_z <- bind_rows(ct_results) %>%
        select(Gene, Celltype, Z, FDR, N) %>%
        mutate(Z_sig = ifelse(FDR < 0.05, Z, 0)) %>%
        select(Gene, Celltype, Z_sig) %>%
        pivot_wider(names_from = Celltype, values_from = Z_sig, values_fill = 0)
      
      merged_z <- merged_z %>%
        filter(if_any(-Gene, ~ abs(.) > 0.3))
      
      write.csv(
        merged_z,
        file.path(table_dir, paste0(anchor, "_merged_zscores.csv")),
        row.names = FALSE
      )
      
      if ("Malignant cells" %in% names(ct_results)) {
        mal_df <- ct_results[["Malignant cells"]] %>%
          select(Gene, Z_malignant = Z, N_malignant = N, FDR_malignant = FDR)
        
        comp_list <- list()
        
        for (ct in common_nonmalignant) {
          if (!ct %in% names(ct_results)) next
          other_df <- ct_results[[ct]] %>%
            select(Gene, Z_other = Z, N_other = N, FDR_other = FDR)
          
          comp <- mal_df %>%
            inner_join(other_df, by = "Gene") %>%
            mutate(
              Z_diff_stat = compare_z_between_groups(Z_malignant, N_malignant, Z_other, N_other),
              P_diff = 2 * pnorm(abs(Z_diff_stat), lower.tail = FALSE),
              Comparison = paste("Malignant_vs", ct, sep = "_")
            )
          
          comp_list[[ct]] <- comp
        }
        
        if (length(comp_list) > 0) {
          comp_all <- bind_rows(comp_list) %>%
            group_by(Comparison) %>%
            mutate(FDR_diff = p.adjust(P_diff, method = "BH")) %>%
            ungroup()
          
          write.csv(
            comp_all,
            file.path(table_dir, paste0(anchor, "_malignant_vs_nonmalignant_comparisons.csv")),
            row.names = FALSE
          )
        }
      }
    }
  }
}

# ---------------------------
# 16) Optional marker analysis
# ---------------------------
tumor_markers <- FindAllMarkers(tumor_obj, min.pct = 0.25, verbose = FALSE)
normal_markers <- FindAllMarkers(normal_obj, min.pct = 0.25, verbose = FALSE)

write_xlsx(
  list(
    Tumor_All = tumor_markers,
    Tumor_HLA = tumor_markers %>% filter(gene %in% hla_genes),
    Normal_All = normal_markers,
    Normal_HLA = normal_markers %>% filter(gene %in% hla_genes)
  ),
  file.path(table_dir, "Markers_HLA_and_all.xlsx")
)

# ---------------------------
# 17) Log and session info
# ---------------------------
capture.output(
  sessionInfo(),
  file = file.path(log_dir, "sessionInfo.txt")
)

writeLines(
  c(
    "12_scRNAseq_hla_workflow.R completed",
    paste("Data directory:", data_dir),
    paste("Results directory:", results_root),
    paste("Tables directory:", table_dir),
    paste("Figures directory:", fig_dir),
    paste("Objects directory:", object_dir),
    paste("Logs directory:", log_dir)
  ),
  con = file.path(log_dir, "12_scRNAseq_hla_workflow_log.txt")
)

message("Single-cell workflow completed successfully.")