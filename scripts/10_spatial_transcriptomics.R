# ==========================================================
# Script: 10_spatial_transcriptomics.R
# Purpose:
#   Perform HLA-I spatial transcriptomics analysis for the
#   HRA000437 cohort, generate pooled and balanced violin
#   plots, spatial feature maps, QC summaries, and export
#   supplementary tables.
#
# Input:
#   data/spatial_transcriptomics/HRA000437/
#     - HCC-*-expr.RDS
#
# Output:
#   results/allele_specific_neoantigen/tables/
#     - HLA_Analysis_Supplement.xlsx
#
#   results/allele_specific_neoantigen/figures/
#     - Main_Figure_Violin_pooled.pdf
#     - Main_Figure_Violin_pooled.png
#     - Supp_Figure_QC_SpotCounts.pdf
#     - Supp_Figure_QC_SpotCounts.png
#     - Supp_Figure_Violin_balanced.pdf
#     - Supp_Figure_Violin_balanced.png
#     - Supp_Figure_FeatureMaps.pdf
#     - Supp_Figure_FeatureMaps.png
#
#   results/allele_specific_neoantigen/logs/
#     - 04_spatial_transcriptomics_log.txt
# ==========================================================

# ---------------------------
# 0) Packages
# ---------------------------
required_pkgs <- c(
  "Seurat", "ggplot2", "dplyr", "tidyr", "openxlsx",
  "stringr", "forcats", "patchwork", "ggpubr"
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
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(ggpubr)
})

# ---------------------------
# 1) Paths
# ---------------------------
project_root <- "."

data_dir <- file.path(
  project_root, "data", "spatial_transcriptomics", "HRA000437"
)

results_root <- file.path(project_root, "results", "allele_specific_neoantigen")
tab_dir <- file.path(results_root, "tables")
fig_dir <- file.path(results_root, "figures")
log_dir <- file.path(results_root, "logs")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(data_dir)) {
  stop("Spatial transcriptomics input directory not found: ", data_dir)
}

rds_files <- list.files(
  path = data_dir,
  pattern = "^HCC-\\d+-expr\\.RDS$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No RDS files found in: ", data_dir)
}

rds_files <- rds_files[order(as.integer(sub(".*HCC-(\\d+).*", "\\1", rds_files)))]

# ---------------------------
# 2) Analysis config
# ---------------------------
genes_hla        <- c("HLA-A", "HLA-B", "HLA-C")
region_lvls_all  <- c("Normal", "Immune", "Stromal", "Tumor")
regions_plot     <- c("Normal", "Stromal", "Tumor")
comparators      <- c("Normal", "Stromal")
do_balanced_pool <- TRUE
k_per_sample_reg <- 200
set.seed(7)

# ---------------------------
# 3) Helpers
# ---------------------------
safe_update <- function(obj) {
  tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
}

safe_features <- function(obj) {
  tryCatch(rownames(obj[["RNA"]]), error = function(e) character(0))
}

# Prefer normalized 'data'; fallback to log1p(counts)
get_expr_vec <- function(obj, gene) {
  feats <- safe_features(obj)
  if (!(gene %in% feats)) {
    return(NULL)
  }
  
  v <- tryCatch(
    GetAssayData(obj, assay = "RNA", slot = "data")[gene, ],
    error = function(e) NULL
  )
  if (!is.null(v)) {
    return(as.numeric(v))
  }
  
  v <- tryCatch(
    GetAssayData(obj, assay = "RNA", slot = "counts")[gene, ],
    error = function(e) NULL
  )
  if (!is.null(v)) {
    return(as.numeric(log1p(v)))
  }
  
  NULL
}

# Harmonize region labels -> Normal/Immune/Stromal/Tumor
get_region_vector <- function(meta) {
  r <- if ("Type" %in% names(meta)) as.character(meta$Type) else NA
  r <- r |>
    stringr::str_trim() |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_to_title()
  
  recode_map <- c(
    "Tumour" = "Tumor",
    "Stroma" = "Stromal",
    "Immune Region" = "Immune",
    "Normal Tissue" = "Normal",
    "Tumoral" = "Tumor"
  )
  
  r <- dplyr::recode(r, !!!recode_map, .default = r)
  r[!r %in% region_lvls_all] <- NA_character_
  
  factor(r, levels = region_lvls_all)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Helper for spatial feature maps (requires imagecol/imagerow in meta)
make_spatial_plot <- function(meta, expr, title = "", pt = 0.35) {
  stopifnot(all(c("imagecol", "imagerow") %in% colnames(meta)))
  
  df <- cbind(meta[, c("imagecol", "imagerow")], expr = expr)
  colnames(df)[1:2] <- c("x", "y")
  
  ggplot(df, aes(x = x, y = y, color = expr)) +
    geom_point(size = pt) +
    scale_color_viridis_c(option = "magma") +
    coord_fixed() +
    scale_y_reverse() +
    theme_void(base_size = 9) +
    labs(title = title, color = "Expr")
}

# ---------------------------
# 4) Build spot-level dataframe (pooled) + collect feature maps
# ---------------------------
df_all <- list()
feature_panels <- list()

for (f in rds_files) {
  obj <- readRDS(f) |> safe_update()
  sid <- sub("-expr$", "", tools::file_path_sans_ext(basename(f)))
  meta <- obj@meta.data
  meta$region <- get_region_vector(meta)
  feats <- safe_features(obj)
  
  # Feature maps per sample × gene
  for (g in genes_hla) {
    expr <- if (g %in% feats) get_expr_vec(obj, g) else NULL
    
    p <- if (is.null(expr)) {
      ggplot() + theme_void() + labs(title = paste0(sid, " – ", g, " (missing)"))
    } else {
      make_spatial_plot(meta, expr, title = paste0(sid, " – ", g), pt = 0.35)
    }
    
    feature_panels[[length(feature_panels) + 1]] <- p
  }
  
  # Assemble pooled df_all
  for (g in genes_hla) {
    if (!(g %in% feats)) next
    
    expr <- get_expr_vec(obj, g)
    if (is.null(expr)) next
    
    df_all[[length(df_all) + 1]] <- data.frame(
      sample = sid,
      gene = g,
      region = meta$region,
      expr = expr
    )
  }
}

df_all <- dplyr::bind_rows(df_all) |>
  dplyr::filter(!is.na(region)) |>
  dplyr::mutate(region = factor(as.character(region), levels = region_lvls_all))

# ---------------------------
# 5) QC counts (for S10 + QC figure)
# ---------------------------
spot_counts <- df_all |>
  dplyr::count(sample, region, name = "n_spots")

spot_totals <- spot_counts |>
  dplyr::group_by(region) |>
  dplyr::summarise(
    total_spots = sum(n_spots),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

region_pal <- c(
  "Normal" = "#8dbfca",
  "Immune" = "#a5d296",
  "Stromal" = "#d4b483",
  "Tumor" = "#dd8a8a"
)

p_counts <- ggplot(spot_counts, aes(x = sample, y = n_spots, fill = region)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = region_pal) +
  theme_minimal(base_size = 11) +
  labs(y = "Spot count", x = NULL, fill = "Region", title = "Spot counts per sample × region")

ggsave(
  file.path(fig_dir, "Supp_Figure_QC_SpotCounts.pdf"),
  p_counts, width = 8, height = 4.6, bg = "white"
)
ggsave(
  file.path(fig_dir, "Supp_Figure_QC_SpotCounts.png"),
  p_counts, width = 8, height = 4.6, dpi = 600, bg = "white"
)

# ---------------------------
# 6) Supplementary Figure – Feature maps
# ---------------------------
Supp_FeatureMaps <- patchwork::wrap_plots(
  feature_panels,
  ncol = length(genes_hla),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Spatial feature maps – HLA-A / HLA-B / HLA-C across HCC sections",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

ggsave(
  file.path(fig_dir, "Supp_Figure_FeatureMaps.pdf"),
  Supp_FeatureMaps, width = 10, height = 16, bg = "white"
)
ggsave(
  file.path(fig_dir, "Supp_Figure_FeatureMaps.png"),
  Supp_FeatureMaps, width = 10, height = 16, dpi = 600, bg = "white"
)

# ---------------------------
# 7) Main plots: pooled violin
# ---------------------------
df_plot <- df_all |>
  dplyr::filter(region %in% regions_plot) |>
  dplyr::mutate(region = factor(as.character(region), levels = regions_plot))

labs_map <- c("HLA-A" = "A. HLA-A", "HLA-B" = "B. HLA-B", "HLA-C" = "C. HLA-C")
comparisons_pooled <- list(c("Normal", "Tumor"), c("Stromal", "Tumor"))

plot_pooled <- ggplot(df_plot, aes(x = region, y = expr, fill = region)) +
  geom_violin(trim = TRUE, alpha = 0.75, width = 0.9) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  stat_summary(fun = median, geom = "point", size = 1.0, color = "black") +
  facet_wrap(~gene, nrow = 1, scales = "fixed", labeller = labeller(gene = labs_map)) +
  scale_fill_manual(values = c("Normal" = "#8dbfca", "Stromal" = "#d4b483", "Tumor" = "#dd8a8a")) +
  scale_y_continuous(breaks = 0:6) +
  coord_cartesian(ylim = c(0, 6.5)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(y = "Expression") +
  ggpubr::stat_compare_means(
    comparisons = comparisons_pooled,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = FALSE,
    step.increase = 0.06,
    size = 4
  )

ggsave(
  file.path(fig_dir, "Main_Figure_Violin_pooled.pdf"),
  plot_pooled, width = 11.5, height = 4.5, bg = "white"
)
ggsave(
  file.path(fig_dir, "Main_Figure_Violin_pooled.png"),
  plot_pooled, width = 11.5, height = 4.5, dpi = 600, bg = "white"
)

# ---------------------------
# 8) Supplementary plots: balanced pooled violin
# ---------------------------
if (do_balanced_pool) {
  set.seed(7)
  
  df_balanced_plot <- df_plot |>
    dplyr::group_by(sample, region, gene) |>
    dplyr::mutate(.rand = runif(dplyr::n())) |>
    dplyr::arrange(sample, region, gene, .rand, .by_group = TRUE) |>
    dplyr::slice_head(n = k_per_sample_reg) |>
    dplyr::ungroup() |>
    dplyr::select(-.rand)
  
  plot_bal <- ggplot(df_balanced_plot, aes(x = region, y = expr, fill = region)) +
    geom_violin(trim = TRUE, alpha = 0.75, width = 0.9) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
    stat_summary(fun = median, geom = "point", size = 1.0, color = "black") +
    facet_wrap(~gene, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("Normal" = "#8dbfca", "Stromal" = "#d4b483", "Tumor" = "#dd8a8a")) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(y = "Expression (balanced)") +
    ggpubr::stat_compare_means(
      comparisons = comparisons_pooled,
      method = "wilcox.test",
      p.adjust.method = "BH",
      label = "p.signif",
      hide.ns = TRUE,
      size = 4
    )
  
  ggsave(
    file.path(fig_dir, "Supp_Figure_Violin_balanced.pdf"),
    plot_bal, width = 11.5, height = 4.5, bg = "white"
  )
  ggsave(
    file.path(fig_dir, "Supp_Figure_Violin_balanced.png"),
    plot_bal, width = 11.5, height = 4.5, dpi = 600, bg = "white"
  )
}

# ---------------------------
# 9) Spot-level statistics (S8)
# ---------------------------
spot_stats_fun <- function(data_src, label) {
  out <- list()
  
  for (s in unique(data_src$sample)) {
    for (g in unique(data_src$gene)) {
      dd <- subset(data_src, sample == s & gene == g)
      if (!nrow(dd)) next
      
      meds <- dd |>
        dplyr::group_by(region) |>
        dplyr::summarise(median = median(expr, na.rm = TRUE), .groups = "drop")
      
      for (comp in comparators) {
        d2 <- subset(dd, region %in% c(comp, "Tumor"))
        
        pval <- if (nrow(d2) > 2 && length(unique(d2$region)) == 2) {
          tryCatch(
            wilcox.test(expr ~ region, data = d2, exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
        
        out[[length(out) + 1]] <- data.frame(
          source = label,
          sample = s,
          gene = g,
          comparison = paste(comp, "vs Tumor"),
          med_comp = meds$median[meds$region == comp] %||% NA_real_,
          med_tumor = meds$median[meds$region == "Tumor"] %||% NA_real_,
          delta = (meds$median[meds$region == "Tumor"] %||% NA_real_) -
            (meds$median[meds$region == comp] %||% NA_real_),
          p_raw = pval
        )
      }
    }
  }
  
  df <- dplyr::bind_rows(out)
  if (nrow(df)) {
    df$p_adj_BH <- p.adjust(df$p_raw, method = "BH")
  }
  
  df
}

spot_stats_pooled <- spot_stats_fun(df_plot, "pooled")
spot_stats_bal <- if (do_balanced_pool) spot_stats_fun(df_balanced_plot, "balanced") else NULL
spot_stats_all <- dplyr::bind_rows(spot_stats_pooled, spot_stats_bal)

# ---------------------------
# 10) Patient-level medians (S9) + paired Wilcoxon (S7)
# ---------------------------
patient_medians <- df_all |>
  dplyr::group_by(sample, gene, region) |>
  dplyr::summarise(median_expr = median(expr, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(region %in% regions_plot) |>
  tidyr::pivot_wider(names_from = region, values_from = median_expr)

paired_summary <- function(wide, region_ref) {
  wd <- wide |>
    dplyr::filter(!is.na(Tumor) & !is.na(.data[[region_ref]]))
  
  wd |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      n_patients = dplyr::n(),
      med_Tumor = stats::median(Tumor, na.rm = TRUE),
      med_Region = stats::median(.data[[region_ref]], na.rm = TRUE),
      delta_med = stats::median(Tumor - .data[[region_ref]], na.rm = TRUE),
      p_paired = if (dplyr::n() >= 3) {
        stats::wilcox.test(
          Tumor, .data[[region_ref]],
          paired = TRUE,
          exact = FALSE
        )$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(comparison = paste("Tumor vs", region_ref))
}

res_patient <- dplyr::bind_rows(
  paired_summary(patient_medians, "Normal"),
  paired_summary(patient_medians, "Stromal")
)

res_patient$p_adj_BH <- p.adjust(res_patient$p_paired, method = "BH")

# ---------------------------
# 11) Export supplementary tables
# ---------------------------
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "PatientLevel_Stats")
openxlsx::writeData(wb, "PatientLevel_Stats", res_patient)

openxlsx::addWorksheet(wb, "SpotLevel_Stats")
openxlsx::writeData(wb, "SpotLevel_Stats", spot_stats_all)

openxlsx::addWorksheet(wb, "PatientLevel_Medians")
openxlsx::writeData(wb, "PatientLevel_Medians", patient_medians)

openxlsx::addWorksheet(wb, "QC_SpotCounts")
openxlsx::writeData(wb, "QC_SpotCounts", spot_counts)

openxlsx::addWorksheet(wb, "QC_SpotTotals")
openxlsx::writeData(wb, "QC_SpotTotals", spot_totals)

for (sh in c(
  "PatientLevel_Stats", "SpotLevel_Stats",
  "PatientLevel_Medians", "QC_SpotCounts", "QC_SpotTotals"
)) {
  openxlsx::setColWidths(wb, sh, cols = 1:50, widths = "auto")
}

excel_path <- file.path(tab_dir, "HLA_Analysis_Supplement.xlsx")
openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)

# ---------------------------
# 12) Log
# ---------------------------
sink(file.path(log_dir, "04_spatial_transcriptomics_log.txt"))
cat("04_spatial_transcriptomics.R completed\n")
cat("Input directory: ", data_dir, "\n", sep = "")
cat("Total RDS files: ", length(rds_files), "\n", sep = "")
cat("Tables directory: ", tab_dir, "\n", sep = "")
cat("Figures directory: ", fig_dir, "\n", sep = "")
cat("Log directory: ", log_dir, "\n", sep = "")
cat("Balanced pooled analysis enabled: ", do_balanced_pool, "\n", sep = "")
sink()

cat("Done: 04_spatial_transcriptomics.R\n")