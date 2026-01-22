#!/usr/bin/env Rscript
# ==============================================================================
# Interactive Plotting Functions for DESeq2 RNA-seq Analysis
# Generates HTML files with hover tooltips for detailed information
# ==============================================================================
#
# Functions:
#   - plot_volcano_interactive()    : Volcano plots with gene hover info
#   - plot_nes_heatmap_interactive(): NES heatmaps with p-value, FDR, and * in cells
#   - plot_eigencorplot_interactive(): PC-metadata correlations with FDR
#   - plot_deg_heatmap_interactive(): log2FC heatmap for genes
#   - plot_expression_heatmap_interactive(): Z-scored expression heatmap
#   - generate_interactive_plots()  : Generate all plots from analysis RDS
#
# Hover Information:
#   - Volcano: gene symbol, log2FC, p-value, padj, baseMean
#   - NES heatmap: signature, comparison, NES, p-value, FDR
#   - Eigencorplot: variable, PC, correlation (r), p-value, FDR
#   - DEG heatmap: gene, comparison, log2FC, p-value, FDR
#
# Cell Annotations (* significance stars displayed in cells):
#   - NES heatmap: *, **, ***, **** based on adjusted p-value
#
# Usage:
#   source("interactive_plots.R")
#   generate_interactive_plots("analysis_data.rds", "output_dir", "deg_dir")
#
# ==============================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(heatmaply)
  library(htmlwidgets)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RColorBrewer)
})

# ==============================================================================
# Interactive Volcano Plot
# ==============================================================================

#' Create interactive volcano plot with hover information
#'
#' @param deg_data Data frame with columns: symbol, log2FoldChange, padj, baseMean
#' @param comparison Name of the comparison
#' @param output_file Output HTML file path
#' @param lfc_threshold Log2FC threshold for significance (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param top_n Number of top genes to label (default: 20)
#' @return plotly object (also saves HTML)
plot_volcano_interactive <- function(deg_data,
                                     comparison = "Comparison",
                                     output_file = "volcano_interactive.html",
                                     lfc_threshold = 1,
                                     padj_threshold = 0.05,
                                     top_n = 20) {

  # Prepare data
  df <- deg_data %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      significance = case_when(
        padj < padj_threshold & log2FoldChange > lfc_threshold ~ "Up",
        padj < padj_threshold & log2FoldChange < -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      # Create hover text
      hover_text = paste0(
        "<b>", symbol, "</b><br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "padj: ", formatC(padj, format = "e", digits = 2), "<br>",
        "baseMean: ", round(baseMean, 1)
      )
    )

  # Color mapping
  colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "#BDC3C7")

  # Get top genes to label
  top_genes <- df %>%
    filter(significance != "NS") %>%
    arrange(padj) %>%
    head(top_n) %>%
    pull(symbol)

  df <- df %>%
    mutate(label = ifelse(symbol %in% top_genes, symbol, ""))

  # Create plotly
  p <- plot_ly(df,
               x = ~log2FoldChange,
               y = ~neg_log10_padj,
               color = ~significance,
               colors = colors,
               type = "scatter",
               mode = "markers",
               text = ~hover_text,
               hoverinfo = "text",
               marker = list(size = 6, opacity = 0.7)) %>%
    layout(
      title = list(text = paste0("<b>Volcano Plot: ", comparison, "</b>"),
                   font = list(size = 16)),
      xaxis = list(title = "log2(Fold Change)", zeroline = TRUE,
                   zerolinecolor = "#999999", zerolinewidth = 1),
      yaxis = list(title = "-log10(adjusted p-value)"),
      shapes = list(
        # Horizontal line at padj threshold
        list(type = "line", x0 = min(df$log2FoldChange, na.rm = TRUE),
             x1 = max(df$log2FoldChange, na.rm = TRUE),
             y0 = -log10(padj_threshold), y1 = -log10(padj_threshold),
             line = list(dash = "dash", color = "#999999", width = 1)),
        # Vertical lines at lfc thresholds
        list(type = "line", x0 = lfc_threshold, x1 = lfc_threshold,
             y0 = 0, y1 = max(df$neg_log10_padj, na.rm = TRUE),
             line = list(dash = "dash", color = "#999999", width = 1)),
        list(type = "line", x0 = -lfc_threshold, x1 = -lfc_threshold,
             y0 = 0, y1 = max(df$neg_log10_padj, na.rm = TRUE),
             line = list(dash = "dash", color = "#999999", width = 1))
      ),
      legend = list(title = list(text = "Significance"))
    ) %>%
    # Add gene labels for top genes
    add_annotations(
      data = df %>% filter(label != ""),
      x = ~log2FoldChange,
      y = ~neg_log10_padj,
      text = ~label,
      showarrow = FALSE,
      font = list(size = 9),
      yshift = 10
    )

  # Add summary annotation
  n_up <- sum(df$significance == "Up")
  n_down <- sum(df$significance == "Down")
  p <- p %>%
    add_annotations(
      x = 0.02, y = 0.98, xref = "paper", yref = "paper",
      text = paste0("Up: ", n_up, " | Down: ", n_down),
      showarrow = FALSE,
      font = list(size = 12, color = "#333333"),
      bgcolor = "rgba(255,255,255,0.8)",
      bordercolor = "#cccccc"
    )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Interactive Heatmap (Expression or DEG)
# ==============================================================================

#' Create interactive heatmap with hover information
#'
#' @param mat Matrix with genes as rows, samples/comparisons as columns
#' @param title Plot title
#' @param output_file Output HTML file path
#' @param row_annotation Optional data frame for row annotations
#' @param col_annotation Optional data frame for column annotations
#' @param scale Scale rows ("row"), columns ("column"), or none ("none")
#' @param colors Color palette
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_cols Cluster columns (default: TRUE)
#' @return heatmaply object (also saves HTML)
plot_heatmap_interactive <- function(mat,
                                     title = "Heatmap",
                                     output_file = "heatmap_interactive.html",
                                     row_annotation = NULL,
                                     col_annotation = NULL,
                                     scale = "row",
                                     colors = c("blue", "white", "red"),
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     show_rownames = TRUE,
                                     fontsize_row = 8,
                                     fontsize_col = 10) {

  # Create hover text matrix
  hover_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      hover_mat[i, j] <- paste0(
        "<b>", rownames(mat)[i], "</b><br>",
        "Sample: ", colnames(mat)[j], "<br>",
        "Value: ", round(mat[i, j], 3)
      )
    }
  }

  # Create heatmaply
  p <- heatmaply(
    mat,
    main = title,
    scale = scale,
    colors = colorRampPalette(colors)(100),
    dendrogram = ifelse(cluster_rows && cluster_cols, "both",
                        ifelse(cluster_rows, "row",
                               ifelse(cluster_cols, "column", "none"))),
    Rowv = cluster_rows,
    Colv = cluster_cols,
    showticklabels = c(TRUE, show_rownames),
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    custom_hovertext = hover_mat,
    row_side_colors = row_annotation,
    col_side_colors = col_annotation,
    plot_method = "plotly"
  )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Interactive NES Heatmap (for GSEA results)
# ==============================================================================

#' Create interactive NES heatmap for GSEA results
#'
#' @param gsea_data Data frame with columns: ID, comparison, NES, p.adjust
#' @param output_file Output HTML file path
#' @param signatures Optional vector of signature names to include
#' @param source_filter Optional source to filter (e.g., "CUSTOM", "GO:BP")
#' @return heatmaply object (also saves HTML)
plot_nes_heatmap_interactive <- function(gsea_data,
                                         output_file = "nes_heatmap_interactive.html",
                                         signatures = NULL,
                                         source_filter = NULL,
                                         title = "GSEA NES Heatmap") {

  df <- gsea_data

  # Filter by source if specified
  if (!is.null(source_filter) && "source" %in% colnames(df)) {
    df <- df %>% filter(source == source_filter)
  }

  # Filter by signatures if specified
  if (!is.null(signatures)) {
    df <- df %>% filter(ID %in% signatures)
  }

  # Create NES matrix
  nes_mat <- df %>%
    select(ID, comparison, NES) %>%
    pivot_wider(names_from = comparison, values_from = NES) %>%
    column_to_rownames("ID") %>%
    as.matrix()

  # Create pvalue matrix for hover
  pval_mat <- df %>%
    select(ID, comparison, pvalue) %>%
    pivot_wider(names_from = comparison, values_from = pvalue) %>%
    column_to_rownames("ID") %>%
    as.matrix()

  # Create padj/FDR matrix for hover
  padj_mat <- df %>%
    select(ID, comparison, p.adjust) %>%
    pivot_wider(names_from = comparison, values_from = p.adjust) %>%
    column_to_rownames("ID") %>%
    as.matrix()

  # Create significance annotation matrix for cell notes
  sig_mat <- matrix("", nrow = nrow(nes_mat), ncol = ncol(nes_mat))
  rownames(sig_mat) <- rownames(nes_mat)
  colnames(sig_mat) <- colnames(nes_mat)

  # Create hover text matrix with both pvalue and FDR
  hover_mat <- matrix("", nrow = nrow(nes_mat), ncol = ncol(nes_mat))
  for (i in 1:nrow(nes_mat)) {
    for (j in 1:ncol(nes_mat)) {
      sig <- ""
      if (!is.na(padj_mat[i, j])) {
        if (padj_mat[i, j] < 0.05) sig <- "*"
        if (padj_mat[i, j] < 0.01) sig <- "**"
        if (padj_mat[i, j] < 0.001) sig <- "***"
        if (padj_mat[i, j] < 0.0001) sig <- "****"
      }
      sig_mat[i, j] <- sig

      hover_mat[i, j] <- paste0(
        "<b>", rownames(nes_mat)[i], "</b><br>",
        "Comparison: ", colnames(nes_mat)[j], "<br>",
        "NES: ", round(nes_mat[i, j], 3), "<br>",
        "p-value: ", formatC(pval_mat[i, j], format = "e", digits = 2), "<br>",
        "FDR: ", formatC(padj_mat[i, j], format = "e", digits = 2)
      )
    }
  }

  # Replace NA with 0 for visualization
  nes_mat[is.na(nes_mat)] <- 0

  # Create heatmaply with cell annotations
  p <- heatmaply(
    nes_mat,
    main = title,
    scale = "none",
    colors = colorRampPalette(c("blue3", "white", "red3"))(100),
    limits = c(-max(abs(nes_mat), na.rm = TRUE), max(abs(nes_mat), na.rm = TRUE)),
    dendrogram = "both",
    showticklabels = c(TRUE, TRUE),
    fontsize_row = 9,
    fontsize_col = 9,
    custom_hovertext = hover_mat,
    cellnote = sig_mat,
    cellnote_textposition = "middle center",
    plot_method = "plotly"
  )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Interactive Eigencorplot
# ==============================================================================

#' Create interactive eigencorplot (PC-metadata correlations)
#'
#' @param pca_obj PCA object from PCAtools
#' @param metadata Data frame with metadata (rows = samples)
#' @param output_file Output HTML file path
#' @param n_pcs Number of PCs to include (default: 10)
#' @return plotly object (also saves HTML)
plot_eigencorplot_interactive <- function(pca_obj,
                                          metadata = NULL,
                                          output_file = "eigencorplot_interactive.html",
                                          n_pcs = 10,
                                          title = "PC-Metadata Correlations") {

  # Use metadata from PCA object if not provided
  if (is.null(metadata)) {
    metadata <- pca_obj$metadata
  }

  # Get PC scores
  pcs <- pca_obj$rotated[, 1:min(n_pcs, ncol(pca_obj$rotated))]

  # Ensure sample order matches
  common_samples <- intersect(rownames(pcs), rownames(metadata))
  pcs <- pcs[common_samples, ]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Calculate correlations and p-values
  n_vars <- ncol(metadata)
  n_pcs_actual <- ncol(pcs)

  cor_mat <- matrix(NA, nrow = n_vars, ncol = n_pcs_actual)
  pval_mat <- matrix(NA, nrow = n_vars, ncol = n_pcs_actual)
  rownames(cor_mat) <- colnames(metadata)
  colnames(cor_mat) <- colnames(pcs)
  rownames(pval_mat) <- colnames(metadata)
  colnames(pval_mat) <- colnames(pcs)

  for (i in 1:n_vars) {
    for (j in 1:n_pcs_actual) {
      tryCatch({
        test <- cor.test(as.numeric(metadata[, i]), pcs[, j], method = "pearson")
        cor_mat[i, j] <- test$estimate
        pval_mat[i, j] <- test$p.value
      }, error = function(e) {
        cor_mat[i, j] <<- NA
        pval_mat[i, j] <<- NA
      })
    }
  }

  # Calculate FDR-adjusted p-values (BH correction across all tests)
  pval_vec <- as.vector(pval_mat)
  fdr_vec <- p.adjust(pval_vec, method = "BH")
  fdr_mat <- matrix(fdr_vec, nrow = n_vars, ncol = n_pcs_actual)
  rownames(fdr_mat) <- rownames(pval_mat)
  colnames(fdr_mat) <- colnames(pval_mat)

  # Create hover text matrix with both p-value and FDR
  hover_mat <- matrix("", nrow = n_vars, ncol = n_pcs_actual)
  for (i in 1:n_vars) {
    for (j in 1:n_pcs_actual) {
      sig <- ""
      if (!is.na(fdr_mat[i, j])) {
        if (fdr_mat[i, j] < 0.05) sig <- " *"
        if (fdr_mat[i, j] < 0.01) sig <- " **"
        if (fdr_mat[i, j] < 0.001) sig <- " ***"
        if (fdr_mat[i, j] < 0.0001) sig <- " ****"
      }

      hover_mat[i, j] <- paste0(
        "<b>", rownames(cor_mat)[i], "</b> vs <b>", colnames(cor_mat)[j], "</b><br>",
        "Correlation (r): ", round(cor_mat[i, j], 3), "<br>",
        "p-value: ", formatC(pval_mat[i, j], format = "e", digits = 2), "<br>",
        "FDR: ", formatC(fdr_mat[i, j], format = "e", digits = 2)
      )
    }
  }

  # Create heatmaply
  p <- heatmaply(
    cor_mat,
    main = title,
    scale = "none",
    colors = colorRampPalette(c("blue3", "white", "red3"))(100),
    limits = c(-1, 1),
    dendrogram = "row",
    Colv = FALSE,  # Don't reorder PCs
    showticklabels = c(TRUE, TRUE),
    fontsize_row = 10,
    fontsize_col = 10,
    custom_hovertext = hover_mat,
    plot_method = "plotly",
    cellnote = round(cor_mat, 2),
    cellnote_textposition = "middle center"
  )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Interactive DEG Heatmap (log2FC across comparisons)
# ==============================================================================

#' Create interactive DEG heatmap showing log2FC across comparisons
#'
#' @param deg_data Data frame with columns: symbol, comparison, log2FoldChange, padj
#' @param genes Vector of gene symbols to include
#' @param output_file Output HTML file path
#' @return heatmaply object (also saves HTML)
plot_deg_heatmap_interactive <- function(deg_data,
                                         genes,
                                         output_file = "deg_heatmap_interactive.html",
                                         title = "DEG Heatmap (log2FC)") {

  # Filter for specified genes
  df <- deg_data %>%
    filter(symbol %in% genes)

  # Create log2FC matrix
  lfc_mat <- df %>%
    select(symbol, comparison, log2FoldChange) %>%
    pivot_wider(names_from = comparison, values_from = log2FoldChange) %>%
    column_to_rownames("symbol") %>%
    as.matrix()

  # Create padj matrix for hover
  padj_mat <- df %>%
    select(symbol, comparison, padj) %>%
    pivot_wider(names_from = comparison, values_from = padj) %>%
    column_to_rownames("symbol") %>%
    as.matrix()

  # Create hover text matrix
  hover_mat <- matrix("", nrow = nrow(lfc_mat), ncol = ncol(lfc_mat))
  for (i in 1:nrow(lfc_mat)) {
    for (j in 1:ncol(lfc_mat)) {
      sig <- ""
      if (!is.na(padj_mat[i, j])) {
        if (padj_mat[i, j] < 0.05) sig <- " *"
        if (padj_mat[i, j] < 0.01) sig <- " **"
        if (padj_mat[i, j] < 0.001) sig <- " ***"
      }

      hover_mat[i, j] <- paste0(
        "<b>", rownames(lfc_mat)[i], "</b><br>",
        "Comparison: ", colnames(lfc_mat)[j], "<br>",
        "log2FC: ", round(lfc_mat[i, j], 3), "<br>",
        "padj: ", formatC(padj_mat[i, j], format = "e", digits = 2)
      )
    }
  }

  # Replace NA with 0
  lfc_mat[is.na(lfc_mat)] <- 0

  # Create heatmaply
  max_val <- max(abs(lfc_mat), na.rm = TRUE)

  p <- heatmaply(
    lfc_mat,
    main = title,
    scale = "none",
    colors = colorRampPalette(c("blue3", "white", "red3"))(100),
    limits = c(-max_val, max_val),
    dendrogram = "both",
    showticklabels = c(TRUE, TRUE),
    fontsize_row = 9,
    fontsize_col = 9,
    custom_hovertext = hover_mat,
    plot_method = "plotly"
  )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Interactive Expression Heatmap
# ==============================================================================

#' Create interactive expression heatmap
#'
#' @param expr_mat Expression matrix (genes x samples)
#' @param metadata Sample metadata with grouping column
#' @param genes Vector of genes to include
#' @param group_col Column name for grouping samples
#' @param output_file Output HTML file path
#' @param scale Scale by row (default: TRUE)
#' @return heatmaply object (also saves HTML)
plot_expression_heatmap_interactive <- function(expr_mat,
                                                metadata,
                                                genes,
                                                group_col = "Treatment",
                                                output_file = "expression_heatmap_interactive.html",
                                                scale = TRUE,
                                                title = "Expression Heatmap") {

  # Filter for genes
  genes_found <- genes[genes %in% rownames(expr_mat)]
  mat <- expr_mat[genes_found, , drop = FALSE]

  # Z-score if requested
  if (scale) {
    mat <- t(scale(t(mat)))
  }

  # Create hover text
  hover_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      sample_name <- colnames(mat)[j]
      group <- ""
      if (!is.null(metadata) && sample_name %in% rownames(metadata)) {
        group <- paste0("<br>Group: ", metadata[sample_name, group_col])
      }

      hover_mat[i, j] <- paste0(
        "<b>", rownames(mat)[i], "</b><br>",
        "Sample: ", sample_name, group, "<br>",
        "Value: ", round(mat[i, j], 3)
      )
    }
  }

  # Create column annotation
  col_annot <- NULL
  if (!is.null(metadata) && group_col %in% colnames(metadata)) {
    col_annot <- data.frame(
      Group = metadata[colnames(mat), group_col],
      row.names = colnames(mat)
    )
  }

  # Create heatmaply
  p <- heatmaply(
    mat,
    main = title,
    scale = "none",  # Already scaled above
    colors = colorRampPalette(c("blue3", "white", "red3"))(100),
    dendrogram = "both",
    showticklabels = c(TRUE, TRUE),
    fontsize_row = 9,
    fontsize_col = 8,
    custom_hovertext = hover_mat,
    col_side_colors = col_annot,
    plot_method = "plotly"
  )

  # Save HTML
  saveWidget(p, output_file, selfcontained = TRUE)
  cat("Saved:", output_file, "\n")

  return(p)
}


# ==============================================================================
# Batch Generation Function
# ==============================================================================

#' Generate all interactive plots for a DESeq2 analysis
#'
#' @param analysis_rds Path to analysis RDS file
#' @param output_dir Output directory for HTML files
#' @param deg_dir Directory containing DEG files
#' @param comparisons Vector of comparison names (NULL = all)
generate_interactive_plots <- function(analysis_rds,
                                       output_dir,
                                       deg_dir = NULL,
                                       comparisons = NULL) {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load analysis data
  data <- readRDS(analysis_rds)

  cat("=== Generating Interactive Plots ===\n")

  # 1. Eigencorplot (with GSVA scores merged into metadata, matching static version)
  cat("\n1. Eigencorplot...\n")
  tryCatch({
    # Get PCA metadata
    pca_obj <- data$pca
    meta <- pca_obj$metadata

    # Combine with GSVA scores if available (same as static eigencorplot)
    if (!is.null(data$gsva_scores)) {
      gsva <- data$gsva_scores
      gsva_t <- as.data.frame(t(gsva))
      common <- intersect(rownames(meta), rownames(gsva_t))
      meta_combined <- cbind(meta[common, , drop = FALSE], gsva_t[common, , drop = FALSE])
    } else {
      meta_combined <- meta
    }

    # Filter metadata - only exclude columns with no variance (same as static)
    meta_filtered <- data.frame(row.names = rownames(meta_combined))

    for (col in colnames(meta_combined)) {
      vals <- meta_combined[[col]]
      n_unique <- length(unique(vals[!is.na(vals)]))

      if (n_unique <= 1) {
        cat("  Excluding (no variance):", col, "\n")
        next
      }

      if (is.numeric(vals)) {
        meta_filtered[[col]] <- vals
      } else {
        meta_filtered[[paste0(col, "_num")]] <- as.numeric(as.factor(vals))
      }
    }

    cat("  Variables for eigencorplot:", ncol(meta_filtered), "\n")

    # Update PCA object with filtered metadata
    pca_obj_eigen <- pca_obj
    pca_obj_eigen$metadata <- meta_filtered

    plot_eigencorplot_interactive(
      pca_obj_eigen,
      output_file = file.path(output_dir, "eigencorplot_interactive.html"),
      title = "PC-Metadata Correlations"
    )
  }, error = function(e) warning("Eigencorplot failed: ", e$message))

  # 2. Volcano plots for each comparison
  cat("\n2. Volcano plots...\n")
  if (!is.null(deg_dir)) {
    deg_files <- list.files(deg_dir, pattern = "_summstats_all\\.txt$", full.names = TRUE)
    if (length(deg_files) > 0) {
      deg_all <- read.delim(deg_files[1])

      comps <- unique(deg_all$comparison)
      if (!is.null(comparisons)) {
        comps <- comps[comps %in% comparisons]
      }

      for (comp in comps) {
        cat("  -", comp, "\n")
        deg_comp <- deg_all %>% filter(comparison == comp)
        tryCatch({
          plot_volcano_interactive(
            deg_comp,
            comparison = comp,
            output_file = file.path(output_dir, paste0("volcano_", gsub("[^A-Za-z0-9]", "_", comp), "_interactive.html"))
          )
        }, error = function(e) warning("Volcano plot failed for ", comp, ": ", e$message))
      }
    }
  }

  # 3. NES heatmap for custom signatures
  cat("\n3. NES heatmap...\n")
  if (!is.null(deg_dir)) {
    gsea_files <- list.files(deg_dir, pattern = "_gsea_all\\.txt$", full.names = TRUE)
    if (length(gsea_files) > 0) {
      gsea_all <- read.delim(gsea_files[1])

      # Custom signatures only
      if ("source" %in% colnames(gsea_all)) {
        gsea_custom <- gsea_all %>% filter(source == "CUSTOM")
        if (nrow(gsea_custom) > 0) {
          tryCatch({
            plot_nes_heatmap_interactive(
              gsea_custom,
              output_file = file.path(output_dir, "nes_heatmap_custom_interactive.html"),
              title = "GSEA NES Heatmap (Custom Signatures)"
            )
          }, error = function(e) warning("NES heatmap failed: ", e$message))
        }
      }
    }
  }

  cat("\n=== Interactive plots saved to:", output_dir, "===\n")
}


# ==============================================================================
# Main (if run as script directly, not sourced)
# ==============================================================================

# Only run CLI when script is executed directly (not when source()'d)
# sys.nframe() == 0 means top-level execution, not sourced from another script
if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) >= 2) {
    analysis_rds <- args[1]
    output_dir <- args[2]
    deg_dir <- if (length(args) >= 3) args[3] else dirname(analysis_rds)

    generate_interactive_plots(analysis_rds, output_dir, deg_dir)
  } else {
    cat("Usage: Rscript interactive_plots.R <analysis.rds> <output_dir> [deg_dir]\n")
  }
}
