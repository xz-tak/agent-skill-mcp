#!/usr/bin/env Rscript
# =============================================================================
# RNA-seq Visualization Functions
# =============================================================================
# Contains 4 plotting functions for RNA-seq analysis visualization:
# 1. plot_signature_heatmap - GSEA/GSVA NES heatmap by comparison
# 2. plot_deg_heatmap - DEG log2FC heatmap by comparison
# 3. plot_expression_heatmap - Expression heatmap grouped by treatment
# 4. plot_gsva_boxplot - GSVA score boxplot by treatment group
#
# Interactive Output (default: interactive = TRUE):
#   - Generates HTML files alongside PNG outputs
#   - Hover tooltips show: gene/signature, comparison, values, p-value, FDR
#   - Cell annotations: *, **, ***, **** significance stars in heatmap cells
#     (signature and DEG heatmaps only)
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggpubr)
  library(RColorBrewer)
  library(grid)
  library(rstatix)  # For Games-Howell test
  library(car)      # For Levene's test
  library(ggsignif) # For significance brackets
  library(heatmaply) # For interactive heatmaps
  library(plotly)    # For interactive plots
  library(htmlwidgets) # For saving HTML widgets
})

# =============================================================================
# Helper Functions
# =============================================================================

#' Convert p-value to significance symbols
#' @param p numeric p-value
#' @return character significance symbol
pval_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

#' Read genes from signature file(s)
#' @param signature_files character vector of file paths
#' @return character vector of gene symbols
read_signature_genes <- function(signature_files) {
  genes <- c()
  for (file in signature_files) {
    if (file.exists(file)) {
      file_genes <- readLines(file, warn = FALSE)
      file_genes <- trimws(file_genes)
      file_genes <- file_genes[file_genes != "" & !startsWith(file_genes, "#")]
      genes <- c(genes, file_genes)
    } else {
      warning(paste("Signature file not found:", file))
    }
  }
  unique(genes)
}

#' Create Blue-White-Red color function centered at 0
#' @param values numeric vector of values to map
#' @param symmetric logical, if TRUE make symmetric around 0
#' @return colorRamp2 function
create_bwr_colormap <- function(values, symmetric = TRUE) {
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
  }

  if (symmetric) {
    max_abs <- max(abs(values), na.rm = TRUE)
    if (max_abs == 0) max_abs <- 1
    breaks <- c(-max_abs, 0, max_abs)
  } else {
    min_val <- min(values, na.rm = TRUE)
    max_val <- max(values, na.rm = TRUE)
    breaks <- c(min_val, 0, max_val)
  }

  colorRamp2(breaks, c("blue", "white", "red"))
}

#' Get static PNG output path in expr_comp_heatmap_boxplot subfolder
#' @param output_file Original PNG output file path
#' @return Path for static PNG file in figures/expr_comp_heatmap_boxplot subfolder
get_static_path <- function(output_file) {
  output_dir <- dirname(output_file)
  base_name <- basename(output_file)

  # Navigate to find figures directory, then create expr_comp_heatmap_boxplot subfolder
  # Supports: figures/file.png -> figures/expr_comp_heatmap_boxplot/file.png

  # Find parent of figures or use output_dir directly
  parent_dir <- dirname(output_dir)
  if (basename(parent_dir) == "figures" || basename(output_dir) == "figures") {
    if (basename(output_dir) == "figures") {
      static_dir <- file.path(output_dir, "expr_comp_heatmap_boxplot")
    } else {
      # Already in a subfolder of figures
      static_dir <- output_dir
    }
  } else {
    # Fallback: create expr_comp_heatmap_boxplot subfolder in output directory
    static_dir <- file.path(output_dir, "expr_comp_heatmap_boxplot")
  }

  dir.create(static_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(static_dir, base_name)
}

#' Get interactive HTML output path
#' @param output_file Original PNG output file path
#' @return Path for interactive HTML file in interactive/expr_comp_heatmap_boxplot subfolder
get_interactive_path <- function(output_file) {
  output_dir <- dirname(output_file)
  base_name <- tools::file_path_sans_ext(basename(output_file))

  # Navigate up to find figures directory, then create interactive subfolder
  # Supports both: figures/xxx/file.png -> figures/interactive/expr_comp_heatmap_boxplot/
  #               output_dir/file.png -> output_dir/interactive/expr_comp_heatmap_boxplot/

  # Find parent of figures or use output_dir directly
  parent_dir <- dirname(output_dir)
  if (basename(parent_dir) == "figures" || basename(output_dir) == "figures") {
    if (basename(output_dir) == "figures") {
      interactive_dir <- file.path(output_dir, "interactive", "expr_comp_heatmap_boxplot")
    } else {
      interactive_dir <- file.path(parent_dir, "interactive", "expr_comp_heatmap_boxplot")
    }
  } else {
    # Fallback: create interactive subfolder in output directory
    interactive_dir <- file.path(output_dir, "interactive", "expr_comp_heatmap_boxplot")
  }

  dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(interactive_dir, paste0(base_name, "_interactive.html"))
}

#' Get statistics TSV output path in expr_comp_heatmap_boxplot subfolder
#' @param output_file Original PNG output file path
#' @param suffix Optional suffix to append before .tsv (e.g., "_pairwise_stats")
#' @return Path for statistics TSV file in figures/expr_comp_heatmap_boxplot subfolder
get_stats_path <- function(output_file, suffix = "_pairwise_stats") {
  output_dir <- dirname(output_file)
  base_name <- tools::file_path_sans_ext(basename(output_file))

  # Use same logic as get_static_path to find the correct subfolder
  parent_dir <- dirname(output_dir)
  if (basename(parent_dir) == "figures" || basename(output_dir) == "figures") {
    if (basename(output_dir) == "figures") {
      stats_dir <- file.path(output_dir, "expr_comp_heatmap_boxplot")
    } else {
      # Already in a subfolder of figures
      stats_dir <- output_dir
    }
  } else {
    # Fallback: create expr_comp_heatmap_boxplot subfolder in output directory
    stats_dir <- file.path(output_dir, "expr_comp_heatmap_boxplot")
  }

  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(stats_dir, paste0(base_name, suffix, ".tsv"))
}

# =============================================================================
# Function 1: GSEA/GSVA Signature Heatmap
# =============================================================================

#' Plot GSEA/GSVA signature heatmap
#'
#' Creates a clustered heatmap with signatures as rows and comparisons as columns.
#' Cell colors represent NES values (blue-white-red centered at 0).
#' Significance annotations (*,**,***,****) based on adjusted p-values.
#'
#' @param gsea_data Either a file path to GSEA results table or data.frame/RDS object
#' @param signatures Character vector of signature names to include (NULL = all)
#' @param output_file Output file path (PNG)
#' @param nes_col Column name for NES values (default: "NES")
#' @param pval_col Column name for adjusted p-values (default: "p.adjust")
#' @param sig_col Column name for signature/term IDs (default: "ID")
#' @param comp_col Column name for comparison (default: "comparison")
#' @param source_col Column name for source filter (default: "source")
#' @param source_filter Filter for source column (default: "CUSTOM")
#' @param title Plot title (default: NULL, auto-generates "Custom Signatures NES" or based on source_filter)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @param fontsize Base font size (default: 10)
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_cols Cluster columns (default: TRUE)
#' @param show_row_names Show row names (default: TRUE)
#' @param show_col_names Show column names (default: TRUE)
#' @param row_names_max_width Max width for row names (default: 20)
#' @param col_names_max_height Max height for column names (default: 15)
#' @return Invisibly returns the Heatmap object
#' @export
plot_signature_heatmap <- function(
    gsea_data,
    signatures = NULL,
    output_file = "figures/signature_NES_heatmap.png",
    nes_col = "NES",
    pval_col = "p.adjust",
    sig_col = "ID",
    comp_col = "comparison",
    source_col = "source",
    source_filter = "CUSTOM",
    title = NULL,
    width = 12,
    height = 8,
    fontsize = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_row_names = TRUE,
    show_col_names = TRUE,
    row_names_max_width = 20,
    col_names_max_height = 15,
    interactive = TRUE
) {
  # Load data
  if (is.character(gsea_data)) {
    if (grepl("\\.rds$", gsea_data, ignore.case = TRUE)) {
      gsea_df <- readRDS(gsea_data)
    } else {
      gsea_df <- read.delim(gsea_data, stringsAsFactors = FALSE)
    }
  } else {
    gsea_df <- as.data.frame(gsea_data)
  }

  # Filter by source if specified
  if (!is.null(source_filter) && source_col %in% colnames(gsea_df)) {
    gsea_df <- gsea_df[gsea_df[[source_col]] == source_filter, ]
  }

  # Filter by signatures if specified
  if (!is.null(signatures)) {
    gsea_df <- gsea_df[gsea_df[[sig_col]] %in% signatures, ]
  }

  if (nrow(gsea_df) == 0) {
    stop("No data remaining after filtering. Check your signature names and source filter.")
  }

  # Pivot to matrices
  nes_wide <- gsea_df %>%
    select(all_of(c(sig_col, comp_col, nes_col))) %>%
    pivot_wider(names_from = all_of(comp_col), values_from = all_of(nes_col)) %>%
    column_to_rownames(sig_col)

  # Extract raw p-value if available
  if ("pvalue" %in% colnames(gsea_df)) {
    rawpval_wide <- gsea_df %>%
      select(all_of(c(sig_col, comp_col, "pvalue"))) %>%
      pivot_wider(names_from = all_of(comp_col), values_from = "pvalue") %>%
      column_to_rownames(sig_col)
    rawpval_mat <- as.matrix(rawpval_wide)
  } else {
    rawpval_mat <- NULL
  }

  pval_wide <- gsea_df %>%
    select(all_of(c(sig_col, comp_col, pval_col))) %>%
    pivot_wider(names_from = all_of(comp_col), values_from = all_of(pval_col)) %>%
    column_to_rownames(sig_col)

  # Convert to matrices
  nes_mat <- as.matrix(nes_wide)
  pval_mat <- as.matrix(pval_wide)  # This is FDR/padj

  # Replace NA with 0 for NES
  nes_mat[is.na(nes_mat)] <- 0

  # Create significance annotation matrix
  sig_mat <- matrix(
    sapply(as.vector(pval_mat), pval_to_stars),
    nrow = nrow(pval_mat),
    dimnames = dimnames(pval_mat)
  )

  # Create color function
  col_fun <- create_bwr_colormap(nes_mat)

  # Auto-generate title if not provided
  if (is.null(title)) {
    if (!is.null(source_filter)) {
      title <- paste0(source_filter, " Signatures NES")
    } else {
      title <- "Signature NES Heatmap"
    }
  }

  # Create heatmap with cell annotations
  ht <- Heatmap(
    nes_mat,
    name = "NES",
    col = col_fun,
    column_title = title,
    column_title_gp = gpar(fontsize = fontsize + 2, fontface = "bold"),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize - 1),
    row_names_max_width = unit(row_names_max_width, "cm"),
    column_names_max_height = unit(col_names_max_height, "cm"),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = "NES",
      title_gp = gpar(fontsize = fontsize),
      labels_gp = gpar(fontsize = fontsize - 1)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(sig_mat[i, j]) && sig_mat[i, j] != "") {
        grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = fontsize - 2, col = "black"))
      }
    }
  )

  # Get output path in expr_comp_heatmap_boxplot subfolder
  static_output <- get_static_path(output_file)

  # Ensure output directory exists
  dir.create(dirname(static_output), showWarnings = FALSE, recursive = TRUE)

  # Save plot
  png(static_output, width = width, height = height, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right")
  dev.off()

  message(paste("Saved signature heatmap to:", static_output))

  # Generate interactive HTML version
  if (interactive) {
    # Create hover text matrix with both p-value and FDR
    hover_mat <- matrix("", nrow = nrow(nes_mat), ncol = ncol(nes_mat))
    for (i in 1:nrow(nes_mat)) {
      for (j in 1:ncol(nes_mat)) {
        sig_text <- sig_mat[i, j]
        pval_text <- if (!is.null(rawpval_mat)) paste0("p-value: ", formatC(rawpval_mat[i, j], format = "e", digits = 2), "<br>") else ""
        hover_mat[i, j] <- paste0(
          "<b>", rownames(nes_mat)[i], "</b><br>",
          "Comparison: ", colnames(nes_mat)[j], "<br>",
          "NES: ", round(nes_mat[i, j], 3), "<br>",
          pval_text,
          "FDR: ", formatC(pval_mat[i, j], format = "e", digits = 2)
        )
      }
    }

    max_val <- max(abs(nes_mat), na.rm = TRUE)

    p_interactive <- heatmaply(
      nes_mat,
      main = title,
      scale = "none",
      colors = colorRampPalette(c("blue3", "white", "red3"))(100),
      limits = c(-max_val, max_val),
      dendrogram = if (cluster_rows && cluster_cols) "both" else if (cluster_rows) "row" else if (cluster_cols) "column" else "none",
      showticklabels = c(show_row_names, show_col_names),
      fontsize_row = fontsize,
      fontsize_col = fontsize - 1,
      custom_hovertext = hover_mat,
      plot_method = "plotly",
      cellnote = sig_mat,
      cellnote_textposition = "middle center"
    )

    html_file <- get_interactive_path(output_file)
    saveWidget(p_interactive, html_file, selfcontained = TRUE)
    message(paste("Saved interactive heatmap to:", html_file))
  }

  invisible(ht)
}

# =============================================================================
# Function 2: DEG Gene Heatmap (log2FC by comparison)
# =============================================================================

#' Plot DEG gene heatmap
#'
#' Creates a clustered heatmap with genes as rows and comparisons as columns.
#' Cell colors represent log2FC values (blue-white-red centered at 0).
#' Significance annotations (*,**,***,****) based on adjusted p-values.
#'
#' @param deg_data Either a file path to DEG results table or data.frame/RDS object
#' @param signature_files Character vector of signature file paths (NULL = use genes param)
#' @param genes Character vector of gene symbols (used if signature_files is NULL)
#' @param output_file Output file path (PNG)
#' @param title Plot title (default: NULL, auto-generates from signature name)
#' @param export_individual Export individual plots per signature file (default: TRUE when multiple signature_files)
#' @param lfc_col Column name for log2FC values (default: "log2FoldChange")
#' @param pval_col Column name for adjusted p-values (default: "padj")
#' @param gene_col Column name for gene symbols (default: "symbol")
#' @param comp_col Column name for comparison (default: "comparison")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @param fontsize Base font size (default: 8)
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_cols Cluster columns (default: TRUE)
#' @param show_row_names Show row names (default: TRUE)
#' @param show_col_names Show column names (default: TRUE)
#' @param row_names_max_width Max width for row names (default: 15)
#' @param col_names_max_height Max height for column names (default: 15)
#' @return Invisibly returns the Heatmap object (or list of Heatmap objects if export_individual=TRUE)
#' @export
plot_deg_heatmap <- function(
    deg_data,
    signature_files = NULL,
    genes = NULL,
    output_file = "figures/deg_gene_heatmap.png",
    title = NULL,
    export_individual = TRUE,
    lfc_col = "log2FoldChange",
    pval_col = "padj",
    gene_col = "symbol",
    comp_col = "comparison",
    width = 12,
    height = 10,
    fontsize = 8,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_row_names = TRUE,
    show_col_names = TRUE,
    row_names_max_width = 15,
    col_names_max_height = 15,
    interactive = TRUE
) {
  # Load data
  if (is.character(deg_data)) {
    if (grepl("\\.rds$", deg_data, ignore.case = TRUE)) {
      deg_df <- readRDS(deg_data)
    } else {
      deg_df <- read.delim(deg_data, stringsAsFactors = FALSE)
    }
  } else {
    deg_df <- as.data.frame(deg_data)
  }

  # Helper function to create single heatmap
  create_deg_heatmap <- function(deg_data_filtered, plot_title, out_file, base_fontsize, base_width, base_height) {
    # Pivot to matrices
    lfc_wide <- deg_data_filtered %>%
      select(all_of(c(gene_col, comp_col, lfc_col))) %>%
      pivot_wider(names_from = all_of(comp_col), values_from = all_of(lfc_col)) %>%
      column_to_rownames(gene_col)

    # Extract raw p-value if available
    if ("pvalue" %in% colnames(deg_data_filtered)) {
      rawpval_wide <- deg_data_filtered %>%
        select(all_of(c(gene_col, comp_col, "pvalue"))) %>%
        pivot_wider(names_from = all_of(comp_col), values_from = "pvalue") %>%
        column_to_rownames(gene_col)
      rawpval_mat <- as.matrix(rawpval_wide)
    } else {
      rawpval_mat <- NULL
    }

    pval_wide <- deg_data_filtered %>%
      select(all_of(c(gene_col, comp_col, pval_col))) %>%
      pivot_wider(names_from = all_of(comp_col), values_from = all_of(pval_col)) %>%
      column_to_rownames(gene_col)

    # Convert to matrices
    lfc_mat <- as.matrix(lfc_wide)
    pval_mat <- as.matrix(pval_wide)  # This is FDR/padj

    # Replace NA with 0 for log2FC
    lfc_mat[is.na(lfc_mat)] <- 0

    # Create significance annotation matrix
    sig_mat <- matrix(
      sapply(as.vector(pval_mat), pval_to_stars),
      nrow = nrow(pval_mat),
      dimnames = dimnames(pval_mat)
    )

    # Create color function
    col_fun <- create_bwr_colormap(lfc_mat)

    # Adjust height based on number of genes
    n_genes <- nrow(lfc_mat)
    adj_height <- base_height
    adj_fontsize <- base_fontsize
    if (n_genes > 50) {
      adj_height <- max(base_height, n_genes * 0.15)
      if (n_genes > 100) adj_fontsize <- max(5, base_fontsize - 2)
    }

    # Create heatmap with cell annotations
    ht <- Heatmap(
      lfc_mat,
      name = "log2FC",
      col = col_fun,
      column_title = plot_title,
      column_title_gp = gpar(fontsize = adj_fontsize + 2, fontface = "bold"),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      show_row_names = show_row_names,
      show_column_names = show_col_names,
      row_names_gp = gpar(fontsize = adj_fontsize),
      column_names_gp = gpar(fontsize = adj_fontsize),
      row_names_max_width = unit(row_names_max_width, "cm"),
      column_names_max_height = unit(col_names_max_height, "cm"),
      column_names_rot = 90,
      heatmap_legend_param = list(
        title = "log2FC",
        title_gp = gpar(fontsize = adj_fontsize + 2),
        labels_gp = gpar(fontsize = adj_fontsize)
      ),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (!is.na(sig_mat[i, j]) && sig_mat[i, j] != "") {
          grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = max(4, adj_fontsize - 3), col = "black"))
        }
      }
    )

    # Get output path in expr_comp_heatmap_boxplot subfolder
    static_out_file <- get_static_path(out_file)

    # Ensure output directory exists
    dir.create(dirname(static_out_file), showWarnings = FALSE, recursive = TRUE)

    # Save plot
    png(static_out_file, width = base_width, height = adj_height, units = "in", res = 300)
    draw(ht, heatmap_legend_side = "right")
    dev.off()

    message(paste("Saved DEG heatmap to:", static_out_file))
    message(paste("Genes plotted:", nrow(lfc_mat)))

    # Generate interactive HTML version
    if (interactive) {
      hover_mat <- matrix("", nrow = nrow(lfc_mat), ncol = ncol(lfc_mat))
      for (i in 1:nrow(lfc_mat)) {
        for (j in 1:ncol(lfc_mat)) {
          sig_text <- sig_mat[i, j]
          pval_text <- if (!is.null(rawpval_mat)) paste0("p-value: ", formatC(rawpval_mat[i, j], format = "e", digits = 2), "<br>") else ""
          hover_mat[i, j] <- paste0(
            "<b>", rownames(lfc_mat)[i], "</b><br>",
            "Comparison: ", colnames(lfc_mat)[j], "<br>",
            "log2FC: ", round(lfc_mat[i, j], 3), "<br>",
            pval_text,
            "FDR: ", formatC(pval_mat[i, j], format = "e", digits = 2)
          )
        }
      }

      max_val <- max(abs(lfc_mat), na.rm = TRUE)

      p_interactive <- heatmaply(
        lfc_mat,
        main = plot_title,
        scale = "none",
        colors = colorRampPalette(c("blue3", "white", "red3"))(100),
        limits = c(-max_val, max_val),
        dendrogram = if (cluster_rows && cluster_cols) "both" else if (cluster_rows) "row" else if (cluster_cols) "column" else "none",
        showticklabels = c(show_row_names, show_col_names),
        fontsize_row = adj_fontsize,
        fontsize_col = adj_fontsize,
        custom_hovertext = hover_mat,
        cellnote = sig_mat,
        cellnote_textposition = "middle center",
        plot_method = "plotly"
      )

      html_file <- get_interactive_path(out_file)
      saveWidget(p_interactive, html_file, selfcontained = TRUE)
      message(paste("Saved interactive heatmap to:", html_file))
    }

    return(ht)
  }

  # If signature_files provided and export_individual is TRUE, export individual plots only
  if (!is.null(signature_files) && export_individual && length(signature_files) > 0) {
    heatmap_list <- list()
    output_dir <- dirname(output_file)
    output_prefix <- tools::file_path_sans_ext(basename(output_file))

    for (sig_file in signature_files) {
      sig_name <- tools::file_path_sans_ext(basename(sig_file))
      sig_genes <- read_signature_genes(sig_file)

      # Filter to signature genes
      deg_filtered <- deg_df[deg_df[[gene_col]] %in% sig_genes, ]

      if (nrow(deg_filtered) == 0) {
        warning(paste("No genes found for signature:", sig_name))
        next
      }

      # Generate title and output path
      sig_title <- if (!is.null(title)) paste0(title, " - ", sig_name) else paste0("DEG Heatmap: ", sig_name)
      sig_output <- file.path(output_dir, paste0(output_prefix, "_", sig_name, ".png"))

      heatmap_list[[sig_name]] <- create_deg_heatmap(
        deg_filtered, sig_title, sig_output, fontsize, width, height
      )
    }

    # No combined plot - only individual exports
    invisible(heatmap_list)

  } else {
    # Single combined plot (original behavior)
    if (!is.null(signature_files)) {
      target_genes <- read_signature_genes(signature_files)
    } else if (!is.null(genes)) {
      target_genes <- genes
    } else {
      stop("Either signature_files or genes must be provided")
    }

    # Filter to target genes
    deg_filtered <- deg_df[deg_df[[gene_col]] %in% target_genes, ]

    if (nrow(deg_filtered) == 0) {
      stop("No genes found in DEG data. Check gene names match between DEG data and signature files.")
    }

    # Auto-generate title if not provided
    plot_title <- if (!is.null(title)) title else "DEG Gene Heatmap"

    ht <- create_deg_heatmap(deg_filtered, plot_title, output_file, fontsize, width, height)
    invisible(ht)
  }
}

# =============================================================================
# Function 3: Expression Heatmap (grouped by sample treatment)
# =============================================================================

#' Plot expression heatmap with samples grouped by treatment
#'
#' Creates a clustered heatmap with genes as rows and samples as columns.
#' Samples are grouped by treatment with clustering within each group.
#' Cell colors represent Z-scored expression (blue-white-red).
#' Column colorbar shows treatment group annotation.
#'
#' @param expr_data Either file path, RDS path, or matrix of expression values
#' @param metadata Either file path or data.frame of sample metadata
#' @param signature_files Character vector of signature file paths (NULL = use genes param)
#' @param genes Character vector of gene symbols (used if signature_files is NULL)
#' @param group_col Column name in metadata for sample grouping (default: "Treatment")
#' @param sample_col Column name in metadata for sample IDs (default: "Sample.ID")
#' @param output_file Output file path (PNG)
#' @param title Plot title (default: NULL, auto-generates from signature name)
#' @param export_individual Export individual plots per signature file (default: TRUE when multiple signature_files)
#' @param zscore Perform Z-score normalization by row (default: TRUE)
#' @param width Plot width in inches (default: 14)
#' @param height Plot height in inches (default: 10)
#' @param fontsize Base font size (default: 8)
#' @param cluster_rows Cluster rows (default: TRUE)
#' @param cluster_within_groups Cluster within treatment groups (default: TRUE)
#' @param show_row_names Show row names (default: TRUE)
#' @param show_col_names Show column names (default: TRUE)
#' @param row_names_max_width Max width for row names (default: 15)
#' @return Invisibly returns the Heatmap object (or list of Heatmap objects if export_individual=TRUE)
#' @export
plot_expression_heatmap <- function(
    expr_data,
    metadata,
    signature_files = NULL,
    genes = NULL,
    group_col = "Treatment",
    sample_col = "Sample.ID",
    output_file = "figures/expression_heatmap.png",
    title = NULL,
    export_individual = TRUE,
    zscore = TRUE,
    width = 14,
    height = 10,
    fontsize = 8,
    cluster_rows = TRUE,
    cluster_within_groups = TRUE,
    show_row_names = TRUE,
    show_col_names = TRUE,
    row_names_max_width = 15,
    interactive = TRUE
) {
  # Load expression data
  if (is.character(expr_data)) {
    if (grepl("\\.rds$", expr_data, ignore.case = TRUE)) {
      rds_data <- readRDS(expr_data)
      # Try to extract expression matrix from RDS
      if ("vsd" %in% names(rds_data)) {
        expr_mat_full <- SummarizedExperiment::assay(rds_data$vsd)
      } else if ("normalized_counts" %in% names(rds_data)) {
        expr_mat_full <- rds_data$normalized_counts
      } else if (is.matrix(rds_data) || is.data.frame(rds_data)) {
        expr_mat_full <- as.matrix(rds_data)
      } else {
        stop("Could not find expression matrix in RDS file")
      }
    } else {
      expr_mat_full <- as.matrix(read.delim(expr_data, row.names = 1))
    }
  } else if (inherits(expr_data, "DESeqTransform")) {
    expr_mat_full <- SummarizedExperiment::assay(expr_data)
  } else {
    expr_mat_full <- as.matrix(expr_data)
  }

  # Load metadata
  if (is.character(metadata)) {
    if (grepl("\\.rds$", metadata, ignore.case = TRUE)) {
      meta_df <- readRDS(metadata)
    } else {
      meta_df <- read.delim(metadata, stringsAsFactors = FALSE)
    }
  } else {
    meta_df <- as.data.frame(metadata)
  }

  # Helper function to create single expression heatmap
  create_expr_heatmap <- function(target_genes, plot_title, out_file, base_fontsize, base_width, base_height) {
    # Filter to target genes that exist in expression matrix
    available_genes <- intersect(target_genes, rownames(expr_mat_full))
    if (length(available_genes) == 0) {
      warning(paste("No target genes found in expression matrix for:", plot_title))
      return(NULL)
    }
    message(paste("Found", length(available_genes), "of", length(target_genes), "genes in expression data"))

    expr_mat <- expr_mat_full[available_genes, , drop = FALSE]

    # Match samples between metadata and expression
    if (sample_col %in% colnames(meta_df)) {
      sample_ids <- meta_df[[sample_col]]
    } else {
      sample_ids <- rownames(meta_df)
    }

    common_samples <- intersect(colnames(expr_mat), sample_ids)
    if (length(common_samples) == 0) {
      stop("No common samples between expression data and metadata")
    }

    # Filter and order
    expr_mat <- expr_mat[, common_samples, drop = FALSE]
    if (sample_col %in% colnames(meta_df)) {
      meta_filtered <- meta_df[match(common_samples, meta_df[[sample_col]]), ]
      rownames(meta_filtered) <- common_samples
    } else {
      meta_filtered <- meta_df[common_samples, ]
    }

    # Z-score normalize if requested
    if (zscore) {
      expr_mat <- t(scale(t(expr_mat)))
      # Replace Inf/-Inf with NA, then NA with 0
      expr_mat[!is.finite(expr_mat)] <- 0
    }

    # Create color function
    col_fun <- create_bwr_colormap(expr_mat)

    # Create group annotation colors
    groups <- factor(meta_filtered[[group_col]])
    n_groups <- length(levels(groups))
    group_colors <- setNames(
      colorRampPalette(brewer.pal(min(n_groups, 12), "Set3"))(n_groups),
      levels(groups)
    )

    # Create column annotation
    col_annotation <- HeatmapAnnotation(
      Treatment = groups,
      col = list(Treatment = group_colors),
      annotation_name_gp = gpar(fontsize = base_fontsize),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = base_fontsize),
        labels_gp = gpar(fontsize = base_fontsize - 1)
      )
    )

    # Adjust height based on number of genes
    n_genes <- nrow(expr_mat)
    adj_height <- base_height
    adj_fontsize <- base_fontsize
    if (n_genes > 50) {
      adj_height <- max(base_height, n_genes * 0.15)
      if (n_genes > 100) adj_fontsize <- max(5, base_fontsize - 2)
    }

    # Create heatmap with column split by treatment (hide column_title for splits)
    ht <- Heatmap(
      expr_mat,
      name = ifelse(zscore, "Z-score", "Expression"),
      col = col_fun,
      column_title = plot_title,
      column_title_gp = gpar(fontsize = adj_fontsize + 2, fontface = "bold"),
      top_annotation = col_annotation,
      column_split = groups,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_within_groups,
      cluster_column_slices = FALSE,
      show_row_names = show_row_names,
      show_column_names = show_col_names,
      row_names_gp = gpar(fontsize = adj_fontsize),
      column_names_gp = gpar(fontsize = adj_fontsize - 1),
      row_names_max_width = unit(row_names_max_width, "cm"),
      column_names_rot = 90,
      # Hide split group labels (keep only color bar annotation)
      show_column_dend = FALSE,
      column_title_side = "top",
      column_gap = unit(1, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = adj_fontsize + 2),
        labels_gp = gpar(fontsize = adj_fontsize)
      )
    )

    # Get output path in expr_comp_heatmap_boxplot subfolder
    static_out_file <- get_static_path(out_file)

    # Ensure output directory exists
    dir.create(dirname(static_out_file), showWarnings = FALSE, recursive = TRUE)

    # Save plot
    png(static_out_file, width = base_width, height = adj_height, units = "in", res = 300)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         column_title = NULL)  # This hides the split labels
    dev.off()

    message(paste("Saved expression heatmap to:", static_out_file))
    message(paste("Genes plotted:", nrow(expr_mat), "| Samples:", ncol(expr_mat)))

    # Generate interactive HTML version
    if (interactive) {
      # Create hover text
      hover_mat <- matrix("", nrow = nrow(expr_mat), ncol = ncol(expr_mat))
      for (i in 1:nrow(expr_mat)) {
        for (j in 1:ncol(expr_mat)) {
          sample_name <- colnames(expr_mat)[j]
          group <- meta_filtered[sample_name, group_col]
          hover_mat[i, j] <- paste0(
            "<b>", rownames(expr_mat)[i], "</b><br>",
            "Sample: ", sample_name, "<br>",
            "Group: ", group, "<br>",
            "Value: ", round(expr_mat[i, j], 3)
          )
        }
      }

      # Create column annotation
      col_annot <- data.frame(
        Group = meta_filtered[colnames(expr_mat), group_col],
        row.names = colnames(expr_mat)
      )

      # Handle NA/Inf values for color scale
      expr_mat_clean <- expr_mat
      expr_mat_clean[!is.finite(expr_mat_clean)] <- 0
      max_val <- max(abs(expr_mat_clean), na.rm = TRUE)
      if (!is.finite(max_val) || max_val == 0) max_val <- 1

      p_interactive <- heatmaply(
        expr_mat_clean,
        main = plot_title,
        scale = "none",
        colors = colorRampPalette(c("blue3", "white", "red3"))(100),
        limits = c(-max_val, max_val),
        dendrogram = "row",
        Colv = FALSE,
        showticklabels = c(show_row_names, show_col_names),
        fontsize_row = adj_fontsize,
        fontsize_col = adj_fontsize - 2,
        custom_hovertext = hover_mat,
        col_side_colors = col_annot,
        plot_method = "plotly"
      )

      html_file <- get_interactive_path(out_file)
      saveWidget(p_interactive, html_file, selfcontained = TRUE)
      message(paste("Saved interactive heatmap to:", html_file))
    }

    return(ht)
  }

  # If signature_files provided and export_individual is TRUE, export individual plots only
  if (!is.null(signature_files) && export_individual && length(signature_files) > 0) {
    heatmap_list <- list()
    output_dir <- dirname(output_file)
    output_prefix <- tools::file_path_sans_ext(basename(output_file))

    for (sig_file in signature_files) {
      sig_name <- tools::file_path_sans_ext(basename(sig_file))
      sig_genes <- read_signature_genes(sig_file)

      # Generate title and output path
      sig_title <- if (!is.null(title)) paste0(title, " - ", sig_name) else paste0("Expression: ", sig_name)
      sig_output <- file.path(output_dir, paste0(output_prefix, "_", sig_name, ".png"))

      ht <- create_expr_heatmap(sig_genes, sig_title, sig_output, fontsize, width, height)
      if (!is.null(ht)) {
        heatmap_list[[sig_name]] <- ht
      }
    }

    # No combined plot - only individual exports
    invisible(heatmap_list)

  } else {
    # Single combined plot (original behavior)
    if (!is.null(signature_files)) {
      target_genes <- read_signature_genes(signature_files)
    } else if (!is.null(genes)) {
      target_genes <- genes
    } else {
      stop("Either signature_files or genes must be provided")
    }

    # Auto-generate title if not provided
    plot_title <- if (!is.null(title)) title else "Expression Heatmap"

    ht <- create_expr_heatmap(target_genes, plot_title, output_file, fontsize, width, height)
    invisible(ht)
  }
}

# =============================================================================
# Function 4: GSVA Boxplot by Treatment Group
# =============================================================================

#' Plot GSVA/ssGSEA boxplot by treatment group
#'
#' Creates boxplots with treatment groups on x-axis and GSVA scores on y-axis.
#' For multiple signatures, creates grouped/faceted boxplots.
#' Includes adaptive statistical comparisons: Levene's test for equal variance,
#' then ANOVA + Tukey HSD (if equal) or Welch's ANOVA + Games-Howell (if unequal).
#'
#' @param gsva_data Either file path, RDS path, or matrix of GSVA scores
#' @param metadata Either file path or data.frame of sample metadata
#' @param signatures Character vector of signature names to plot
#' @param group_col Column name in metadata for sample grouping (default: "Treatment")
#' @param sample_col Column name in metadata for sample IDs (default: "Sample.ID")
#' @param output_file Output file path (PNG)
#' @param title Plot title (default: NULL, auto-generates "GSVA Scores by Treatment")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @param fontsize Base font size (default: 12)
#' @param palette Color palette name (default: "Set2")
#' @param add_points Add individual data points (default: TRUE)
#' @param rotate_x_labels Rotate x-axis labels (default: TRUE, at 90 degrees)
#' @param facet_ncol Number of columns for faceting multiple signatures (default: NULL = auto)
#' @param stat_compare Statistical comparison method: "adaptive" (Levene → Tukey/Games-Howell),
#'        "kruskal" (Kruskal-Wallis global), "wilcox_ref" (Wilcoxon vs reference),
#'        "wilcox_all" (all pairwise), or "none" (default: "adaptive")
#' @param ref_group Reference group for Wilcoxon comparisons (default: first level)
#' @param p_adjust_method P-value adjustment method (default: "BH" for Benjamini-Hochberg)
#' @param hide_ns Hide non-significant comparisons (default: TRUE)
#' @param max_pairs Maximum number of significant pairs to display on plot (default: 10)
#' @param export_stats Export full pairwise statistics to TSV file (default: TRUE)
#' @return Invisibly returns the ggplot object
#' @export
plot_gsva_boxplot <- function(
    gsva_data,
    metadata,
    signatures,
    group_col = "Treatment",
    sample_col = "Sample.ID",
    output_file = "figures/gsva_boxplot.png",
    title = NULL,
    width = 12,
    height = 8,
    fontsize = 12,
    palette = "Set2",
    add_points = TRUE,
    rotate_x_labels = TRUE,
    facet_ncol = NULL,
    stat_compare = "adaptive",
    ref_group = NULL,
    p_adjust_method = "BH",
    hide_ns = TRUE,
    max_pairs = 10,
    export_stats = TRUE,
    interactive = TRUE
) {
  # Load GSVA data
  if (is.character(gsva_data)) {
    if (grepl("\\.rds$", gsva_data, ignore.case = TRUE)) {
      rds_data <- readRDS(gsva_data)
      # Try to extract GSVA scores from RDS
      if ("gsva_scores" %in% names(rds_data)) {
        gsva_mat <- rds_data$gsva_scores
      } else if (is.matrix(rds_data) || is.data.frame(rds_data)) {
        gsva_mat <- as.matrix(rds_data)
      } else {
        stop("Could not find GSVA scores in RDS file")
      }
    } else {
      gsva_mat <- as.matrix(read.delim(gsva_data, row.names = 1))
    }
  } else {
    gsva_mat <- as.matrix(gsva_data)
  }

  # Load metadata
  if (is.character(metadata)) {
    if (grepl("\\.rds$", metadata, ignore.case = TRUE)) {
      meta_df <- readRDS(metadata)
    } else {
      meta_df <- read.delim(metadata, stringsAsFactors = FALSE)
    }
  } else {
    meta_df <- as.data.frame(metadata)
  }

  # Filter to requested signatures
  available_sigs <- intersect(signatures, rownames(gsva_mat))
  if (length(available_sigs) == 0) {
    stop("None of the requested signatures found in GSVA data")
  }
  message(paste("Found", length(available_sigs), "of", length(signatures), "signatures"))

  gsva_mat <- gsva_mat[available_sigs, , drop = FALSE]

  # Match samples between metadata and GSVA
  if (sample_col %in% colnames(meta_df)) {
    sample_ids <- meta_df[[sample_col]]
  } else {
    sample_ids <- rownames(meta_df)
  }

  common_samples <- intersect(colnames(gsva_mat), sample_ids)
  if (length(common_samples) == 0) {
    stop("No common samples between GSVA data and metadata")
  }

  # Filter and match order
  gsva_mat <- gsva_mat[, common_samples, drop = FALSE]
  if (sample_col %in% colnames(meta_df)) {
    meta_df <- meta_df[match(common_samples, meta_df[[sample_col]]), ]
  } else {
    meta_df <- meta_df[common_samples, ]
  }

  # Reshape to long format
  gsva_long <- gsva_mat %>%
    as.data.frame() %>%
    rownames_to_column("Signature") %>%
    pivot_longer(
      cols = -Signature,
      names_to = "Sample",
      values_to = "GSVA_Score"
    ) %>%
    left_join(
      meta_df %>%
        mutate(Sample = if(sample_col %in% colnames(.)) .[[sample_col]] else rownames(.)) %>%
        select(Sample, all_of(group_col)),
      by = "Sample"
    )

  # Ensure group column is a factor
  gsva_long[[group_col]] <- factor(gsva_long[[group_col]])

  # Determine plot layout based on number of signatures
  n_sigs <- length(available_sigs)
  n_treatments <- length(levels(gsva_long[[group_col]]))

  # Auto-generate title if not provided
  plot_title <- if (!is.null(title)) title else "GSVA Scores by Treatment"

  # Create extended color palette for many groups
  if (n_treatments > 8) {
    group_colors <- colorRampPalette(brewer.pal(8, palette))(n_treatments)
  } else {
    group_colors <- brewer.pal(max(3, n_treatments), palette)[1:n_treatments]
  }
  names(group_colors) <- levels(gsva_long[[group_col]])

  # Helper function for adaptive statistical testing per signature
  perform_adaptive_stats <- function(data, sig_name, grp_col) {
    sig_data <- data[data$Signature == sig_name, ]

    # Levene's test for equal variances
    levene_p <- tryCatch({
      formula_str <- paste("GSVA_Score ~", grp_col)
      levene_result <- leveneTest(as.formula(formula_str), data = sig_data)
      levene_result$`Pr(>F)`[1]
    }, error = function(e) {
      1  # Assume equal variance if test fails
    })

    if (is.na(levene_p)) levene_p <- 1

    # Perform appropriate post-hoc test
    tryCatch({
      if (levene_p > 0.05) {
        # Equal variance: ANOVA + Tukey HSD
        aov_result <- aov(as.formula(paste("GSVA_Score ~", grp_col)), data = sig_data)
        tukey_result <- TukeyHSD(aov_result)
        pairwise <- as.data.frame(tukey_result[[1]])
        pairwise$comparison <- rownames(pairwise)
        pairwise$p.adj <- pairwise$`p adj`
        method_used <- "Tukey HSD"
      } else {
        # Unequal variance: Games-Howell
        games_result <- sig_data %>%
          games_howell_test(as.formula(paste("GSVA_Score ~", grp_col)))
        pairwise <- as.data.frame(games_result)
        pairwise$comparison <- paste(pairwise$group1, pairwise$group2, sep = "-")
        method_used <- "Games-Howell"
      }

      # Add significance annotations to all pairs
      all_pairs <- pairwise %>%
        mutate(
          p.signif = case_when(
            p.adj < 0.0001 ~ "****",
            p.adj < 0.001 ~ "***",
            p.adj < 0.01 ~ "**",
            p.adj < 0.05 ~ "*",
            TRUE ~ "ns"
          ),
          significant = p.adj < 0.05,
          Signature = sig_name
        )

      # Filter to significant pairs for plotting
      sig_pairs <- all_pairs %>% filter(significant)

      return(list(
        pairs = sig_pairs,
        all_pairs = all_pairs,
        method = method_used,
        levene_p = levene_p
      ))
    }, error = function(e) {
      message(paste("Statistical test failed for", sig_name, ":", e$message))
      return(list(pairs = data.frame(), all_pairs = data.frame(), method = "failed", levene_p = levene_p))
    })
  }

  # Helper function to create single signature boxplot with stats
  # Returns list with plot and full stats for export
  create_single_boxplot <- function(sig_data, sig_name, grp_col, sig_title, show_stats = TRUE,
                                    max_pairs_display = max_pairs) {
    p <- ggplot(sig_data, aes(x = .data[[grp_col]], y = GSVA_Score, fill = .data[[grp_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      labs(
        title = sig_title,
        x = grp_col,
        y = "GSVA Score",
        fill = grp_col
      )

    # Add individual points if requested
    if (add_points) {
      p <- p + geom_jitter(width = 0.2, alpha = 0.5, size = 1)
    }

    # Initialize stats_result for return
    full_stats <- NULL
    all_stats <- NULL
    n_sig_total <- 0
    n_shown <- 0
    method_used <- ""

    # Add statistical comparisons
    if (show_stats && stat_compare == "adaptive") {
      stats_result <- perform_adaptive_stats(gsva_long, sig_name, grp_col)
      full_stats <- stats_result$pairs  # significant pairs only (for backwards compatibility)
      all_stats <- stats_result$all_pairs  # all pairs (for export)
      method_used <- stats_result$method
      n_sig_total <- nrow(stats_result$pairs)

      message(paste("Signature:", sig_name, "| Method:", stats_result$method,
                    "| Levene p:", round(stats_result$levene_p, 4),
                    "| Significant pairs:", n_sig_total))

      if (n_sig_total > 0) {
        # Sort by p.adj and get top significant pairs
        stats_result$pairs <- stats_result$pairs[order(stats_result$pairs$p.adj), ]

        # Limit to max_pairs_display most significant pairs
        n_shown <- min(max_pairs_display, n_sig_total)
        top_pairs <- stats_result$pairs[1:n_shown, ]

        # Parse comparisons and annotations for ggsignif
        comparisons_list <- list()
        annotations_vec <- c()

        for (i in 1:nrow(top_pairs)) {
          pair <- top_pairs$comparison[i]
          if (grepl("-", pair)) {
            parts <- strsplit(pair, "-")[[1]]
            if (length(parts) >= 2) {
              comparisons_list[[length(comparisons_list) + 1]] <- c(parts[1], parts[2])
              annotations_vec <- c(annotations_vec, top_pairs$p.signif[i])
            }
          }
        }

        if (length(comparisons_list) > 0) {
          # Use ggsignif with pre-computed annotations
          p <- p + geom_signif(
            comparisons = comparisons_list,
            annotations = annotations_vec,
            step_increase = 0.08,
            tip_length = 0.01,
            textsize = 3,
            vjust = 0.2
          )
        }

        # Add subtitle showing X of Y significant pairs if some are hidden
        if (n_sig_total > n_shown) {
          p <- p + labs(subtitle = paste0("Showing ", n_shown, " of ", n_sig_total,
                                          " significant pairs (", method_used, ")"))
        } else {
          p <- p + labs(subtitle = paste0(n_sig_total, " significant pairs (", method_used, ")"))
        }
      } else {
        p <- p + labs(subtitle = paste0("No significant pairs (", method_used, ")"))
      }
    } else if (show_stats && stat_compare == "kruskal") {
      p <- p + stat_compare_means(method = "kruskal.test", label.y.npc = 0.95)
    } else if (show_stats && stat_compare == "wilcox_ref") {
      ref <- if (is.null(ref_group)) levels(sig_data[[grp_col]])[1] else ref_group
      p <- p + stat_compare_means(
        method = "wilcox.test",
        ref.group = ref,
        label = "p.signif",
        hide.ns = hide_ns,
        p.adjust.method = p_adjust_method
      )
    } else if (show_stats && stat_compare == "wilcox_all") {
      comparisons <- combn(levels(sig_data[[grp_col]]), 2, simplify = FALSE)
      p <- p + stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = hide_ns,
        p.adjust.method = p_adjust_method
      )
    }

    # Apply theme and styling with manual color scale
    p <- p +
      scale_fill_manual(values = group_colors) +
      theme_bw(base_size = fontsize) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = fontsize + 2)
      )

    # Rotate x-axis labels to 90 degrees (vertical)
    if (rotate_x_labels) {
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    # Return both plot and stats (all_stats for TSV export, full_stats for significant only)
    return(list(plot = p, stats = full_stats, all_stats = all_stats, method = method_used,
                n_sig_total = n_sig_total, n_shown = n_shown))
  }

  # For adaptive stats with multiple signatures, generate individual plots
  if (stat_compare == "adaptive" && n_sigs > 1) {
    output_dir <- dirname(output_file)
    output_prefix <- tools::file_path_sans_ext(basename(output_file))
    plot_list <- list()
    all_stats <- list()

    for (sig in available_sigs) {
      sig_data <- gsva_long[gsva_long$Signature == sig, ]
      sig_title <- paste0(plot_title, ": ", sig)
      sig_output <- file.path(output_dir, paste0(output_prefix, "_", sig, ".png"))

      result <- create_single_boxplot(sig_data, sig, group_col, sig_title, show_stats = TRUE)
      p <- result$plot

      # Adjust dimensions
      adj_width <- if (n_treatments > 6) max(width, n_treatments * 1.2) else width

      # Get output path in expr_comp_heatmap_boxplot subfolder
      static_sig_output <- get_static_path(sig_output)

      # Ensure output directory exists
      dir.create(dirname(static_sig_output), showWarnings = FALSE, recursive = TRUE)

      # Save individual plot
      ggsave(static_sig_output, p, width = adj_width, height = height, dpi = 300)
      message(paste("Saved GSVA boxplot to:", static_sig_output))

      # Generate interactive HTML version
      if (interactive) {
        p_interactive <- ggplotly(p, tooltip = c("x", "y", "fill"))
        html_file <- get_interactive_path(sig_output)
        saveWidget(p_interactive, html_file, selfcontained = TRUE)
        message(paste("Saved interactive boxplot to:", html_file))
      }

      # Export full pairwise statistics to TSV if requested (all pairs, not just significant)
      if (export_stats && !is.null(result$all_stats) && nrow(result$all_stats) > 0) {
        stats_output <- get_stats_path(sig_output)
        write.table(result$all_stats, stats_output, sep = "\t", row.names = FALSE, quote = FALSE)
        message(paste("Saved pairwise statistics to:", stats_output,
                      "(", sum(result$all_stats$significant), "of", nrow(result$all_stats), "significant)"))
      }

      plot_list[[sig]] <- p
      all_stats[[sig]] <- result$all_stats
    }

    message(paste("Generated", length(available_sigs), "individual GSVA boxplots"))
    message(paste("Statistical comparison:", stat_compare))
    invisible(list(plots = plot_list, stats = all_stats))

  } else if (n_sigs == 1) {
    # Single signature
    sig_data <- gsva_long
    sig_title <- paste0(plot_title, ": ", available_sigs[1])

    result <- create_single_boxplot(sig_data, available_sigs[1], group_col, sig_title, show_stats = TRUE)
    p <- result$plot

    # Adjust dimensions
    if (n_treatments > 6) {
      width <- max(width, n_treatments * 1.2)
    }

    # Get output path in expr_comp_heatmap_boxplot subfolder
    static_output <- get_static_path(output_file)

    # Ensure output directory exists
    dir.create(dirname(static_output), showWarnings = FALSE, recursive = TRUE)

    # Save plot
    ggsave(static_output, p, width = width, height = height, dpi = 300)
    message(paste("Saved GSVA boxplot to:", static_output))

    # Generate interactive HTML version
    if (interactive) {
      p_interactive <- ggplotly(p, tooltip = c("x", "y", "fill"))
      html_file <- get_interactive_path(output_file)
      saveWidget(p_interactive, html_file, selfcontained = TRUE)
      message(paste("Saved interactive boxplot to:", html_file))
    }

    # Export full pairwise statistics to TSV if requested (all pairs, not just significant)
    if (export_stats && !is.null(result$all_stats) && nrow(result$all_stats) > 0) {
      stats_output <- get_stats_path(output_file)
      write.table(result$all_stats, stats_output, sep = "\t", row.names = FALSE, quote = FALSE)
      message(paste("Saved pairwise statistics to:", stats_output,
                    "(", sum(result$all_stats$significant), "of", nrow(result$all_stats), "significant)"))
    }

    message(paste("Saved GSVA boxplot to:", output_file))
    message(paste("Signatures plotted:", 1, "| Treatment groups:", n_treatments))
    message(paste("Statistical comparison:", stat_compare))
    invisible(list(plot = p, stats = result$stats, all_stats = result$all_stats))

  } else {
    # Multiple signatures with non-adaptive stats - use faceted plot
    p <- ggplot(gsva_long, aes(x = .data[[group_col]], y = GSVA_Score, fill = .data[[group_col]])) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      facet_wrap(~ Signature, ncol = facet_ncol, scales = "free_y") +
      labs(
        title = plot_title,
        x = group_col,
        y = "GSVA Score",
        fill = group_col
      )

    # Add individual points if requested
    if (add_points) {
      p <- p + geom_jitter(width = 0.2, alpha = 0.5, size = 1)
    }

    # Add statistical comparisons for faceted plots
    if (stat_compare == "kruskal") {
      p <- p + stat_compare_means(method = "kruskal.test", label.y.npc = 0.95, size = 3)
    } else if (stat_compare == "wilcox_ref") {
      ref <- if (is.null(ref_group)) levels(gsva_long[[group_col]])[1] else ref_group
      p <- p + stat_compare_means(
        method = "wilcox.test",
        ref.group = ref,
        label = "p.signif",
        hide.ns = hide_ns,
        p.adjust.method = p_adjust_method,
        size = 3
      )
    }

    # Apply theme and styling with manual color scale
    p <- p +
      scale_fill_manual(values = group_colors) +
      theme_bw(base_size = fontsize) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = fontsize - 1),
        plot.title = element_text(face = "bold", size = fontsize + 2)
      )

    # Rotate x-axis labels to 90 degrees (vertical)
    if (rotate_x_labels) {
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    # Adjust dimensions for many treatments
    if (n_treatments > 6) {
      width <- max(width, n_treatments * 1.2)
    }
    if (n_sigs > 3) {
      height <- max(height, ceiling(n_sigs / (facet_ncol %||% 3)) * 4)
    }

    # Get output path in expr_comp_heatmap_boxplot subfolder
    static_output <- get_static_path(output_file)

    # Ensure output directory exists
    dir.create(dirname(static_output), showWarnings = FALSE, recursive = TRUE)

    # Save plot
    ggsave(static_output, p, width = width, height = height, dpi = 300)

    message(paste("Saved GSVA boxplot to:", static_output))
    message(paste("Signatures plotted:", length(available_sigs), "| Treatment groups:", n_treatments))
    message(paste("Statistical comparison:", stat_compare))

    # Generate interactive HTML version
    if (interactive) {
      p_interactive <- ggplotly(p, tooltip = c("x", "y", "fill"))
      html_file <- get_interactive_path(output_file)
      saveWidget(p_interactive, html_file, selfcontained = TRUE)
      message(paste("Saved interactive boxplot to:", html_file))
    }

    invisible(p)
  }
}

# =============================================================================
# Null-coalescing operator (if not available)
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# Print usage information when sourced
# =============================================================================
message("RNA-seq Visualization Functions loaded successfully!")
message("Available functions:")
message("  1. plot_signature_heatmap() - GSEA/GSVA NES heatmap by comparison")
message("  2. plot_deg_heatmap() - DEG log2FC heatmap by comparison")
message("  3. plot_expression_heatmap() - Expression heatmap grouped by treatment")
message("  4. plot_gsva_boxplot() - GSVA score boxplot by treatment group")
message("")
message("Use ?function_name or help(function_name) for detailed usage.")
