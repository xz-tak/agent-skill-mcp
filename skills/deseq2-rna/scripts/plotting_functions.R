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

# Set global random seed for reproducibility
set.seed(42)

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

# =============================================================================
# Internal functions for rerunning GSEA/GSVA (matching DEG pipeline methods)
# =============================================================================

#' Run GSEA for plotting (mirrors run_gsea_custom from deseq2_analysis.R)
#' @param deg_df Data frame with DEG results (must have padj, log2FoldChange, symbol columns)
#' @param gene_sets Named list of gene sets (character vectors of gene symbols)
#' @param pval_cutoff P-value cutoff for GSEA (default: 1 to include all results)
#' @return Data frame with GSEA results
run_gsea_for_plot <- function(deg_df, gene_sets, pval_cutoff = 1) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package required for rerunning GSEA")
  }

  # Create ranked gene list (same method as DEG pipeline)
  res_clean <- deg_df[!is.na(deg_df$padj) & !is.na(deg_df$log2FoldChange), ]

  # Handle symbol column (might be 'symbol' or 'gene')
  symbol_col <- if ("symbol" %in% colnames(res_clean)) "symbol" else "gene"

  gene_list <- -log10(res_clean$padj + 1e-300) * sign(res_clean$log2FoldChange)
  names(gene_list) <- res_clean[[symbol_col]]
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # Create TERM2GENE data frame
  term2gene <- do.call(rbind, lapply(names(gene_sets), function(term) {
    data.frame(term = term, gene = gene_sets[[term]], stringsAsFactors = FALSE)
  }))

  tryCatch({
    gsea_res <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = pval_cutoff,
      minGSSize = 3,
      maxGSSize = Inf,
      verbose = FALSE
    )

    if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
      warning("No GSEA results returned")
      return(NULL)
    }

    result_df <- gsea_res@result
    result_df$source <- "CUSTOM"  # Mark as custom
    return(result_df)
  }, error = function(e) {
    warning(paste("GSEA failed:", e$message))
    return(NULL)
  })
}

#' Run GSVA for plotting (mirrors run_gsva from deseq2_analysis.R)
#' @param counts_mat Expression/counts matrix (genes as rows, samples as columns)
#' @param gene_sets Named list of gene sets (character vectors of gene symbols)
#' @param kcdf Kernel cumulative distribution function: "auto" (default), "Poisson" (integer counts), or "Gaussian" (normalized)
#' @return Matrix of GSVA scores (signatures as rows, samples as columns)
run_gsva_for_plot <- function(counts_mat, gene_sets, kcdf = "auto") {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package required for rerunning GSVA")
  }

  # Auto-detect kcdf if not specified
  if (kcdf == "auto") {
    # Check if matrix contains only integers
    is_integer_counts <- all(counts_mat == floor(counts_mat), na.rm = TRUE)
    kcdf <- if (is_integer_counts) "Poisson" else "Gaussian"
    message(paste("Auto-detected kcdf:", kcdf,
                  "(", if(is_integer_counts) "integer counts" else "normalized values", ")"))
  }

  tryCatch({
    # Use GSVA 2.0+ API
    gsva_param <- GSVA::gsvaParam(
      exprData = counts_mat,
      geneSets = gene_sets,
      kcdf = kcdf,
      maxDiff = TRUE
    )
    gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)
    return(gsva_scores)
  }, error = function(e) {
    warning(paste("GSVA failed:", e$message))
    return(NULL)
  })
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
#' @param gsea_data Either a file path to GSEA results table, data.frame/RDS object, or NULL if rerun_gsea=TRUE
#' @param signatures Character vector of signature names to include (NULL = all)
#' @param comparisons Character vector of comparison names to include, in desired order (NULL = all)
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
#' @param interactive Generate interactive HTML version (default: TRUE)
#' @param rerun_gsea If TRUE, run GSEA instead of using existing data (default: FALSE)
#' @param deg_data Required if rerun_gsea=TRUE: file path to DEG results (summstats file)
#' @param gene_sets Required if rerun_gsea=TRUE: named list of gene sets (character vectors)
#' @param pval_cutoff P-value cutoff for GSEA when rerunning (default: 1 to include all results)
#' @param export_table Export GSEA results table with leading edge genes (default: TRUE)
#' @return Invisibly returns the Heatmap object
#' @export
plot_signature_heatmap <- function(
    gsea_data = NULL,
    signatures = NULL,
    comparisons = NULL,
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
    interactive = TRUE,
    rerun_gsea = FALSE,
    deg_data = NULL,
    gene_sets = NULL,
    pval_cutoff = 1,
    export_table = TRUE
) {
  # Handle rerun_gsea case
  if (rerun_gsea) {
    # Save user's requested comparisons before overwriting
    requested_comparisons <- comparisons

    if (is.null(deg_data)) {
      stop("deg_data is required when rerun_gsea=TRUE")
    }
    if (is.null(gene_sets) || length(gene_sets) == 0) {
      stop("gene_sets is required when rerun_gsea=TRUE (named list of character vectors)")
    }

    message("Running GSEA with provided gene sets...")

    # Load DEG data
    if (is.character(deg_data)) {
      if (grepl("\\.rds$", deg_data, ignore.case = TRUE)) {
        deg_df <- readRDS(deg_data)
      } else {
        deg_df <- read.delim(deg_data, stringsAsFactors = FALSE)
      }
    } else {
      deg_df <- as.data.frame(deg_data)
    }

    # Run GSEA for each comparison if comparison column exists
    if (comp_col %in% colnames(deg_df)) {
      all_comparisons <- unique(deg_df[[comp_col]])
      all_results <- list()

      for (comp in all_comparisons) {
        comp_deg <- deg_df[deg_df[[comp_col]] == comp, ]
        comp_gsea <- run_gsea_for_plot(comp_deg, gene_sets, pval_cutoff)
        if (!is.null(comp_gsea)) {
          comp_gsea[[comp_col]] <- comp
          all_results[[comp]] <- comp_gsea
        }
      }

      if (length(all_results) == 0) {
        stop("No GSEA results returned for any comparison")
      }
      gsea_df <- do.call(rbind, all_results)
    } else {
      # Single comparison
      gsea_df <- run_gsea_for_plot(deg_df, gene_sets, pval_cutoff)
      if (is.null(gsea_df)) {
        stop("No GSEA results returned")
      }
      gsea_df[[comp_col]] <- "comparison"
    }

    message(paste("GSEA completed:", nrow(gsea_df), "results"))

  } else {
    # Original behavior: load existing data
    if (is.null(gsea_data)) {
      stop("gsea_data is required when rerun_gsea=FALSE")
    }

    if (is.character(gsea_data)) {
      if (grepl("\\.rds$", gsea_data, ignore.case = TRUE)) {
        gsea_df <- readRDS(gsea_data)
      } else {
        gsea_df <- read.delim(gsea_data, stringsAsFactors = FALSE)
      }
    } else {
      gsea_df <- as.data.frame(gsea_data)
    }
  }

  # Filter by source if specified
  if (!is.null(source_filter) && source_col %in% colnames(gsea_df)) {
    gsea_df <- gsea_df[gsea_df[[source_col]] == source_filter, ]
  }

  # Filter by signatures if specified
  if (!is.null(signatures)) {
    gsea_df <- gsea_df[gsea_df[[sig_col]] %in% signatures, ]
  }

  # Filter and order comparisons if specified (use saved value if rerun_gsea)
  filter_comparisons <- if (exists("requested_comparisons")) requested_comparisons else comparisons
  if (!is.null(filter_comparisons)) {
    gsea_df <- gsea_df[gsea_df[[comp_col]] %in% filter_comparisons, ]
    # Reorder to match specified order
    gsea_df[[comp_col]] <- factor(gsea_df[[comp_col]], levels = filter_comparisons)
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
  png(static_output, width = width, height = height, units = "in", res = 720)
  draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 0.5, 0.5, 0.5), "cm"))
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

  # Export GSEA results table with leading edge genes
  if (export_table) {
    table_output <- get_stats_path(output_file, suffix = "_gsea_results")

    # Build list of columns to export
    export_cols <- c(sig_col, comp_col, nes_col, pval_col)

    # Add optional columns if they exist
    optional_cols <- c("pvalue", "setSize", "enrichmentScore", "leading_edge", "core_enrichment")
    for (col in optional_cols) {
      if (col %in% colnames(gsea_df)) {
        export_cols <- c(export_cols, col)
      }
    }

    export_df <- gsea_df[, intersect(export_cols, colnames(gsea_df)), drop = FALSE]

    # Convert core_enrichment from /-separated to semicolon-separated for easier parsing
    if ("core_enrichment" %in% colnames(export_df)) {
      export_df$core_enrichment <- gsub("/", ";", export_df$core_enrichment)
    }

    write.table(export_df, table_output, sep = "\t", row.names = FALSE, quote = FALSE)
    message(paste("Saved GSEA results table to:", table_output))
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
    png(static_out_file, width = base_width, height = adj_height, units = "in", res = 720)
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
    png(static_out_file, width = base_width, height = adj_height, units = "in", res = 720)
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
#' @param gsva_data Either file path, RDS path, matrix of GSVA scores, or NULL if rerun_gsva=TRUE
#' @param metadata Either file path or data.frame of sample metadata
#' @param signatures Character vector of signature names to plot (used as names for gene_sets if rerun_gsva=TRUE)
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
#' @param interactive Generate interactive HTML version (default: TRUE)
#' @param rerun_gsva If TRUE, run GSVA instead of using existing data (default: FALSE)
#' @param counts_data Required if rerun_gsva=TRUE: expression/counts matrix file path or matrix
#' @param gene_sets Required if rerun_gsva=TRUE: named list of gene sets (character vectors)
#' @param kcdf Kernel cumulative distribution function for GSVA: "auto" (default), "Poisson" (integer counts), or "Gaussian" (normalized)
#' @return Invisibly returns the ggplot object
#' @export
plot_gsva_boxplot <- function(
    gsva_data = NULL,
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
    interactive = TRUE,
    rerun_gsva = FALSE,
    counts_data = NULL,
    gene_sets = NULL,
    kcdf = "auto"
) {
  # Handle rerun_gsva case
  if (rerun_gsva) {
    if (is.null(counts_data)) {
      stop("counts_data is required when rerun_gsva=TRUE")
    }
    if (is.null(gene_sets) || length(gene_sets) == 0) {
      stop("gene_sets is required when rerun_gsva=TRUE (named list of character vectors)")
    }

    message("Running GSVA with provided gene sets...")

    # Load counts data
    if (is.character(counts_data)) {
      if (grepl("\\.rds$", counts_data, ignore.case = TRUE)) {
        rds_data <- readRDS(counts_data)
        # Try to extract counts matrix from RDS
        if ("counts" %in% names(rds_data)) {
          counts_mat <- rds_data$counts
        } else if ("normalized_counts" %in% names(rds_data)) {
          counts_mat <- rds_data$normalized_counts
        } else if ("dds" %in% names(rds_data)) {
          # Extract counts from DESeqDataSet
          counts_mat <- DESeq2::counts(rds_data$dds, normalized = FALSE)
          mode(counts_mat) <- "integer"
          message("Extracted raw counts from DESeqDataSet")
        } else if ("vsd" %in% names(rds_data)) {
          # Use VST-transformed data (normalized)
          counts_mat <- SummarizedExperiment::assay(rds_data$vsd)
          message("Using VST-transformed data (normalized)")
        } else if (is.matrix(rds_data) || is.data.frame(rds_data)) {
          counts_mat <- as.matrix(rds_data)
        } else {
          stop("Could not find counts matrix in RDS file. Expected: counts, normalized_counts, dds, or vsd")
        }
      } else {
        counts_mat <- as.matrix(read.delim(counts_data, row.names = 1))
      }
    } else {
      counts_mat <- as.matrix(counts_data)
    }

    # Run GSVA
    gsva_mat <- run_gsva_for_plot(counts_mat, gene_sets, kcdf)
    if (is.null(gsva_mat)) {
      stop("GSVA failed to return results")
    }

    message(paste("GSVA completed:", nrow(gsva_mat), "signatures x", ncol(gsva_mat), "samples"))

    # Update signatures to use gene_sets names (what we actually computed)
    signatures <- names(gene_sets)

  } else {
    # Original behavior: load existing data
    if (is.null(gsva_data)) {
      stop("gsva_data is required when rerun_gsva=FALSE")
    }

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
      ggsave(static_sig_output, p, width = adj_width, height = height, dpi = 720)
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
    ggsave(static_output, p, width = width, height = height, dpi = 720)
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
    ggsave(static_output, p, width = width, height = height, dpi = 720)

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
# Function 5: GSVA + Limma Differential Expression Analysis
# =============================================================================

#' Parse GMT file to named list of gene vectors
#' @param gmt_file Path to GMT file
#' @return Named list where names are set names and values are gene vectors
parse_gmt_file <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop(paste("GMT file not found:", gmt_file))
  }

  gmt_lines <- readLines(gmt_file)
  gene_sets <- list()

  for (line in gmt_lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next

    set_name <- fields[1]
    # Skip description (field 2), genes start at field 3
    genes <- fields[3:length(fields)]
    genes <- genes[genes != ""]
    gene_sets[[set_name]] <- genes
  }

  gene_sets
}

#' Run GSVA + Limma Differential Expression Analysis
#'
#' Computes GSVA scores for custom gene signatures, runs limma-based differential
#' expression across treatment comparisons, and generates volcano plots and summary heatmap.
#'
#' @param counts_file Path to raw counts file (CSV/TSV with Gene.name column or gene symbols as rownames)
#' @param metadata_file Path to metadata file (TSV with sample IDs as rownames)
#' @param gmt_file Path to GMT file with gene signatures (or named list of gene vectors)
#' @param comparisons List of comparisons, each with name/numerator/denominator elements
#' @param group_col Column in metadata for treatment groups (default "Treatment_std")
#' @param output_dir Base output directory (default ".")
#' @param padj_cutoff Adjusted p-value threshold for significance markers (default 0.05)
#' @param pval_cutoff Nominal p-value threshold for volcano/heatmap significance (default 0.05, used for asterisks)
#' @param cluster_rows Cluster rows in heatmap (default FALSE for reproducible ordering)
#' @param cluster_columns Cluster columns in heatmap (default FALSE for reproducible ordering)
#' @param interactive Generate interactive HTML plots (default TRUE)
#' @param width Heatmap width in inches (default 14)
#' @param height Heatmap height in inches (default 10)
#' @param volcano_width Volcano plot width in inches (default 10)
#' @param volcano_height Volcano plot height in inches (default 8)
#' @return List with gsva_matrix, limma_results, output_files
#' @export
run_gsva_limma_de <- function(
    counts_file,
    metadata_file,
    gmt_file,
    comparisons,
    group_col = "Treatment_std",
    output_dir = ".",
    padj_cutoff = 0.05,
    pval_cutoff = 0.05,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_order = NULL,
    comparisons_filter = NULL,
    row_order = NULL,
    column_labels = NULL,
    limma_results_file = NULL,
    interactive = TRUE,
    width = 14,
    height = 10,
    volcano_width = 10,
    volcano_height = 8
) {
  # Load additional required packages
  suppressPackageStartupMessages({
    library(GSVA)
    library(limma)
    library(ggrepel)
  })

  # Create output directories
  gsva_deg_dir <- file.path(output_dir, "deg/gsva")
  gsva_fig_dir <- file.path(output_dir, "figures/gsva")
  dir.create(gsva_deg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gsva_fig_dir, recursive = TRUE, showWarnings = FALSE)

  # ============================================================
  # REPLOT MODE: Load existing limma results and skip to heatmap
  # ============================================================
  if (!is.null(limma_results_file)) {
    message("========================================")
    message("REPLOT MODE: Loading existing limma results")
    message("========================================\n")

    combined_results <- read.delim(limma_results_file, stringsAsFactors = FALSE)
    message(sprintf("Loaded %d rows from: %s", nrow(combined_results), limma_results_file))
    message(sprintf("Signatures: %d", length(unique(combined_results$signature))))
    message(sprintf("Comparisons: %d", length(unique(combined_results$comparison))))

    # Jump directly to heatmap generation
    # (combined_results is already available, skip GSVA and limma)
  } else {
    message("========================================")
    message("GSVA + Limma Differential Expression Analysis")
    message("========================================\n")

  # ============================================================
  # 1. Load gene signatures (GMT file or list)
  # ============================================================
  message("Loading gene signatures...")

  if (is.list(gmt_file) && !is.character(gmt_file)) {
    # Already a named list of gene vectors
    gene_sets <- gmt_file
  } else if (is.character(gmt_file)) {
    gene_sets <- parse_gmt_file(gmt_file)
  } else {
    stop("gmt_file must be a file path or named list of gene vectors")
  }

  message(sprintf("Loaded %d gene signatures:", length(gene_sets)))
  for (name in names(gene_sets)) {
    message(sprintf("  - %s: %d genes", name, length(gene_sets[[name]])))
  }

  # ============================================================
  # 2. Load raw counts and metadata
  # ============================================================
  message("\nLoading raw counts...")

  # Detect delimiter by reading first line
  first_line <- readLines(counts_file, n = 1)
  if (grepl(",", first_line)) {
    counts <- read.csv(counts_file, row.names = NULL, check.names = FALSE)
  } else {
    counts <- read.delim(counts_file, row.names = NULL, check.names = FALSE)
  }
  message(sprintf("Raw counts dimensions: %d genes x %d columns", nrow(counts), ncol(counts)))

  # Handle Gene.name column if present (may be any position)
  if ("Gene.name" %in% colnames(counts)) {
    # Handle NA values
    if ("ID" %in% colnames(counts)) {
      counts$Gene.name[is.na(counts$Gene.name)] <- counts$ID[is.na(counts$Gene.name)]
    }
    rownames(counts) <- make.unique(counts$Gene.name)
    counts <- counts[, !colnames(counts) %in% c("ID", "Gene.name")]
  } else if ("ID" %in% colnames(counts)) {
    # Use ID column as rownames if no Gene.name
    rownames(counts) <- make.unique(counts$ID)
    counts <- counts[, colnames(counts) != "ID"]
  }

  counts <- as.matrix(counts)
  message(sprintf("Processed counts dimensions: %d genes x %d samples", nrow(counts), ncol(counts)))

  # Load metadata
  message("\nLoading metadata...")
  metadata <- read.delim(metadata_file, row.names = 1, stringsAsFactors = FALSE)
  message(sprintf("Metadata: %d samples", nrow(metadata)))

  # Verify group column exists
  if (!group_col %in% colnames(metadata)) {
    stop(sprintf("Group column '%s' not found in metadata. Available columns: %s",
                 group_col, paste(colnames(metadata), collapse = ", ")))
  }

  # Match samples between counts and metadata
  common_samples <- intersect(colnames(counts), rownames(metadata))
  if (length(common_samples) == 0) {
    stop("No matching samples between counts columns and metadata rows")
  }
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, , drop = FALSE]
  message(sprintf("Matched samples: %d", length(common_samples)))

  # Check gene overlap
  all_sig_genes <- unique(unlist(gene_sets))
  overlap <- sum(all_sig_genes %in% rownames(counts))
  message(sprintf("\nGene overlap: %d/%d signature genes found in counts (%.1f%%)",
                  overlap, length(all_sig_genes), 100*overlap/length(all_sig_genes)))

  # ============================================================
  # 3. Compute GSVA (v2.0 API with Poisson kernel)
  # ============================================================
  message("\n========================================")
  message("Computing GSVA scores...")
  message("========================================")
  message("Parameters: kcdf='Poisson' (for count data)")

  gsva_param <- gsvaParam(
    exprData = counts,
    geneSets = gene_sets,
    kcdf = "Poisson"
  )

  gsva_matrix <- gsva(gsva_param, verbose = TRUE)
  message(sprintf("\nGSVA results dimensions: %d signatures x %d samples",
                  nrow(gsva_matrix), ncol(gsva_matrix)))

  # Save GSVA matrix
  gsva_file <- file.path(gsva_deg_dir, "gsva_scores.txt")
  write.table(gsva_matrix, gsva_file, sep = "\t", quote = FALSE, row.names = TRUE)
  message(sprintf("GSVA matrix saved to: %s", gsva_file))

  # ============================================================
  # 4. Validate comparisons
  # ============================================================
  message(sprintf("\nValidating %d comparisons...", length(comparisons)))

  # Make Treatment column a factor
  metadata[[group_col]] <- factor(metadata[[group_col]])
  available_groups <- levels(metadata[[group_col]])

  # Validate each comparison
  valid_comparisons <- list()
  for (comp in comparisons) {
    if (!all(c("name", "numerator", "denominator") %in% names(comp))) {
      warning(sprintf("Comparison missing required fields (name, numerator, denominator): %s",
                      paste(names(comp), collapse = ", ")))
      next
    }

    if (!comp$numerator %in% available_groups) {
      warning(sprintf("Numerator '%s' not found in %s. Available: %s",
                      comp$numerator, group_col, paste(available_groups, collapse = ", ")))
      next
    }

    if (!comp$denominator %in% available_groups) {
      warning(sprintf("Denominator '%s' not found in %s. Available: %s",
                      comp$denominator, group_col, paste(available_groups, collapse = ", ")))
      next
    }

    valid_comparisons <- append(valid_comparisons, list(comp))
  }

  if (length(valid_comparisons) == 0) {
    stop("No valid comparisons found")
  }
  message(sprintf("Valid comparisons: %d", length(valid_comparisons)))

  # ============================================================
  # 5. Run limma DE for each comparison
  # ============================================================
  message("\n========================================")
  message("Running limma differential expression...")
  message("========================================")

  # Cell-means model: ~0 + group_col
  design <- model.matrix(as.formula(paste0("~0 + ", group_col)), data = metadata)

  # Get original column names (with special characters)
  original_colnames <- gsub(group_col, "", colnames(design))

  # Convert to valid R names for makeContrasts
  valid_colnames <- make.names(original_colnames)
  colnames(design) <- valid_colnames

  # Create mapping from original names to valid names
  name_map <- setNames(valid_colnames, original_colnames)

  message("\nDesign matrix groups (original -> valid R names):")
  for (i in seq_along(original_colnames)) {
    message(sprintf("  %s -> %s", original_colnames[i], valid_colnames[i]))
  }

  # Fit initial model on GSVA matrix
  fit <- lmFit(gsva_matrix, design)

  # Store all limma results
  all_results <- list()

  for (comp in valid_comparisons) {
    message(sprintf("\nProcessing: %s", comp$name))

    # Convert comparison names to valid R names
    num_valid <- name_map[comp$numerator]
    den_valid <- name_map[comp$denominator]

    contrast_string <- paste0(num_valid, " - ", den_valid)
    message(sprintf("  Contrast: %s", contrast_string))

    contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)

    # Fit contrast and apply eBayes
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    # Extract results
    res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    res$signature <- rownames(res)
    res$comparison <- comp$name

    # Save individual limma results
    res_file <- file.path(gsva_deg_dir, paste0(comp$name, "_limma.txt"))
    write.table(res, res_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("  Results saved: %s", basename(res_file)))
    message(sprintf("  Significant (padj < %.2f): %d signatures", padj_cutoff, sum(res$adj.P.Val < padj_cutoff)))

    all_results[[comp$name]] <- res
  }

  # ============================================================
  # 6. Combine all results
  # ============================================================
  combined_results <- bind_rows(all_results)
  combined_file <- file.path(gsva_deg_dir, "all_limma_combined.txt")
  write.table(combined_results, combined_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("\nCombined results saved: %s", combined_file))
  message(sprintf("Total rows: %d (%d signatures x %d comparisons)",
                  nrow(combined_results), nrow(gsva_matrix), length(all_results)))

  # ============================================================
  # 7. Generate volcano plots (using nominal P.Value for y-axis)
  # ============================================================
  message("\n========================================")
  message("Generating volcano plots...")
  message("========================================")

  volcano_colors <- c("down" = "#3366CC", "up" = "#CC3333", "NS" = "grey50")
  volcano_files <- character()

  for (comp_name in names(all_results)) {
    res <- all_results[[comp_name]]

    # Classify significance using nominal p-value (better for small signature panels)
    res$significance <- ifelse(res$P.Value >= pval_cutoff, "NS",
                               ifelse(res$logFC > 0, "up", "down"))

    # Label only significant signatures
    res$delabel <- ifelse(res$P.Value < pval_cutoff, res$signature, NA)

    # Create volcano plot
    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value), color = significance)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = volcano_colors,
                         labels = c("down" = "Down", "up" = "Up", "NS" = "Not Sig.")) +
      geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey60") +
      geom_text_repel(aes(label = delabel),
                      size = 3,
                      max.overlaps = 20,
                      box.padding = 0.5,
                      segment.color = "grey50",
                      na.rm = TRUE) +
      labs(title = comp_name,
           subtitle = sprintf("p < %.2f: %d signatures", pval_cutoff, sum(res$P.Value < pval_cutoff)),
           x = "Log2 Fold Change (GSVA score)",
           y = "-log10(p-value)",
           color = "Regulation") +
      theme_bw() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = "right",
            panel.grid.minor = element_blank())

    volcano_file <- file.path(gsva_fig_dir, paste0("volcano_", comp_name, ".png"))
    ggsave(volcano_file, p, width = volcano_width, height = volcano_height, dpi = 150)
    message(sprintf("  Saved: %s", basename(volcano_file)))
    volcano_files <- c(volcano_files, volcano_file)

    # Generate interactive version
    if (interactive) {
      interactive_dir <- file.path(gsva_fig_dir, "interactive")
      dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

      p_interactive <- ggplotly(p, tooltip = c("x", "y", "color"))
      html_file <- file.path(interactive_dir, paste0("volcano_", comp_name, "_interactive.html"))
      saveWidget(p_interactive, html_file, selfcontained = TRUE)
      message(sprintf("  Saved interactive: %s", basename(html_file)))
    }
  }
  } # End of else block (full pipeline mode)

  # ============================================================
  # 8. Generate summary heatmap (using nominal P.Value for asterisks)
  # ============================================================
  message("\n========================================")
  message("Generating summary heatmap...")
  message("========================================")

  # Filter and order comparisons for heatmap
  combined_for_heatmap <- combined_results

  # Filter comparisons if specified
  if (!is.null(comparisons_filter)) {
    combined_for_heatmap <- combined_for_heatmap %>%
      filter(comparison %in% comparisons_filter)
    message(sprintf("  Filtered to %d comparisons for heatmap", length(unique(combined_for_heatmap$comparison))))
  }

  # Apply custom column order if specified
  if (!is.null(column_order)) {
    # Validate column_order contains valid comparisons
    valid_cols <- intersect(column_order, unique(combined_for_heatmap$comparison))
    if (length(valid_cols) == 0) {
      warning("No valid comparisons in column_order, using default order")
    } else {
      combined_for_heatmap$comparison <- factor(
        combined_for_heatmap$comparison,
        levels = valid_cols
      )
      combined_for_heatmap <- combined_for_heatmap %>%
        filter(!is.na(comparison)) %>%
        arrange(comparison)
      message(sprintf("  Applied custom column order: %d comparisons", length(valid_cols)))
    }
  }

  # Apply custom row order if specified
  if (!is.null(row_order)) {
    valid_rows <- intersect(row_order, unique(combined_for_heatmap$signature))
    if (length(valid_rows) == 0) {
      warning("No valid signatures in row_order, using default order")
    } else {
      combined_for_heatmap$signature <- factor(
        combined_for_heatmap$signature,
        levels = valid_rows
      )
      combined_for_heatmap <- combined_for_heatmap %>%
        filter(!is.na(signature)) %>%
        arrange(signature)
      message(sprintf("  Applied custom row order: %d signatures", length(valid_rows)))
    }
  }

  # Pivot to wide format for logFC
  logfc_matrix <- combined_for_heatmap %>%
    select(signature, comparison, logFC) %>%
    pivot_wider(names_from = comparison, values_from = logFC) %>%
    column_to_rownames("signature") %>%
    as.matrix()

  # Pivot to wide format for nominal p-value (better for small signature panels)
  pval_matrix <- combined_for_heatmap %>%
    select(signature, comparison, P.Value) %>%
    pivot_wider(names_from = comparison, values_from = P.Value) %>%
    column_to_rownames("signature") %>%
    as.matrix()

  # Apply asterisks using existing pval_to_stars helper
  asterisk_matrix <- matrix(sapply(pval_matrix, pval_to_stars),
                            nrow = nrow(pval_matrix),
                            ncol = ncol(pval_matrix),
                            dimnames = dimnames(pval_matrix))

  # Determine color scale limits (centered at 0)
  max_abs <- max(abs(logfc_matrix), na.rm = TRUE)
  col_limit <- ceiling(max_abs * 10) / 10  # Round up to 1 decimal
  if (col_limit == 0) col_limit <- 1

  # BWR colormap centered at 0
  col_fun <- colorRamp2(c(-col_limit, 0, col_limit), c("blue", "white", "red"))

  # Apply custom column labels or use comparison names directly
  if (!is.null(column_labels)) {
    short_names <- column_labels[colnames(logfc_matrix)]
    # Fall back to original names for any missing labels
    short_names[is.na(short_names)] <- colnames(logfc_matrix)[is.na(short_names)]
    message(sprintf("  Applied custom column labels"))
  } else {
    # Default: use comparison names directly from matrix
    short_names <- colnames(logfc_matrix)
  }

  # Create heatmap
  ht <- Heatmap(logfc_matrix,
                name = "logFC",
                col = col_fun,
                cluster_rows = cluster_rows,
                cluster_columns = cluster_columns,
                show_row_dend = cluster_rows,
                show_column_dend = cluster_columns,
                row_names_side = "left",
                column_names_side = "bottom",
                column_names_rot = 45,
                column_labels = short_names,
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 8),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(asterisk_matrix[i, j], x, y,
                            gp = gpar(fontsize = 8, col = "black"))
                },
                heatmap_legend_param = list(
                  title = "logFC",
                  at = c(-col_limit, 0, col_limit),
                  labels = c(sprintf("-%.1f", col_limit), "0", sprintf("%.1f", col_limit))
                ),
                column_title = "GSVA Signature Differential Expression",
                column_title_gp = gpar(fontsize = 14, fontface = "bold"))

  # Save heatmap
  heatmap_file <- file.path(gsva_fig_dir, "signatures_heatmap.png")
  png(heatmap_file, width = width, height = height, units = "in", res = 150)
  draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
  dev.off()
  message(sprintf("Heatmap saved: %s", heatmap_file))

  # Generate interactive heatmap
  if (interactive) {
    interactive_dir <- file.path(gsva_fig_dir, "interactive")
    dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

    # Create hover text matrix
    hover_matrix <- matrix(
      paste0("Signature: ", rownames(logfc_matrix)[row(logfc_matrix)],
             "\nComparison: ", colnames(logfc_matrix)[col(logfc_matrix)],
             "\nlogFC: ", round(logfc_matrix, 3),
             "\np-value: ", format(pval_matrix, digits = 3, scientific = TRUE),
             "\nSignificance: ", asterisk_matrix),
      nrow = nrow(logfc_matrix),
      ncol = ncol(logfc_matrix),
      dimnames = dimnames(logfc_matrix)
    )

    # Create interactive heatmap with heatmaply
    p_heatmap <- heatmaply(
      logfc_matrix,
      custom_hovertext = hover_matrix,
      colors = colorRampPalette(c("blue", "white", "red"))(100),
      dendrogram = if (cluster_rows || cluster_columns) "both" else "none",
      showticklabels = c(TRUE, TRUE),
      main = "GSVA Signature Differential Expression",
      xlab = "Comparison",
      ylab = "Signature",
      cellnote = asterisk_matrix,
      cellnote_textposition = "middle center",
      cellnote_size = 10
    )

    html_heatmap_file <- file.path(interactive_dir, "signatures_heatmap_interactive.html")
    saveWidget(p_heatmap, html_heatmap_file, selfcontained = TRUE)
    message(sprintf("Interactive heatmap saved: %s", basename(html_heatmap_file)))
  }

  # ============================================================
  # Summary
  # ============================================================
  message("\n========================================")
  message("ANALYSIS COMPLETE")
  message("========================================")

  if (!is.null(limma_results_file)) {
    # Replot mode - simpler summary
    message(sprintf("Loaded from: %s", limma_results_file))
    message(sprintf("Summary heatmap: %s", heatmap_file))
    message("\nSignificance summary (p < ", pval_cutoff, "):")
    print(table(combined_results$P.Value < pval_cutoff))

    output_files <- list(heatmap = heatmap_file)
    invisible(list(
      combined_results = combined_results,
      output_files = output_files
    ))
  } else {
    # Full pipeline mode
    message(sprintf("GSVA matrix: %s", gsva_file))
    message(sprintf("Combined limma results: %s", combined_file))
    message(sprintf("Volcano plots: %s/volcano_*.png", gsva_fig_dir))
    message(sprintf("Summary heatmap: %s", heatmap_file))
    message("\nSignificance summary (p < ", pval_cutoff, "):")
    print(table(combined_results$P.Value < pval_cutoff))

    output_files <- list(
      gsva_matrix = gsva_file,
      combined_results = combined_file,
      volcano_plots = volcano_files,
      heatmap = heatmap_file
    )

    invisible(list(
      gsva_matrix = gsva_matrix,
      limma_results = all_results,
      combined_results = combined_results,
      output_files = output_files
    ))
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
message("  5. run_gsva_limma_de() - GSVA + limma DE analysis (signatures x comparisons)")
message("")
message("Use ?function_name or help(function_name) for detailed usage.")
