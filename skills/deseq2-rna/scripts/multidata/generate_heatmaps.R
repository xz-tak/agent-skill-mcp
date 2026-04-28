#!/usr/bin/env Rscript
# DEG Multi-Dataset Analysis - Heatmap Generator
#
# This script generates heatmaps from TSV tables with log2FC values and
# significance annotations.
#
# Usage:
#   Rscript generate_heatmaps.R --config config.json
#   Rscript generate_heatmaps.R --up_table UP.tsv --down_table DOWN.tsv --output_dir ./output
#
# The script supports:
# - Column grouping with labeled separators
# - Score row annotation
# - Significance star overlay
# - Customizable color scales

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
  library(jsonlite)
  library(heatmaply)
  library(htmlwidgets)
})

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Default output directory name (created in working directory)
DEFAULT_OUTPUT_DIR <- "deg_multidata_output"

# Default column label mappings (customize per analysis)
# Format: "column_name" = "Display\nLabel"
DEFAULT_COL_MAPPING <- c()

# Default column groups (customize per analysis)
# Will be populated from config or inferred from data
DEFAULT_COLUMN_GROUPS <- list()

# Default heatmap settings
DEFAULT_HEATMAP_CONFIG <- list(
  color_scale = c("blue", "white", "red"),
  row_annotation = "Score",
  width = 14,
  height = 10,
  fontsize = 9,
  star_fontsize = 6,
  column_gap = 3,            # mm
  cluster_rows = FALSE,      # Default: score-sorted (no clustering)
  clustering_method = "ward.D2",  # Method when clustering enabled
  html_fontsize_col = 12     # Column label font size for HTML
)

# Pathway mode flag (set by --pathway_mode argument)
PATHWAY_MODE <- FALSE


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

extract_lfc <- function(x) {
  #' Extract numeric log2FC from formatted string.
  #'
  #' @param x Formatted string (e.g., "-1.45***")
  #' @return Numeric log2FC value
  if (length(x) == 0 || is.null(x)) return(NA)
  if (is.na(x)) return(NA)
  x_str <- as.character(x)
  if (x_str == "NA" || x_str == "" || x_str == "0") return(0)
  # Remove trailing stars first, then trailing dot (significance marker)
  x_str <- gsub("\\*+$", "", x_str)  # Remove stars
  x_str <- gsub("\\.$", "", x_str)   # Remove trailing dot (not decimal point)
  as.numeric(x_str)
}


extract_sig <- function(x) {
  #' Extract significance stars from formatted string.
  #'
  #' @param x Formatted string (e.g., "-1.45***")
  #' @return String of stars only
  if (length(x) == 0 || is.null(x)) return("")
  if (is.na(x)) return("")
  x_str <- as.character(x)
  if (x_str == "NA" || x_str == "" || x_str == "0") return("")
  # Only keep stars, remove dots
  gsub("[^*]", "", x_str)
}


create_heatmap <- function(data_file, output_file, title, config, pathway_mode = FALSE) {
  #' Create heatmap from TSV data file.
  #'
  #' Generates three outputs with consistent row ordering:
  #' 1. PNG heatmap (ComplexHeatmap)
  #' 2. TSV table with rows in same order as PNG
  #' 3. Interactive HTML heatmap (heatmaply) with same row order
  #'
  #' @param data_file Path to TSV file with formatted values
  #' @param output_file Path for output PNG
  #' @param title Heatmap title
  #' @param config Configuration list with column_mapping, column_groups, heatmap settings
  #' @param pathway_mode If TRUE, use "Pathway" column and NES values instead of Gene/log2FC

  # Load data
  df <- read.delim(data_file, stringsAsFactors = FALSE, check.names = FALSE)

  # Determine row identifier column (Gene vs Pathway)
  if (pathway_mode || "Pathway" %in% names(df)) {
    row_col <- "Pathway"
    value_label <- "NES"
    pathway_mode <- TRUE
  } else {
    row_col <- "Gene"
    value_label <- "log2FC"
  }

  if (!row_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in data file.", row_col))
  }

  # Get configuration
  col_mapping <- config$column_mapping
  column_groups <- config$column_groups
  heatmap_config <- modifyList(DEFAULT_HEATMAP_CONFIG, config$heatmap %||% list())

  # Determine data columns (all columns that are in column_groups)
  all_group_cols <- unlist(column_groups)
  data_cols <- all_group_cols[all_group_cols %in% names(df)]

  if (length(data_cols) == 0) {
    stop("No data columns found in file. Check column_groups configuration.")
  }

  # Extract matrices
  n_rows <- nrow(df)
  n_cols <- length(data_cols)

  value_mat <- matrix(NA, nrow = n_rows, ncol = n_cols)
  sig_mat <- matrix("", nrow = n_rows, ncol = n_cols)

  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      val <- df[i, data_cols[j]]
      value_mat[i, j] <- extract_lfc(val)  # Works for both log2FC and NES
      sig_mat[i, j] <- extract_sig(val)
    }
  }

  # Set row and column names
  rownames(value_mat) <- df[[row_col]]
  if (!is.null(col_mapping) && length(col_mapping) > 0) {
    colnames(value_mat) <- sapply(data_cols, function(x) {
      if (x %in% names(col_mapping)) col_mapping[[x]] else x
    })
  } else {
    colnames(value_mat) <- data_cols
  }
  rownames(sig_mat) <- df[[row_col]]
  colnames(sig_mat) <- colnames(value_mat)

  # Get scores
  scores <- df$Score

  # Set color scale based on data range
  max_abs <- max(abs(value_mat), na.rm = TRUE)
  if (is.na(max_abs) || max_abs == 0) max_abs <- 1

  colors <- heatmap_config$color_scale
  col_fun <- colorRamp2(c(-max_abs, 0, max_abs), colors)

  # =============================================================================
  # DETERMINE ROW ORDER (score-sorted or clustered)
  # =============================================================================
  cluster_rows <- heatmap_config$cluster_rows %||% FALSE

  if (cluster_rows && n_rows > 1) {
    # Hierarchical clustering
    row_dist <- dist(value_mat)
    row_hclust <- hclust(row_dist, method = heatmap_config$clustering_method)
    row_dend <- as.dendrogram(row_hclust)
    row_order_vec <- row_hclust$order
  } else {
    # Keep original score order (default)
    row_dend <- FALSE
    row_order_vec <- 1:n_rows
  }

  # Score color scale
  score_range <- range(scores, na.rm = TRUE)
  if (score_range[1] == score_range[2]) {
    score_col <- colorRamp2(c(0, max(1, score_range[2])), c("grey90", "darkgreen"))
  } else {
    score_col <- colorRamp2(c(score_range[1], score_range[2]), c("grey90", "darkgreen"))
  }

  # Row annotation (use ordered scores for consistency)
  row_ha <- rowAnnotation(
    Score = scores,
    col = list(Score = score_col),
    annotation_name_gp = gpar(fontsize = 10)
  )

  # Create column split factor based on groups
  col_split_labels <- character(length(data_cols))
  for (group_name in names(column_groups)) {
    group_cols <- column_groups[[group_name]]
    col_split_labels[data_cols %in% group_cols] <- group_name
  }
  col_split <- factor(col_split_labels, levels = names(column_groups))

  # Calculate dimensions
  plot_height <- max(heatmap_config$height, n_rows * 0.25)
  fontsize <- heatmap_config$fontsize

  # Adjust row width for pathway names (usually longer than gene names)
  row_width <- if (pathway_mode) unit(12, "cm") else unit(8, "cm")

  # Bottom annotation with column labels
  bottom_ha <- HeatmapAnnotation(
    empty = anno_empty(border = FALSE, height = unit(0.5, "cm")),
    labels = anno_text(
      colnames(value_mat),
      rot = 0,
      just = "center",
      gp = gpar(fontsize = fontsize),
      location = unit(0.5, "npc")
    ),
    which = "column"
  )

  # =============================================================================
  # 1. GENERATE PNG HEATMAP
  # =============================================================================
  ht <- Heatmap(
    value_mat,
    name = value_label,
    col = col_fun,
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = row_dend,  # FALSE or dendrogram
    cluster_columns = FALSE,
    column_split = col_split,
    column_gap = unit(heatmap_config$column_gap, "mm"),
    show_row_names = TRUE,
    show_column_names = FALSE,
    bottom_annotation = bottom_ha,
    row_names_gp = gpar(fontsize = fontsize),
    row_names_max_width = row_width,
    right_annotation = row_ha,
    heatmap_legend_param = list(title = value_label),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(sig_mat[i, j]) && sig_mat[i, j] != "") {
        grid.text(sig_mat[i, j], x, y,
                  gp = gpar(fontsize = heatmap_config$star_fontsize, col = "black"))
      }
    }
  )

  png(output_file, width = heatmap_config$width, height = plot_height + 1.5,
      units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  item_label <- if (pathway_mode) "pathways" else "genes"
  message(paste("Saved PNG:", output_file, "-", n_rows, item_label))

  # =============================================================================
  # 2. EXPORT TABLE TSV (same row order as PNG)
  # =============================================================================
  table_output <- gsub("\\.png$", "_table.tsv", output_file)
  ordered_df <- df[row_order_vec, ]
  write.table(ordered_df, table_output, sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste("Saved TSV:", table_output))

  # =============================================================================
  # 3. GENERATE INTERACTIVE HTML HEATMAP (same row order and colors as PNG)
  # =============================================================================
  html_output <- gsub("\\.png$", ".html", output_file)

  # Create BWR color palette matching PNG (centered at 0)
  bwr_colors <- colorRampPalette(heatmap_config$color_scale)(100)

  # Order matrix to match PNG row order
  mat_ordered <- value_mat[row_order_vec, , drop = FALSE]
  scores_ordered <- scores[row_order_vec]

  # Create column side colors for group annotation
  col_side_df <- data.frame(
    Group = col_split_labels,
    row.names = colnames(mat_ordered)
  )

  # Create row side colors for scores
  row_side_df <- data.frame(
    Score = scores_ordered,
    row.names = rownames(mat_ordered)
  )

  # Determine whether to show row labels (hide if too many)
  show_row_labels <- n_rows <= 100

  # Generate interactive heatmap with heatmaply
  p <- heatmaply(
    mat_ordered,
    Rowv = if (cluster_rows && n_rows > 1) row_dend else FALSE,
    dendrogram = if (cluster_rows && n_rows > 1) "row" else "none",
    Colv = FALSE,
    colors = bwr_colors,
    limits = c(-max_abs, max_abs),
    col_side_colors = col_side_df,
    row_side_colors = row_side_df,
    fontsize_col = heatmap_config$html_fontsize_col,
    fontsize_row = if (show_row_labels) heatmap_config$fontsize else 0,
    showticklabels = c(TRUE, show_row_labels),
    plot_method = "plotly",
    main = title,
    xlab = "",
    ylab = ""
  )

  # Save as self-contained HTML
  saveWidget(p, html_output, selfcontained = TRUE)
  message(paste("Saved HTML:", html_output))
}


load_config <- function(config_path) {
  #' Load configuration from JSON file.
  #'
  #' @param config_path Path to JSON config file
  #' @return Configuration list
  config <- fromJSON(config_path, simplifyVector = FALSE)
  return(config)
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  config_path <- NULL
  up_table <- NULL
  down_table <- NULL
  output_dir <- NULL
  prefix <- NULL
  pathway_mode <- FALSE

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--config") {
      config_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--up_table") {
      up_table <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--down_table") {
      down_table <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output_dir") {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--prefix") {
      prefix <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--pathway_mode") {
      pathway_mode <- TRUE
      i <- i + 1
    } else {
      i <- i + 1
    }
  }

  # Load config or use defaults
  if (!is.null(config_path)) {
    config <- load_config(config_path)

    # Get paths from config if not specified
    if (is.null(output_dir)) output_dir <- config$output$directory
    if (is.null(prefix)) prefix <- config$output$prefix
  } else {
    config <- list(
      column_mapping = DEFAULT_COL_MAPPING,
      column_groups = DEFAULT_COLUMN_GROUPS,
      heatmap = DEFAULT_HEATMAP_CONFIG
    )
  }

  # Validate required paths
  if (is.null(up_table) && !is.null(output_dir) && !is.null(prefix)) {
    up_table <- file.path(output_dir, paste0(prefix, "_UP_by_stim_DOWN_by_treatment.tsv"))
    down_table <- file.path(output_dir, paste0(prefix, "_DOWN_by_stim_UP_by_treatment.tsv"))
  }

  if (is.null(up_table) || is.null(down_table)) {
    stop("Must specify --config or (--up_table and --down_table)")
  }

  if (is.null(output_dir)) {
    output_dir <- DEFAULT_OUTPUT_DIR
  }
  if (is.null(prefix)) {
    prefix <- "analysis"
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("DEG Multi-Dataset Analysis - Heatmap Generator")
  message("=" |> rep(50) |> paste(collapse = ""))

  # Check if pathway mode from config
  if (!pathway_mode && !is.null(config$mode)) {
    pathway_mode <- config$mode == "pathway_analysis"
  }

  # Adjust titles for pathway mode
  if (pathway_mode) {
    up_title <- "Pathways Enriched (UP) in Stimulation & Reversed (DOWN) by Treatment"
    down_title <- "Pathways Enriched (DOWN) in Stimulation & Reversed (UP) by Treatment"
  } else {
    up_title <- "Upregulated by Stimulation & Reversed (Downregulated) by Treatment"
    down_title <- "Downregulated by Stimulation & Reversed (Upregulated) by Treatment"
  }

  # Generate UP heatmap
  if (file.exists(up_table)) {
    message(paste("\nProcessing UP table:", up_table))
    create_heatmap(
      up_table,
      file.path(output_dir, paste0(prefix, "_UP_stim_DOWN_treatment_heatmap.png")),
      up_title,
      config,
      pathway_mode
    )
  } else {
    message(paste("Warning: UP table not found:", up_table))
  }

  # Generate DOWN heatmap
  if (file.exists(down_table)) {
    message(paste("\nProcessing DOWN table:", down_table))
    create_heatmap(
      down_table,
      file.path(output_dir, paste0(prefix, "_DOWN_stim_UP_treatment_heatmap.png")),
      down_title,
      config,
      pathway_mode
    )
  } else {
    message(paste("Warning: DOWN table not found:", down_table))
  }

  message("\nDone!")
}


# Run main if script is executed directly
if (!interactive()) {
  main()
}
