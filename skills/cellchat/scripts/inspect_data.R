#!/usr/bin/env Rscript
# =============================================================================
# Inspect Seurat RDS File for CellChat Analysis
# =============================================================================
# Inspects a Seurat RDS file and reports metadata for CellChat configuration
#
# Usage:
#   Rscript inspect_data.R <file.rds> [--json]
#
# Arguments:
#   file.rds  - Path to Seurat RDS file (required)
#   --json    - Output in JSON format (optional)
#
# Note: For h5ad files, use the anndata-seurat-conversion skill first:
#   python validate_h5ad.py <file.h5ad>
#   Rscript convert_h5ad_to_seurat.R <file.h5ad> <output.rds>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
json_output <- FALSE
input_file <- NULL

for (arg in args) {
  if (arg == "--json") {
    json_output <- TRUE
  } else if (is.null(input_file)) {
    input_file <- arg
  }
}

# Validate input
if (is.null(input_file)) {
  cat("ERROR: No input file specified.\n\n")
  cat("Usage: Rscript inspect_data.R <file.rds> [--json]\n\n")
  cat("Arguments:\n")
  cat("  file.rds   - Path to Seurat RDS file (required)\n")
  cat("  --json     - Output in JSON format (optional)\n\n")
  cat("Note: For h5ad files, convert first using:\n")
  cat("  Rscript convert_h5ad_to_seurat.R <file.h5ad> <output.rds>\n")
  quit(status = 1)
}

if (!file.exists(input_file)) {
  cat(sprintf("ERROR: File not found: %s\n", input_file))
  quit(status = 1)
}

# Check file extension
file_ext <- tolower(tools::file_ext(input_file))
if (file_ext == "h5ad") {
  cat("ERROR: Input is h5ad format. Please convert to Seurat RDS first:\n\n")
  cat(sprintf("  Rscript /home/sagemaker-user/.claude/skills/anndata-seurat-conversion/scripts/convert_h5ad_to_seurat.R %s\n\n", input_file))
  quit(status = 1)
}

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

# Helper function to get unique values with truncation
get_unique_values <- function(x, max_values = 10) {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals)]
  n_total <- length(vals)

  if (n_total > max_values) {
    return(list(
      values = head(vals, max_values),
      truncated = TRUE,
      total = n_total
    ))
  } else {
    return(list(
      values = vals,
      truncated = FALSE,
      total = n_total
    ))
  }
}

# Helper to detect column type suitability
detect_column_type <- function(x) {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals)]
  n_unique <- length(vals)
  n_total <- length(x)

  if (n_unique <= 1) {
    return("constant")
  } else if (n_unique <= 30) {
    return("categorical")  # Good for cell_type or condition
  } else if (n_unique <= 100) {
    return("high_cardinality")  # Could be patient_id
  } else {
    return("continuous")  # Likely numeric or too many unique values
  }
}

# Load Seurat object
cat("Loading Seurat object...\n")
tryCatch({
  seurat_obj <- readRDS(input_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to load RDS file: %s\n", e$message))
  quit(status = 1)
})

# Validate it's a Seurat object
if (!inherits(seurat_obj, "Seurat")) {
  cat("ERROR: File is not a Seurat object\n")
  quit(status = 1)
}

# Build inspection results
results <- list(
  valid = TRUE,
  file = input_file,
  file_size_mb = round(file.info(input_file)$size / 1024^2, 2),
  n_cells = ncol(seurat_obj),
  n_genes = nrow(seurat_obj),
  assays = Assays(seurat_obj),
  default_assay = DefaultAssay(seurat_obj),
  reductions = names(seurat_obj@reductions),
  metadata_columns = list(),
  suggested_cell_type_cols = character(),
  suggested_condition_cols = character()
)

# Analyze metadata columns
meta_cols <- colnames(seurat_obj@meta.data)
results$n_metadata_cols <- length(meta_cols)

for (col in meta_cols) {
  col_data <- seurat_obj@meta.data[[col]]
  col_type <- detect_column_type(col_data)
  unique_info <- get_unique_values(col_data)

  results$metadata_columns[[col]] <- list(
    type = col_type,
    n_unique = unique_info$total,
    values = unique_info$values,
    truncated = unique_info$truncated
  )

  # Suggest cell type columns (categorical with 2-30 values)
  if (col_type == "categorical" && unique_info$total >= 2 && unique_info$total <= 30) {
    # Prefer columns with "cell", "type", "cluster" in name
    if (grepl("cell|type|cluster|label|annot", col, ignore.case = TRUE)) {
      results$suggested_cell_type_cols <- c(results$suggested_cell_type_cols, col)
    }
  }

  # Suggest condition columns (categorical with 2-10 values)
  if (col_type == "categorical" && unique_info$total >= 2 && unique_info$total <= 10) {
    if (grepl("cond|group|treat|sample|batch|disease|status|stage", col, ignore.case = TRUE)) {
      results$suggested_condition_cols <- c(results$suggested_condition_cols, col)
    }
  }
}

# Check RNA assay layers
if ("RNA" %in% Assays(seurat_obj)) {
  results$rna_layers <- Layers(seurat_obj[["RNA"]])

  # Check for counts and data
  has_counts <- tryCatch({
    !is.null(seurat_obj[["RNA"]]$counts)
  }, error = function(e) FALSE)

  has_data <- tryCatch({
    !is.null(seurat_obj[["RNA"]]$data)
  }, error = function(e) FALSE)

  results$has_counts <- has_counts
  results$has_data <- has_data
}

# Helper function to generate rationale for column
generate_rationale <- function(col_name, n_unique, values) {
  col_lower <- tolower(col_name)

  # Cell type rationale
  if (grepl("cluster|cell|type|label|annot", col_lower)) {
    if (n_unique <= 5) {
      return("Coarse cell type grouping")
    } else if (n_unique <= 15) {
      return("Detailed cell subtype annotation")
    } else {
      return("Fine-grained clusters (consider merging)")
    }
  }

  # Condition rationale
  if (grepl("cond|group|disease|health|status|treat", col_lower)) {
    if (n_unique == 2) {
      return("Binary comparison (case/control)")
    } else if (n_unique <= 5) {
      return("Multi-group comparison")
    } else {
      return("Many groups (may need filtering)")
    }
  }

  # Patient/sample rationale
  if (grepl("patient|subject|sample|donor|id", col_lower)) {
    return("Sample/patient identifier")
  }

  # Tissue rationale
  if (grepl("tissue|organ|region|location", col_lower)) {
    return("Tissue/region annotation")
  }

  # Ontology/popv rationale
  if (grepl("popv|ontology", col_lower)) {
    return("Automated cell type prediction")
  }

  # Default
  return("")
}

# Output results
if (json_output) {
  cat(toJSON(results, pretty = TRUE, auto_unbox = TRUE))
  cat("\n")
} else {
  # Human-readable output
  cat("\n")
  cat("=============================================================================\n")
  cat("Seurat Object Inspection\n")
  cat("=============================================================================\n")
  cat(sprintf("File:       %s\n", results$file))
  cat(sprintf("Size:       %.2f MB\n", results$file_size_mb))
  cat(sprintf("Cells:      %d\n", results$n_cells))
  cat(sprintf("Genes:      %d\n", results$n_genes))
  cat(sprintf("Assays:     %s\n", paste(results$assays, collapse = ", ")))
  cat(sprintf("Reductions: %s\n", ifelse(length(results$reductions) > 0,
                                         paste(results$reductions, collapse = ", "),
                                         "None")))

  if ("RNA" %in% Assays(seurat_obj)) {
    cat(sprintf("RNA Layers: %s\n", paste(results$rna_layers, collapse = ", ")))
    cat(sprintf("Has counts: %s\n", ifelse(results$has_counts, "Yes", "No")))
    cat(sprintf("Has data:   %s\n", ifelse(results$has_data, "Yes", "No")))
  }

  cat("\n")
  cat("=============================================================================\n")
  cat(sprintf("Metadata Columns (%d total)\n", results$n_metadata_cols))
  cat("=============================================================================\n")

  # Get categorical columns
  categorical_cols <- names(results$metadata_columns)[
    sapply(results$metadata_columns, function(x) x$type == "categorical")
  ]

  # Separate into cell type candidates and condition candidates
  cell_type_cols <- character()
  condition_cols <- character()
  other_categorical <- character()

  for (col in categorical_cols) {
    col_lower <- tolower(col)
    info <- results$metadata_columns[[col]]

    if (grepl("cell|type|cluster|label|annot|popv", col_lower) && info$n_unique >= 2) {
      cell_type_cols <- c(cell_type_cols, col)
    } else if (grepl("cond|group|disease|health|status|treat|strictur|fibrosis", col_lower) && info$n_unique >= 2 && info$n_unique <= 10) {
      condition_cols <- c(condition_cols, col)
    } else {
      other_categorical <- c(other_categorical, col)
    }
  }

  # Print Cell Type Candidates
  if (length(cell_type_cols) > 0) {
    cat("\n--- Cell Type Candidates ---\n")
    cat(sprintf("%-30s | %6s | %-40s | %s\n", "Column", "Values", "Sample Values", "Rationale"))
    cat(paste(rep("-", 110), collapse = ""), "\n")
    for (col in cell_type_cols) {
      info <- results$metadata_columns[[col]]
      values_str <- paste(head(info$values, 3), collapse = ", ")
      if (info$n_unique > 3) values_str <- paste0(values_str, ", ...")
      rationale <- generate_rationale(col, info$n_unique, info$values)
      cat(sprintf("%-30s | %6d | %-40s | %s\n",
                  substr(col, 1, 30), info$n_unique, substr(values_str, 1, 40), rationale))
    }
  }

  # Print Condition Candidates
  if (length(condition_cols) > 0) {
    cat("\n--- Condition Candidates ---\n")
    cat(sprintf("%-30s | %6s | %-40s | %s\n", "Column", "Values", "Sample Values", "Rationale"))
    cat(paste(rep("-", 110), collapse = ""), "\n")
    for (col in condition_cols) {
      info <- results$metadata_columns[[col]]
      values_str <- paste(head(info$values, 4), collapse = ", ")
      if (info$n_unique > 4) values_str <- paste0(values_str, ", ...")
      rationale <- generate_rationale(col, info$n_unique, info$values)
      cat(sprintf("%-30s | %6d | %-40s | %s\n",
                  substr(col, 1, 30), info$n_unique, substr(values_str, 1, 40), rationale))
    }
  }

  # Print Other Categorical
  if (length(other_categorical) > 0) {
    cat("\n--- Other Categorical Columns ---\n")
    cat(sprintf("%-30s | %6s | %-40s\n", "Column", "Values", "Sample Values"))
    cat(paste(rep("-", 85), collapse = ""), "\n")
    for (col in other_categorical) {
      info <- results$metadata_columns[[col]]
      values_str <- paste(head(info$values, 3), collapse = ", ")
      if (info$n_unique > 3) values_str <- paste0(values_str, ", ...")
      cat(sprintf("%-30s | %6d | %-40s\n",
                  substr(col, 1, 30), info$n_unique, substr(values_str, 1, 40)))
    }
  }

  # Print high cardinality columns
  high_card_cols <- names(results$metadata_columns)[
    sapply(results$metadata_columns, function(x) x$type == "high_cardinality")
  ]

  if (length(high_card_cols) > 0) {
    cat("\n--- High Cardinality Columns (patient/sample IDs) ---\n")
    for (col in high_card_cols) {
      info <- results$metadata_columns[[col]]
      cat(sprintf("  %s (%d unique values)\n", col, info$n_unique))
    }
  }

  # Recommendations
  cat("\n")
  cat("=============================================================================\n")
  cat("Recommendations for CellChat Analysis\n")
  cat("=============================================================================\n")

  if (length(cell_type_cols) > 0) {
    # Recommend based on number of values
    best_cell_type <- cell_type_cols[which.min(sapply(cell_type_cols, function(x) {
      abs(results$metadata_columns[[x]]$n_unique - 5)  # Prefer ~5 cell types
    }))]
    cat(sprintf("\nRecommended cell_type: %s (%d types)\n",
                best_cell_type, results$metadata_columns[[best_cell_type]]$n_unique))
    cat(sprintf("  Values: %s\n", paste(results$metadata_columns[[best_cell_type]]$values, collapse = ", ")))
  }

  if (length(condition_cols) > 0) {
    cat(sprintf("\nRecommended condition: %s (%d groups)\n",
                condition_cols[1], results$metadata_columns[[condition_cols[1]]]$n_unique))
    cat(sprintf("  Values: %s\n", paste(results$metadata_columns[[condition_cols[1]]]$values, collapse = ", ")))

    # Suggest reference condition
    cond_values <- results$metadata_columns[[condition_cols[1]]]$values
    normal_idx <- grep("normal|control|healthy|wt|baseline", cond_values, ignore.case = TRUE)
    if (length(normal_idx) > 0) {
      cat(sprintf("  Suggested reference: %s\n", cond_values[normal_idx[1]]))
    }
  }

  cat("\n")
  cat("=============================================================================\n")
}
