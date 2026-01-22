#!/usr/bin/env Rscript

# ==============================================================================
# Metadata Inspection Script
#
# This script inspects metadata columns in Seurat RDS or h5ad files and provides
# a summary to help users identify appropriate columns for analysis.
#
# Usage:
#   Rscript inspect_metadata.R <input_file>
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript inspect_metadata.R <input_file>")
}

input_file <- args[1]

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Determine file type and load
file_ext <- tolower(tools::file_ext(input_file))

if (file_ext %in% c("rds", "rdata")) {
  message("Loading Seurat RDS file...")
  seurat_obj <- readRDS(input_file)
  metadata <- seurat_obj@meta.data
} else if (file_ext == "h5ad") {
  stop("h5ad files must be converted to Seurat RDS first. Use the anndata-seurat-conversion skill.")
} else {
  stop("Unsupported file format: ", file_ext, ". Supported: .rds, .h5ad")
}

# Print basic info
cat("\n", strrep("=", 60), "\n")
cat("METADATA INSPECTION REPORT\n")
cat(strrep("=", 60), "\n\n")

cat("Input file:", input_file, "\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Total genes:", nrow(seurat_obj), "\n")
cat("Metadata columns:", ncol(metadata), "\n\n")

cat(strrep("-", 60), "\n")
cat("=== Metadata Column Summary ===\n")
cat(strrep("-", 60), "\n\n")

# Helper function to interpret column type
interpret_column <- function(col_name, values, n_unique) {
  col_lower <- tolower(col_name)

  # Check for condition/disease columns
  if (grepl("condition|disease|treatment|group|status|diagnosis", col_lower)) {
    return("Condition/disease status (likely comparison column)")
  }

  # Check for sample columns
  if (grepl("sample|orig.ident|patient|subject|donor|individual", col_lower)) {
    return("Sample identifier")
  }

  # Check for cell type columns
  if (grepl("cell.?type|celltype|annotation|cluster|cell.?identity|label", col_lower)) {
    return("Cell type annotation")
  }

  # Check for batch columns
  if (grepl("batch|experiment|date|run|lane|plate|library", col_lower)) {
    return("Batch/technical variable")
  }

  # Check for QC columns
  if (grepl("ncount|nfeature|percent|mt|ribo|qc", col_lower)) {
    return("QC metric (numeric)")
  }

  # Infer from data characteristics
  if (n_unique == 1) {
    return("Constant value (not useful for analysis)")
  } else if (n_unique == 2) {
    return("Binary variable (potential condition or batch)")
  } else if (n_unique >= 3 && n_unique <= 10) {
    return("Categorical variable (potential condition, batch, or covariate)")
  } else if (n_unique > 10 && n_unique < ncol(seurat_obj) * 0.5) {
    # Check if values look like cell types
    sample_values <- head(unique(values), 5)
    if (any(grepl("cell|cyte|blast|phage|neuron|fibro|endo|epi|macro|mono|lymph|T cell|B cell|NK", sample_values, ignore.case = TRUE))) {
      return("Cell type annotation (many cell types)")
    }
    return("Sample identifier or grouping variable")
  } else {
    return("High cardinality (likely sample-level or continuous)")
  }
}

# Analyze each column
for (col_name in colnames(metadata)) {
  col_values <- metadata[[col_name]]

  # Get unique values
  if (is.factor(col_values)) {
    unique_vals <- levels(col_values)
  } else {
    unique_vals <- unique(col_values)
  }

  n_unique <- length(unique_vals)

  # Get examples
  if (n_unique <= 5) {
    examples <- paste(unique_vals, collapse = ", ")
  } else {
    examples <- paste(c(head(unique_vals, 3), "..."), collapse = ", ")
  }

  # Get interpretation
  interpretation <- interpret_column(col_name, unique_vals, n_unique)

  cat("Column:", col_name, "\n")
  cat("  - Data type:", class(col_values)[1], "\n")
  cat("  - Unique values:", n_unique, "\n")
  cat("  - Examples:", examples, "\n")
  cat("  - Interpretation:", interpretation, "\n\n")
}

# Provide recommendations
cat(strrep("-", 60), "\n")
cat("=== Recommended Column Assignments ===\n")
cat(strrep("-", 60), "\n\n")

# Find best candidates for each role
condition_candidates <- c()
sample_candidates <- c()
celltype_candidates <- c()
batch_candidates <- c()

for (col_name in colnames(metadata)) {
  col_values <- metadata[[col_name]]
  n_unique <- length(unique(col_values))
  col_lower <- tolower(col_name)

  # Condition candidates (2-10 unique values)
  if (n_unique >= 2 && n_unique <= 10) {
    if (grepl("condition|disease|treatment|group|status|diagnosis", col_lower)) {
      condition_candidates <- c(condition_candidates, paste0(col_name, " (", n_unique, " values)"))
    }
  }

  # Sample candidates
  if (grepl("sample|orig.ident|patient|subject|donor|individual", col_lower)) {
    sample_candidates <- c(sample_candidates, paste0(col_name, " (", n_unique, " values)"))
  }

  # Cell type candidates
  if (grepl("cell.?type|celltype|annotation|cluster|cell.?identity|label", col_lower)) {
    celltype_candidates <- c(celltype_candidates, paste0(col_name, " (", n_unique, " values)"))
  }

  # Batch candidates
  if (grepl("batch|experiment|date|run|lane|plate|library", col_lower)) {
    batch_candidates <- c(batch_candidates, paste0(col_name, " (", n_unique, " values)"))
  }
}

cat("Potential CONDITION columns:\n")
if (length(condition_candidates) > 0) {
  for (c in condition_candidates) cat("  -", c, "\n")
} else {
  cat("  (No obvious candidates - look for columns with 2-10 unique values)\n")
}

cat("\nPotential SAMPLE columns:\n")
if (length(sample_candidates) > 0) {
  for (c in sample_candidates) cat("  -", c, "\n")
} else {
  cat("  (No obvious candidates - look for columns like 'orig.ident')\n")
}

cat("\nPotential CELLTYPE columns:\n")
if (length(celltype_candidates) > 0) {
  for (c in celltype_candidates) cat("  -", c, "\n")
} else {
  cat("  (No obvious candidates - look for annotation columns)\n")
}

cat("\nPotential BATCH columns:\n")
if (length(batch_candidates) > 0) {
  for (c in batch_candidates) cat("  -", c, "\n")
} else {
  cat("  (No obvious candidates - batch correction may not be needed)\n")
}

cat("\n")
