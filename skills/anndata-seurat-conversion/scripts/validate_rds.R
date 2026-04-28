#!/usr/bin/env Rscript
# =============================================================================
# Validate that a file is a valid Seurat RDS object and return structure info.
#
# Usage:
#   Rscript validate_rds.R <file.rds>
#
# Output: JSON with validation results and object structure
# Exit codes: 0 = valid Seurat, 1 = invalid/error, 2 = bad arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat('{"valid": false, "error": "No file path provided. Usage: Rscript validate_rds.R <file.rds>"}\n')
  quit(status = 2)
}

input_file <- args[1]

if (!file.exists(input_file)) {
  cat(sprintf('{"valid": false, "error": "File not found: %s"}\n', input_file))
  quit(status = 1)
}

# Try loading
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

obj <- tryCatch(readRDS(input_file), error = function(e) {
  cat(sprintf('{"valid": false, "error": "Failed to read RDS: %s"}\n', gsub('"', '\\\\"', e$message)))
  quit(status = 1)
})

# Check if Seurat object
if (!inherits(obj, "Seurat")) {
  cat(sprintf('{"valid": false, "error": "Object is class %s, not Seurat"}\n', paste(class(obj), collapse = "/")))
  quit(status = 1)
}

# Extract structure
assay_names <- names(obj@assays)
default_assay <- DefaultAssay(obj)
reduction_names <- names(obj@reductions)
meta_cols <- colnames(obj@meta.data)

# Per-assay info
assay_info <- list()
for (a in assay_names) {
  layers <- Layers(obj[[a]])
  n_features <- nrow(obj[[a]])
  var_cols <- colnames(obj[[a]][[]])
  assay_info[[a]] <- list(
    n_features = n_features,
    layers = layers,
    var_columns = var_cols
  )
}

# Per-reduction info
reduc_info <- list()
for (r in reduction_names) {
  emb <- Embeddings(obj, r)
  load_mat <- tryCatch(Loadings(obj, r), silent = TRUE, error = function(e) matrix(nrow = 0, ncol = 0))
  reduc_info[[r]] <- list(
    dims = ncol(emb),
    cells = nrow(emb),
    has_loadings = nrow(load_mat) > 0 && ncol(load_mat) > 0,
    loadings_genes = nrow(load_mat)
  )
}

# Factor columns in metadata
factor_cols <- names(which(sapply(obj@meta.data, is.factor)))

# Build JSON manually (avoids jsonlite dependency)
to_json_array <- function(x) paste0("[", paste0('"', x, '"', collapse = ", "), "]")
to_json_num <- function(x) paste0("[", paste0(x, collapse = ", "), "]")

assay_json <- paste0(sapply(names(assay_info), function(a) {
  ai <- assay_info[[a]]
  sprintf('"%s": {"n_features": %d, "layers": %s, "var_columns": %s}',
          a, ai$n_features, to_json_array(ai$layers), to_json_array(ai$var_columns))
}), collapse = ", ")

reduc_json <- paste0(sapply(names(reduc_info), function(r) {
  ri <- reduc_info[[r]]
  sprintf('"%s": {"dims": %d, "cells": %d, "has_loadings": %s, "loadings_genes": %d}',
          r, ri$dims, ri$cells, tolower(as.character(ri$has_loadings)), ri$loadings_genes)
}), collapse = ", ")

cat(sprintf(paste0(
  '{\n',
  '  "valid": true,\n',
  '  "file_path": "%s",\n',
  '  "class": "Seurat",\n',
  '  "n_cells": %d,\n',
  '  "default_assay": "%s",\n',
  '  "assays": {%s},\n',
  '  "reductions": {%s},\n',
  '  "metadata_columns": %s,\n',
  '  "factor_columns": %s,\n',
  '  "n_metadata_cols": %d\n',
  '}\n'),
  input_file,
  ncol(obj),
  default_assay,
  assay_json,
  reduc_json,
  to_json_array(meta_cols),
  to_json_array(factor_cols),
  length(meta_cols)
))

quit(status = 0)
