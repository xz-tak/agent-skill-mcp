#!/usr/bin/env Rscript
# =============================================================================
# Convert Seurat (v5) RDS to H5AD (AnnData) with full preservation of:
#   - X: log-normalized expression from RNA$data
#   - layers: raw counts from RNA$counts (+ additional layers)
#   - obs: cell metadata (factors/logicals converted for h5ad)
#   - var: gene metadata from the source assay
#   - obsm: all dimensional reductions (PCA, UMAP, t-SNE, etc.)
#   - varm: gene loadings from PCA reductions (zero-padded to full gene set)
#
# Usage:
#   Rscript convert_seurat_to_h5ad.R <input.rds> [output.h5ad] [output.log] [assay]
#
# Arguments:
#   input.rds    - Path to input Seurat RDS file (required)
#   output.h5ad  - Path to output h5ad file (optional, defaults to same name in working dir)
#   output.log   - Path to log file (optional, defaults to output with .log extension)
#   assay        - Assay to convert (optional, defaults to "RNA")
#
# Environment:
#   RETICULATE_PYTHON - Path to Python with anndata installed (required)
#                       Set before running or export in shell
# =============================================================================

# -----------------------------------------------------------------------------
# Parse Command Line Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("ERROR: No input file specified.\n\n")
  cat("Usage: Rscript convert_seurat_to_h5ad.R <input.rds> [output.h5ad] [output.log] [assay]\n\n")
  cat("Arguments:\n")
  cat("  input.rds    - Path to input Seurat RDS file (required)\n")
  cat("  output.h5ad  - Path to output h5ad file (optional)\n")
  cat("  output.log   - Path to log file (optional)\n")
  cat("  assay        - Assay to convert (optional, default: RNA)\n\n")
  cat("Environment:\n")
  cat("  RETICULATE_PYTHON - Path to Python with anndata installed\n")
  quit(status = 1)
}

input_file <- args[1]

if (!file.exists(input_file)) {
  cat(sprintf("ERROR: Input file not found: %s\n", input_file))
  quit(status = 1)
}

# Output defaults to working directory with same basename, .h5ad extension
if (length(args) >= 2) {
  output_file <- args[2]
} else {
  base_name <- sub("\\.[Rr][Dd][Ss]$", "", basename(input_file))
  output_file <- file.path(getwd(), paste0(base_name, ".h5ad"))
}

if (length(args) >= 3) {
  log_file <- args[3]
} else {
  log_file <- sub("\\.h5ad$", ".log", output_file, ignore.case = TRUE)
}

source_assay <- if (length(args) >= 4) args[4] else "RNA"

# Create output directory if needed
output_dir <- dirname(output_file)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check RETICULATE_PYTHON
reticulate_python <- Sys.getenv("RETICULATE_PYTHON", "")
if (reticulate_python == "") {
  cat("ERROR: RETICULATE_PYTHON environment variable not set.\n")
  cat("Set it to a Python with anndata installed, e.g.:\n")
  cat("  export RETICULATE_PYTHON=/path/to/conda/envs/myenv/bin/python3\n")
  quit(status = 1)
}

if (!file.exists(reticulate_python)) {
  cat(sprintf("ERROR: RETICULATE_PYTHON path not found: %s\n", reticulate_python))
  quit(status = 1)
}

# -----------------------------------------------------------------------------
# Setup logging
# -----------------------------------------------------------------------------
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message", append = TRUE)

warnings_list <- character()
log_warning <- function(msg) {
  warnings_list <<- c(warnings_list, msg)
  cat(sprintf("WARNING: %s\n", msg))
}

cat("=============================================================================\n")
cat("Seurat RDS to H5AD (AnnData) Conversion\n")
cat("=============================================================================\n")
cat(sprintf("Start time: %s\n\n", Sys.time()))
cat(sprintf("Input:  %s\n", input_file))
cat(sprintf("Output: %s\n", output_file))
cat(sprintf("Log:    %s\n", log_file))
cat(sprintf("Assay:  %s\n", source_assay))
cat(sprintf("Python: %s\n\n", reticulate_python))

# -----------------------------------------------------------------------------
# 1. Libraries
# NOTE: Load Seurat FIRST, then anndata LAST (anndata masks SeuratObject::Layers)
# All Seurat extraction must happen before library(anndata) is called.
# -----------------------------------------------------------------------------
cat("Step 1: Loading libraries (Seurat + Matrix)...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(reticulate)
})

cat(sprintf("  - Seurat version: %s\n", packageVersion("Seurat")))

if (packageVersion("Seurat") < "5.0.0") {
  sink(type = "output"); sink(type = "message"); close(log_conn)
  stop("This script requires Seurat v5+.")
}

# -----------------------------------------------------------------------------
# 2. Load Seurat object
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading Seurat object...\n")
seu <- tryCatch(readRDS(input_file), error = function(e) {
  cat(sprintf("ERROR: Failed to read RDS: %s\n", e$message))
  sink(type = "output"); sink(type = "message"); close(log_conn)
  quit(status = 1)
})

if (!inherits(seu, "Seurat")) {
  cat(sprintf("ERROR: Object is class '%s', not Seurat.\n", paste(class(seu), collapse = "/")))
  sink(type = "output"); sink(type = "message"); close(log_conn)
  quit(status = 1)
}

if (!(source_assay %in% names(seu@assays))) {
  cat(sprintf("ERROR: Assay '%s' not found. Available: %s\n",
              source_assay, paste(names(seu@assays), collapse = ", ")))
  sink(type = "output"); sink(type = "message"); close(log_conn)
  quit(status = 1)
}

n_cells <- ncol(seu)
n_genes <- nrow(seu[[source_assay]])
cat(sprintf("  - Loaded: %d cells x %d genes (assay: %s)\n", n_cells, n_genes, source_assay))
cat(sprintf("  - Available assays: %s\n", paste(names(seu@assays), collapse = ", ")))
cat(sprintf("  - Available layers in %s: %s\n", source_assay,
            paste(Layers(seu[[source_assay]]), collapse = ", ")))

# -----------------------------------------------------------------------------
# 3. Extract X (normalized data, cells x genes)
# -----------------------------------------------------------------------------
cat("\nStep 3: Extracting X (normalized data)...\n")

assay_layers <- Layers(seu[[source_assay]])

if ("data" %in% assay_layers) {
  X <- t(LayerData(seu, assay = source_assay, layer = "data"))
  cat(sprintf("  - X from '%s$data': %d x %d\n", source_assay, nrow(X), ncol(X)))
} else if ("counts" %in% assay_layers) {
  log_warning(sprintf("No 'data' layer in %s; using 'counts' as X.", source_assay))
  X <- t(LayerData(seu, assay = source_assay, layer = "counts"))
  cat(sprintf("  - X from '%s$counts' (fallback): %d x %d\n", source_assay, nrow(X), ncol(X)))
} else {
  cat(sprintf("ERROR: No 'data' or 'counts' layer in assay '%s'.\n", source_assay))
  sink(type = "output"); sink(type = "message"); close(log_conn)
  quit(status = 1)
}

# Report X range
if (inherits(X, "sparseMatrix") && length(X@x) > 0) {
  cat(sprintf("  - X range: [%.4f, %.4f]\n", min(X@x), max(X@x)))
} else if (!inherits(X, "sparseMatrix")) {
  cat(sprintf("  - X range: [%.4f, %.4f]\n", min(X), max(X)))
}

# -----------------------------------------------------------------------------
# 4. Extract layers (counts + any additional)
# -----------------------------------------------------------------------------
cat("\nStep 4: Extracting layers...\n")

layers <- list()

# Raw counts
if ("counts" %in% assay_layers) {
  counts_mat <- t(LayerData(seu, assay = source_assay, layer = "counts"))
  layers[["counts"]] <- counts_mat
  cat(sprintf("  - counts: %d x %d\n", nrow(counts_mat), ncol(counts_mat)))
  if (inherits(counts_mat, "sparseMatrix") && length(counts_mat@x) > 0) {
    cat(sprintf("    range: [%.0f, %.0f]\n", min(counts_mat@x), max(counts_mat@x)))
  }
} else {
  log_warning("No 'counts' layer found; layers will not include raw counts.")
}

# Additional layers (skip 'data', 'counts', 'scale.data' which are standard Seurat)
skip_layers <- c("data", "counts", "scale.data")
extra_layers <- setdiff(assay_layers, skip_layers)
if (length(extra_layers) > 0) {
  for (ln in extra_layers) {
    tryCatch({
      layer_mat <- t(LayerData(seu, assay = source_assay, layer = ln))
      layers[[ln]] <- layer_mat
      cat(sprintf("  - %s: %d x %d\n", ln, nrow(layer_mat), ncol(layer_mat)))
    }, error = function(e) {
      log_warning(sprintf("Failed to extract layer '%s': %s", ln, e$message))
    })
  }
}

if (length(layers) == 0) {
  cat("  - No layers extracted.\n")
}

# -----------------------------------------------------------------------------
# 5. Extract obs (cell metadata)
# -----------------------------------------------------------------------------
cat("\nStep 5: Extracting obs (cell metadata)...\n")

obs <- seu@meta.data

# Convert factors to character for h5ad compatibility
factor_cols <- names(which(sapply(obs, is.factor)))
if (length(factor_cols) > 0) {
  cat(sprintf("  - Converting factor columns: %s\n", paste(factor_cols, collapse = ", ")))
  obs[factor_cols] <- lapply(obs[factor_cols], as.character)
}

# Convert logical to character (avoids h5ad boolean edge cases with NAs)
logical_cols <- names(which(sapply(obs, is.logical)))
if (length(logical_cols) > 0) {
  cat(sprintf("  - Converting logical columns: %s\n", paste(logical_cols, collapse = ", ")))
  obs[logical_cols] <- lapply(obs[logical_cols], function(x) {
    ifelse(is.na(x), NA_character_, as.character(x))
  })
}

cat(sprintf("  - obs: %d cells x %d columns\n", nrow(obs), ncol(obs)))

# -----------------------------------------------------------------------------
# 6. Extract var (gene metadata)
# -----------------------------------------------------------------------------
cat("\nStep 6: Extracting var (gene metadata)...\n")

gene_names <- rownames(seu[[source_assay]])
var_df <- data.frame(row.names = gene_names)

# Add feature metadata from the source assay
assay_var <- seu[[source_assay]][[]]
if (ncol(assay_var) > 0) {
  for (col in colnames(assay_var)) {
    var_df[[col]] <- assay_var[[col]]
  }
  cat(sprintf("  - Added %d var columns from %s: %s\n",
              ncol(assay_var), source_assay, paste(colnames(assay_var), collapse = ", ")))
}

# If source assay differs from default, also check other assays for var metadata
other_assays <- setdiff(names(seu@assays), source_assay)
for (other_a in other_assays) {
  other_var <- seu[[other_a]][[]]
  if (ncol(other_var) > 0) {
    other_genes <- rownames(other_var)
    new_cols <- setdiff(colnames(other_var), colnames(var_df))
    if (length(new_cols) > 0) {
      for (col in new_cols) {
        var_df[[col]] <- NA
        matching <- intersect(other_genes, gene_names)
        if (length(matching) > 0) {
          var_df[matching, col] <- other_var[matching, col]
        }
      }
      cat(sprintf("  - Added %d var columns from %s: %s\n",
                  length(new_cols), other_a, paste(new_cols, collapse = ", ")))
    }
  }
}

cat(sprintf("  - var: %d genes x %d columns\n", nrow(var_df), ncol(var_df)))

# -----------------------------------------------------------------------------
# 7. Extract obsm (dimensional reductions -> embeddings)
# -----------------------------------------------------------------------------
cat("\nStep 7: Extracting obsm (reductions)...\n")

obsm <- list()
reduction_names <- names(seu@reductions)

if (length(reduction_names) == 0) {
  cat("  - No reductions found.\n")
} else {
  for (rname in reduction_names) {
    emb <- Embeddings(seu, rname)
    # AnnData convention: prefix with X_
    ad_key <- paste0("X_", rname)
    # Strip row/col names for clean numpy arrays
    colnames(emb) <- NULL
    rownames(emb) <- NULL
    obsm[[ad_key]] <- emb
    cat(sprintf("  - %s: %d x %d\n", ad_key, nrow(emb), ncol(emb)))
  }
}

# -----------------------------------------------------------------------------
# 8. Extract varm (gene loadings from PCA-like reductions)
# -----------------------------------------------------------------------------
cat("\nStep 8: Extracting varm (gene loadings)...\n")

varm <- list()

for (rname in reduction_names) {
  load_mat <- tryCatch(Loadings(seu, rname), silent = TRUE, error = function(e) {
    matrix(nrow = 0, ncol = 0)
  })
  if (nrow(load_mat) > 0 && ncol(load_mat) > 0) {
    n_pcs <- ncol(load_mat)
    # Zero-pad to full gene set for AnnData compatibility
    full_load <- matrix(0, nrow = length(gene_names), ncol = n_pcs)
    rownames(full_load) <- gene_names
    matching <- intersect(rownames(load_mat), gene_names)
    full_load[matching, ] <- load_mat[matching, ]
    # Clean names
    colnames(full_load) <- NULL
    rownames(full_load) <- NULL

    # Name: use reduction name with _loadings suffix for clarity
    varm_key <- if (rname == "pca") "PCs" else paste0(rname, "_loadings")
    varm[[varm_key]] <- full_load
    cat(sprintf("  - %s: %d x %d (%d non-zero genes)\n",
                varm_key, nrow(full_load), ncol(full_load), length(matching)))
  }
}

if (length(varm) == 0) {
  cat("  - No gene loadings found.\n")
}

# -----------------------------------------------------------------------------
# 9. Create AnnData and write h5ad
# NOTE: Load anndata NOW (after all Seurat extraction) to avoid masking Layers()
# -----------------------------------------------------------------------------
cat("\nStep 9: Loading anndata R package and creating AnnData object...\n")
suppressPackageStartupMessages(library(anndata))
cat(sprintf("  - anndata (R) version: %s\n", packageVersion("anndata")))

adata_args <- list(X = X, obs = obs, var = var_df)
if (length(layers) > 0) adata_args$layers <- layers
if (length(obsm) > 0)   adata_args$obsm <- obsm
if (length(varm) > 0)   adata_args$varm <- varm

adata <- do.call(AnnData, adata_args)

cat(sprintf("  - AnnData: %d cells x %d genes\n", adata$n_obs, adata$n_vars))

cat("\nStep 10: Writing h5ad...\n")
adata$write_h5ad(output_file)

file_size_mb <- file.info(output_file)$size / 1e6
cat(sprintf("  - Saved to: %s\n", output_file))
cat(sprintf("  - File size: %.1f MB\n", file_size_mb))

# -----------------------------------------------------------------------------
# 10. Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("CONVERSION SUMMARY\n")
cat("=============================================================================\n")
cat(sprintf("Status:      SUCCESS\n"))
cat(sprintf("Input:       %s\n", input_file))
cat(sprintf("Output:      %s (%.1f MB)\n", output_file, file_size_mb))
cat(sprintf("Log:         %s\n", log_file))
cat(sprintf("Assay:       %s\n", source_assay))
cat(sprintf("Cells:       %d\n", adata$n_obs))
cat(sprintf("Genes:       %d\n", adata$n_vars))
cat(sprintf("X:           %s$data (normalized)\n", source_assay))
cat(sprintf("Layers:      %s\n", if (length(layers) > 0) paste(names(layers), collapse = ", ") else "None"))
cat(sprintf("Obs:         %d columns\n", ncol(obs)))
cat(sprintf("Var:         %d columns\n", ncol(var_df)))
cat(sprintf("Obsm:        %s\n", if (length(obsm) > 0) paste(names(obsm), collapse = ", ") else "None"))
cat(sprintf("Varm:        %s\n", if (length(varm) > 0) paste(names(varm), collapse = ", ") else "None"))

if (length(warnings_list) > 0) {
  cat(sprintf("\nWarnings (%d):\n", length(warnings_list)))
  for (w in warnings_list) cat(sprintf("  - %s\n", w))
}

cat(sprintf("\nEnd time: %s\n", Sys.time()))
cat("=============================================================================\n")

sink(type = "output")
sink(type = "message")
close(log_conn)

quit(status = 0)
