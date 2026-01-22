#!/usr/bin/env Rscript
# =============================================================================
# Convert H5AD (AnnData) to Seurat (v5) with:
#   - counts: from ad$raw$X (preferred) or counts-like layers if valid integer/nonneg
#   - data: from ad$X (assumed log1p-normalized)
#   - layers: all ad$layers copied into Seurat assay layers
#   - obsm: all ad$obsm copied into Seurat reductions (DimReduc)
#
# Usage:
#   Rscript convert_h5ad_to_seurat.R <input.h5ad> [output.rds] [output.log]
#
# Arguments:
#   input.h5ad   - Path to input h5ad file (required)
#   output.rds   - Path to output RDS file (optional, defaults to input with .rds extension)
#   output.log   - Path to log file (optional, defaults to output with .log extension)
# =============================================================================

# -----------------------------------------------------------------------------
# Parse Command Line Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("ERROR: No input file specified.\n\n")
  cat("Usage: Rscript convert_h5ad_to_seurat.R <input.h5ad> [output.rds] [output.log]\n\n")
  cat("Arguments:\n")
  cat("  input.h5ad   - Path to input h5ad file (required)\n")
  cat("  output.rds   - Path to output RDS file (optional)\n")
  cat("  output.log   - Path to log file (optional)\n\n")
  cat("Examples:\n")
  cat("  Rscript convert_h5ad_to_seurat.R my_data.h5ad\n")
  cat("  Rscript convert_h5ad_to_seurat.R my_data.h5ad converted.rds\n")
  cat("  Rscript convert_h5ad_to_seurat.R my_data.h5ad converted.rds conversion.log\n")
  quit(status = 1)
}

input_file <- args[1]

# Check if input file exists
if (!file.exists(input_file)) {
  cat(sprintf("ERROR: Input file not found: %s\n", input_file))
  quit(status = 1)
}

# Set output file
if (length(args) >= 2) {
  output_file <- args[2]
} else {
  output_file <- sub("\\.h5ad$", ".rds", input_file, ignore.case = TRUE)
}

# Set log file
if (length(args) >= 3) {
  log_file <- args[3]
} else {
  log_file <- sub("\\.rds$", ".log", output_file, ignore.case = TRUE)
}

# Create output directory if needed
output_dir <- dirname(output_file)
if (output_dir != "." && !dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Setup logging - capture all output to both console and log file
# -----------------------------------------------------------------------------
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message", append = TRUE)

# Track warnings for summary
warnings_list <- character()
log_warning <- function(msg) {
  warnings_list <<- c(warnings_list, msg)
  cat(sprintf("WARNING: %s\n", msg))
}

# Track conversion summary
conversion_summary <- list(
  status = "SUCCESS",
  input_file = input_file,
  output_file = output_file,
  log_file = log_file,
  cells = 0,
  genes = 0,
  metadata_cols = 0,
  counts_source = "",
  layers = character(),
  reductions = character(),
  warnings = character()
)

cat("=============================================================================\n")
cat("H5AD to Seurat v5 Conversion (X + layers + raw + obsm)\n")
cat("=============================================================================\n")
cat(sprintf("Start time: %s\n\n", Sys.time()))

cat(sprintf("Input:  %s\n", input_file))
cat(sprintf("Output: %s\n", output_file))
cat(sprintf("Log:    %s\n\n", log_file))

# -----------------------------------------------------------------------------
# 1. Libraries
# -----------------------------------------------------------------------------
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(anndata)
  library(Seurat)
  library(Matrix)
  library(reticulate)
})

cat(sprintf("  - Seurat version: %s\n", packageVersion("Seurat")))
if (packageVersion("Seurat") < "5.0.0") {
  sink(type = "output")
  sink(type = "message")
  close(log_conn)
  stop("This script is intended for Seurat v5+. Please upgrade Seurat.")
}
cat(sprintf("  - anndata version: %s\n", packageVersion("anndata")))
cat(sprintf("  - Matrix version: %s\n", packageVersion("Matrix")))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
sanitize_layer_name <- function(x) {
  y <- gsub("[^A-Za-z0-9_]", "_", x)
  if (grepl("^[0-9]", y)) y <- paste0("L_", y)
  y
}

make_reduc_name <- function(nm) {
  nm2 <- sub("^X_", "", nm)
  nm2 <- gsub("[^A-Za-z0-9_]", "_", nm2)
  tolower(nm2)
}

make_reduc_key <- function(rname) {
  if (rname == "pca")  return("PC_")
  if (rname == "umap") return("UMAP_")
  if (rname == "tsne") return("tSNE_")
  key <- toupper(gsub("[^A-Za-z0-9]", "", rname))
  if (nchar(key) == 0) key <- "OBSM"
  paste0(substr(key, 1, 10), "_")
}

# Convert AnnData matrix (cells x genes) -> Seurat (genes x cells) dgCMatrix
to_gene_by_cell_dgC <- function(X, var_names, obs_names) {
  # Use reticulate to convert Python objects first (handles memory better)
  Xr <- tryCatch(reticulate::py_to_r(X), error = function(e) X)

  if (inherits(Xr, "dgRMatrix")) {
    Xt <- as(Xr, "TsparseMatrix")
    m  <- t(Xt)
    m  <- as(m, "CsparseMatrix")
  } else if (inherits(Xr, "dgCMatrix")) {
    m <- as(t(Xr), "CsparseMatrix")
  } else if (inherits(Xr, "sparseMatrix")) {
    m <- as(t(Xr), "CsparseMatrix")
  } else {
    m <- as(t(as.matrix(Xr)), "dgCMatrix")
  }

  rownames(m) <- var_names
  colnames(m) <- obs_names
  m
}

# Check matrix is non-negative and integer-like (within tol).
# For speed on huge matrices, checks at most max_check non-zeros.
is_nonneg_integerish <- function(m, tol = 1e-6, max_check = 2e6) {
  sampled <- FALSE
  if (inherits(m, "sparseMatrix")) {
    x <- m@x
    if (length(x) == 0) return(list(ok = TRUE, sampled = FALSE))
    if (length(x) > max_check) {
      x <- x[sample.int(length(x), max_check)]
      sampled <- TRUE
    }
    ok <- all(x >= -tol) && all(abs(x - round(x)) <= tol)
    return(list(ok = ok, sampled = sampled))
  } else {
    x <- as.numeric(m)
    if (length(x) == 0) return(list(ok = TRUE, sampled = FALSE))
    if (length(x) > max_check) {
      x <- x[sample.int(length(x), max_check)]
      sampled <- TRUE
    }
    ok <- all(x >= -tol) && all(abs(x - round(x)) <= tol)
    return(list(ok = ok, sampled = sampled))
  }
}

# -----------------------------------------------------------------------------
# 2. Load H5AD
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading h5ad...\n")

ad <- tryCatch({
  read_h5ad(input_file)
}, error = function(e) {
  conversion_summary$status <<- "FAILED"
  conversion_summary$warnings <<- c(conversion_summary$warnings, sprintf("Failed to load h5ad: %s", e$message))
  cat(sprintf("ERROR: Failed to load h5ad file: %s\n", e$message))
  cat("\nThis file may not be a valid AnnData h5ad object.\n")
  sink(type = "output")
  sink(type = "message")
  close(log_conn)
  quit(status = 1)
})

cat(sprintf("  - Loaded AnnData: %d cells x %d genes\n", nrow(ad$obs), length(ad$var_names)))
conversion_summary$cells <- nrow(ad$obs)
conversion_summary$genes <- length(ad$var_names)

# -----------------------------------------------------------------------------
# 3. Extract metadata, X, layers, raw.X and pick counts
# -----------------------------------------------------------------------------
cat("\nStep 3: Extracting X / layers / raw and selecting counts...\n")

# Metadata
meta_data <- as.data.frame(ad$obs)
rownames(meta_data) <- ad$obs_names
cat(sprintf("  - Metadata columns: %d\n", ncol(meta_data)))
conversion_summary$metadata_cols <- ncol(meta_data)

# data from X (assumed log1p normalized)
cat("  - Converting ad$X -> data_matrix...\n")
data_matrix <- to_gene_by_cell_dgC(ad$X, ad$var_names, ad$obs_names)
cat(sprintf("    data_matrix dims: %d genes x %d cells\n", nrow(data_matrix), ncol(data_matrix)))

# Data range check
data_vals <- data_matrix@x
if (length(data_vals) > 0) {
  cat(sprintf("    data range: [%.4f, %.4f], mean: %.4f\n",
              min(data_vals), max(data_vals), mean(data_vals)))
}

# layers
layer_names <- tryCatch(names(ad$layers), error = function(e) NULL)
layer_mats  <- list()
if (!is.null(layer_names) && length(layer_names) > 0) {
  cat(sprintf("  - Found layers (%d): %s\n", length(layer_names), paste(layer_names, collapse = ", ")))
  for (ln in layer_names) {
    cat(sprintf("    * Converting layer '%s'...\n", ln))
    tryCatch({
      layer_mats[[ln]] <- to_gene_by_cell_dgC(ad$layers[[ln]], ad$var_names, ad$obs_names)
      cat(sprintf("      dims: %d x %d\n", nrow(layer_mats[[ln]]), ncol(layer_mats[[ln]])))
    }, error = function(e) {
      log_warning(sprintf("Failed to convert layer '%s': %s", ln, e$message))
    })
  }
} else {
  cat("  - No layers found.\n")
}

# raw.X
raw_X_mat <- NULL
if (!is.null(ad$raw)) {
  raw_X <- tryCatch(ad$raw$X, error = function(e) NULL)
  if (!is.null(raw_X)) {
    cat("  - Found ad$raw$X; converting...\n")
    tryCatch({
      raw_X_mat <- to_gene_by_cell_dgC(raw_X, ad$var_names, ad$obs_names)
      cat(sprintf("    raw_X_mat dims: %d x %d\n", nrow(raw_X_mat), ncol(raw_X_mat)))
    }, error = function(e) {
      log_warning(sprintf("Failed to convert raw$X: %s", e$message))
    })
  } else {
    cat("  - ad$raw exists but raw$X not accessible.\n")
  }
} else {
  cat("  - No ad$raw found.\n")
}

# Pick counts: raw.X first if valid, else counts-like layers if valid, else fallback to X
counts_candidates <- c("count", "counts", "raw_count", "raw_counts")
counts_matrix <- NULL
counts_source <- NULL

if (!is.null(raw_X_mat)) {
  chk <- is_nonneg_integerish(raw_X_mat)
  cat(sprintf("  - raw.X integer/nonneg check: %s%s\n",
              if (chk$ok) "OK" else "FAIL",
              if (chk$sampled) " (sampled)" else ""))
  if (chk$ok) {
    counts_matrix <- raw_X_mat
    counts_source <- "raw.X"
  }
}

if (is.null(counts_matrix) && length(layer_mats) > 0) {
  layer_lc <- tolower(names(layer_mats))
  for (nm in counts_candidates) {
    hit <- which(layer_lc == nm)
    if (length(hit) == 0) next
    ln <- names(layer_mats)[hit[1]]
    chk <- is_nonneg_integerish(layer_mats[[ln]])
    cat(sprintf("  - layer '%s' integer/nonneg check: %s%s\n",
                ln,
                if (chk$ok) "OK" else "FAIL",
                if (chk$sampled) " (sampled)" else ""))
    if (chk$ok) {
      counts_matrix <- layer_mats[[ln]]
      counts_source <- paste0("layer:", ln)
      break
    }
  }
}

if (is.null(counts_matrix)) {
  log_warning("No valid integer/non-negative counts found in raw.X or counts-like layers. Using X as fallback.")
  counts_matrix <- data_matrix
  counts_source <- "X_fallback"
}

cat(sprintf("  - Selected counts source: %s\n", counts_source))
cat(sprintf("    counts_matrix dims: %d genes x %d cells\n", nrow(counts_matrix), ncol(counts_matrix)))
conversion_summary$counts_source <- counts_source

# -----------------------------------------------------------------------------
# 4. Create Seurat object using chosen counts, then set data = X
# -----------------------------------------------------------------------------
cat("\nStep 4: Creating Seurat object (counts) + setting data (X)...\n")

seurat_obj <- CreateSeuratObject(
  counts    = counts_matrix,
  meta.data = meta_data,
  project   = "h5ad_converted"
)

# Set normalized data from X
seurat_obj[["RNA"]]$data <- data_matrix

cat(sprintf("  - Created Seurat: %d features x %d cells\n", nrow(seurat_obj), ncol(seurat_obj)))

# -----------------------------------------------------------------------------
# 5. Add layers to Seurat assay (skip the layer used as counts to avoid duplication)
# -----------------------------------------------------------------------------
cat("\nStep 5: Adding AnnData layers into Seurat assay layers...\n")

# Determine which layer was used as counts (to skip it)
counts_layer_used <- NULL
if (grepl("^layer:", counts_source)) {
  counts_layer_used <- sub("^layer:", "", counts_source)
  cat(sprintf("  - Skipping layer '%s' (already used as counts)\n", counts_layer_used))
}

if (length(layer_mats) == 0) {
  cat("  - No additional layers to add.\n")
} else {
  for (ln in names(layer_mats)) {
    # Skip the layer that was used as counts (avoid storing twice -> C stack overflow)
    if (!is.null(counts_layer_used) && ln == counts_layer_used) {
      next
    }

    target <- sanitize_layer_name(ln)
    # Avoid overwriting counts/data
    if (target %in% c("counts", "data")) target <- paste0("layer_", target)

    tryCatch({
      # Use proper Seurat v5 method to add layer
      seurat_obj[["RNA"]][[target]] <- layer_mats[[ln]]
      cat(sprintf("  - Added layer: %s (from '%s')\n", target, ln))
    }, error = function(e) {
      log_warning(sprintf("Failed to add layer '%s': %s", target, e$message))
    })
  }
}

# Only add raw.X as separate layer if it wasn't used as counts
if (!is.null(raw_X_mat) && counts_source != "raw.X") {
  cat("  - Skipping rawX layer (would duplicate data or cause memory issues)\n")
}

# List all layers
all_layers <- Layers(seurat_obj[["RNA"]])
cat(sprintf("  - All layers in RNA assay: %s\n", paste(all_layers, collapse = ", ")))
conversion_summary$layers <- all_layers

# -----------------------------------------------------------------------------
# 6. Add obsm -> Seurat reductions
# -----------------------------------------------------------------------------
cat("\nStep 6: Adding AnnData obsm as Seurat reductions...\n")

obsm_names <- tryCatch(names(ad$obsm), error = function(e) NULL)
if (is.null(obsm_names) || length(obsm_names) == 0) {
  cat("  - No obsm found.\n")
} else {
  cat(sprintf("  - Found obsm (%d): %s\n", length(obsm_names), paste(obsm_names, collapse = ", ")))

  # Known dimensional reduction prefixes (will be converted to Seurat reductions)
  dim_reduc_patterns <- c("^X_pca", "^X_umap", "^X_tsne", "^X_draw_graph", "^X_diffmap", "^X_phate")

  for (nm in obsm_names) {
    cat(sprintf("    * Processing obsm['%s']...\n", nm))

    # Check if this looks like a dimensional reduction (PCA, UMAP, etc.)
    is_dim_reduc <- any(sapply(dim_reduc_patterns, function(p) grepl(p, nm, ignore.case = TRUE)))

    tryCatch({
      emb_py <- ad$obsm[[nm]]
      emb_r  <- tryCatch(reticulate::py_to_r(emb_py), error = function(e) emb_py)

      # DimReduc embeddings are usually dense; convert to numeric matrix
      if (inherits(emb_r, "data.frame")) emb_r <- as.matrix(emb_r)
      if (inherits(emb_r, "dgRMatrix")) {
        emb_r <- as(emb_r, "TsparseMatrix")
        emb_r <- as(emb_r, "CsparseMatrix")
        emb_r <- as.matrix(emb_r)
      } else if (inherits(emb_r, "sparseMatrix")) {
        emb_r <- as(emb_r, "CsparseMatrix")
        emb_r <- as.matrix(emb_r)
      } else {
        emb_r <- as.matrix(emb_r)
      }

      # Check if column names are non-numeric (indicates non-dimensional data like pathway scores)
      col_names <- colnames(emb_r)
      if (!is.null(col_names) && length(col_names) > 0) {
        # If column names are not numeric-like (e.g., "Pathway1", "GREM1", etc.), skip
        numeric_like <- suppressWarnings(!any(is.na(as.numeric(gsub("^[A-Za-z_]+", "", col_names)))))
        if (!numeric_like && !is_dim_reduc) {
          cat(sprintf("      - Skipping '%s': non-dimensional data (named columns: %s...)\n",
                      nm, paste(head(col_names, 3), collapse = ", ")))
          next
        }
      }

      # Skip known non-reduction obsm types
      if (grepl("estimate|pvals|scores|activity", nm, ignore.case = TRUE) && !is_dim_reduc) {
        cat(sprintf("      - Skipping '%s': appears to be scores/estimates, not dimensional reduction\n", nm))
        next
      }

      # Expect cells x dims; fix if transposed
      if (nrow(emb_r) != ncol(seurat_obj) && ncol(emb_r) == ncol(seurat_obj)) {
        cat("      - Detected transposed embedding; transposing...\n")
        emb_r <- t(emb_r)
      }

      if (nrow(emb_r) != ncol(seurat_obj)) {
        cat(sprintf("      - Skipping '%s': rows (%d) != #cells (%d)\n",
                    nm, nrow(emb_r), ncol(seurat_obj)))
        next
      }

      rownames(emb_r) <- colnames(seurat_obj)

      reduc_name <- make_reduc_name(nm)
      reduc_key  <- make_reduc_key(reduc_name)

      # Avoid overwriting existing reductions
      if (reduc_name %in% names(seurat_obj@reductions)) {
        reduc_name <- paste0(reduc_name, "_2")
      }

      seurat_obj[[reduc_name]] <- CreateDimReducObject(
        embeddings = emb_r,
        key        = reduc_key,
        assay      = DefaultAssay(seurat_obj)
      )

      cat(sprintf("      - Added reduction: %s (%d dims)\n", reduc_name, ncol(emb_r)))

    }, error = function(e) {
      log_warning(sprintf("Failed to add obsm '%s': %s", nm, e$message))
    })
  }
}

conversion_summary$reductions <- names(seurat_obj@reductions)

# -----------------------------------------------------------------------------
# 7. Validation
# -----------------------------------------------------------------------------
cat("\nStep 7: Validation...\n")
cat(sprintf("  - counts source: %s\n", counts_source))
cat(sprintf("  - RNA counts dims: %d x %d\n", nrow(seurat_obj[["RNA"]]$counts), ncol(seurat_obj[["RNA"]]$counts)))
cat(sprintf("  - RNA data dims:   %d x %d\n", nrow(seurat_obj[["RNA"]]$data),   ncol(seurat_obj[["RNA"]]$data)))

# Verify counts vs data differ (if not fallback)
if (counts_source != "X_fallback") {
  counts_sum <- sum(seurat_obj[["RNA"]]$counts@x[1:min(1000, length(seurat_obj[["RNA"]]$counts@x))])
  data_sum <- sum(seurat_obj[["RNA"]]$data@x[1:min(1000, length(seurat_obj[["RNA"]]$data@x))])
  cat(sprintf("  - Sample sums: counts=%.2f, data=%.2f (should differ)\n", counts_sum, data_sum))
}

if (length(Layers(seurat_obj[["RNA"]])) > 0) {
  cat(sprintf("  - Layers: %s\n", paste(Layers(seurat_obj[["RNA"]]), collapse = ", ")))
}

if (length(seurat_obj@reductions) > 0) {
  cat(sprintf("  - Reductions: %s\n", paste(names(seurat_obj@reductions), collapse = ", ")))
} else {
  cat("  - Reductions: none\n")
}

# Check metadata
cat(sprintf("  - Metadata columns: %d\n", ncol(seurat_obj@meta.data)))

# -----------------------------------------------------------------------------
# 8. Save
# -----------------------------------------------------------------------------
cat("\nStep 8: Saving Seurat object...\n")
saveRDS(seurat_obj, file = output_file)

cat(sprintf("  - Saved to: %s\n", output_file))
cat(sprintf("  - File size: %.2f MB\n", file.info(output_file)$size / 1024^2))

# -----------------------------------------------------------------------------
# 9. Final Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("CONVERSION SUMMARY\n")
cat("=============================================================================\n")
cat(sprintf("Status:      %s\n", conversion_summary$status))
cat(sprintf("Input:       %s\n", input_file))
cat(sprintf("Output:      %s\n", output_file))
cat(sprintf("Log:         %s\n", log_file))
cat(sprintf("Cells:       %d\n", ncol(seurat_obj)))
cat(sprintf("Genes:       %d\n", nrow(seurat_obj)))
cat(sprintf("Metadata:    %d columns\n", ncol(seurat_obj@meta.data)))
cat(sprintf("Counts from: %s\n", counts_source))
cat(sprintf("Layers:      %s\n", paste(Layers(seurat_obj[["RNA"]]), collapse = ", ")))
cat(sprintf("Reductions:  %s\n",
            ifelse(length(seurat_obj@reductions) > 0,
                   paste(names(seurat_obj@reductions), collapse = ", "),
                   "None")))

if (length(warnings_list) > 0) {
  cat(sprintf("\nWarnings (%d):\n", length(warnings_list)))
  for (w in warnings_list) {
    cat(sprintf("  - %s\n", w))
  }
  conversion_summary$warnings <- warnings_list
}

cat(sprintf("\nEnd time: %s\n", Sys.time()))
cat("=============================================================================\n")

# Close log file
sink(type = "output")
sink(type = "message")
close(log_conn)

# Exit with appropriate status
if (conversion_summary$status == "SUCCESS") {
  quit(status = 0)
} else {
  quit(status = 1)
}
