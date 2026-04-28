#!/usr/bin/env Rscript
# =============================================================================
# BPCells Helper Functions for Large scRNA-seq Datasets
#
# Reusable utilities for downstream tools that can't handle BPCells
# IterableMatrix objects directly (scran, presto, etc).
#
# Usage:
#   source("bpcells_helpers.R")
#   is_bpcells_matrix(counts(sce))
#   sce_agg <- aggregate_bpcells_per_cluster(sce, celltype_col, agg_cols)
#   markers <- bpcells_find_all_markers(seurat_obj, celltype_col)
# =============================================================================

#' Check if a matrix is a BPCells IterableMatrix
is_bpcells_matrix <- function(mat) {
  requireNamespace("BPCells", quietly = TRUE) && inherits(mat, "IterableMatrix")
}

#' Estimate nnz for an h5ad file without loading the full matrix
#' Returns list(n_obs, n_vars, estimated_nnz, exceeds_limit)
estimate_h5ad_size <- function(h5ad_path, nnz_limit = 1.5e9, cell_limit = 1e6) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package required for h5ad size estimation")
  }

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  on.exit(h5$close_all())

  obs_idx <- h5[["obs"]][["_index"]]$read()
  n_obs <- length(obs_idx)

  var_idx <- h5[["var"]][["_index"]]$read()
  n_vars <- length(var_idx)

  estimated_nnz <- NA
  x_group <- NULL
  if ("X" %in% names(h5)) {
    x_item <- h5[["X"]]
    if (inherits(x_item, "H5Group") && "data" %in% names(x_item)) {
      estimated_nnz <- x_item[["data"]]$dims
    }
  }
  if (is.na(estimated_nnz)) {
    # Rough estimate: ~1500 nnz per cell (typical scRNA-seq sparsity)
    estimated_nnz <- n_obs * 1500
  }

  list(
    n_obs = n_obs,
    n_vars = n_vars,
    estimated_nnz = estimated_nnz,
    exceeds_cell_limit = n_obs > cell_limit,
    exceeds_nnz_limit = estimated_nnz > nnz_limit,
    needs_bpcells = (n_obs > cell_limit) || (estimated_nnz > nnz_limit)
  )
}

#' Aggregate BPCells-backed SCE per cluster
#'
#' scran::aggregateAcrossCells can't handle BPCells IterableMatrix.
#' This function materializes per-cluster subsets (each fits in memory)
#' and aggregates them individually.
#'
#' @param sce SingleCellExperiment with BPCells counts
#' @param celltype_col Column name in colData for cell type grouping
#' @param agg_cols Character vector of columns to aggregate by
#' @param use_assay_type Assay to aggregate (default "counts")
#' @return Aggregated SCE with in-memory dgCMatrix counts
aggregate_bpcells_per_cluster <- function(sce, celltype_col, agg_cols,
                                           use_assay_type = "counts") {
  if (!requireNamespace("scuttle", quietly = TRUE)) {
    stop("scuttle package required for aggregation")
  }

  clusters <- unique(SummarizedExperiment::colData(sce)[[celltype_col]])
  clusters <- clusters[!is.na(clusters) & clusters != "" & clusters != "NaN"]
  message("  Aggregating ", length(clusters), " clusters via BPCells per-cluster materialization...")

  agg_list <- list()
  for (cl in clusters) {
    idx <- which(SummarizedExperiment::colData(sce)[[celltype_col]] == cl)
    message("    ", cl, ": ", length(idx), " cells")
    sce_sub <- sce[, idx]
    SummarizedExperiment::assay(sce_sub, use_assay_type) <-
      as(SummarizedExperiment::assay(sce_sub, use_assay_type), "dgCMatrix")
    agg_list[[cl]] <- scuttle::aggregateAcrossCells(
      sce_sub,
      id = SummarizedExperiment::colData(sce_sub)[, agg_cols],
      use.assay.type = use_assay_type
    )
    rm(sce_sub)
    gc(verbose = FALSE)
  }

  result <- do.call(SingleCellExperiment::cbind, agg_list)
  rm(agg_list)
  gc(verbose = FALSE)
  message("  Aggregated to ", ncol(result), " pseudobulk samples")
  result
}

#' Find all markers using Seurat (BPCells-compatible)
#'
#' Replacement for presto::wilcoxauc and scran::findMarkers
#' when the underlying matrix is BPCells on-disk.
#' Seurat v5 FindAllMarkers supports BPCells natively.
#'
#' @param seurat_obj Seurat object (BPCells or in-memory)
#' @param celltype_col Metadata column for cell type identity
#' @param min_pct Minimum fraction of cells expressing (default 0.1)
#' @param logfc_threshold Log2FC filter threshold (default 0.1)
#' @param test_use Statistical test (default "wilcox")
#' @return data.frame in presto-like format (feature, group, auc, pval, padj, logFC)
bpcells_find_all_markers <- function(seurat_obj, celltype_col,
                                      min_pct = 0.1, logfc_threshold = 0.1,
                                      test_use = "wilcox") {
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[celltype_col]]

  valid <- !is.na(Seurat::Idents(seurat_obj)) &
           Seurat::Idents(seurat_obj) != "NaN" &
           Seurat::Idents(seurat_obj) != ""
  if (sum(!valid) > 0) {
    message("  Filtering ", sum(!valid), " cells with NA/empty cell type")
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[valid])
  }

  message("  Running FindAllMarkers on ", ncol(seurat_obj), " cells, ",
          length(unique(Seurat::Idents(seurat_obj))), " cell types...")

  markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = FALSE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = test_use
  )

  # Convert to presto-like format
  result <- data.frame(
    feature = markers$gene,
    group = as.character(markers$cluster),
    auc = (markers$pct.1 - markers$pct.2 + 1) / 2,
    pval = markers$p_val,
    padj = markers$p_val_adj,
    logFC = markers$avg_log2FC,
    stringsAsFactors = FALSE
  )

  result <- result[!is.na(result$group) & result$group != "NaN", ]
  rownames(result) <- NULL
  message("  Found ", nrow(result), " markers across ", length(unique(result$group)), " groups")
  result
}

#' Read h5ad metadata using hdf5r (handles categoricals)
#'
#' Reads obs from h5ad without loading expression data.
#' Handles AnnData categorical encoding (codes + categories arrays).
#'
#' @param h5ad_path Path to h5ad file
#' @return data.frame with cell barcodes as rownames
read_h5ad_metadata <- function(h5ad_path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package required")
  }

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  on.exit(h5$close_all())

  obs <- h5[["obs"]]
  obs_names <- obs[["_index"]]$read()
  n_cells <- length(obs_names)

  col_names <- setdiff(names(obs), c("_index", "__categories"))
  meta_list <- list()
  failed_cols <- character(0)

  for (col in col_names) {
    tryCatch({
      item <- obs[[col]]
      if (inherits(item, "H5Group")) {
        codes <- item[["codes"]]$read()
        cats <- item[["categories"]]$read()
        vals <- rep(NA_character_, length(codes))
        valid <- codes >= 0
        vals[valid] <- cats[codes[valid] + 1L]
        meta_list[[col]] <- vals
      } else {
        vals <- item$read()
        if (length(vals) == n_cells) {
          meta_list[[col]] <- vals
        } else {
          failed_cols <- c(failed_cols, col)
        }
      }
    }, error = function(e) {
      failed_cols <<- c(failed_cols, col)
    })
  }

  if (length(failed_cols) > 0) {
    message("  Warning: Could not read columns via hdf5r: ",
            paste(failed_cols, collapse = ", "))
  }

  meta <- as.data.frame(meta_list, stringsAsFactors = FALSE)
  rownames(meta) <- obs_names
  meta
}

#' Validate BPCells Seurat RDS
#'
#' Checks that the BPCells on-disk directory exists alongside the RDS.
#' @param rds_path Path to Seurat RDS file
#' @return TRUE if valid, stops with error if not
validate_bpcells_rds <- function(rds_path) {
  if (!file.exists(rds_path)) stop("RDS file not found: ", rds_path)

  obj <- readRDS(rds_path)
  counts_mat <- obj[["RNA"]]$counts

  if (is_bpcells_matrix(counts_mat)) {
    # Check that the referenced directory exists
    bpcells_dir <- BPCells::matrix_dir(counts_mat)
    if (!dir.exists(bpcells_dir)) {
      stop("BPCells on-disk directory not found: ", bpcells_dir,
           "\nThe RDS references this directory. Both must be kept together.")
    }
    message("BPCells validation passed: ", bpcells_dir)
  }

  invisible(TRUE)
}
