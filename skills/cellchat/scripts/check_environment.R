#!/usr/bin/env Rscript
# =============================================================================
# CellChat Environment Check
# =============================================================================
# Validates R environment for CellChat analysis
# Usage: Rscript check_environment.R
# =============================================================================

cat("=============================================================================\n")
cat("CellChat Environment Check\n")
cat("=============================================================================\n")
cat(sprintf("Check time: %s\n\n", Sys.time()))

# Initialize status tracking
all_ok <- TRUE
warnings_list <- character()
errors_list <- character()

# -----------------------------------------------------------------------------
# 1. Check R version
# -----------------------------------------------------------------------------
cat("1. R Version\n")
cat(sprintf("   %s\n", R.version.string))
r_major <- as.numeric(R.version$major)
r_minor <- as.numeric(R.version$minor)
if (r_major < 4) {
  errors_list <- c(errors_list, "R version < 4.0 (CellChat requires R >= 4.0)")
  all_ok <- FALSE
}
cat("\n")

# -----------------------------------------------------------------------------
# 2. Check required packages
# -----------------------------------------------------------------------------
cat("2. Required Packages\n")

required_packages <- list(
  CellChat = list(min_ver = "2.0.0", critical = TRUE),
  Seurat = list(min_ver = "5.0.0", critical = TRUE),
  Matrix = list(min_ver = "1.5.0", critical = TRUE),
  reticulate = list(min_ver = "1.30", critical = FALSE),
  patchwork = list(min_ver = "1.0.0", critical = FALSE),
  ggplot2 = list(min_ver = "3.0.0", critical = FALSE),
  dplyr = list(min_ver = "1.0.0", critical = FALSE),
  png = list(min_ver = "0.1", critical = FALSE)
)

for (pkg in names(required_packages)) {
  info <- required_packages[[pkg]]
  installed <- requireNamespace(pkg, quietly = TRUE)

  if (installed) {
    ver <- as.character(packageVersion(pkg))
    status <- "[OK]"

    # Version check
    if (compareVersion(ver, info$min_ver) < 0) {
      status <- "[WARN]"
      warnings_list <- c(warnings_list, sprintf("%s version %s < recommended %s", pkg, ver, info$min_ver))
    }

    cat(sprintf("   %s %s: v%s\n", status, pkg, ver))
  } else {
    if (info$critical) {
      status <- "[ERROR]"
      errors_list <- c(errors_list, sprintf("%s: NOT INSTALLED (required)", pkg))
      all_ok <- FALSE
    } else {
      status <- "[WARN]"
      warnings_list <- c(warnings_list, sprintf("%s: NOT INSTALLED (optional)", pkg))
    }
    cat(sprintf("   %s %s: NOT INSTALLED\n", status, pkg))
  }
}
cat("\n")

# -----------------------------------------------------------------------------
# 3. Check CellChat specific features
# -----------------------------------------------------------------------------
cat("3. CellChat Features\n")

if (requireNamespace("CellChat", quietly = TRUE)) {
  library(CellChat, quietly = TRUE)

  # Check CellChatDB
  tryCatch({
    db_human <- CellChatDB.human
    db_mouse <- CellChatDB.mouse
    cat(sprintf("   [OK] CellChatDB.human: %d interactions\n", nrow(db_human$interaction)))
    cat(sprintf("   [OK] CellChatDB.mouse: %d interactions\n", nrow(db_mouse$interaction)))
  }, error = function(e) {
    cat(sprintf("   [ERROR] CellChatDB: %s\n", e$message))
    errors_list <<- c(errors_list, sprintf("CellChatDB error: %s", e$message))
    all_ok <<- FALSE
  })
} else {
  cat("   [SKIP] CellChat not installed\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 4. Check Python/reticulate for UMAP
# -----------------------------------------------------------------------------
cat("4. Python/UMAP (for similarity embedding)\n")

if (requireNamespace("reticulate", quietly = TRUE)) {
  library(reticulate, quietly = TRUE)

  # Check Python availability
  py_avail <- tryCatch(py_available(), error = function(e) FALSE)

  if (py_avail) {
    py_path <- tryCatch(py_config()$python, error = function(e) "unknown")
    cat(sprintf("   [OK] Python: %s\n", py_path))

    # Check UMAP
    umap_ok <- tryCatch({
      umap <- import("umap", delay_load = TRUE)
      TRUE
    }, error = function(e) FALSE)

    if (umap_ok) {
      cat("   [OK] umap-learn: available\n")
    } else {
      cat("   [WARN] umap-learn: NOT AVAILABLE\n")
      warnings_list <- c(warnings_list, "Python umap-learn not installed (needed for similarity embedding)")
    }
  } else {
    cat("   [WARN] Python: NOT AVAILABLE via reticulate\n")
    warnings_list <- c(warnings_list, "Python not available (needed for UMAP similarity plots)")
  }
} else {
  cat("   [SKIP] reticulate not installed\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 5. Check H5AD Conversion Support
# -----------------------------------------------------------------------------
cat("5. H5AD Conversion Support\n")

# Get script directory from command args (works with Rscript)
get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(dirname(script_path)))
  }
  # Fallback: try to get from source() context
  if (sys.nframe() > 0) {
    for (i in seq_len(sys.nframe())) {
      ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(ofile)) {
        return(normalizePath(dirname(ofile)))
      }
    }
  }
  return(getwd())
}

script_dir <- get_script_dir()
convert_script <- file.path(script_dir, "convert_h5ad_to_seurat.R")

if (file.exists(convert_script)) {
  cat(sprintf("   [OK] convert_h5ad_to_seurat.R: found (%s)\n", convert_script))
} else {
  cat("   [WARN] convert_h5ad_to_seurat.R: NOT FOUND\n")
  warnings_list <- c(warnings_list, "H5AD conversion script not found")
}

# Check anndata R package (needed for conversion)
if (requireNamespace("anndata", quietly = TRUE)) {
  cat(sprintf("   [OK] anndata R package: v%s\n", packageVersion("anndata")))
} else {
  cat("   [WARN] anndata R package: NOT INSTALLED\n")
  warnings_list <- c(warnings_list, "anndata R package not installed (needed for h5ad conversion)")
}
cat("\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("=============================================================================\n")
cat("SUMMARY\n")
cat("=============================================================================\n")

if (length(errors_list) > 0) {
  cat(sprintf("\nERRORS (%d):\n", length(errors_list)))
  for (err in errors_list) {
    cat(sprintf("  - %s\n", err))
  }
}

if (length(warnings_list) > 0) {
  cat(sprintf("\nWARNINGS (%d):\n", length(warnings_list)))
  for (warn in warnings_list) {
    cat(sprintf("  - %s\n", warn))
  }
}

if (all_ok && length(warnings_list) == 0) {
  cat("\n[OK] Environment ready for CellChat analysis!\n")
} else if (all_ok) {
  cat("\n[OK] Environment ready (with warnings)\n")
} else {
  cat("\n[ERROR] Environment NOT ready - please fix errors above\n")
  cat("\n")
  cat("=============================================================================\n")
  cat("INSTALLATION INSTRUCTIONS\n")
  cat("=============================================================================\n")

  # Check what needs to be installed
  cellchat_missing <- !requireNamespace("CellChat", quietly = TRUE)
  seurat_missing <- !requireNamespace("Seurat", quietly = TRUE)

  if (cellchat_missing || seurat_missing) {
    cat("\n--- Required R Packages ---\n\n")

    if (seurat_missing) {
      cat("Install Seurat v5:\n")
      cat("  install.packages('Seurat')\n\n")
    }

    if (cellchat_missing) {
      cat("Install CellChat:\n")
      cat("  # Install dependencies first\n")
      cat("  install.packages(c('NMF', 'circlize', 'ComplexHeatmap'))\n")
      cat("  BiocManager::install(c('Biobase', 'BiocGenerics', 'AnnotationDbi', 'GO.db'))\n")
      cat("  \n")
      cat("  # Install CellChat\n")
      cat("  devtools::install_github('jinworks/CellChat')\n\n")
    }
  }

  # Python UMAP
  umap_missing <- tryCatch({
    library(reticulate, quietly = TRUE)
    !py_module_available("umap")
  }, error = function(e) TRUE)

  if (umap_missing) {
    cat("--- Python UMAP (for similarity analysis) ---\n\n")
    cat("Install umap-learn:\n")
    cat("  pip install umap-learn\n")
    cat("  # or with conda:\n")
    cat("  conda install -c conda-forge umap-learn\n\n")
  }

  cat("After installation, run this check again:\n")
  cat("  Rscript check_environment.R\n")
}

cat("\n=============================================================================\n")

# Exit with appropriate code
if (all_ok) {
  quit(status = 0)
} else {
  quit(status = 1)
}
