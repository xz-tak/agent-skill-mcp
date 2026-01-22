#!/usr/bin/env Rscript
# DESeq2 RNA-seq Analysis - Environment Validation Script
# Checks R environment, packages, and input files

# ==============================================================================
# Parse Arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Default values
counts_file <- NULL
metadata_file <- NULL
species <- "hsapiens"

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--counts" && i < length(args)) {
    counts_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--metadata" && i < length(args)) {
    metadata_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--species" && i < length(args)) {
    species <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# ==============================================================================
# Package Check Functions
# ==============================================================================

check_package <- function(pkg, verbose = TRUE) {
  result <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    version <- as.character(packageVersion(pkg))
    if (verbose) cat(sprintf("  [OK] %s (%s)\n", pkg, version))
    list(status = "OK", version = version)
  }, error = function(e) {
    if (verbose) cat(sprintf("  [MISSING] %s\n", pkg))
    list(status = "MISSING", version = NA)
  })
  return(result)
}

# ==============================================================================
# Main Validation
# ==============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("DESeq2 RNA-seq Environment Validation\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

# Track overall status
all_ok <- TRUE
issues <- c()

# ------------------------------------------------------------------------------
# 1. Check R Version
# ------------------------------------------------------------------------------

cat("1. R Version\n")
cat("   ", R.version.string, "\n\n")

# ------------------------------------------------------------------------------
# 2. Check Core Packages
# ------------------------------------------------------------------------------

cat("2. Core Packages\n")

core_packages <- c(
  "DESeq2",
  "clusterProfiler",
  "msigdbr",
  "gprofiler2",
  "enrichplot",
  "PCAtools",
  "ggplot2",
  "pheatmap",
  "glue",
  "dplyr",
  "GSVA",
  "BiocParallel",
  "limma",
  "tidyr",
  "tibble",
  "ggrepel",
  "jsonlite"
)

core_results <- sapply(core_packages, check_package, simplify = FALSE)
missing_core <- names(core_results)[sapply(core_results, function(x) x$status == "MISSING")]

if (length(missing_core) > 0) {
  all_ok <- FALSE
  issues <- c(issues, paste("Missing core packages:", paste(missing_core, collapse = ", ")))
}

cat("\n")

# ------------------------------------------------------------------------------
# 2b. Check Plotting Function Packages (for ad-hoc visualization)
# ------------------------------------------------------------------------------

cat("2b. Plotting Function Packages (optional, for ad-hoc visualization)\n")

plotting_packages <- c(
  "ComplexHeatmap",
  "ggsignif",
  "rstatix",
  "car",
  "RColorBrewer",
  "circlize",
  "ggpubr"
)

plotting_results <- sapply(plotting_packages, check_package, simplify = FALSE)
missing_plotting <- names(plotting_results)[sapply(plotting_results, function(x) x$status == "MISSING")]

if (length(missing_plotting) > 0) {
  cat("  [WARN] Some plotting packages missing - ad-hoc visualization functions may not work\n")
  issues <- c(issues, paste("Missing plotting packages (optional):", paste(missing_plotting, collapse = ", ")))
}

cat("\n")

# ------------------------------------------------------------------------------
# 3. Check Species-Specific Packages
# ------------------------------------------------------------------------------

cat("3. Species-Specific Packages\n")

species_packages <- list(
  hsapiens = "org.Hs.eg.db",
  mmusculus = "org.Mm.eg.db",
  rnorvegicus = "org.Rn.eg.db"
)

# Check requested species
if (species %in% names(species_packages)) {
  pkg <- species_packages[[species]]
  result <- check_package(pkg)
  if (result$status == "MISSING") {
    all_ok <- FALSE
    issues <- c(issues, paste("Missing species package for", species, ":", pkg))
  }
} else {
  cat("  [WARN] Unknown species:", species, "\n")
  issues <- c(issues, paste("Unknown species:", species))
}

# Also check other species packages availability
cat("  Other species:\n")
for (sp in names(species_packages)) {
  if (sp != species) {
    pkg <- species_packages[[sp]]
    check_package(pkg, verbose = TRUE)
  }
}

cat("\n")

# ------------------------------------------------------------------------------
# 4. Validate Counts File
# ------------------------------------------------------------------------------

if (!is.null(counts_file)) {
  cat("4. Counts File Validation\n")
  cat("   File:", counts_file, "\n")

  if (!file.exists(counts_file)) {
    cat("   [ERROR] File does not exist\n")
    all_ok <- FALSE
    issues <- c(issues, "Counts file does not exist")
  } else {
    # Detect delimiter
    first_line <- readLines(counts_file, n = 1)
    if (grepl("\t", first_line)) {
      delim <- "\t"
      cat("   Format: Tab-separated\n")
    } else if (grepl(",", first_line)) {
      delim <- ","
      cat("   Format: Comma-separated\n")
    } else {
      delim <- "\t"
      cat("   Format: Assuming tab-separated\n")
    }

    # Read file
    counts <- tryCatch({
      if (delim == ",") {
        read.csv(counts_file, row.names = 1, check.names = FALSE, nrows = 100)
      } else {
        read.table(counts_file, sep = delim, header = TRUE, row.names = 1,
                   check.names = FALSE, nrows = 100)
      }
    }, error = function(e) {
      cat("   [ERROR] Failed to read file:", e$message, "\n")
      NULL
    })

    if (!is.null(counts)) {
      # Check for Gene.name column and remove if present
      if ("Gene.name" %in% colnames(counts)) {
        counts <- counts[, colnames(counts) != "Gene.name"]
        cat("   Note: Gene.name column detected (will be used for gene symbols)\n")
      }

      # Check dimensions
      cat("   Genes (preview):", nrow(counts), "\n")
      cat("   Samples:", ncol(counts), "\n")

      # Check for non-negative integers
      numeric_cols <- sapply(counts, is.numeric)
      if (!all(numeric_cols)) {
        cat("   [WARN] Some columns are not numeric\n")
      }

      # Check for negative values
      if (any(counts[, numeric_cols] < 0, na.rm = TRUE)) {
        cat("   [ERROR] Negative values detected - counts must be non-negative\n")
        all_ok <- FALSE
        issues <- c(issues, "Counts file contains negative values")
      } else {
        cat("   [OK] All values are non-negative\n")
      }

      # Check for non-integer values
      counts_numeric <- as.matrix(counts[, numeric_cols])
      if (any(counts_numeric != round(counts_numeric), na.rm = TRUE)) {
        cat("   [WARN] Non-integer values detected - will be rounded\n")
      }

      # Store sample names for metadata check
      count_samples <- colnames(counts)
    }
  }
  cat("\n")
} else {
  cat("4. Counts File: Not provided\n\n")
}

# ------------------------------------------------------------------------------
# 5. Validate Metadata File
# ------------------------------------------------------------------------------

if (!is.null(metadata_file)) {
  cat("5. Metadata File Validation\n")
  cat("   File:", metadata_file, "\n")

  if (!file.exists(metadata_file)) {
    cat("   [ERROR] File does not exist\n")
    all_ok <- FALSE
    issues <- c(issues, "Metadata file does not exist")
  } else {
    # Read metadata
    metadata <- tryCatch({
      read.table(metadata_file, sep = "\t", header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      cat("   [ERROR] Failed to read file:", e$message, "\n")
      NULL
    })

    if (!is.null(metadata)) {
      cat("   Samples:", nrow(metadata), "\n")
      cat("   Columns:", ncol(metadata), "\n")
      cat("   Column names:", paste(colnames(metadata), collapse = ", "), "\n")

      # Check sample matching if counts were loaded
      if (exists("count_samples")) {
        # Try to find sample column
        sample_col <- NULL
        for (col in colnames(metadata)) {
          if (all(count_samples %in% gsub(" ", "", metadata[[col]]))) {
            sample_col <- col
            break
          }
        }

        if (is.null(sample_col)) {
          # Check first column or row names
          if (all(count_samples %in% gsub(" ", "", metadata[[1]]))) {
            sample_col <- colnames(metadata)[1]
          }
        }

        if (!is.null(sample_col)) {
          cat("   [OK] Sample column identified:", sample_col, "\n")
        } else {
          cat("   [WARN] Could not match samples to counts file\n")
          cat("   Counts samples:", paste(head(count_samples, 5), collapse = ", "), "...\n")
          cat("   Metadata first col:", paste(head(metadata[[1]], 5), collapse = ", "), "...\n")
        }
      }

      # Summarize columns
      cat("\n   Column Summary:\n")
      for (col in colnames(metadata)) {
        vals <- metadata[[col]]
        if (is.numeric(vals)) {
          cat(sprintf("   - %s: numeric, range [%.2f, %.2f]\n",
                      col, min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
        } else {
          unique_vals <- unique(vals)
          n_unique <- length(unique_vals)
          if (n_unique <= 5) {
            cat(sprintf("   - %s: categorical (%d levels: %s)\n",
                        col, n_unique, paste(unique_vals, collapse = ", ")))
          } else {
            cat(sprintf("   - %s: categorical (%d levels)\n", col, n_unique))
          }
        }
      }
    }
  }
  cat("\n")
} else {
  cat("5. Metadata File: Not provided\n\n")
}

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("=" , rep("=", 60), "\n", sep = "")
cat("Validation Summary\n")
cat("=" , rep("=", 60), "\n\n", sep = "")

if (all_ok && length(issues) == 0) {
  cat("[SUCCESS] All checks passed. Environment is ready for DESeq2 analysis.\n\n")
  quit(status = 0)
} else {
  cat("[ISSUES FOUND]\n\n")
  for (issue in issues) {
    cat("  - ", issue, "\n")
  }
  cat("\n")

  if (length(missing_core) > 0) {
    cat("To install missing packages, run:\n")
    cat("  conda install -n r_env -c conda-forge -c bioconda \\\n")
    for (pkg in missing_core) {
      pkg_conda <- tolower(gsub("\\.", "-", pkg))
      if (grepl("^org\\.", pkg)) {
        pkg_conda <- paste0("bioconductor-", pkg_conda)
      } else if (pkg %in% c("DESeq2", "clusterProfiler", "enrichplot", "PCAtools",
                            "GSVA", "limma", "BiocParallel")) {
        pkg_conda <- paste0("bioconductor-", tolower(pkg))
      } else {
        pkg_conda <- paste0("r-", tolower(pkg))
      }
      cat("    ", pkg_conda, " \\\n")
    }
    cat("\n")
  }

  quit(status = 1)
}
