#!/usr/bin/env Rscript
# =============================================================================
# Sample-Level Expression Analysis
# =============================================================================
# Main entry point for sample-level expression analysis across internal and
# external studies. Runs comparison boxplots, GSVA scoring, and correlations.
#
# Usage:
#   Rscript sample_level_analysis.R \
#     --expr-uri <S3_SOMA_URI> \
#     --deg-uri <S3_SOMA_URI> \
#     --target-name <name> \
#     --targets <GENE1,GENE2,...> \
#     --signatures <SigName:Gene1,Gene2;Sig2:GeneA,GeneB> \
#     --output-dir <path>/sample_level \
#     --modules comparison,gsva,corr \
#     --config-json <path_to_internal_study_configs.json> \
#     --extract-script <path_to_soma_expr_extract.py> \
#     --conda-env <r_env>
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(digest)
})

# --- Resolve script directory for sourcing helper files ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  return(getwd())
}

SCRIPT_DIR <- get_script_dir()

# Source helper modules
source(file.path(SCRIPT_DIR, "sla_plot_helpers.R"))
source(file.path(SCRIPT_DIR, "sla_study_dispatch.R"))
source(file.path(SCRIPT_DIR, "correlation_modules.R"))

# =============================================================================
# Argument Parsing
# =============================================================================
parse_args <- function() {
  parser <- ArgumentParser(description = "Sample-Level Expression Analysis")

  parser$add_argument("--expr-uri", required = TRUE,
                      help = "SOMA experiment URI for expression data")
  parser$add_argument("--deg-uri", required = TRUE,
                      help = "SOMA experiment URI for DEG statistics")
  parser$add_argument("--target-name", required = TRUE,
                      help = "Target name for labeling")
  parser$add_argument("--targets", required = TRUE,
                      help = "Comma-separated target gene symbols")
  parser$add_argument("--signatures", default = NULL,
                      help = "Signatures: SigName:Gene1,Gene2;Sig2:GeneA,GeneB")
  parser$add_argument("--output-dir", required = TRUE,
                      help = "Output directory path")
  parser$add_argument("--sources", default = "internal,curated,omicsoft",
                      help = "Comma-separated source categories to include")
  parser$add_argument("--per-sample-studies", default = NULL,
                      help = "Comma-separated external study IDs to include in per-sample analysis")
  parser$add_argument("--modules", default = "comparison,gsva,corr",
                      help = "Comma-separated modules to run")
  parser$add_argument("--config-json", default = NULL,
                      help = "Path to internal_study_configs.json (default: SCRIPT_DIR/../references/internal_study_configs.json)")
  parser$add_argument("--extract-script", default = NULL,
                      help = "Path to soma_expr_extract.py (default: SCRIPT_DIR/soma_expr_extract.py)")
  parser$add_argument("--slice-script", default = NULL,
                      help = "Path to slice_h5ad.py (file-based extraction)")
  parser$add_argument("--backend", default = "soma",
                      help = "Data backend: soma, h5ad, or auto (default: soma)")
  parser$add_argument("--conda-env", default = "r_env",
                      help = "Conda environment name (default: r_env)")
  parser$add_argument("--no-cache", action = "store_true", default = FALSE,
                      help = "Force re-extraction (skip cached data)")

  parser$parse_args()
}

# =============================================================================
# Parse Signatures
# =============================================================================
parse_signatures <- function(sig_str) {
  if (is.null(sig_str) || sig_str == "") return(list())

  sigs <- list()
  entries <- strsplit(sig_str, ";")[[1]]
  for (entry in entries) {
    entry <- trimws(entry)
    if (!grepl(":", entry)) next
    parts <- strsplit(entry, ":", fixed = TRUE)[[1]]
    sig_name <- trimws(parts[1])
    genes <- trimws(strsplit(parts[2], ",")[[1]])
    genes <- genes[genes != ""]
    if (length(genes) > 0) sigs[[sig_name]] <- genes
  }
  sigs
}

# =============================================================================
# Data Loading via Python Extraction Script
# Supports both pipe-based (soma_expr_extract.py) and file-based (slice_h5ad.py)
# =============================================================================
load_expr_data <- function(uri, genes, studies, extract_script, conda_env,
                           slice_script = NULL, backend = "soma") {
  gene_str <- paste(genes, collapse = ",")
  study_str <- if (!is.null(studies) && length(studies) > 0) {
    paste(studies, collapse = ",")
  } else NULL

  if (!is.null(slice_script) && file.exists(slice_script)) {
    output_file <- tempfile(fileext = ".tsv")
    cmd_parts <- c(
      "conda", "run", "-n", conda_env,
      "python3", slice_script,
      "--path", shQuote(uri),
      "--mode", "expr",
      "--backend", backend
    )
    if (!is.null(study_str)) {
      cmd_parts <- c(cmd_parts, "--study", shQuote(study_str))
    }
    cmd_parts <- c(cmd_parts, "--output", shQuote(output_file))
    cmd <- paste(cmd_parts, collapse = " ")
    message(paste("[Data] Loading EXPR via slice_h5ad:", substr(cmd, 1, 200), "..."))
    status <- system(cmd, intern = FALSE)
    if (status != 0 || !file.exists(output_file)) {
      message("ERROR: slice_h5ad.py extraction failed")
      return(NULL)
    }
    df <- read.delim(output_file, stringsAsFactors = FALSE, row.names = 1)
    unlink(output_file)
  } else {
    output_file <- tempfile(fileext = ".tsv")
    cmd_parts <- c(
      "conda", "run", "-n", conda_env,
      "python3", extract_script,
      "--uri", shQuote(uri),
      "--mode", "expr",
      "--genes", "ALL"
    )
    if (!is.null(study_str)) {
      cmd_parts <- c(cmd_parts, "--per-sample-studies", shQuote(study_str))
    }
    cmd <- paste(c(cmd_parts, ">", shQuote(output_file)), collapse = " ")
    message(paste("[Data] Loading EXPR via:", substr(cmd, 1, 200), "..."))
    status <- system(cmd, intern = FALSE)
    df <- NULL
    if (status == 0 && file.exists(output_file) && file.size(output_file) > 0) {
      df <- tryCatch(
        read.delim(output_file, stringsAsFactors = FALSE),
        error = function(e) {
          message(paste("ERROR loading EXPR data:", e$message))
          NULL
        }
      )
    } else {
      message("ERROR: EXPR extraction failed or produced empty output")
    }
    unlink(output_file)
  }

  if (is.null(df) || nrow(df) == 0) {
    message("ERROR: No expression data loaded")
    return(NULL)
  }

  message(paste("[Data] EXPR loaded:", nrow(df), "samples x", ncol(df), "columns"))
  df
}

load_deg_data <- function(uri, genes, studies, extract_script, conda_env,
                          slice_script = NULL, backend = "soma") {
  gene_str <- paste(genes, collapse = ",")
  study_str <- if (!is.null(studies) && length(studies) > 0) {
    paste(studies, collapse = ",")
  } else NULL

  if (!is.null(slice_script) && file.exists(slice_script)) {
    output_file <- tempfile(fileext = ".tsv")
    cmd_parts <- c(
      "conda", "run", "-n", conda_env,
      "python3", slice_script,
      "--path", shQuote(uri),
      "--mode", "deg",
      "--backend", backend,
      "--genes", shQuote(gene_str)
    )
    if (!is.null(study_str)) {
      cmd_parts <- c(cmd_parts, "--study", shQuote(study_str))
    }
    cmd_parts <- c(cmd_parts, "--output", shQuote(output_file))
    cmd <- paste(cmd_parts, collapse = " ")
    message(paste("[Data] Loading DEG via slice_h5ad:", substr(cmd, 1, 200), "..."))
    status <- system(cmd, intern = FALSE)
    if (status != 0 || !file.exists(output_file)) {
      message("WARNING: DEG slice_h5ad extraction failed")
      return(NULL)
    }
    df <- read.delim(output_file, stringsAsFactors = FALSE)
    unlink(output_file)
  } else {
    cmd_parts <- c(
      "conda", "run", "-n", conda_env,
      "python3", extract_script,
      "--uri", shQuote(uri),
      "--mode", "deg",
      "--genes", shQuote(gene_str)
    )
    if (!is.null(study_str)) {
      cmd_parts <- c(cmd_parts, "--per-sample-studies", shQuote(study_str))
    }
    cmd <- paste(cmd_parts, collapse = " ")
    message(paste("[Data] Loading DEG via:", substr(cmd, 1, 200), "..."))
    df <- tryCatch(
      read.delim(pipe(cmd), stringsAsFactors = FALSE),
      error = function(e) {
        message(paste("WARNING: DEG loading failed:", e$message))
        NULL
      }
    )
  }

  if (!is.null(df)) {
    message(paste("[Data] DEG loaded:", nrow(df), "comparisons x", ncol(df), "columns"))
  }
  df
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  args <- parse_args()

  # Auto-resolve config-json and extract-script from SCRIPT_DIR if not provided

  if (is.null(args$config_json)) {
    args$config_json <- file.path(dirname(SCRIPT_DIR), "references", "internal_study_configs.json")
  }
  if (is.null(args$extract_script)) {
    args$extract_script <- file.path(SCRIPT_DIR, "soma_expr_extract.py")
  }

  # Parse parameters
  target_genes <- trimws(strsplit(args$targets, ",")[[1]])
  target_genes <- target_genes[target_genes != ""]
  signatures <- parse_signatures(args$signatures)
  modules <- trimws(strsplit(args$modules, ",")[[1]])
  sources <- trimws(strsplit(args$sources, ",")[[1]])
  external_studies <- if (!is.null(args$per_sample_studies) && args$per_sample_studies != "") {
    trimws(strsplit(args$per_sample_studies, ",")[[1]])
  } else NULL

  message("==========================================================")
  message(paste("SAMPLE-LEVEL EXPRESSION ANALYSIS:", args$target_name))
  message("==========================================================")
  message(paste("  Targets:", paste(target_genes, collapse = ", ")))
  message(paste("  Signatures:", length(signatures)))
  message(paste("  Modules:", paste(modules, collapse = ", ")))
  message(paste("  Sources:", paste(sources, collapse = ", ")))
  message(paste("  Output:", args$output_dir))

  # Create output directory
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  no_cache <- isTRUE(args$no_cache)

  # Load internal study configs
  config <- fromJSON(args$config_json, simplifyVector = FALSE)
  internal_studies <- names(config)
  message(paste("  Internal studies:", paste(internal_studies, collapse = ", ")))

  # Collect all genes needed (targets + signature genes)
  all_genes <- unique(c(target_genes, unlist(signatures)))
  message(paste("  Total unique genes to query:", length(all_genes)))

  # =========================================================================
  # Cache strategy: digest-keyed in skill cache dir (persistent across runs)
  # Same key formula as target_query_plot.R for consistency
  # =========================================================================
  skill_cache_dir <- file.path(dirname(SCRIPT_DIR), "cache")
  dir.create(skill_cache_dir, recursive = TRUE, showWarnings = FALSE)

  # --no-cache: wipe entire skill cache dir and force fresh extraction
  if (no_cache) {
    existing_cache <- list.files(skill_cache_dir, full.names = TRUE)
    if (length(existing_cache) > 0) {
      message(paste("[Cache] --no-cache: removing", length(existing_cache), "cached files"))
      unlink(existing_cache)
    }
  }

  expr_cache_key <- digest(list(internal_studies, "ALL_GENES", NULL, "expr", args$backend))
  expr_cache_path <- file.path(skill_cache_dir, paste0("expr_", expr_cache_key, ".tsv"))

  deg_cache_key <- digest(list(internal_studies, NULL, NULL, target_genes, "deg", args$backend))
  deg_cache_path <- file.path(skill_cache_dir, paste0("deg_", deg_cache_key, ".tsv"))

  # =========================================================================
  # Load INTERNAL EXPR (with digest-keyed caching)
  # =========================================================================
  internal_expr_df <- NULL

  if (!no_cache && file.exists(expr_cache_path)) {
    message(paste("[Cache] Using skill cache:", expr_cache_path))
    internal_expr_df <- read.delim(expr_cache_path, stringsAsFactors = FALSE, row.names = 1)
    has_metadata <- any(grepl("^metadata_", colnames(internal_expr_df)))
    if (!has_metadata) {
      message("[Cache] Stale cache (no metadata_* columns). Re-extracting.")
      internal_expr_df <- NULL
    }
  }

  if (is.null(internal_expr_df)) {
    internal_expr_df <- load_expr_data(
      args$expr_uri, all_genes, internal_studies,
      args$extract_script, args$conda_env,
      slice_script = args$slice_script, backend = args$backend
    )
    if (!is.null(internal_expr_df) && nrow(internal_expr_df) > 0) {
      write.table(internal_expr_df, expr_cache_path,
                  sep = "\t", row.names = TRUE, quote = FALSE)
      message(paste("[Cache] Saved EXPR:", expr_cache_path))
    }
  }

  # =========================================================================
  # Load INTERNAL DEG (with digest-keyed caching)
  # =========================================================================
  internal_deg_df <- NULL

  if (!no_cache && file.exists(deg_cache_path)) {
    message(paste("[Cache] Using skill DEG cache:", deg_cache_path))
    internal_deg_df <- read.delim(deg_cache_path, stringsAsFactors = FALSE)
  }

  if (is.null(internal_deg_df)) {
    internal_deg_df <- load_deg_data(
      args$deg_uri, target_genes, internal_studies,
      args$extract_script, args$conda_env,
      slice_script = args$slice_script, backend = args$backend
    )
    if (!is.null(internal_deg_df) && nrow(internal_deg_df) > 0) {
      write.table(internal_deg_df, deg_cache_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      message(paste("[Cache] Saved DEG:", deg_cache_path))
    }
  }

  # =========================================================================
  # Load EXTERNAL studies (curated/omicsoft — NEVER cached)
  # =========================================================================
  external_expr_df <- NULL
  external_deg_df <- NULL
  if (!is.null(external_studies) && length(external_studies) > 0) {
    external_expr_df <- load_expr_data(
      args$expr_uri, all_genes, external_studies,
      args$extract_script, args$conda_env,
      slice_script = args$slice_script, backend = args$backend
    )
    external_deg_df <- load_deg_data(
      args$deg_uri, target_genes, external_studies,
      args$extract_script, args$conda_env,
      slice_script = args$slice_script, backend = args$backend
    )
  }

  # =========================================================================
  # Combine data and build expression matrix
  # Use bind_rows to preserve all metadata columns (fills missing with NA)
  # =========================================================================
  expr_df <- NULL
  if (!is.null(internal_expr_df) && nrow(internal_expr_df) > 0) {
    expr_df <- internal_expr_df
  }
  if (!is.null(external_expr_df) && nrow(external_expr_df) > 0) {
    if (!is.null(expr_df)) {
      expr_df <- dplyr::bind_rows(expr_df, external_expr_df)
    } else {
      expr_df <- external_expr_df
    }
  }

  if (is.null(expr_df) || nrow(expr_df) == 0) {
    message("FATAL: Could not load expression data. Exiting.")
    quit(status = 1)
  }

  deg_df <- NULL
  if (!is.null(internal_deg_df) && nrow(internal_deg_df) > 0) {
    deg_df <- internal_deg_df
  }
  if (!is.null(external_deg_df) && nrow(external_deg_df) > 0) {
    if (!is.null(deg_df)) {
      common_cols <- intersect(colnames(deg_df), colnames(external_deg_df))
      deg_df <- rbind(deg_df[, common_cols, drop = FALSE],
                      external_deg_df[, common_cols, drop = FALSE])
    } else {
      deg_df <- external_deg_df
    }
  }

  # Identify gene vs metadata columns
  # SOMA obs uses "metadata_" prefix; slice_h5ad uses "meta_" prefix — catch both
  meta_cols <- colnames(expr_df)[grepl("^meta_|^metadata_", colnames(expr_df))]
  gene_cols <- setdiff(colnames(expr_df), meta_cols)

  if (length(gene_cols) == 0) {
    message("FATAL: No gene expression columns found in data. Exiting.")
    quit(status = 1)
  }

  # Set row identifiers (rownames may already be set from row.names=1 read)
  if ("sample_id" %in% colnames(expr_df)) {
    rownames(expr_df) <- make.unique(expr_df$sample_id)
  } else if (all(rownames(expr_df) == as.character(seq_len(nrow(expr_df))))) {
    rownames(expr_df) <- paste0("S", seq_len(nrow(expr_df)))
  }

  expr_mat <- t(as.matrix(expr_df[, gene_cols, drop = FALSE]))

  message(paste("[Matrix] Expression:", nrow(expr_mat), "genes x",
                ncol(expr_mat), "samples"))

  # Determine project_id column
  pid_col <- if ("meta_project_id" %in% colnames(expr_df)) "meta_project_id"
             else if ("project_id" %in% colnames(expr_df)) "project_id"
             else if ("study" %in% colnames(expr_df)) "study"
             else NULL

  if (is.null(pid_col)) {
    message("FATAL: No project_id/study column found. Exiting.")
    quit(status = 1)
  }

  # Initialize manifest
  manifest <- list(
    status = "success",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    studies_processed = character(0),
    studies_skipped = list(),
    modules_run = modules,
    warnings = character(0),
    package_versions = list(
      R = paste0(R.version$major, ".", R.version$minor),
      GSVA = if (requireNamespace("GSVA", quietly = TRUE))
        as.character(packageVersion("GSVA")) else "not installed",
      limma = if (requireNamespace("limma", quietly = TRUE))
        as.character(packageVersion("limma")) else "not installed"
    ),
    soma_uris = list(expr = args$expr_uri, deg = args$deg_uri)
  )

  # ==========================================================================
  # Process Internal Studies (ALWAYS run)
  # ==========================================================================
  message("\n\n=== INTERNAL STUDIES ===")

  for (study_name in internal_studies) {
    study_mask <- expr_df[[pid_col]] == study_name
    if (sum(study_mask, na.rm = TRUE) == 0) {
      msg <- paste(study_name, ": 0 samples in data")
      message(paste("  SKIP:", msg))
      manifest$studies_skipped[[length(manifest$studies_skipped) + 1]] <-
        list(name = study_name, reason = "0 samples after filter")
      manifest$warnings <- c(manifest$warnings, msg)
      next
    }

    study_df <- expr_df[study_mask, , drop = FALSE]
    study_config <- config[[study_name]]

    # Get DEG data for this study
    study_deg <- NULL
    if (!is.null(deg_df)) {
      deg_pid_col <- if ("project_id" %in% colnames(deg_df)) "project_id" else "study"
      if (deg_pid_col %in% colnames(deg_df)) {
        study_deg <- deg_df[deg_df[[deg_pid_col]] == study_name, , drop = FALSE]
      }
    }

    tryCatch({
      run_internal_study(
        study_name, study_df, expr_mat, study_deg,
        study_config, target_genes, signatures,
        gene_cols, args$output_dir
      )
      manifest$studies_processed <- c(manifest$studies_processed, study_name)
    }, error = function(e) {
      msg <- paste(study_name, ":", e$message)
      message(paste("  ERROR:", msg))
      manifest$studies_skipped[[length(manifest$studies_skipped) + 1]] <<-
        list(name = study_name, reason = e$message)
      manifest$warnings <<- c(manifest$warnings, msg)
    })
  }

  # ==========================================================================
  # Process External Studies (user-selected)
  # ==========================================================================
  if (!is.null(external_studies) && length(external_studies) > 0) {
    message("\n\n=== EXTERNAL STUDIES ===")

    for (study_id in external_studies) {
      study_mask <- expr_df[[pid_col]] == study_id
      if (sum(study_mask, na.rm = TRUE) == 0) {
        msg <- paste(study_id, ": 0 samples in data")
        message(paste("  SKIP:", msg))
        manifest$studies_skipped[[length(manifest$studies_skipped) + 1]] <-
          list(name = study_id, reason = "0 samples after filter")
        next
      }

      study_df <- expr_df[study_mask, , drop = FALSE]

      # Get DEG data for this study
      study_deg <- NULL
      if (!is.null(deg_df)) {
        deg_pid_col <- if ("project_id" %in% colnames(deg_df)) "project_id" else "study"
        if (deg_pid_col %in% colnames(deg_df)) {
          study_deg <- deg_df[deg_df[[deg_pid_col]] == study_id, , drop = FALSE]
        }
      }

      tryCatch({
        run_external_study(
          study_id, study_df, expr_mat, study_deg,
          target_genes, signatures, gene_cols, args$output_dir
        )
        manifest$studies_processed <- c(manifest$studies_processed, study_id)
      }, error = function(e) {
        msg <- paste(study_id, ":", e$message)
        message(paste("  ERROR:", msg))
        manifest$studies_skipped[[length(manifest$studies_skipped) + 1]] <<-
          list(name = study_id, reason = e$message)
        manifest$warnings <<- c(manifest$warnings, msg)
      })
    }
  }

  # ==========================================================================
  # Write Manifest
  # ==========================================================================
  manifest_path <- file.path(args$output_dir, "manifest.json")
  write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  message(paste("\n[Manifest] Saved:", manifest_path))

  # Summary
  message("\n==========================================================")
  message("SAMPLE-LEVEL ANALYSIS COMPLETE")
  message("==========================================================")
  message(paste("  Studies processed:", length(manifest$studies_processed)))
  message(paste("  Studies skipped:", length(manifest$studies_skipped)))
  message(paste("  Warnings:", length(manifest$warnings)))
  message(paste("  Output:", args$output_dir))
}

# Run
main()
