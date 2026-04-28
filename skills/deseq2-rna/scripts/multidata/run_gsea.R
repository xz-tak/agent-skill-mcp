#!/usr/bin/env Rscript
# DEG Multi-Dataset Analysis - GSEA Runner
#
# This script runs Gene Set Enrichment Analysis (GSEA) on multiple DEG data
# sources and combines results into a single output file.
#
# Usage:
#   Rscript run_gsea.R --config config.json
#   Rscript run_gsea.R --deg_file deg.xlsx --label "TGFb_stim" --output_dir ./output
#
# The script supports:
# - MSigDB Hallmark gene sets
# - GO Biological Process gene sets
# - Custom GMT files
# - Multiple data sources with different formats

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(jsonlite)
  library(dplyr)
  library(readxl)
})

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Default output directory name (created in working directory)
DEFAULT_OUTPUT_DIR <- "deg_multidata_output"

# Default GSEA settings
DEFAULT_GSEA_SETTINGS <- list(
  seed = 42,
  pvalueCutoff = 1,
  minGSSize = 3,
  maxGSSize = Inf,
  nPermSimple = 10000
)


# =============================================================================
# GENE SET LOADING FUNCTIONS
# =============================================================================

load_hallmark_gene_sets <- function() {
  #' Load MSigDB Hallmark gene sets for Homo sapiens.
  #'
  #' @return Data frame with term (gene set name) and gene (gene symbol)
  message("  Loading MSigDB Hallmark gene sets...")

  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  gene_sets <- hallmark %>%
    select(gs_name, gene_symbol) %>%
    rename(term = gs_name, gene = gene_symbol)

  message(sprintf("    Loaded %d gene sets", length(unique(gene_sets$term))))
  return(gene_sets)
}


load_gobp_gene_sets <- function() {
  #' Load GO Biological Process gene sets for Homo sapiens.
  #'
  #' @return Data frame with term (GO term name) and gene (gene symbol)
  message("  Loading GO:BP gene sets...")

  gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  gene_sets <- gobp %>%
    select(gs_name, gene_symbol) %>%
    rename(term = gs_name, gene = gene_symbol)

  message(sprintf("    Loaded %d gene sets", length(unique(gene_sets$term))))
  return(gene_sets)
}


load_custom_gmt <- function(gmt_path) {
  #' Load custom gene sets from GMT file.
  #'
  #' @param gmt_path Path to GMT file
  #' @return Data frame with term and gene columns
  message(sprintf("  Loading custom GMT: %s", gmt_path))

  if (!file.exists(gmt_path)) {
    stop(sprintf("GMT file not found: %s", gmt_path))
  }

  gene_sets <- read.gmt(gmt_path)

  message(sprintf("    Loaded %d gene sets", length(unique(gene_sets$term))))
  return(gene_sets)
}


load_gene_sets <- function(config) {
  #' Load all configured gene sets.
  #'
  #' @param config Configuration list with gene_sets section
  #' @return Combined data frame with term and gene columns
  gene_sets_config <- config$gene_sets
  all_gene_sets <- data.frame()

  # Load Hallmark if requested
  if (isTRUE(gene_sets_config$hallmark)) {
    hallmark <- load_hallmark_gene_sets()
    all_gene_sets <- rbind(all_gene_sets, hallmark)
  }

  # Load GO:BP if requested
  if (isTRUE(gene_sets_config$gobp)) {
    gobp <- load_gobp_gene_sets()
    all_gene_sets <- rbind(all_gene_sets, gobp)
  }

  # Load custom GMT if provided
  if (!is.null(gene_sets_config$custom_gmt) && gene_sets_config$custom_gmt != "") {
    custom <- load_custom_gmt(gene_sets_config$custom_gmt)
    all_gene_sets <- rbind(all_gene_sets, custom)
  }

  if (nrow(all_gene_sets) == 0) {
    stop("No gene sets loaded. Enable hallmark, gobp, or provide custom_gmt.")
  }

  message(sprintf("\nTotal gene sets loaded: %d", length(unique(all_gene_sets$term))))
  return(all_gene_sets)
}


# =============================================================================
# DEG DATA LOADING FUNCTIONS
# =============================================================================

load_deg_data <- function(source_config) {
  #' Load DEG data from a single source.
  #'
  #' @param source_config Data source configuration
  #' @return Data frame with gene, log2fc, padj columns
  file_path <- source_config$file

  # Determine file type and load
  if (grepl("\\.xlsx$|\\.xls$", file_path, ignore.case = TRUE)) {
    sheet_name <- source_config$sheet
    if (is.null(sheet_name)) {
      df <- read_excel(file_path)
    } else {
      df <- read_excel(file_path, sheet = sheet_name)
    }
  } else if (grepl("\\.tsv$|\\.txt$", file_path, ignore.case = TRUE)) {
    df <- read.delim(file_path, stringsAsFactors = FALSE)
  } else {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
  }

  # Apply contrast filter if specified
  if (!is.null(source_config$contrast_filter)) {
    filter_col <- source_config$contrast_filter$column
    filter_val <- source_config$contrast_filter$value
    df <- df[grepl(filter_val, df[[filter_col]], ignore.case = TRUE), ]
  }

  # Extract relevant columns
  gene_col <- source_config$gene_col
  log2fc_col <- source_config$log2fc_col
  padj_col <- source_config$padj_col

  result <- data.frame(
    gene = df[[gene_col]],
    log2fc = as.numeric(df[[log2fc_col]]),
    padj = as.numeric(df[[padj_col]]),
    stringsAsFactors = FALSE
  )

  # Remove NA values
  result <- result[!is.na(result$gene) & !is.na(result$log2fc), ]

  return(result)
}


create_gene_ranking <- function(deg_data) {
  #' Create gene ranking for GSEA.
  #'
  #' Ranking formula: -log10(padj + 1e-300) * sign(log2fc)
  #'
  #' @param deg_data Data frame with gene, log2fc, padj columns
  #' @return Named numeric vector (genes as names, scores as values)

  # Handle missing padj values
  deg_data$padj[is.na(deg_data$padj)] <- 1

  # Calculate ranking score
  # Using -log10(padj) * sign(log2fc) to incorporate both significance and direction
  deg_data$score <- -log10(deg_data$padj + 1e-300) * sign(deg_data$log2fc)

  # Create named vector, sorted by score (descending)
  deg_data <- deg_data[order(deg_data$score, decreasing = TRUE), ]
  gene_list <- deg_data$score
  names(gene_list) <- deg_data$gene

  # Remove duplicates (keep first = highest score)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  return(gene_list)
}


# =============================================================================
# GSEA FUNCTIONS
# =============================================================================

run_gsea_single <- function(gene_list, gene_sets, gsea_settings, source_label) {
  #' Run GSEA on a single gene list.
  #'
  #' @param gene_list Named numeric vector (gene ranking)
  #' @param gene_sets Data frame with term and gene columns
  #' @param gsea_settings GSEA parameter settings
  #' @param source_label Label for this data source
  #' @return Data frame with GSEA results

  message(sprintf("\nRunning GSEA for: %s", source_label))
  message(sprintf("  Input genes: %d", length(gene_list)))

  # Set seed for reproducibility
  set.seed(gsea_settings$seed)

  # Convert gene sets to TERM2GENE format
  term2gene <- gene_sets[, c("term", "gene")]

  # Run GSEA
  gsea_result <- tryCatch({
    GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = gsea_settings$pvalueCutoff,
      minGSSize = gsea_settings$minGSSize,
      maxGSSize = gsea_settings$maxGSSize,
      nPermSimple = gsea_settings$nPermSimple,
      verbose = FALSE
    )
  }, error = function(e) {
    message(sprintf("  Warning: GSEA failed - %s", e$message))
    return(NULL)
  })

  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    message("  No significant enrichments found")
    return(data.frame())
  }

  # Extract results
  result_df <- gsea_result@result
  result_df$source <- source_label

  message(sprintf("  Found %d enriched pathways", nrow(result_df)))

  return(result_df)
}


run_gsea_all_sources <- function(config, gene_sets) {
  #' Run GSEA on all configured data sources.
  #'
  #' @param config Full configuration
  #' @param gene_sets Combined gene sets data frame
  #' @return Combined data frame with all GSEA results

  gsea_settings <- modifyList(DEFAULT_GSEA_SETTINGS, config$gsea_settings %||% list())

  all_results <- data.frame()

  for (source_config in config$data_sources) {
    label <- source_config$label

    # Load DEG data
    message(sprintf("\nLoading DEG data: %s", label))
    deg_data <- load_deg_data(source_config)
    message(sprintf("  Loaded %d genes", nrow(deg_data)))

    # Create gene ranking
    gene_list <- create_gene_ranking(deg_data)

    # Run GSEA
    gsea_result <- run_gsea_single(gene_list, gene_sets, gsea_settings, label)

    if (nrow(gsea_result) > 0) {
      all_results <- rbind(all_results, gsea_result)
    }
  }

  return(all_results)
}


# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

save_gsea_results <- function(results, output_dir, prefix) {
  #' Save GSEA results to file.
  #'
  #' @param results Combined GSEA results data frame
  #' @param output_dir Output directory
  #' @param prefix Output file prefix

  # Create gsea subdirectory
  gsea_dir <- file.path(output_dir, "gsea")
  if (!dir.exists(gsea_dir)) {
    dir.create(gsea_dir, recursive = TRUE)
  }

  # Save combined results
  output_file <- file.path(gsea_dir, paste0(prefix, "_gsea_all.txt"))
  write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("\nSaved GSEA results: %s", output_file))
  message(sprintf("  Total rows: %d", nrow(results)))
  message(sprintf("  Unique pathways: %d", length(unique(results$ID))))
  message(sprintf("  Data sources: %d", length(unique(results$source))))

  return(output_file)
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

load_config <- function(config_path) {
  #' Load configuration from JSON file.
  #'
  #' @param config_path Path to JSON config file
  #' @return Configuration list
  config <- fromJSON(config_path, simplifyVector = FALSE)
  return(config)
}


main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  config_path <- NULL
  output_dir <- NULL
  prefix <- NULL

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--config") {
      config_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output_dir") {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--prefix") {
      prefix <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # Validate config path
  if (is.null(config_path)) {
    stop("Must specify --config")
  }

  message("DEG Multi-Dataset Analysis - GSEA Runner")
  message("=" |> rep(50) |> paste(collapse = ""))

  # Load configuration
  config <- load_config(config_path)

  # Override output settings if specified
  if (is.null(output_dir)) output_dir <- config$output$directory
  if (is.null(prefix)) prefix <- config$output$prefix
  if (is.null(output_dir)) output_dir <- DEFAULT_OUTPUT_DIR
  if (is.null(prefix)) prefix <- "analysis"

  # Validate mode
  if (config$mode != "pathway_analysis") {
    stop("This script is for pathway_analysis mode only. Use generate_tables.py for gene-level analysis.")
  }

  # Load gene sets
  message("\nLoading gene sets...")
  gene_sets <- load_gene_sets(config)

  # Run GSEA on all sources
  message("\nRunning GSEA on all data sources...")
  all_results <- run_gsea_all_sources(config, gene_sets)

  if (nrow(all_results) == 0) {
    stop("No GSEA results generated. Check input data and gene sets.")
  }

  # Save results
  save_gsea_results(all_results, output_dir, prefix)

  message("\nDone!")
}


# Run main if script is executed directly
if (!interactive()) {
  main()
}
