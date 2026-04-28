#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 RNA-seq Analysis Script
# Parameterized script for differential expression analysis
# Supports command-line arguments or JSON config file
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(msigdbr)
  library(gprofiler2)
  library(enrichplot)
  library(PCAtools)
  library(ggplot2)
  library(pheatmap)
  library(glue)
  library(dplyr)
  library(GSVA)
  library(BiocParallel)
  library(limma)
  library(tidyr)
  library(tibble)
  library(ggrepel)
  library(jsonlite)
})

# Set global random seed for reproducibility
set.seed(42)

# ==============================================================================
# Parameter Parsing
# ==============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default values
  params <- list(
    counts = NULL,
    metadata = NULL,
    design = NULL,
    comparison_col = NULL,
    reference = NULL,
    comparisons = "all_pairwise",
    msigdb = "H",
    custom_sigs = NULL,
    species = "hsapiens",
    padj = 0.05,
    lfc = 1,
    output_dir = ".",
    prefix = "RNAseq",
    conda_env = "r_env",
    config = NULL
  )

  # Parse arguments
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--config" && i < length(args)) {
      params$config <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--counts" && i < length(args)) {
      params$counts <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--metadata" && i < length(args)) {
      params$metadata <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--design" && i < length(args)) {
      params$design <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--comparison_col" && i < length(args)) {
      params$comparison_col <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--reference" && i < length(args)) {
      params$reference <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--comparisons" && i < length(args)) {
      params$comparisons <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--advanced_comparisons" && i < length(args)) {
      params$advanced_comparisons <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--msigdb" && i < length(args)) {
      params$msigdb <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--custom_sigs" && i < length(args)) {
      params$custom_sigs <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--species" && i < length(args)) {
      params$species <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--padj" && i < length(args)) {
      params$padj <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--lfc" && i < length(args)) {
      params$lfc <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--output_dir" && i < length(args)) {
      params$output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--prefix" && i < length(args)) {
      params$prefix <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--conda_env" && i < length(args)) {
      params$conda_env <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # If config file provided, load and merge
  if (!is.null(params$config)) {
    config <- fromJSON(params$config)
    for (name in names(config)) {
      if (!is.null(config[[name]])) {
        params[[name]] <- config[[name]]
      }
    }
  }

  return(params)
}

# ==============================================================================
# Species Configuration
# ==============================================================================

get_species_config <- function(species) {
  configs <- list(
    hsapiens = list(
      org_db = "org.Hs.eg.db",
      msigdb_species = "Homo sapiens",
      gprofiler_organism = "hsapiens",
      keytypes = c("SYMBOL", "ENSEMBL", "ENTREZID")
    ),
    mmusculus = list(
      org_db = "org.Mm.eg.db",
      msigdb_species = "Mus musculus",
      gprofiler_organism = "mmusculus",
      keytypes = c("SYMBOL", "ENSEMBL", "ENTREZID")
    ),
    rnorvegicus = list(
      org_db = "org.Rn.eg.db",
      msigdb_species = "Rattus norvegicus",
      gprofiler_organism = "rnorvegicus",
      keytypes = c("SYMBOL", "ENSEMBL", "ENTREZID")
    )
  )

  if (!species %in% names(configs)) {
    stop(glue("Unknown species: {species}. Supported: hsapiens, mmusculus, rnorvegicus"))
  }

  config <- configs[[species]]
  suppressPackageStartupMessages(library(config$org_db, character.only = TRUE))
  config$org_db_obj <- get(config$org_db)

  return(config)
}

# ==============================================================================
# Signature Loading Functions
# ==============================================================================

read_gmt <- function(file_path) {
  lines <- readLines(file_path)
  gene_sets <- list()

  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      set_name <- parts[1]
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]
      gene_sets[[set_name]] <- genes
    }
  }

  return(gene_sets)
}

read_text_signatures <- function(file_path) {
  lines <- readLines(file_path)
  gene_sets <- list()
  current_set <- NULL
  current_genes <- c()

  # Check if file has any > headers (multi-set format)
  has_headers <- any(grepl("^>", lines))

  if (!has_headers) {
    # Single signature format: just genes, one per line
    genes <- c()
    for (line in lines) {
      line <- trimws(line)
      if (line != "" && !startsWith(line, "#")) {
        genes <- c(genes, line)
      }
    }
    if (length(genes) > 0) {
      gene_sets[["CUSTOM_SIGNATURE"]] <- genes
      cat(glue("    Single signature: {length(genes)} genes\n"))
    }
    return(gene_sets)
  }

  # Multi-set format with > headers
  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) next

    if (startsWith(line, ">")) {
      # Save previous set
      if (!is.null(current_set) && length(current_genes) > 0) {
        gene_sets[[current_set]] <- current_genes
      }
      # Start new set
      current_set <- gsub("^>", "", line)
      current_genes <- c()
    } else {
      current_genes <- c(current_genes, line)
    }
  }

  # Save last set
  if (!is.null(current_set) && length(current_genes) > 0) {
    gene_sets[[current_set]] <- current_genes
  }

  return(gene_sets)
}

load_custom_signatures <- function(file_path) {
  if (is.null(file_path) || file_path == "" || file_path == "none") {
    return(list(gene_sets = list(), sources = list()))
  }

  if (!file.exists(file_path)) {
    warning(glue("Custom signature file not found: {file_path}"))
    return(list(gene_sets = list(), sources = list()))
  }

  # Auto-detect format
  first_line <- readLines(file_path, n = 1)

  if (grepl("\t", first_line)) {
    cat("  Detected GMT format\n")
    gene_sets <- read_gmt(file_path)
  } else {
    cat("  Detected text format\n")
    gene_sets <- read_text_signatures(file_path)
  }

  # Add source tracking - all custom signatures get "CUSTOM" source
  sources <- setNames(rep("CUSTOM", length(gene_sets)), names(gene_sets))

  return(list(gene_sets = gene_sets, sources = as.list(sources)))
}

# ==============================================================================
# Data Loading
# ==============================================================================

load_counts <- function(counts_path) {
  # Detect delimiter
  first_line <- readLines(counts_path, n = 1)
  if (grepl("\t", first_line)) {
    counts <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1,
                         check.names = FALSE)
  } else {
    counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  }

  # Handle Gene.name column if present
  gene_symbols <- NULL
  if ("Gene.name" %in% colnames(counts)) {
    gene_symbols <- counts$Gene.name
    names(gene_symbols) <- rownames(counts)
    counts <- counts[, colnames(counts) != "Gene.name"]
  }

  # Convert to integer matrix
  counts <- as.matrix(counts)
  mode(counts) <- "integer"

  return(list(counts = counts, gene_symbols = gene_symbols))
}

load_metadata <- function(metadata_path, counts) {
  metadata <- read.table(metadata_path, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)

  # Find sample column
  sample_col <- NULL
  for (col in colnames(metadata)) {
    if (all(colnames(counts) %in% gsub(" ", "", metadata[[col]]))) {
      sample_col <- col
      break
    }
  }

  if (is.null(sample_col)) {
    # Try first column
    if (all(colnames(counts) %in% gsub(" ", "", metadata[[1]]))) {
      sample_col <- colnames(metadata)[1]
    } else {
      stop("Could not match metadata samples to counts columns")
    }
  }

  # Set row names and order
  rownames(metadata) <- gsub(" ", "", metadata[[sample_col]])
  metadata <- metadata[colnames(counts), ]

  return(metadata)
}

# ==============================================================================
# MSigDB Loading
# ==============================================================================

load_msigdb_collections <- function(collections, species_config) {
  gene_sets <- list()
  gene_set_sources <- list()  # Track source per gene set

  collections <- strsplit(collections, ",")[[1]]
  collections <- trimws(collections)

  for (coll in collections) {
    cat(glue("  Loading MSigDB {coll}...\n"))

    tryCatch({
      if (coll == "H") {
        msig <- msigdbr(species = species_config$msigdb_species, category = "H")
      } else if (startsWith(coll, "C")) {
        msig <- msigdbr(species = species_config$msigdb_species, category = coll)
      } else {
        warning(glue("Unknown MSigDB collection: {coll}"))
        next
      }

      # Convert to list format
      sets <- split(msig$gene_symbol, msig$gs_name)
      gene_sets <- c(gene_sets, sets)

      # Track source for each gene set
      for (set_name in names(sets)) {
        gene_set_sources[[set_name]] <- coll
      }

      cat(glue("    Loaded {length(sets)} gene sets\n"))
    }, error = function(e) {
      warning(glue("Failed to load MSigDB {coll}: {e$message}"))
    })
  }

  return(list(gene_sets = gene_sets, sources = gene_set_sources))
}

# ==============================================================================
# GSVA Analysis
# ==============================================================================

run_gsva <- function(counts, gene_sets, method = "gsva") {
  cat("Running GSVA...\n")

  # Use GSVA 2.0+ API
  gsva_param <- gsvaParam(
    exprData = counts,
    geneSets = gene_sets,
    kcdf = "Poisson",
    maxDiff = TRUE
  )

  gsva_scores <- gsva(gsva_param, verbose = FALSE)

  return(gsva_scores)
}

# ==============================================================================
# DESeq2 Analysis
# ==============================================================================

run_deseq2 <- function(counts, coldata, design_formula, comparison_col, reference) {
  cat("Running DESeq2...\n")

  # Set reference level
  coldata[[comparison_col]] <- factor(coldata[[comparison_col]])
  coldata[[comparison_col]] <- relevel(coldata[[comparison_col]], ref = reference)

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = as.formula(design_formula)
  )

  # Filter low counts
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Run DESeq2
  dds <- DESeq(dds)

  # Get VST for visualization
  vsd <- vst(dds, blind = FALSE)

  return(list(dds = dds, vsd = vsd))
}

# ==============================================================================
# Extract Results
# ==============================================================================

parse_comparisons <- function(comparisons_param, comparison_col, coldata) {
  # If already a list/matrix (from JSON config), convert to list of pairs
  if (is.list(comparisons_param) || is.matrix(comparisons_param)) {
    if (is.matrix(comparisons_param)) {
      # Convert matrix to list of row vectors
      pairs <- lapply(1:nrow(comparisons_param), function(i) comparisons_param[i, ])
    } else {
      pairs <- comparisons_param
    }
    return(pairs)
  }

  # String input (from command line)
  if (is.character(comparisons_param) && length(comparisons_param) == 1) {
    if (comparisons_param == "all_pairwise") {
      levels <- levels(coldata[[comparison_col]])
      pairs <- list()
      for (i in 1:(length(levels) - 1)) {
        for (j in (i + 1):length(levels)) {
          pairs[[length(pairs) + 1]] <- c(levels[j], levels[i])
        }
      }
      return(pairs)
    }

    # Parse JSON string
    parsed <- fromJSON(comparisons_param)
    if (is.matrix(parsed)) {
      return(lapply(1:nrow(parsed), function(i) parsed[i, ]))
    }
    return(parsed)
  }

  stop("Invalid comparisons format. Expected 'all_pairwise', JSON string, or list of pairs.")
}

#' Parse advanced comparisons (list contrasts) from config
#'
#' @param advanced_param Advanced comparisons from config (list or JSON string)
#' @return List of advanced comparison specifications, each with: name, contrast, listValues
parse_advanced_comparisons <- function(advanced_param) {
  if (is.null(advanced_param) || length(advanced_param) == 0) {
    return(list())
  }

  # If JSON string, parse it
  # CRITICAL: simplifyVector=FALSE is required to preserve nested list structure

  # Without this, fromJSON() converts listValues arrays into vectors/matrices,
  # which breaks DESeq2's list contrast mechanism for interaction comparisons.
  # Example: {"listValues": [1, -1]} must remain as list(1, -1), not c(1, -1)
  # This fix enables interaction contrasts like: (TGFb+ALK5i vs TGFb) vs (TGFb vs Non-treatment)
  if (is.character(advanced_param) && length(advanced_param) == 1) {
    advanced_param <- fromJSON(advanced_param, simplifyVector = FALSE)
  }

  # Ensure it's a list of comparisons (handle single comparison case)
  if (!is.null(names(advanced_param)) && "name" %in% names(advanced_param)) {
    # Single comparison passed as object, wrap in list
    advanced_param <- list(advanced_param)
  }

  # Validate each comparison has required fields
  validated <- list()
  for (i in seq_along(advanced_param)) {
    comp <- advanced_param[[i]]
    if (is.null(comp$name) || is.null(comp$contrast) || is.null(comp$listValues)) {
      warning(paste("Advanced comparison", i, "missing required fields (name, contrast, listValues). Skipping."))
      next
    }
    validated[[length(validated) + 1]] <- list(
      name = comp$name,
      contrast = unlist(comp$contrast),
      listValues = as.numeric(unlist(comp$listValues))
    )
  }

  return(validated)
}

#' Print available coefficient names from DESeq2 model
#'
#' @param dds DESeqDataSet object after running DESeq()
print_available_coefficients <- function(dds) {
  coefs <- resultsNames(dds)
  cat("\n=== Available Model Coefficients ===\n")
  cat("Use these names for list contrast specifications:\n\n")
  for (i in seq_along(coefs)) {
    cat(sprintf("  [%d] %s\n", i, coefs[i]))
  }
  cat("\n")
  return(invisible(coefs))
}

extract_results <- function(dds, comparison_col, numerator, denominator, gene_symbols = NULL) {
  contrast <- c(comparison_col, numerator, denominator)

  res <- results(dds, contrast = contrast)
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  # Add gene symbols if available
  if (!is.null(gene_symbols)) {
    res$symbol <- gene_symbols[res$gene]
  } else {
    res$symbol <- res$gene
  }

  # Create comparison name
  comp_name <- paste(
    gsub("[^A-Za-z0-9]", "_", numerator),
    "vs",
    gsub("[^A-Za-z0-9]", "_", denominator),
    sep = "_"
  )

  res$comparison <- comp_name

  return(list(results = res, name = comp_name))
}

#' Extract DESeq2 results using list contrast (for interaction/complex comparisons)
#'
#' @param dds DESeqDataSet object
#' @param contrast_spec List specifying the contrast coefficients
#' @param list_values Numeric vector of weights for the contrast
#' @param comp_name Name for this comparison
#' @param gene_symbols Optional named vector of gene symbols
#' @return List with results data frame and comparison name
extract_results_list <- function(dds, contrast_spec, list_values, comp_name, gene_symbols = NULL) {
  # Convert contrast_spec to proper R list format
  # contrast_spec can be:
  #   - ["coef1", "coef2"] -> list("coef1", "coef2")
  #   - [["coef1", "coef2"], "coef3"] -> list(c("coef1", "coef2"), "coef3")

  # Build the contrast list
  contrast_list <- lapply(contrast_spec, function(x) {
    if (is.list(x) || (is.character(x) && length(x) > 1)) {
      return(as.character(x))  # Vector of coefficient names
    } else {
      return(as.character(x))  # Single coefficient name
    }
  })

  # Validate coefficient names exist
  available_coefs <- resultsNames(dds)
  all_coefs <- unlist(contrast_spec)
  missing_coefs <- setdiff(all_coefs, available_coefs)
  if (length(missing_coefs) > 0) {
    stop(paste0("Coefficient(s) not found in model: ", paste(missing_coefs, collapse = ", "),
                "\nAvailable coefficients: ", paste(available_coefs, collapse = ", ")))
  }

  # Run DESeq2 results with list contrast
  res <- results(dds, contrast = contrast_list, listValues = list_values)
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  # Add gene symbols if available
  if (!is.null(gene_symbols)) {
    res$symbol <- gene_symbols[res$gene]
  } else {
    res$symbol <- res$gene
  }

  res$comparison <- comp_name

  return(list(results = res, name = comp_name))
}

# ==============================================================================
# Enrichment Analysis
# ==============================================================================

run_gsea_gobp <- function(res, species_config, pval_cutoff = 1) {
  # Create ranked gene list
  res_clean <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

  gene_list <- -log10(res_clean$padj + 1e-300) * sign(res_clean$log2FoldChange)
  names(gene_list) <- res_clean$symbol
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Remove duplicates
  gene_list <- gene_list[!duplicated(names(gene_list))]

  tryCatch({
    gsea_res <- gseGO(
      geneList = gene_list,
      OrgDb = species_config$org_db_obj,
      ont = "BP",
      keyType = "SYMBOL",
      pvalueCutoff = pval_cutoff,
      minGSSize = 3,
      maxGSSize = Inf,
      verbose = FALSE
    )
    return(gsea_res)
  }, error = function(e) {
    warning(glue("GSEA failed: {e$message}"))
    return(NULL)
  })
}

run_gprofiler <- function(genes_up, genes_dn, species_config) {
  results <- list()

  if (length(genes_up) > 0) {
    tryCatch({
      gp_up <- gost(
        query = genes_up,
        organism = species_config$gprofiler_organism,
        sources = c("GO:BP", "KEGG", "REAC"),
        significant = TRUE,
        user_threshold = 0.05
      )
      if (!is.null(gp_up$result)) {
        gp_up$result$DE <- "up"
        results$up <- gp_up$result
      }
    }, error = function(e) {
      warning(glue("g:Profiler (up) failed: {e$message}"))
    })
  }

  if (length(genes_dn) > 0) {
    tryCatch({
      gp_dn <- gost(
        query = genes_dn,
        organism = species_config$gprofiler_organism,
        sources = c("GO:BP", "KEGG", "REAC"),
        significant = TRUE,
        user_threshold = 0.05
      )
      if (!is.null(gp_dn$result)) {
        gp_dn$result$DE <- "dn"
        results$dn <- gp_dn$result
      }
    }, error = function(e) {
      warning(glue("g:Profiler (dn) failed: {e$message}"))
    })
  }

  return(results)
}

# Helper function to flatten list columns for write.table
flatten_list_columns <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        paste(x, collapse = ";")
      })
    }
  }
  return(df)
}

# ==============================================================================
# GSEA with Custom Gene Sets (MSigDB + Custom signatures)
# ==============================================================================

run_gsea_custom <- function(res, gene_sets, gene_set_sources, pval_cutoff = 1) {
  if (length(gene_sets) == 0) {
    return(NULL)
  }

  # Create ranked gene list (same as run_gsea_gobp)
  res_clean <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  gene_list <- -log10(res_clean$padj + 1e-300) * sign(res_clean$log2FoldChange)
  names(gene_list) <- res_clean$symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # Create TERM2GENE data frame from gene sets list
  term2gene <- do.call(rbind, lapply(names(gene_sets), function(term) {
    data.frame(term = term, gene = gene_sets[[term]], stringsAsFactors = FALSE)
  }))

  tryCatch({
    gsea_res <- GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = pval_cutoff,
      minGSSize = 3,
      maxGSSize = Inf,
      verbose = FALSE
    )

    # Add source column from gene_set_sources
    if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
      gsea_res@result$source <- sapply(gsea_res@result$ID, function(id) {
        if (id %in% names(gene_set_sources)) gene_set_sources[[id]] else "UNKNOWN"
      })
    }

    return(gsea_res)
  }, error = function(e) {
    warning(glue("GSEA (custom) failed: {e$message}"))
    return(NULL)
  })
}

# ==============================================================================
# ORA with Custom Gene Sets (MSigDB + Custom signatures)
# ==============================================================================

run_ora_custom <- function(genes_up, genes_dn, gene_sets, gene_set_sources) {
  if (length(gene_sets) == 0) {
    return(list())
  }

  results <- list()

  # Create TERM2GENE data frame from gene sets list
  term2gene <- do.call(rbind, lapply(names(gene_sets), function(term) {
    data.frame(term = term, gene = gene_sets[[term]], stringsAsFactors = FALSE)
  }))

  # Upregulated genes
  if (length(genes_up) > 0) {
    tryCatch({
      ora_up <- enricher(
        gene = genes_up,
        TERM2GENE = term2gene,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
        ora_up@result$DE <- "up"
        ora_up@result$source <- sapply(ora_up@result$ID, function(id) {
          if (id %in% names(gene_set_sources)) gene_set_sources[[id]] else "UNKNOWN"
        })
        results$up <- ora_up@result
      }
    }, error = function(e) {
      warning(glue("ORA (up) failed: {e$message}"))
    })
  }

  # Downregulated genes
  if (length(genes_dn) > 0) {
    tryCatch({
      ora_dn <- enricher(
        gene = genes_dn,
        TERM2GENE = term2gene,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
      if (!is.null(ora_dn) && nrow(ora_dn@result) > 0) {
        ora_dn@result$DE <- "dn"
        ora_dn@result$source <- sapply(ora_dn@result$ID, function(id) {
          if (id %in% names(gene_set_sources)) gene_set_sources[[id]] else "UNKNOWN"
        })
        results$dn <- ora_dn@result
      }
    }, error = function(e) {
      warning(glue("ORA (dn) failed: {e$message}"))
    })
  }

  return(results)
}

# ==============================================================================
# Harmonize g:Profiler and enricher Output Columns
# ==============================================================================

harmonize_ora_results <- function(gost_results, enricher_results) {
  harmonized_list <- list()

  # Process g:Profiler (gost) results
  if (length(gost_results) > 0) {
    for (name in names(gost_results)) {
      df <- gost_results[[name]]
      if (!is.null(df) && nrow(df) > 0) {
        harmonized <- data.frame(
          term_id = df$term_id,
          term_name = df$term_name,
          p_value = df$p_value,
          p_adjust = df$p_value,  # gost p_value is already adjusted
          source = df$source,
          DE = df$DE,
          genes = sapply(df$intersection, paste, collapse = ";"),
          gene_count = df$intersection_size,
          term_size = df$term_size,
          stringsAsFactors = FALSE
        )
        harmonized_list[[paste0("gost_", name)]] <- harmonized
      }
    }
  }

  # Process enricher results
  if (length(enricher_results) > 0) {
    for (name in names(enricher_results)) {
      df <- enricher_results[[name]]
      if (!is.null(df) && nrow(df) > 0) {
        # Parse GeneRatio to get term_size (format: "count/total")
        gene_ratio_parts <- strsplit(df$GeneRatio, "/")
        term_sizes <- sapply(gene_ratio_parts, function(x) as.numeric(x[2]))

        harmonized <- data.frame(
          term_id = df$ID,
          term_name = df$Description,
          p_value = df$pvalue,
          p_adjust = df$p.adjust,
          source = df$source,
          DE = df$DE,
          genes = gsub("/", ";", df$geneID),  # Convert / separator to ;
          gene_count = df$Count,
          term_size = term_sizes,
          stringsAsFactors = FALSE
        )
        harmonized_list[[paste0("enricher_", name)]] <- harmonized
      }
    }
  }

  # Combine all
  if (length(harmonized_list) > 0) {
    return(do.call(rbind, harmonized_list))
  } else {
    return(NULL)
  }
}

# ==============================================================================
# Plotting Functions
# ==============================================================================

plot_volcano <- function(res, padj_thresh, lfc_thresh, output_path) {
  res_plot <- res[!is.na(res$padj), ]

  res_plot$significance <- "NS"
  res_plot$significance[res_plot$padj < padj_thresh & res_plot$log2FoldChange > lfc_thresh] <- "Up"
  res_plot$significance[res_plot$padj < padj_thresh & res_plot$log2FoldChange < -lfc_thresh] <- "Down"

  # Label top genes
  res_plot$label <- ""
  top_up <- res_plot %>% filter(significance == "Up") %>% arrange(padj) %>% head(10)
  top_dn <- res_plot %>% filter(significance == "Down") %>% arrange(padj) %>% head(10)
  res_plot$label[res_plot$gene %in% c(top_up$gene, top_dn$gene)] <-
    res_plot$symbol[res_plot$gene %in% c(top_up$gene, top_dn$gene)]

  p <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey60")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
    labs(
      title = res$comparison[1],
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(output_path, p, width = 10, height = 8, dpi = 720)
}

plot_gsea_dotplot <- function(gsea_res, output_path, n_show = 10) {
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    return(NULL)
  }

  tryCatch({
    p <- dotplot(gsea_res, showCategory = n_show, split = ".sign") +
      facet_grid(~.sign) +
      theme_bw()

    ggsave(output_path, p, width = 12, height = 8, dpi = 720)
  }, error = function(e) {
    warning(glue("Failed to create GSEA dotplot: {e$message}"))
  })
}

plot_pca <- function(vsd, coldata, color_by, output_path) {
  # Run PCA
  pca_data <- plotPCA(vsd, intgroup = color_by, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = paste0("PCA - ", color_by)
    ) +
    theme_bw() +
    theme(legend.position = "right")

  ggsave(output_path, p, width = 10, height = 8, dpi = 720)
}

plot_pca_advanced <- function(vsd, coldata, custom_gsva_scores, figures_dir, prefix) {
  # PCAtools analysis
  vst_mat <- assay(vsd)

  # Remove zero-variance genes
  vars <- apply(vst_mat, 1, var)
  vst_mat <- vst_mat[vars > 0, ]

  pca_obj <- pca(vst_mat, metadata = coldata, removeVar = 0.1)

  # Scree plot
  png(file.path(figures_dir, paste0(prefix, "_PCA_screeplot.png")),
      width = 10, height = 8, units = "in", res = 720)
  print(screeplot(pca_obj, components = 1:10))
  dev.off()

  # Eigencorplot with metadata + CUSTOM signatures only (not MSigDB)
  # This keeps the plot manageable and focused on user-defined signatures
  if (!is.null(custom_gsva_scores) && nrow(custom_gsva_scores) > 0) {
    meta_combined <- cbind(coldata, t(custom_gsva_scores[, rownames(coldata)]))
  } else {
    meta_combined <- coldata
  }

  # Filter metadata for eigencorplot:
  # - Exclude ONLY columns with single unique value (no variance - causes correlation error)
  # - Convert non-numeric columns to numeric factors for correlation
  meta_filtered <- data.frame(row.names = rownames(meta_combined))

  for (col in colnames(meta_combined)) {
    vals <- meta_combined[[col]]
    n_unique <- length(unique(vals[!is.na(vals)]))

    # Skip columns with single unique value (no variance)
    if (n_unique <= 1) {
      cat(glue("  Eigencorplot: excluding '{col}' (no variance)\n"))
      next
    }

    # Add to filtered metadata
    if (is.numeric(vals)) {
      meta_filtered[[col]] <- vals
    } else {
      # Convert categorical to numeric factor for correlation
      meta_filtered[[paste0(col, "_num")]] <- as.numeric(as.factor(vals))
    }
  }

  # Update PCA object metadata for eigencorplot
  pca_obj_eigen <- pca_obj
  pca_obj_eigen$metadata <- meta_filtered

  tryCatch({
    if (ncol(meta_filtered) > 0) {
      png(file.path(figures_dir, paste0(prefix, "_PCA_eigencorplot.png")),
          width = 12, height = 10, units = "in", res = 720)
      print(eigencorplot(
        pca_obj_eigen,
        metavars = colnames(meta_filtered),
        cexLabX = 0.7,
        cexLabY = 0.7,
        col = c("blue3", "white", "red3"),
        main = "PC-Metadata Correlations"
      ))
      dev.off()
      cat(glue("  Eigencorplot generated with {ncol(meta_filtered)} variables\n"))
    } else {
      warning("No suitable columns for eigencorplot after filtering")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    warning(glue("Eigencorplot failed: {e$message}"))
  })

  return(pca_obj)
}

# ==============================================================================
# PCA Biplots for Categorical Metadata Columns
# ==============================================================================

plot_pca_biplots <- function(pca_obj, coldata, figures_dir, prefix) {
  # Identify categorical/character columns only (exclude numeric like GSVA scores)
  cat_cols <- sapply(coldata, function(x) is.factor(x) || is.character(x))
  cat_col_names <- names(cat_cols)[cat_cols]

  if (length(cat_col_names) == 0) {
    warning("No categorical columns found in metadata for biplots")
    return(NULL)
  }

  cat(glue("  Generating biplots for {length(cat_col_names)} categorical columns...\n"))

  for (col in cat_col_names) {
    tryCatch({
      p <- biplot(pca_obj,
                  showLoadings = TRUE,
                  lab = NULL,
                  pointSize = 3,
                  sizeLoadingsNames = 3,
                  colby = col,
                  hline = NULL,
                  vline = NULL,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendTitleSize = 12,
                  legendIconSize = 10,
                  title = glue("{prefix} - PCA biplot colored by {col}"))

      output_path <- file.path(figures_dir, paste0(prefix, "_PCA_biplot_", col, ".png"))
      ggsave(output_path, p, width = 12, height = 10, dpi = 720)
      cat(glue("    Saved: {col}\n"))
    }, error = function(e) {
      warning(glue("Biplot for {col} failed: {e$message}"))
    })
  }
}

# ==============================================================================
# PC-Specific GSEA Analysis
# ==============================================================================

run_pc_gsea <- function(pca_obj, species_config, figures_dir, prefix, pcs = c("PC1", "PC2", "PC3")) {
  pca_loadings <- pca_obj$loadings
  pc_gsea_results <- list()

  for (pc in pcs) {
    if (!pc %in% colnames(pca_loadings)) {
      warning(glue("PC {pc} not found in PCA loadings"))
      next
    }

    cat(glue("  Running GSEA for {pc}...\n"))

    # Extract loadings for this PC, sorted high to low
    pc_contribute <- pca_loadings[order(pca_loadings[, pc], decreasing = TRUE), pc, drop = FALSE]

    # Create gene list with loadings as values
    gene_list <- pc_contribute[, 1]
    names(gene_list) <- rownames(pc_contribute)
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Remove duplicates
    gene_list <- gene_list[!duplicated(names(gene_list))]

    if (length(gene_list) > 0) {
      tryCatch({
        gse <- gseGO(
          geneList = gene_list,
          ont = "BP",
          keyType = "SYMBOL",
          pvalueCutoff = 0.05,
          verbose = FALSE,
          OrgDb = species_config$org_db_obj,
          pAdjustMethod = "BH",
          eps = 0
        )

        if (!is.null(gse) && nrow(gse@result) > 0) {
          # Add PC identifier to results
          gse@result$PC <- pc
          pc_gsea_results[[pc]] <- gse@result

          # Generate dotplot
          p <- dotplot(gse, showCategory = 20, split = ".sign",
                       title = glue("{prefix} - {pc} Top Loadings GSEA")) +
            facet_grid(. ~ .sign)

          ggsave(file.path(figures_dir, paste0("GSEA_", prefix, "_", pc, "_loadings.png")),
                 p, width = 14, height = 10, dpi = 720)

          cat(glue("    {pc}: {nrow(gse@result)} enriched pathways\n"))
        } else {
          cat(glue("    {pc}: No significant enrichment\n"))
        }
      }, error = function(e) {
        warning(glue("GSEA for {pc} failed: {e$message}"))
      })
    }
  }

  # Combine all PC GSEA results
  if (length(pc_gsea_results) > 0) {
    combined <- do.call(rbind, pc_gsea_results)
    return(combined)
  }
  return(NULL)
}

plot_de_summary <- function(de_counts, output_path) {
  de_long <- de_counts %>%
    pivot_longer(cols = c(up, dn), names_to = "direction", values_to = "count")

  p <- ggplot(de_long, aes(x = comparison, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("up" = "red", "dn" = "blue"),
                      labels = c("up" = "Upregulated", "dn" = "Downregulated")) +
    labs(
      title = "Differentially Expressed Genes by Comparison",
      x = "Comparison",
      y = "Number of DEGs",
      fill = "Direction"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(output_path, p, width = 14, height = 8, dpi = 720)
}

# ==============================================================================
# GSEA Plots by Source
# ==============================================================================

plot_gsea_by_source <- function(gsea_df, figures_dir, comp_name, n_show = 15) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(NULL)
  }

  sources <- unique(gsea_df$source)

  for (src in sources) {
    src_data <- gsea_df[gsea_df$source == src, ]
    if (nrow(src_data) == 0) next

    # Sort by absolute NES and take top entries
    src_data <- src_data[order(-abs(src_data$NES)), ]
    src_data <- head(src_data, n_show)

    # Create dotplot
    p <- ggplot(src_data, aes(x = NES, y = reorder(Description, NES))) +
      geom_point(aes(size = -log10(pvalue), color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                            name = "NES") +
      scale_size_continuous(name = "-log10(p)") +
      labs(
        title = paste0(comp_name, " - GSEA: ", src),
        x = "Normalized Enrichment Score",
        y = ""
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "right"
      )

    # Sanitize source name for filename
    src_clean <- gsub("[^A-Za-z0-9]", "_", src)
    src_clean <- gsub("_+", "_", src_clean)
    output_path <- file.path(figures_dir, paste0("GSEA_", comp_name, "_", src_clean, ".png"))

    tryCatch({
      ggsave(output_path, p, width = 12, height = 8, dpi = 720)
    }, error = function(e) {
      warning(glue("Failed to create GSEA plot for {src}: {e$message}"))
    })
  }
}

# ==============================================================================
# ORA Plots by Source
# ==============================================================================

plot_ora_by_source <- function(ora_df, figures_dir, comp_name, n_show = 15) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(NULL)
  }

  sources <- unique(ora_df$source)

  for (src in sources) {
    src_data <- ora_df[ora_df$source == src, ]
    if (nrow(src_data) == 0) next

    # Sort by p-value and take top entries
    src_data <- src_data[order(src_data$p_adjust), ]
    src_data <- head(src_data, n_show)

    # Compute gene ratio for visualization
    src_data$gene_ratio <- src_data$gene_count / src_data$term_size

    # Create dotplot with up/down split
    p <- ggplot(src_data, aes(x = gene_ratio, y = reorder(term_name, -p_adjust))) +
      geom_point(aes(size = gene_count, color = -log10(p_adjust))) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "-log10(padj)") +
      scale_size_continuous(name = "Gene Count") +
      facet_wrap(~DE, labeller = labeller(DE = c("up" = "Upregulated", "dn" = "Downregulated"))) +
      labs(
        title = paste0(comp_name, " - ORA: ", src),
        x = "Gene Ratio",
        y = ""
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90")
      )

    # Sanitize source name for filename
    src_clean <- gsub("[^A-Za-z0-9]", "_", src)
    src_clean <- gsub("_+", "_", src_clean)
    output_path <- file.path(figures_dir, paste0("ORA_", comp_name, "_", src_clean, ".png"))

    tryCatch({
      ggsave(output_path, p, width = 14, height = 8, dpi = 720)
    }, error = function(e) {
      warning(glue("Failed to create ORA plot for {src}: {e$message}"))
    })
  }
}

# ==============================================================================
# Save Reproducible Script
# ==============================================================================

save_reproducible_script <- function(params, output_dir) {
  # Build comparisons string
  if (is.list(params$comparisons) || is.matrix(params$comparisons)) {
    comparisons_str <- deparse(params$comparisons, width.cutoff = 500)
    comparisons_str <- paste(comparisons_str, collapse = "")
  } else {
    comparisons_str <- paste0('"', params$comparisons, '"')
  }

  # Handle custom signatures
  custom_sigs_str <- if (is.null(params$custom_sigs) || params$custom_sigs == "" || params$custom_sigs == "none") {
    "NULL"
  } else {
    paste0('"', params$custom_sigs, '"')
  }

  script_content <- glue('
#!/usr/bin/env Rscript
# ==============================================================================
# Reproducible DESeq2 Analysis Script
# Generated by deseq2-rna skill
# Date: {Sys.time()}
# ==============================================================================

# To reproduce this analysis, run:
#   conda run -n {params$conda_env} Rscript {params$prefix}_reproducible.R

# --- Parameters ---
counts_file <- "{params$counts}"
metadata_file <- "{params$metadata}"
design_formula <- "{params$design}"
comparison_col <- "{params$comparison_col}"
reference_level <- "{params$reference}"
comparisons <- {comparisons_str}
msigdb_collections <- "{params$msigdb}"
custom_signatures <- {custom_sigs_str}
species <- "{params$species}"
padj_threshold <- {params$padj}
lfc_threshold <- {params$lfc}
output_directory <- "{params$output_dir}"
prefix <- "{params$prefix}"
conda_env <- "{params$conda_env}"

# --- Run Analysis ---
# Option 1: Run via command line
library(jsonlite)
cmd <- paste(
  paste0("conda run -n ", conda_env, " Rscript"),
  "/home/sagemaker-user/.claude/skills/deseq2-rna/scripts/deseq2_analysis.R",
  "--counts", shQuote(counts_file),
  "--metadata", shQuote(metadata_file),
  "--design", shQuote(design_formula),
  "--comparison_col", comparison_col,
  "--reference", reference_level,
  "--comparisons", shQuote(toJSON(comparisons)),
  "--msigdb", shQuote(msigdb_collections),
  if (!is.null(custom_signatures)) paste("--custom_sigs", shQuote(custom_signatures)) else "",
  "--species", species,
  "--padj", padj_threshold,
  "--lfc", lfc_threshold,
  "--output_dir", shQuote(output_directory),
  "--prefix", prefix,
  "--conda_env", conda_env
)
cat("Run command:\\n", cmd, "\\n")
# system(cmd)  # Uncomment to execute
')

  script_path <- file.path(output_dir, paste0(params$prefix, "_reproducible.R"))
  writeLines(script_content, script_path)
  cat(glue("Reproducible script saved: {script_path}\n"))
}

# ==============================================================================
# Main Analysis
# ==============================================================================

main <- function() {
  cat("\n")
  cat("=", rep("=", 60), "\n", sep = "")
  cat("DESeq2 RNA-seq Analysis\n")
  cat("=", rep("=", 60), "\n\n", sep = "")

  # Parse parameters
  params <- parse_args()

  # Validate required parameters
  required <- c("counts", "metadata", "design", "comparison_col", "reference")
  missing <- required[sapply(required, function(x) is.null(params[[x]]))]
  if (length(missing) > 0) {
    stop(glue("Missing required parameters: {paste(missing, collapse = ', ')}"))
  }

  # Print parameters
  cat("Parameters:\n")
  cat(glue("  Counts: {params$counts}\n"))
  cat(glue("  Metadata: {params$metadata}\n"))
  cat(glue("  Design: {params$design}\n"))
  cat(glue("  Comparison column: {params$comparison_col}\n"))
  cat(glue("  Reference: {params$reference}\n"))
  cat(glue("  Comparisons: {params$comparisons}\n"))
  cat(glue("  MSigDB: {params$msigdb}\n"))
  cat(glue("  Custom signatures: {params$custom_sigs}\n"))
  cat(glue("  Species: {params$species}\n"))
  cat(glue("  padj threshold: {params$padj}\n"))
  cat(glue("  log2FC threshold: {params$lfc}\n"))
  cat(glue("  Output directory: {params$output_dir}\n"))
  cat(glue("  Prefix: {params$prefix}\n"))
  cat("\n")

  # Create output directories
  figures_dir <- file.path(params$output_dir, "figures")
  deg_dir <- file.path(params$output_dir, "deg")
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

  # Get species configuration
  cat("Loading species configuration...\n")
  species_config <- get_species_config(params$species)

  # Load data
  cat("Loading counts data...\n")
  counts_data <- load_counts(params$counts)
  counts <- counts_data$counts
  gene_symbols <- counts_data$gene_symbols
  cat(glue("  Loaded {nrow(counts)} genes x {ncol(counts)} samples\n"))

  cat("Loading metadata...\n")
  coldata <- load_metadata(params$metadata, counts)
  cat(glue("  Loaded {nrow(coldata)} samples\n"))

  # Load gene sets (now returns list with gene_sets and sources)
  cat("Loading MSigDB gene sets...\n")
  msigdb_data <- load_msigdb_collections(params$msigdb, species_config)
  msigdb_sets <- msigdb_data$gene_sets
  msigdb_sources <- msigdb_data$sources

  cat("Loading custom signatures...\n")
  custom_data <- load_custom_signatures(params$custom_sigs)
  custom_sets <- custom_data$gene_sets
  custom_sources <- custom_data$sources

  # Combine gene sets and sources for GSVA and enrichment analyses
  all_gene_sets <- c(msigdb_sets, custom_sets)
  all_gene_set_sources <- c(msigdb_sources, custom_sources)
  cat(glue("  Total gene sets: {length(all_gene_sets)}\n"))

  # Run GSVA on all gene sets
  gsva_scores <- NULL
  custom_gsva_scores <- NULL
  if (length(all_gene_sets) > 0) {
    tryCatch({
      gsva_scores <- run_gsva(counts, all_gene_sets)
      cat(glue("  GSVA computed for {nrow(gsva_scores)} gene sets\n"))

      # Extract custom signature scores for eigencorplot (keeps it manageable)
      if (length(custom_sets) > 0) {
        custom_sig_names <- names(custom_sets)
        custom_gsva_scores <- gsva_scores[rownames(gsva_scores) %in% custom_sig_names, , drop = FALSE]
        cat(glue("  Custom signatures for eigencorplot: {nrow(custom_gsva_scores)}\n"))
      }
    }, error = function(e) {
      warning(glue("GSVA failed: {e$message}"))
    })
  }

  # Run DESeq2
  deseq_results <- run_deseq2(counts, coldata, params$design,
                              params$comparison_col, params$reference)
  dds <- deseq_results$dds
  vsd <- deseq_results$vsd
  cat(glue("  DESeq2 complete: {nrow(dds)} genes\n"))

  # Print available coefficients for list contrast reference
  available_coefs <- print_available_coefficients(dds)

  # Parse simple pairwise comparisons
  comparison_pairs <- parse_comparisons(params$comparisons, params$comparison_col, coldata)
  cat(glue("  Simple comparisons to analyze: {length(comparison_pairs)}\n"))

  # Parse advanced comparisons (list contrasts)
  advanced_comparisons <- parse_advanced_comparisons(params$advanced_comparisons)
  if (length(advanced_comparisons) > 0) {
    cat(glue("  Advanced (list) comparisons to analyze: {length(advanced_comparisons)}\n"))
  }

  # Save checkpoint after DESeq2 + GSVA
  cat("\nSaving checkpoint (DESeq2 + GSVA complete)...\n")
  analysis_data <- list(
    dds = dds,
    vsd = vsd,
    gsva_scores = gsva_scores,
    counts = counts,
    coldata = coldata,
    gene_symbols = gene_symbols,
    params = params,
    stage = "deseq2_gsva_complete"
  )
  saveRDS(analysis_data, file.path(params$output_dir, paste0(params$prefix, "_analysis_data.rds")))

  # PCA plots
  cat("\nGenerating PCA plots...\n")
  plot_pca(vsd, coldata, params$comparison_col,
           file.path(figures_dir, paste0(params$prefix, "_PCA_basic_", params$comparison_col, ".png")))

  pca_obj <- plot_pca_advanced(vsd, coldata, custom_gsva_scores, figures_dir, params$prefix)

  # PCA biplots for each categorical metadata column
  cat("\nGenerating PCA biplots...\n")
  plot_pca_biplots(pca_obj, coldata, figures_dir, params$prefix)

  # PC-specific GSEA analysis
  cat("\nRunning PC-specific GSEA...\n")
  pc_gsea_results <- run_pc_gsea(pca_obj, species_config, figures_dir, params$prefix)

  # Save PC GSEA results
  if (!is.null(pc_gsea_results)) {
    write.table(pc_gsea_results,
                file.path(deg_dir, paste0(params$prefix, "_PC_gsea.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Process each comparison
  cat("\nProcessing comparisons...\n")
  all_results <- list()
  all_gsea <- list()
  all_gprofiler <- list()
  de_counts <- data.frame()

  for (pair in comparison_pairs) {
    numerator <- pair[1]
    denominator <- pair[2]
    cat(glue("  {numerator} vs {denominator}...\n"))

    # Extract results
    res_data <- extract_results(dds, params$comparison_col, numerator, denominator, gene_symbols)
    res <- res_data$results
    comp_name <- res_data$name

    all_results[[comp_name]] <- res

    # Count DEGs
    sig_up <- sum(res$padj < params$padj & res$log2FoldChange > params$lfc, na.rm = TRUE)
    sig_dn <- sum(res$padj < params$padj & res$log2FoldChange < -params$lfc, na.rm = TRUE)
    de_counts <- rbind(de_counts, data.frame(
      comparison = comp_name,
      up = sig_up,
      dn = sig_dn
    ))

    # Save DEG results
    write.table(res, file.path(deg_dir, paste0(params$prefix, "_", comp_name, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Volcano plot
    plot_volcano(res, params$padj, params$lfc,
                 file.path(figures_dir, paste0("volcano_", comp_name, ".png")))

    # ============= Combined GSEA Analysis =============
    gsea_combined <- NULL

    # GSEA GO:BP
    gsea_gobp <- run_gsea_gobp(res, species_config)
    if (!is.null(gsea_gobp) && nrow(gsea_gobp@result) > 0) {
      gsea_gobp@result$source <- "GO:BP"
      gsea_combined <- gsea_gobp@result
    }

    # GSEA MSigDB + Custom gene sets
    if (length(all_gene_sets) > 0) {
      gsea_custom <- run_gsea_custom(res, all_gene_sets, all_gene_set_sources)
      if (!is.null(gsea_custom) && nrow(gsea_custom@result) > 0) {
        if (is.null(gsea_combined)) {
          gsea_combined <- gsea_custom@result
        } else {
          gsea_combined <- rbind(gsea_combined, gsea_custom@result)
        }
      }
    }

    # Store and save combined GSEA results
    if (!is.null(gsea_combined) && nrow(gsea_combined) > 0) {
      gsea_combined$comparison <- comp_name
      all_gsea[[comp_name]] <- gsea_combined

      # Save combined GSEA results
      write.table(gsea_combined,
                  file.path(deg_dir, paste0(params$prefix, "_", comp_name, "_gsea.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # Generate plots by source
      plot_gsea_by_source(gsea_combined, figures_dir, comp_name)
    }

    # ============= Combined ORA Analysis =============
    sig_genes <- res[res$padj < params$padj & !is.na(res$padj), ]
    genes_up <- sig_genes$symbol[sig_genes$log2FoldChange > params$lfc]
    genes_dn <- sig_genes$symbol[sig_genes$log2FoldChange < -params$lfc]

    # g:Profiler (GO:BP, KEGG, REAC)
    gp_res <- run_gprofiler(genes_up, genes_dn, species_config)

    # ORA MSigDB + Custom gene sets
    ora_custom <- list()
    if (length(all_gene_sets) > 0) {
      ora_custom <- run_ora_custom(genes_up, genes_dn, all_gene_sets, all_gene_set_sources)
    }

    # Harmonize and combine ORA results
    ora_combined <- harmonize_ora_results(gp_res, ora_custom)

    if (!is.null(ora_combined) && nrow(ora_combined) > 0) {
      ora_combined$comparison <- comp_name
      all_gprofiler[[comp_name]] <- ora_combined

      # Save combined ORA results
      write.table(ora_combined,
                  file.path(deg_dir, paste0(params$prefix, "_", comp_name, "_ora.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # Generate plots by source
      plot_ora_by_source(ora_combined, figures_dir, comp_name)
    }
  }

  # ===========================================================================
  # Process Advanced Comparisons (List Contrasts)
  # ===========================================================================
  if (length(advanced_comparisons) > 0) {
    cat("\nProcessing advanced (list) comparisons...\n")

    for (adv_comp in advanced_comparisons) {
      comp_name <- adv_comp$name
      cat(glue("  {comp_name} (list contrast)...\n"))

      # Extract results using list contrast
      tryCatch({
        res_data <- extract_results_list(
          dds,
          adv_comp$contrast,
          adv_comp$listValues,
          comp_name,
          gene_symbols
        )
        res <- res_data$results

        all_results[[comp_name]] <- res

        # Count DEGs
        sig_up <- sum(res$padj < params$padj & res$log2FoldChange > params$lfc, na.rm = TRUE)
        sig_dn <- sum(res$padj < params$padj & res$log2FoldChange < -params$lfc, na.rm = TRUE)
        de_counts <- rbind(de_counts, data.frame(
          comparison = comp_name,
          up = sig_up,
          dn = sig_dn
        ))
        cat(glue("    DEGs: {sig_up} up, {sig_dn} down\n"))

        # Save DEG results
        write.table(res, file.path(deg_dir, paste0(params$prefix, "_", comp_name, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Volcano plot
        plot_volcano(res, params$padj, params$lfc,
                     file.path(figures_dir, paste0("volcano_", comp_name, ".png")))

        # ============= Combined GSEA Analysis =============
        gsea_combined <- NULL

        # GSEA GO:BP
        gsea_gobp <- run_gsea_gobp(res, species_config)
        if (!is.null(gsea_gobp) && nrow(gsea_gobp@result) > 0) {
          gsea_gobp@result$source <- "GO:BP"
          gsea_combined <- gsea_gobp@result
        }

        # GSEA MSigDB + Custom gene sets
        if (length(all_gene_sets) > 0) {
          gsea_custom <- run_gsea_custom(res, all_gene_sets, all_gene_set_sources)
          if (!is.null(gsea_custom) && nrow(gsea_custom@result) > 0) {
            if (is.null(gsea_combined)) {
              gsea_combined <- gsea_custom@result
            } else {
              gsea_combined <- rbind(gsea_combined, gsea_custom@result)
            }
          }
        }

        if (!is.null(gsea_combined)) {
          gsea_combined$comparison <- comp_name
          all_gsea[[comp_name]] <- gsea_combined

          # Save combined GSEA results
          write.table(gsea_combined,
                      file.path(deg_dir, paste0(params$prefix, "_", comp_name, "_gsea.txt")),
                      sep = "\t", quote = FALSE, row.names = FALSE)

          # Generate plots by source
          plot_gsea_by_source(gsea_combined, figures_dir, comp_name)
        }

        # ============= Combined ORA Analysis =============
        sig_genes <- res[res$padj < params$padj & !is.na(res$padj), ]
        genes_up <- sig_genes$symbol[sig_genes$log2FoldChange > params$lfc]
        genes_dn <- sig_genes$symbol[sig_genes$log2FoldChange < -params$lfc]

        # g:Profiler (GO:BP, KEGG, REAC)
        gp_res <- run_gprofiler(genes_up, genes_dn, species_config)

        # ORA MSigDB + Custom gene sets
        ora_custom <- list()
        if (length(all_gene_sets) > 0) {
          ora_custom <- run_ora_custom(genes_up, genes_dn, all_gene_sets, all_gene_set_sources)
        }

        # Harmonize and combine ORA results
        ora_combined <- harmonize_ora_results(gp_res, ora_custom)

        if (!is.null(ora_combined) && nrow(ora_combined) > 0) {
          ora_combined$comparison <- comp_name
          all_gprofiler[[comp_name]] <- ora_combined

          # Save combined ORA results
          write.table(ora_combined,
                      file.path(deg_dir, paste0(params$prefix, "_", comp_name, "_ora.txt")),
                      sep = "\t", quote = FALSE, row.names = FALSE)

          # Generate plots by source
          plot_ora_by_source(ora_combined, figures_dir, comp_name)
        }

      }, error = function(e) {
        warning(glue("Advanced comparison '{comp_name}' failed: {e$message}"))
      })
    }
  }

  # Save checkpoint after comparisons
  cat("\nSaving checkpoint (comparisons complete)...\n")
  analysis_data <- list(
    dds = dds,
    vsd = vsd,
    gsva_scores = gsva_scores,
    pca = pca_obj,
    results = all_results,
    gsea = all_gsea,
    gprofiler = all_gprofiler,
    de_counts = de_counts,
    params = params,
    stage = "comparisons_complete"
  )
  saveRDS(analysis_data, file.path(params$output_dir, paste0(params$prefix, "_analysis_data.rds")))

  # Save combined results
  cat("\nSaving combined results...\n")

  # Combined DEG statistics
  all_results_df <- do.call(rbind, all_results)
  write.table(all_results_df,
              file.path(deg_dir, paste0(params$prefix, "_summstats_all.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # DE counts summary
  write.table(de_counts,
              file.path(deg_dir, paste0(params$prefix, "_DE_lfc", params$lfc, "padj", params$padj, "_count.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Combined GSEA (GO:BP + MSigDB + Custom)
  if (length(all_gsea) > 0) {
    all_gsea_df <- do.call(rbind, all_gsea)
    write.table(all_gsea_df,
                file.path(deg_dir, paste0(params$prefix, "_gsea_all.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Combined ORA (g:Profiler + MSigDB + Custom)
  if (length(all_gprofiler) > 0) {
    all_ora_df <- do.call(rbind, all_gprofiler)
    write.table(all_ora_df,
                file.path(deg_dir, paste0(params$prefix, "_ora_all.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # DE summary barplot
  plot_de_summary(de_counts, file.path(figures_dir, "DE_genes_summary_barplot.png"))

  # Save final R objects
  cat("Saving final R objects...\n")
  analysis_data <- list(
    dds = dds,
    vsd = vsd,
    results = all_results,
    gsea = all_gsea,
    gprofiler = all_gprofiler,
    de_counts = de_counts,
    gsva_scores = gsva_scores,
    pca = pca_obj,
    params = params,
    stage = "complete"
  )

  saveRDS(analysis_data, file.path(params$output_dir, paste0(params$prefix, "_analysis_data.rds")))
  save.image(file.path(params$output_dir, paste0(params$prefix, "_analysis.RData")))

  cat("\n")
  cat("=", rep("=", 60), "\n", sep = "")
  cat("Analysis Complete\n")
  cat("=", rep("=", 60), "\n\n", sep = "")
  cat(glue("Output directory: {params$output_dir}\n"))
  cat(glue("Figures: {figures_dir}\n"))
  cat(glue("DEG results: {deg_dir}\n"))
  cat(glue("R objects: {params$prefix}_analysis_data.rds\n"))
  cat("\n")

  # Print DE summary
  cat("DE Gene Summary:\n")
  print(de_counts)

  # Generate draft report
  cat("\nGenerating draft report...\n")
  script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
  if (length(script_dir) == 0) script_dir <- "."
  report_script <- file.path(script_dir, "generate_report.py")

  if (file.exists(report_script)) {
    report_output <- file.path(params$output_dir, paste0(params$prefix, "_Draft_Report.md"))

    # Save config for report script
    config_for_report <- list(
      input_dir = params$output_dir,
      output = report_output,
      prefix = params$prefix,
      padj = params$padj,
      lfc = params$lfc,
      title = paste0(params$prefix, " RNA-seq Analysis Report"),
      description = if (!is.null(params$description)) params$description else ""
    )
    config_path <- file.path(params$output_dir, "report_config.json")
    write(toJSON(config_for_report, auto_unbox = TRUE), config_path)

    report_cmd <- sprintf('python3 "%s" --config "%s"', report_script, config_path)

    tryCatch({
      system(report_cmd)
      cat(glue("Draft report generated: {report_output}\n"))
    }, error = function(e) {
      warning(glue("Report generation failed: {e$message}"))
    })
  } else {
    cat(glue("Report script not found at {report_script}, skipping report generation.\n"))
  }

  # Generate interactive plots (always generated alongside static plots)
  cat("\nGenerating interactive HTML plots...\n")
  interactive_script <- file.path(script_dir, "interactive_plots.R")

  if (file.exists(interactive_script)) {
    tryCatch({
      source(interactive_script)
      interactive_dir <- file.path(figures_dir, "interactive")
      dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

      rds_path <- file.path(params$output_dir, paste0(params$prefix, "_analysis_data.rds"))
      generate_interactive_plots(
        analysis_rds = rds_path,
        output_dir = interactive_dir,
        deg_dir = deg_dir
      )
      cat(glue("Interactive plots saved to: {interactive_dir}\n"))
    }, error = function(e) {
      warning(glue("Interactive plot generation failed: {e$message}"))
    })
  } else {
    cat(glue("Interactive plots script not found at {interactive_script}, skipping.\n"))
  }

  # Save reproducible script with all parameters
  save_reproducible_script(params, params$output_dir)
}

# Run main function
main()
