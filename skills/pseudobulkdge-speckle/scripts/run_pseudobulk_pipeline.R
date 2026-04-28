#!/usr/bin/env Rscript

# ==============================================================================
# Pseudobulk DEG Analysis Pipeline
#
# This script performs comprehensive single-cell RNA-seq pseudobulk differential
# expression analysis using R/Bioconductor packages (scran, edgeR, fgsea, speckle).
#
# Usage:
#   Rscript run_pseudobulk_pipeline.R \
#     --input <input_file.rds> \
#     --output <output_dir> \
#     --condition <condition_column> \
#     --sample <sample_column> \
#     --celltype <celltype_column> \
#     --comparisons "group1,group2;group1,group3"
#
# Optional arguments:
#   --batch <batch_column>
#   --covariate <covariate_column>
#   --custom_signatures <path_to_gmt_or_genelist>
#   --min_count <int> (default: 3)
#   --min_total_count <int> (default: 5)
#   --min_prop <float> (default: 0.2)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scran)
  library(scuttle)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(reactome.db)
  library(org.Hs.eg.db)
  library(speckle)
  library(limma)
  library(dplyr)
  library(tibble)
  library(glue)
  library(presto)
  library(optparse)
  library(parallel)
  library(EnhancedVolcano)
  library(gridExtra)
  library(plotly)
  library(htmlwidgets)
})

# ==============================================================================
# Command Line Argument Parsing
# ==============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input Seurat RDS file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "./pseudobulk_results",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-c", "--condition"), type = "character", default = NULL,
              help = "Condition/disease column name", metavar = "COLUMN"),
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "Sample identifier column name", metavar = "COLUMN"),
  make_option(c("-t", "--celltype"), type = "character", default = NULL,
              help = "Cell type annotation column name", metavar = "COLUMN"),
  make_option(c("-b", "--batch"), type = "character", default = NULL,
              help = "Batch column name (optional)", metavar = "COLUMN"),
  make_option(c("-v", "--covariate"), type = "character", default = NULL,
              help = "Covariate column name (optional)", metavar = "COLUMN"),
  make_option(c("-p", "--comparisons"), type = "character", default = NULL,
              help = "Comparison pairs: 'g1,g2;g1,g3'", metavar = "STRING"),
  make_option(c("--custom_signatures"), type = "character", default = NULL,
              help = "Path to custom signatures (GMT or text file) - uses GENE SYMBOLS", metavar = "FILE"),
  make_option(c("--min_count"), type = "integer", default = 3,
              help = "filterByExpr min.count [default: %default]"),
  make_option(c("--min_total_count"), type = "integer", default = 5,
              help = "filterByExpr min.total.count [default: %default]"),
  make_option(c("--min_prop"), type = "double", default = 0.2,
              help = "filterByExpr min.prop [default: %default]"),
  make_option(c("--volcano_fdr"), type = "double", default = 0.05,
              help = "FDR cutoff for volcano plot significance [default: %default]"),
  make_option(c("--volcano_log2fc"), type = "double", default = 1,
              help = "|log2FC| cutoff for volcano plot significance [default: %default]"),

  # Module selection flags
  make_option(c("--modules"), type = "character", default = "all",
              help = "Comma-separated modules to run: deg,fgsea_groups,fgsea_clusters,speckle,markers OR 'all' [default: %default]"),
  make_option(c("--skip_fgsea_clusters"), action = "store_true", default = FALSE,
              help = "Skip cluster-level fgsea analysis (faster)"),
  make_option(c("--deg_only"), action = "store_true", default = FALSE,
              help = "Run only DEG analysis"),
  make_option(c("--fgsea_only"), action = "store_true", default = FALSE,
              help = "Run only fgsea on existing DEG results (requires --deg_input)"),
  make_option(c("--speckle_only"), action = "store_true", default = FALSE,
              help = "Run only speckle cell composition analysis"),
  make_option(c("--deg_input"), type = "character", default = NULL,
              help = "Path to directory with existing DEG results (for --fgsea_only mode)", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ==============================================================================
# Module Selection Logic
# ==============================================================================

#' Parse module selection from CLI arguments
#' @param opt Parsed options from optparse
#' @return Character vector of modules to run
parse_modules <- function(opt) {
  # Handle convenience flags first (they override --modules)
  if (opt$deg_only) {
    return(c("deg"))
  }
  if (opt$fgsea_only) {
    return(c("fgsea_groups"))
  }
  if (opt$speckle_only) {
    return(c("speckle"))
  }

  # Parse --modules argument
  if (opt$modules == "all") {
    modules <- c("deg", "fgsea_groups", "fgsea_clusters", "speckle", "markers")
  } else {
    modules <- trimws(strsplit(opt$modules, ",")[[1]])
  }

  # Apply --skip_fgsea_clusters flag

  if (opt$skip_fgsea_clusters) {
    modules <- modules[modules != "fgsea_clusters"]
  }

  return(modules)
}

# Parse which modules to run
selected_modules <- parse_modules(opt)
message("Selected modules: ", paste(selected_modules, collapse = ", "))

# ==============================================================================
# Argument Validation (mode-dependent)
# ==============================================================================

# For fgsea_only mode, we need --deg_input instead of full input requirements
if (opt$fgsea_only) {
  if (is.null(opt$deg_input)) {
    stop("--deg_input is required when using --fgsea_only mode")
  }
  if (!dir.exists(opt$deg_input)) {
    stop("DEG input directory does not exist: ", opt$deg_input)
  }
  # Check for DEG files in the input directory
  deg_files <- list.files(opt$deg_input, pattern = "_DE\\.tsv$", full.names = TRUE)
  if (length(deg_files) == 0) {
    stop("No DEG result files (*_DE.tsv) found in: ", opt$deg_input)
  }
  message("Found ", length(deg_files), " DEG result file(s) for fgsea rerun")
} else {
  # Standard validation for other modes
  if (is.null(opt$input)) stop("--input is required")
  if (is.null(opt$condition)) stop("--condition is required")
  if (is.null(opt$sample)) stop("--sample is required")
  if (is.null(opt$celltype)) stop("--celltype is required")
  if (is.null(opt$comparisons)) stop("--comparisons is required")
}

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Load custom gene signatures from various formats
#'
#' IMPORTANT: Custom signatures use GENE SYMBOLS (same as MSigDB)
#' They are run with the same ranked gene list as MSigDB analysis
#'
#' @param input Path to GMT file, text file (one gene per row), or comma-separated gene list
#' @return Named list of gene sets (each element is a character vector of gene symbols)
load_custom_signatures <- function(input) {
  if (is.null(input) || input == "" || input == "NULL") {
    return(NULL)
  }

  if (file.exists(input) && grepl("\\.gmt$", input, ignore.case = TRUE)) {
    # GMT file format - standard format with gene symbols
    message("Loading custom signatures from GMT file: ", input)
    signatures <- fgsea::gmtPathways(input)
    message("  Loaded ", length(signatures), " gene sets")
    for (name in names(signatures)[1:min(3, length(signatures))]) {
      message("    - ", name, ": ", length(signatures[[name]]), " genes")
    }
    return(signatures)
  } else if (file.exists(input)) {
    # Text file - one gene symbol per row OR tab-separated: signature_name\tgene1,gene2,...
    message("Loading custom signatures from text file: ", input)
    lines <- readLines(input)
    lines <- lines[lines != ""]  # Remove empty lines

    # Check if tab-separated format (name\tgenes)
    if (any(grepl("\t", lines))) {
      signatures <- list()
      for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) {
          name <- parts[1]
          genes <- unlist(strsplit(parts[2], "[,;]"))
          signatures[[name]] <- trimws(genes)
        }
      }
      message("  Loaded ", length(signatures), " named gene sets")
    } else {
      # Simple list - one gene per row
      signatures <- list("custom_signature" = trimws(lines))
      message("  Loaded 1 gene set with ", length(lines), " genes")
    }
    return(signatures)
  } else {
    # Direct gene list (comma or newline separated)
    message("Parsing custom signatures from gene list string")
    genes <- unlist(strsplit(input, "[,\n;]+"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    message("  Parsed ", length(genes), " genes")
    return(list("custom_signature" = genes))
  }
}

#' Modify DE results dataframe
modify_DE_results <- function(current_de_results, group1, group2, cluster = NULL, clust_res = NULL) {
  current_de_results$group1 <- group1
  current_de_results$group2 <- group2
  if (!is.null(cluster) && !is.null(clust_res)) {
    current_de_results$cluster <- cluster
    current_de_results$clust_res <- clust_res
  }

  # Replace NA values
  if ("logFC" %in% colnames(current_de_results)) {
    current_de_results$logFC[is.na(current_de_results$logFC)] <- 0
  }
  if ("logCPM" %in% colnames(current_de_results)) {
    current_de_results$logCPM[is.na(current_de_results$logCPM)] <- 0
  }
  if ("F" %in% colnames(current_de_results)) {
    current_de_results$F[is.na(current_de_results$F)] <- 0
  }
  if ("PValue" %in% colnames(current_de_results)) {
    current_de_results$PValue[is.na(current_de_results$PValue)] <- 1
  }
  if ("FDR" %in% colnames(current_de_results)) {
    current_de_results$FDR[is.na(current_de_results$FDR)] <- 1
  }

  current_de_results$gene_name <- rownames(current_de_results)
  rownames(current_de_results) <- NULL

  return(current_de_results)
}

#' Update master dataframe with current results
update_DE_results_df <- function(de_results_dataframe, current_de_results) {
  if (dim(de_results_dataframe)[1] == 0) {
    de_results_dataframe <- current_de_results
  } else {
    de_results_dataframe <- rbind(de_results_dataframe, current_de_results)
  }
  return(de_results_dataframe)
}

#' Store edgeR metrics in nested list
writeToNestedList <- function(in_list, cluster_name, group1_name, group2_name, value) {
  group1_name <- toString(group1_name)
  group2_name <- toString(group2_name)

  if (!cluster_name %in% names(in_list)) {
    in_list[[cluster_name]] <- list()
  }
  if (!group1_name %in% names(in_list[[cluster_name]])) {
    in_list[[cluster_name]][[group1_name]] <- list()
  }

  if (!group2_name %in% names(in_list[[cluster_name]][[group1_name]])) {
    in_list[[cluster_name]][[group1_name]][[group2_name]] <- value
  }

  return(in_list)
}

fill_edge_metrics_by_groups <- function(de_metrics_list, cluster_name, group1_name, group2_name, value) {
  de_metrics_list <- writeToNestedList(
    in_list = de_metrics_list,
    cluster_name = cluster_name,
    group1_name = group1_name,
    group2_name = group2_name,
    value = value
  )
  de_metrics_list <- writeToNestedList(
    in_list = de_metrics_list,
    cluster_name = cluster_name,
    group1_name = group2_name,
    group2_name = group1_name,
    value = value
  )
  return(de_metrics_list)
}

#' Calculate adjusted p-value
calc_adjusted_pval <- function(dataset, method = "bonferroni") {
  if ("PValue" %in% names(dataset)) {
    p.val <- as.numeric(dataset$PValue)
    p.adj <- p.adjust(p.val, method = method)
    names(dataset)[names(dataset) == "PValue"] <- "pval"
    dataset$padj <- p.adj
    return(dataset)
  } else {
    stop("Column PValue not found in the dataset.")
  }
}

#' Get MSigDB gene sets (uses GENE SYMBOLS)
get_fgsea_sets_msigdbr <- function() {
  # Note: msigdbr 10.0.0+ uses 'collection' instead of deprecated 'category'
  msig.df <- msigdbr(species = "Homo sapiens", collection = "C8")
  msig.df1 <- msigdbr(species = "Homo sapiens", collection = "C2")
  msig.df2 <- msigdbr(species = "Homo sapiens", collection = "H")

  fgsea.set <- msig.df %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea.set1 <- msig.df1 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea.set2 <- msig.df2 %>% split(x = .$gene_symbol, f = .$gs_name)
  return(list("C8" = fgsea.set, "C2" = fgsea.set1, "H" = fgsea.set2))
}

#' Convert gene symbols to Entrez IDs (for Reactome only)
convert_to_entrez_ids_df <- function(dataset, entrez.db) {
  g.entrezid <- mapIds(entrez.db,
                       keys = dataset$gene_name,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

  dataset$gene_name <- g.entrezid
  dataset <- dataset[complete.cases(dataset), ]
  names(dataset)[names(dataset) == "gene_name"] <- "feature"

  return(dataset)
}

#' Perform fGSEA multilevel analysis
performGSEAmultilevel <- function(markers, fgsea.set, metric, cluster, n.lines, minS, maxS, sSize, gseaParam, scoreType, eps, coll.path = FALSE, sample_name) {
  message("Running fgseaMultilevel for cluster: ", cluster)

  if (metric == "auc") {
    clstr.genes <- markers %>%
      dplyr::filter(group == cluster) %>%
      arrange(desc(auc)) %>%
      dplyr::select(feature, auc)
  } else {
    markers$padj[markers$padj == 0] <- 2e-230
    markers$pval[markers$pval == 0] <- 2e-230
    fcSign <- markers$logFC * (-log10(markers$pval))
    markers$metric2 <- fcSign

    clstr.genes <- markers %>%
      dplyr::filter(cluster == .env$cluster) %>%
      arrange(desc(metric2)) %>%
      dplyr::select(feature, metric2)
  }
  rnks <- deframe(clstr.genes)

  fgsea.res <- fgsea(fgsea.set, rnks, minSize = minS, maxSize = maxS,
                     sampleSize = sSize, gseaParam = gseaParam,
                     scoreType = scoreType, eps = eps, nproc = 8)

  fgsea.res.tidy <- fgsea.res %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::select(-ES, -log2err) %>%
    arrange(padj) %>%
    head(n = n.lines)

  if (coll.path) {
    collapsed.pathways <- collapsePathways(
      fgsea.res[order(pval)][padj < 0.01],
      fgsea.set, rnks
    )
    main.pathways <- fgsea.res[pathway %in% collapsed.pathways$mainPathways] %>%
      arrange(desc(NES))
  } else {
    main.pathways <- fgsea.res %>%
      arrange(desc(NES))
  }
  main.pathways$leadingEdge <- vapply(main.pathways$leadingEdge, function(x) paste(unlist(x), collapse = ","), character(1), USE.NAMES = FALSE)

  top.pathways <- fgsea.res[NES > 0][head(order(-NES), n = 30), pathway]

  tryCatch({
    fig <- plotGseaTable(fgsea.set[top.pathways], rnks, fgsea.res, gseaParam = 1, colwidths = c(10, 0, 1, 1, 1), render = FALSE)
    grid.arrange(top = paste0("Top pathways in cluster ", cluster), fig)
  }, error = function(e) {
    message("Warning: Could not generate GSEA plot for cluster ", cluster, ": ", e$message)
  })

  message("Finished processing cluster ", cluster)

  main.pathways$cluster <- cluster
  tem.str <- strsplit(sample_name, "_")[[1]]
  main.pathways$database <- tem.str[length(tem.str)]
  if ('comparison' %in% colnames(markers)) {
    main.pathways$comparison <- unique(markers$comparison)
  } else {
    main.pathways$comparison <- cluster
  }
  return(main.pathways)
}

#' Perform fGSEA by groups
performGSEAmultilevel_by_groups <- function(de_df_list, fgsea.sets, file_namings_info, analysis_type) {
  if (!is.list(de_df_list)) stop("Dataset must be passed as a list.")
  if (!is.list(fgsea.sets)) stop("fgsea.sets must be passed as a list of sets.")
  if (!is.list(file_namings_info)) stop("Please provide file naming info as a list.")

  df_list <- list()
  for (fgsea.set_i in 1:length(fgsea.sets)) {
    fgsea.set_name <- names(fgsea.sets[fgsea.set_i])
    if (is.null(fgsea.set_name)) {
      set_name <- ""
    } else {
      set_name <- fgsea.set_name
    }

    df_siglist <- list()
    for (i in 1:length(de_df_list)) {
      if (analysis_type == "all") {
        de_df_groups <- unique(de_df_list[[i]][["group"]])
        comparison <- ""
      } else {
        de_df_groups <- unique(de_df_list[[i]][["cluster"]])
        if (!is.null(file_namings_info$comparison)) {
          comparison <- paste0(file_namings_info$comparison, "_")
        } else {
          comparison <- paste0(names(de_df_list)[i], "_")
        }
      }
      pdf_name <- file.path(file_namings_info$folder_name, paste0(comparison, file_namings_info$level, "_", file_namings_info$db_name, set_name, "_pathways.pdf"))
      sample_name <- file.path(file_namings_info$folder_name, paste0(file_namings_info$level, "_", comparison, file_namings_info$db_name, set_name))

      pdf(pdf_name, width = 11, height = 5)
      df.vec <- do.call(rbind,
                        lapply(seq_along(de_df_groups),
                               function(n) {
                                 performGSEAmultilevel(markers = de_df_list[[i]], fgsea.set = fgsea.sets[[fgsea.set_i]], metric = "pval",
                                                       cluster = de_df_groups[n], n.lines = 10, minS = 1, maxS = Inf, sSize = 101, gseaP = 1,
                                                       scoreT = "std", eps = 1 * 10^-10, coll.path = FALSE, sample_name = sample_name)
                               }))
      dev.off()
      df_siglist[[i]] <- df.vec
    }
    df_list[[fgsea.set_i]] <- do.call(rbind, df_siglist)
  }
  df.concat <- do.call(rbind, df_list)
  return(df.concat)
}

#' Render volcano plot
render_volcano_plot <- function(current_de_results, group1, group2, plots_subfolder, comparison_name, cluster = NULL) {
  filter_pvalues <- current_de_results[!is.na(current_de_results$PValue) & current_de_results$PValue <= 0.05, ]
  filter_logFC <- filter_pvalues[!is.na(filter_pvalues$logFC) & abs(filter_pvalues$logFC) > 2, c("gene_name", "logFC", "PValue")]
  as.data.frame(filter_logFC[order(abs(filter_logFC$logFC), decreasing = TRUE), ])

  if (!is.null(cluster)) {
    title_suffix <- paste0(comparison_name, "_", group1, "_vs_", group2, "_", cluster)
  } else {
    title_suffix <- paste0(comparison_name, "_", group1, "_vs_", group2)
  }

  reso <- 300
  length <- 2.25 * reso / 72
  png(file.path(plots_subfolder, paste0("Volcano_plot_for_", title_suffix, ".png")), width = length, height = length, units = "in", res = reso)
  tryCatch({
    plot(EnhancedVolcano(current_de_results,
                         lab = current_de_results[, "gene_name"],
                         x = "logFC",
                         y = "PValue",
                         selectLab = filter_pvalues$gene_name[head(order(abs(filter_pvalues$logFC), decreasing = TRUE), n = 20)],
                         pointSize = 2,
                         labSize = 3,
                         title = paste("Volcano plot for", title_suffix),
                         titleLabSize = 10,
                         subtitle = NULL,
                         pCutoff = 0.05,
                         FCcutoff = 2
    ))
  }, error = function(e) {
    message("Warning: Could not generate volcano plot: ", e$message)
  })
  dev.off()
}

#' Render interactive volcano plot as HTML with hover information
#' @param current_de_results Data frame with DEG results (gene_name, logFC, PValue, FDR)
#' @param group1 First group name (reference)
#' @param group2 Second group name (comparison)
#' @param plots_subfolder Output folder for plots
#' @param comparison_name Name prefix for the comparison
#' @param cluster Optional cluster name for cluster-specific plots
#' @param fdr_cutoff FDR/adjusted p-value cutoff for significance (default: 0.05)
#' @param log2fc_cutoff Absolute log2 fold change cutoff for significance (default: 1)
render_interactive_volcano_plot <- function(current_de_results, group1, group2, plots_subfolder,
                                            comparison_name, cluster = NULL,
                                            fdr_cutoff = 0.05, log2fc_cutoff = 1) {

  if (!is.null(cluster)) {
    title_suffix <- paste0(comparison_name, "_", group1, "_vs_", group2, "_", cluster)
    plot_title <- paste0("Volcano Plot: ", group2, " vs ", group1, "\nCluster: ", cluster,
                        "\n(FDR < ", fdr_cutoff, ", |log2FC| > ", log2fc_cutoff, ")")
  } else {
    title_suffix <- paste0(comparison_name, "_", group1, "_vs_", group2)
    plot_title <- paste0("Volcano Plot: ", group2, " vs ", group1,
                        "\n(FDR < ", fdr_cutoff, ", |log2FC| > ", log2fc_cutoff, ")")
  }

  tryCatch({
    # Prepare data for plotting
    df <- current_de_results
    df$neg_log10_pval <- -log10(df$PValue)

    # Determine significance category for coloring based on FDR and log2FC cutoffs
    df$significance <- "Not Significant"
    df$significance[df$FDR < fdr_cutoff & df$logFC > log2fc_cutoff] <- "Up-regulated"
    df$significance[df$FDR < fdr_cutoff & df$logFC < -log2fc_cutoff] <- "Down-regulated"
    df$significance[df$FDR < fdr_cutoff & abs(df$logFC) <= log2fc_cutoff] <- "Significant (small FC)"

    # Define colors
    color_map <- c(
      "Up-regulated" = "#E41A1C",
      "Down-regulated" = "#377EB8",
      "Significant (small FC)" = "#FF7F00",
      "Not Significant" = "#999999"
    )

    # Create hover text with gene name as title
    df$hover_text <- paste0(
      "<b>", df$gene_name, "</b><br>",
      "Log2 Fold Change: ", round(df$logFC, 3), "<br>",
      "P-value: ", signif(df$PValue, 3), "<br>",
      "FDR: ", signif(df$FDR, 3)
    )

    # Create plotly scatter plot
    p <- plot_ly(
      data = df,
      x = ~logFC,
      y = ~neg_log10_pval,
      type = "scatter",
      mode = "markers",
      color = ~significance,
      colors = color_map,
      text = ~hover_text,
      hoverinfo = "text",
      marker = list(size = 5, opacity = 0.7)
    ) %>%
      layout(
        title = list(text = plot_title, font = list(size = 14)),
        xaxis = list(
          title = "Log2 Fold Change",
          zeroline = TRUE,
          zerolinecolor = "#969696",
          zerolinewidth = 1
        ),
        yaxis = list(
          title = "-Log10(P-value)",
          zeroline = FALSE
        ),
        # Add significance threshold lines
        shapes = list(
          # Vertical line for negative log2FC cutoff
          list(
            type = "line",
            x0 = -log2fc_cutoff,
            x1 = -log2fc_cutoff,
            y0 = 0,
            y1 = max(df$neg_log10_pval, na.rm = TRUE) * 1.1,
            line = list(color = "#969696", dash = "dash", width = 1)
          ),
          # Vertical line for positive log2FC cutoff
          list(
            type = "line",
            x0 = log2fc_cutoff,
            x1 = log2fc_cutoff,
            y0 = 0,
            y1 = max(df$neg_log10_pval, na.rm = TRUE) * 1.1,
            line = list(color = "#969696", dash = "dash", width = 1)
          )
        ),
        legend = list(
          orientation = "h",
          x = 0.5,
          xanchor = "center",
          y = -0.15
        ),
        hoverlabel = list(
          bgcolor = "white",
          font = list(size = 12)
        )
      )

    # Save as HTML
    html_file <- file.path(plots_subfolder, paste0("Volcano_plot_", title_suffix, ".html"))
    saveWidget(p, html_file, selfcontained = TRUE)
    message("Saved interactive volcano plot: ", html_file)

  }, error = function(e) {
    message("Warning: Could not generate interactive volcano plot: ", e$message)
  })
}

#' Generate DEG completion report after pseudobulk analysis
#'
#' Creates three output files:
#' 1. deg_completion_report.tsv - per-cluster detailed stats
#' 2. deg_filtering_stats.tsv - filterByExpr breakdown
#' 3. deg_summary.txt - human-readable summary
#'
#' @param deg_results List of DEG results from custom_pseudoBulkDGE
#' @param sample_counts Table of sample counts per cluster/condition
#' @param comparison_name Name of the comparison (e.g., "Healthy_vs_IPF")
#' @param output_dir Output directory for report files
#' @param filter_params List with min.count, min.total.count, min.prop values used
generate_deg_completion_report <- function(deg_results, sample_counts, comparison_name,
                                           output_dir, filter_params) {

  message("\nGenerating DEG completion report for: ", comparison_name)

  # Initialize data frames for reports
  completion_report <- data.frame(
    comparison = character(),
    cluster = character(),
    status = character(),
    n_samples_group1 = integer(),
    n_samples_group2 = integer(),
    n_genes_input = integer(),
    n_genes_tested = integer(),
    n_genes_filtered = integer(),
    pct_genes_retained = numeric(),
    n_deg_total = integer(),
    n_deg_up = integer(),
    n_deg_down = integer(),
    reason_if_excluded = character(),
    stringsAsFactors = FALSE
  )

  filtering_stats <- data.frame(
    comparison = character(),
    cluster = character(),
    filter_step = character(),
    genes_removed = integer(),
    genes_remaining = integer(),
    stringsAsFactors = FALSE
  )

  # Get all clusters (from sample_counts)
  all_clusters <- rownames(sample_counts)
  failed_clusters <- S4Vectors::metadata(deg_results)$failed
  successful_clusters <- setdiff(all_clusters, failed_clusters)

  # Summary counters
  total_deg_up <- 0
  total_deg_down <- 0
  cluster_deg_counts <- list()

  # Process each cluster
  for (cluster in all_clusters) {
    cluster_counts <- sample_counts[cluster, ]
    group1_n <- as.integer(cluster_counts[1])
    group2_n <- as.integer(cluster_counts[2])

    if (cluster %in% failed_clusters) {
      # Determine failure reason
      if (any(cluster_counts == 0)) {
        reason <- paste0("no_samples_in_", names(cluster_counts)[cluster_counts == 0])
      } else if (any(cluster_counts < 2)) {
        reason <- "insufficient_replicates"
      } else {
        reason <- "design_matrix_singularity"
      }

      completion_report <- rbind(completion_report, data.frame(
        comparison = comparison_name,
        cluster = cluster,
        status = "excluded",
        n_samples_group1 = group1_n,
        n_samples_group2 = group2_n,
        n_genes_input = NA,
        n_genes_tested = NA,
        n_genes_filtered = NA,
        pct_genes_retained = NA,
        n_deg_total = NA,
        n_deg_up = NA,
        n_deg_down = NA,
        reason_if_excluded = reason,
        stringsAsFactors = FALSE
      ))
    } else {
      # Get results for this cluster
      cluster_result <- deg_results@listData[[cluster]]

      if (!is.null(cluster_result)) {
        n_genes_input <- nrow(cluster_result)
        n_genes_tested <- sum(!is.na(cluster_result$PValue))
        n_genes_filtered <- n_genes_input - n_genes_tested

        # Count DEGs (FDR < 0.05)
        sig_genes <- cluster_result[!is.na(cluster_result$FDR) & cluster_result$FDR < 0.05, ]
        n_deg_total <- nrow(sig_genes)
        n_deg_up <- sum(sig_genes$logFC > 0, na.rm = TRUE)
        n_deg_down <- sum(sig_genes$logFC < 0, na.rm = TRUE)

        total_deg_up <- total_deg_up + n_deg_up
        total_deg_down <- total_deg_down + n_deg_down
        cluster_deg_counts[[cluster]] <- n_deg_total

        completion_report <- rbind(completion_report, data.frame(
          comparison = comparison_name,
          cluster = cluster,
          status = "completed",
          n_samples_group1 = group1_n,
          n_samples_group2 = group2_n,
          n_genes_input = n_genes_input,
          n_genes_tested = n_genes_tested,
          n_genes_filtered = n_genes_filtered,
          pct_genes_retained = round(n_genes_tested / n_genes_input * 100, 1),
          n_deg_total = n_deg_total,
          n_deg_up = n_deg_up,
          n_deg_down = n_deg_down,
          reason_if_excluded = "",
          stringsAsFactors = FALSE
        ))

        # Add filtering stats entry
        filtering_stats <- rbind(filtering_stats, data.frame(
          comparison = comparison_name,
          cluster = cluster,
          filter_step = paste0("filterByExpr(min.count=", filter_params$min.count,
                               ",min.total.count=", filter_params$min.total.count,
                               ",min.prop=", filter_params$min.prop, ")"),
          genes_removed = n_genes_filtered,
          genes_remaining = n_genes_tested,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Write completion report
  report_file <- file.path(output_dir, "deg_completion_report.tsv")
  if (file.exists(report_file)) {
    existing_report <- read.table(report_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    completion_report <- rbind(existing_report, completion_report)
  }
  write.table(completion_report, file = report_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Saved: deg_completion_report.tsv")

  # Write filtering stats
  filter_file <- file.path(output_dir, "deg_filtering_stats.tsv")
  if (file.exists(filter_file)) {
    existing_filter <- read.table(filter_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    filtering_stats <- rbind(existing_filter, filtering_stats)
  }
  write.table(filtering_stats, file = filter_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Saved: deg_filtering_stats.tsv")

  # Write human-readable summary
  summary_file <- file.path(output_dir, "deg_summary.txt")
  summary_lines <- c(
    "=== DEG Analysis Summary ===",
    "",
    paste0("Comparison: ", comparison_name),
    paste0("  - Clusters analyzed: ", length(successful_clusters), "/", length(all_clusters)),
    paste0("  - Clusters excluded: ", length(failed_clusters)),
    paste0("  - Total DEGs (FDR < 0.05): ", total_deg_up + total_deg_down,
           " (up: ", total_deg_up, ", down: ", total_deg_down, ")")
  )

  # Find top cluster by DEGs
  if (length(cluster_deg_counts) > 0) {
    top_cluster <- names(which.max(unlist(cluster_deg_counts)))
    summary_lines <- c(summary_lines,
                       paste0("  - Top cluster by DEGs: ", top_cluster,
                              " (", cluster_deg_counts[[top_cluster]], " DEGs)"))
  }

  # Add excluded clusters section
  if (length(failed_clusters) > 0) {
    summary_lines <- c(summary_lines, "", "Excluded clusters:")
    for (fc in failed_clusters) {
      fc_reason <- completion_report$reason_if_excluded[completion_report$cluster == fc &
                                                          completion_report$comparison == comparison_name]
      fc_counts <- sample_counts[fc, ]
      summary_lines <- c(summary_lines,
                         paste0("  - ", fc, ": ", fc_reason,
                                " (", paste(names(fc_counts), "=", fc_counts, collapse = ", "), ")"))
    }
  }

  summary_lines <- c(summary_lines, "", paste0("Generated: ", Sys.time()))

  # Append to existing summary or create new
  if (file.exists(summary_file)) {
    existing_summary <- readLines(summary_file)
    summary_lines <- c(existing_summary, "", strrep("-", 60), "", summary_lines)
  }
  writeLines(summary_lines, summary_file)
  message("  Saved: deg_summary.txt")

  return(invisible(completion_report))
}

#' Custom pseudoBulkDGE implementation with configurable filterByExpr parameters
#'
#' This function implements pseudobulk DEG analysis using edgeR directly,
#' allowing custom filterByExpr thresholds. The scran::pseudoBulkDGE function
#' uses hardcoded filterByExpr defaults (min.count=10, min.total.count=15, min.prop=0.7)
#' which are too stringent for many scRNA-seq datasets.
#'
#' @param sce SingleCellExperiment object (already aggregated to pseudobulk)
#' @param label_col Column name for cell type labels
#' @param design_formula Design formula for the model
#' @param coef Coefficient to test
#' @param min.count Minimum count threshold for filterByExpr (default: 3)
#' @param min.total.count Minimum total count across all samples (default: 5)
#' @param min.prop Minimum proportion of samples with counts above threshold (default: 0.2)
#' @return List of DEG results per cluster, compatible with scran output format
custom_pseudoBulkDGE <- function(sce, label_col, design_formula, coef,
                                  min.count = 3, min.total.count = 5, min.prop = 0.2) {

  # Get unique cell type labels
  labels <- unique(colData(sce)[[label_col]])
  message("Processing ", length(labels), " cell types...")

  results <- list()
  failed <- character(0)

  for (label in labels) {
    message("  Processing: ", label)

    # Subset to current cell type
    idx <- which(colData(sce)[[label_col]] == label)
    if (length(idx) < 2) {
      message("    -> Skipped: fewer than 2 samples")
      failed <- c(failed, label)
      results[[label]] <- NULL
      next
    }

    current_sce <- sce[, idx]

    tryCatch({
      # Get counts matrix
      counts_mat <- assay(current_sce, "counts")

      # Create DGEList
      y <- DGEList(counts = counts_mat)

      # Get design matrix
      col_data <- as.data.frame(colData(current_sce))
      design <- model.matrix(design_formula, data = col_data)

      # Check design matrix validity
      if (ncol(design) > ncol(counts_mat)) {
        message("    -> Skipped: more coefficients than samples")
        failed <- c(failed, label)
        results[[label]] <- NULL
        next
      }

      # Check for rank deficiency and auto-drop redundant columns
      design_rank <- qr(design)$rank
      if (design_rank < ncol(design)) {
        message("    -> WARNING: Design matrix rank-deficient (rank=", design_rank,
                ", ncol=", ncol(design), "). Auto-dropping ", ncol(design) - design_rank, " redundant columns.")
        pivot_cols <- qr(design)$pivot[1:design_rank]
        dropped_cols <- colnames(design)[setdiff(seq_len(ncol(design)), pivot_cols)]
        message("    -> Dropped columns: ", paste(dropped_cols, collapse = ", "))
        design <- design[, pivot_cols, drop = FALSE]
      }

      # Check if coefficient exists in design
      if (!coef %in% colnames(design)) {
        # Try to find matching column
        matching_cols <- grep(gsub("^.*\\.", "", coef), colnames(design), value = TRUE)
        if (length(matching_cols) > 0) {
          coef <- matching_cols[1]
        } else {
          message("    -> Skipped: coefficient '", coef, "' not found in design")
          message("       Available coefficients: ", paste(colnames(design), collapse = ", "))
          failed <- c(failed, label)
          results[[label]] <- NULL
          next
        }
      }

      # Apply filterByExpr with custom thresholds
      # Extract condition for group-aware filtering
      group <- col_data[[names(col_data)[1]]]  # Use first column as group

      keep <- filterByExpr(y, design = design, group = group,
                           min.count = min.count,
                           min.total.count = min.total.count,
                           min.prop = min.prop)

      n_kept <- sum(keep)
      n_total <- length(keep)
      message("    -> Genes passing filter: ", n_kept, "/", n_total,
              " (", round(n_kept/n_total*100, 1), "%)")

      if (n_kept < 10) {
        message("    -> Skipped: fewer than 10 genes pass filter")
        failed <- c(failed, label)
        results[[label]] <- NULL
        next
      }

      y <- y[keep, , keep.lib.sizes = FALSE]

      # Normalize
      y <- calcNormFactors(y)

      # Estimate dispersion
      y <- estimateDisp(y, design, robust = TRUE)

      # Fit model
      fit <- glmQLFit(y, design, robust = TRUE)

      # Test coefficient
      qlf <- glmQLFTest(fit, coef = coef)

      # Extract results
      res <- topTags(qlf, n = Inf, sort.by = "none")$table

      # Create full result matrix (include filtered genes as NA)
      full_res <- data.frame(
        logFC = rep(NA_real_, n_total),
        logCPM = rep(NA_real_, n_total),
        F = rep(NA_real_, n_total),
        PValue = rep(NA_real_, n_total),
        FDR = rep(NA_real_, n_total),
        row.names = rownames(counts_mat)
      )

      # Fill in results for kept genes
      full_res[rownames(res), ] <- res

      # Store as DataFrame with metadata (compatible with scran output)
      result_df <- S4Vectors::DataFrame(full_res)
      S4Vectors::metadata(result_df)$y <- y  # Store DGEList for downstream use

      results[[label]] <- result_df
      message("    -> Completed: ", sum(!is.na(full_res$PValue)), " genes tested")

    }, error = function(e) {
      message("    -> Error: ", e$message)
      failed <<- c(failed, label)
      results[[label]] <<- NULL
    })
  }

  # Create output list with metadata
  output <- S4Vectors::SimpleList(results)
  S4Vectors::metadata(output)$failed <- failed

  message("\nDEG analysis complete:")
  message("  Successful: ", length(labels) - length(failed), " cell types")
  message("  Failed: ", length(failed), " cell types")

  return(output)
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

message("\n", strrep("=", 80))
message("PSEUDOBULK DEG ANALYSIS PIPELINE")
message(strrep("=", 80))

# Parse comparisons
comparison_pairs <- strsplit(opt$comparisons, ";")[[1]]
group_combinations <- matrix(unlist(lapply(comparison_pairs, function(x) strsplit(x, ",")[[1]])), nrow = 2)

# Setup output directory
output_dir <- opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_dir <- normalizePath(output_dir)  # Convert to absolute path before setwd
setwd(output_dir)

# Load custom signatures if provided (uses GENE SYMBOLS)
custom_signatures <- load_custom_signatures(opt$custom_signatures)

# Store column names
condition.obj <- opt$condition
celltype.obj <- opt$celltype
sample.obj <- opt$sample
batch.obj <- opt$batch

# Parse covariate(s) - supports comma-separated list for multiple covariates
covariate_raw <- opt$covariate
if (!is.null(covariate_raw) && covariate_raw != "NULL" && nchar(trimws(covariate_raw)) > 0) {
  covariate.obj_list <- trimws(strsplit(covariate_raw, ",")[[1]])
  covariate.obj_list <- covariate.obj_list[nchar(covariate.obj_list) > 0]
} else {
  covariate.obj_list <- character(0)
}
# Keep single covariate.obj for backward compatibility (first covariate or NULL)
covariate.obj <- if (length(covariate.obj_list) > 0) covariate.obj_list[1] else NULL

# Create safe variable names for formulas
condition <- make.names(condition.obj)
celltype <- make.names(celltype.obj)
sample <- make.names(sample.obj)
batch <- if (!is.null(batch.obj) && batch.obj != "NULL") make.names(batch.obj) else NULL
# Create safe names for all covariates
covariates <- if (length(covariate.obj_list) > 0) make.names(covariate.obj_list) else character(0)
# Legacy single covariate (backward compat)
covariate <- if (length(covariates) > 0) covariates[1] else NULL

# Print configuration
message("\n=== Analysis Configuration ===")
message("Selected modules: ", paste(selected_modules, collapse = ", "))
if (opt$fgsea_only) {
  message("Mode: fgsea_only (rerunning on existing DEG results)")
  message("DEG input directory: ", opt$deg_input)
} else {
  message("Input file: ", opt$input)
}
message("Output directory: ", output_dir)
if (!opt$fgsea_only) {
  message("Condition column: ", condition.obj, " -> ", condition)
  message("Cell type column: ", celltype.obj, " -> ", celltype)
  message("Sample column: ", sample.obj, " -> ", sample)
  message("Batch column: ", ifelse(is.null(batch.obj) || batch.obj == "NULL", "None", batch.obj))
  message("Covariate columns: ", ifelse(length(covariate.obj_list) == 0, "None", paste(covariate.obj_list, collapse = ", ")))
  message("Number of comparisons: ", ncol(group_combinations))
  for (i in 1:ncol(group_combinations)) {
    message("  Comparison ", i, ": ", group_combinations[1, i], " vs ", group_combinations[2, i])
  }
  message("filterByExpr: min.count=", opt$min_count, ", min.total.count=", opt$min_total_count, ", min.prop=", opt$min_prop)
  message("Volcano plot cutoffs: FDR < ", opt$volcano_fdr, ", |log2FC| > ", opt$volcano_log2fc)
}
message("Custom signatures: ", ifelse(is.null(custom_signatures), "None", paste(names(custom_signatures), collapse = ", ")))

# ==============================================================================
# Step 1: Load Data (skip if fgsea_only mode)
# ==============================================================================

if (!opt$fgsea_only) {
  message("\n", strrep("-", 60))
  message("STEP 1: Loading input data...")
  message(strrep("-", 60))

  seurat_obj <- readRDS(opt$input)
  message("Loaded Seurat object with ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " genes")

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
  message("Converted to SingleCellExperiment")
} else {
  message("\n", strrep("-", 60))
  message("STEP 1: SKIPPED (fgsea_only mode - using existing DEG results)")
  message(strrep("-", 60))
}

# ==============================================================================
# Step 2: Pseudobulk Aggregation (needed for DEG and Speckle)
# ==============================================================================

# Only run aggregation if DEG or speckle is selected (and not fgsea_only mode)
needs_aggregation <- !opt$fgsea_only && (("deg" %in% selected_modules) || ("speckle" %in% selected_modules))

if (needs_aggregation) {
  message("\n", strrep("-", 60))
  message("STEP 2: Pseudobulk aggregation...")
  message(strrep("-", 60))

  # Build aggregation columns (include all covariates for metadata preservation)
  agg_cols.obj <- c(condition.obj, celltype.obj, sample.obj)
  if (length(covariate.obj_list) > 0) agg_cols.obj <- c(agg_cols.obj, covariate.obj_list)
  if (!is.null(batch.obj) && batch.obj != "NULL" && !(batch.obj %in% agg_cols.obj)) {
    agg_cols.obj <- c(agg_cols.obj, batch.obj)
  }

  # Aggregate cell counts
  # BPCells on-disk matrices can't be aggregated directly (not a dgCMatrix).
  # Materialize per-cluster subsets on-the-fly — each cluster fits in memory
  # (largest cluster ~435K cells, ~430M nnz, well under 2.1B int32 limit).
  if (requireNamespace("BPCells", quietly = TRUE) && inherits(counts(sce), "IterableMatrix")) {
    message("  BPCells detected: aggregating per-cluster to avoid 2.1B nnz limit...")
    clusters <- unique(colData(sce)[[celltype.obj]])
    agg_list <- list()
    for (cl in clusters) {
      idx <- which(colData(sce)[[celltype.obj]] == cl)
      message("    ", cl, ": ", length(idx), " cells")
      sce_sub <- sce[, idx]
      counts(sce_sub) <- as(counts(sce_sub), "dgCMatrix")
      agg_list[[cl]] <- aggregateAcrossCells(
        sce_sub,
        id = colData(sce_sub)[, agg_cols.obj],
        use.assay.type = "counts"
      )
      rm(sce_sub); gc(verbose = FALSE)
    }
    sce_aggregated <- do.call(cbind, agg_list)
    rm(agg_list); gc(verbose = FALSE)
  } else {
    sce_aggregated <- aggregateAcrossCells(
      sce,
      id = colData(sce)[, agg_cols.obj],
      use.assay.type = "counts"
    )
  }
  message("Aggregated to ", ncol(sce_aggregated), " pseudobulk samples")

  # QC check
  table_clusters <- as.data.frame.matrix(table(sce_aggregated[[condition.obj]], sce_aggregated[[celltype.obj]]))
  write.table(table_clusters, "pseudobulk_aggregation_qc.tsv", quote = FALSE, sep = "\t")
  message("Saved aggregation QC table")
} else {
  message("\n", strrep("-", 60))
  message("STEP 2: SKIPPED (not needed for selected modules)")
  message(strrep("-", 60))
}

# ==============================================================================
# Step 3: Pseudobulk Differential Expression (conditional on 'deg' module)
# ==============================================================================

# Initialize result containers (needed for fgsea even if DEG is skipped)
de_results_dataframe <- data.frame()
de_metrics_list <- list()
result_name <- NULL

if ("deg" %in% selected_modules && !opt$fgsea_only) {
  message("\n", strrep("-", 60))
  message("STEP 3: Pseudobulk differential expression analysis...")
  message(strrep("-", 60))

  # Build design formula with all covariates
  all_formula_terms <- c(condition)
  if (length(covariates) > 0) all_formula_terms <- c(all_formula_terms, covariates)
  if (!is.null(batch) && !(batch %in% all_formula_terms)) {
    all_formula_terms <- c(all_formula_terms, batch)
  }
  formula <- as.formula(paste("~", paste(all_formula_terms, collapse = " + ")))
  message("Design formula: ", deparse(formula))

  # Create volcano plots directory
  if (!dir.exists("volcano_plots")) {
    dir.create("volcano_plots")
  }

  # Process each comparison
  for (i in 1:dim(group_combinations)[2]) {
    group1 <- group_combinations[1, i]
    group2 <- group_combinations[2, i]

    message("\n=== Comparison ", i, ": ", group1, " vs ", group2, " ===")

    # Filter to current groups
    current_sce <- sce_aggregated[, sce_aggregated[[condition.obj]] %in% c(group1, group2)]
    current_sce[[condition]] <- relevel(factor(current_sce[[condition.obj]]), ref = group2)

    # Pre-analysis diagnostics
    cluster_labels <- colData(current_sce)[, celltype.obj]
    condition_values <- current_sce[[condition.obj]]
    sample_counts <- table(cluster_labels, condition_values)

    message("Sample counts per cluster per condition:")
    print(sample_counts)

    # Check for problematic clusters
    clusters_missing_condition <- rownames(sample_counts)[
      apply(sample_counts, 1, function(x) any(x == 0))
    ]
    if (length(clusters_missing_condition) > 0) {
      message("WARNING: Clusters with ZERO samples in one condition:")
      for (cl in clusters_missing_condition) {
        message("  - ", cl, ": ", paste(colnames(sample_counts), "=", sample_counts[cl, ], collapse = ", "))
      }
    }

    clusters_low_replicates <- rownames(sample_counts)[
      apply(sample_counts, 1, function(x) any(x > 0 & x < 2))
    ]
    if (length(clusters_low_replicates) > 0) {
      message("WARNING: Clusters with <2 replicates:")
      for (cl in clusters_low_replicates) {
        message("  - ", cl, ": ", paste(colnames(sample_counts), "=", sample_counts[cl, ], collapse = ", "))
      }
    }

    # Check batch/covariate confounding
    if (!is.null(batch.obj) && batch.obj != "NULL") {
      batch_values <- current_sce[[batch.obj]]
      confound_table <- table(condition_values, batch_values)
      batch_confounded <- apply(confound_table, 2, function(x) sum(x > 0) == 1)
      if (any(batch_confounded)) {
        message("WARNING: BATCH CONFOUNDING DETECTED!")
        for (b in names(batch_confounded)[batch_confounded]) {
          cond_with_batch <- rownames(confound_table)[confound_table[, b] > 0]
          message("  Batch '", b, "' only in condition '", cond_with_batch, "'")
        }
      }
    }

    # Run custom pseudobulk DEG with relaxed filterByExpr thresholds
    # The scran::pseudoBulkDGE doesn't expose filterByExpr params, so we implement directly
    message("Running pseudobulk DEG with custom filterByExpr thresholds...")
    message("  filterByExpr: min.count=", opt$min_count, ", min.total.count=", opt$min_total_count, ", min.prop=", opt$min_prop)

    current_de_result <- custom_pseudoBulkDGE(
      current_sce,
      label_col = celltype.obj,
      design_formula = formula,
      coef = paste0(condition, make.names(group1)),
      min.count = opt$min_count,
      min.total.count = opt$min_total_count,
      min.prop = opt$min_prop
    )

    # Report failed clusters
    failed_clusters <- metadata(current_de_result)$failed
    if (length(failed_clusters) > 0) {
      message("\nFAILED CLUSTERS (", length(failed_clusters), " of ", length(unique(cluster_labels)), "):")
      for (fc in failed_clusters) {
        fc_counts <- sample_counts[fc, ]
        if (any(fc_counts == 0)) {
          reason <- paste0("Cell type only in condition '", names(fc_counts)[fc_counts > 0], "' (n=", fc_counts[fc_counts > 0], ")")
        } else if (any(fc_counts < 2)) {
          reason <- paste0("Insufficient replicates: ", paste(names(fc_counts), "=", fc_counts, collapse = ", "), " (need >=2 each)")
        } else {
          reason <- "Design matrix singularity (likely batch/covariate confounding)"
        }
        message("  - ", fc, ": ", reason)
      }
    }

    # Generate DEG completion report
    comparison_name <- paste0(group1, "_vs_", group2)
    filter_params <- list(
      min.count = opt$min_count,
      min.total.count = opt$min_total_count,
      min.prop = opt$min_prop
    )
    generate_deg_completion_report(current_de_result, sample_counts, comparison_name,
                                   output_dir, filter_params)

    # Extract results per cluster
    de_results_clusters <- current_de_result
    for (cluster in names(de_results_clusters)) {
      message("Extracting results for cluster ", cluster)

      current_de_results <- de_results_clusters@listData[[cluster]]

      if (is.null(current_de_results)) {
        message("  -> NULL result (see FAILED CLUSTERS above)")
        next
      }

      # Report filtering statistics
      n_genes_total <- nrow(current_de_results)
      n_genes_tested <- sum(!is.na(current_de_results$PValue))
      n_genes_filtered <- n_genes_total - n_genes_tested
      if (n_genes_filtered > 0) {
        message("  -> ", n_genes_tested, "/", n_genes_total, " genes tested (", n_genes_filtered, " filtered)")
      }

      cluster <- toString(cluster)

      # Store metrics
      de_metrics_list <- fill_edge_metrics_by_groups(de_metrics_list, cluster, group1, group2, metadata(current_de_results)$y)

      # Modify and store results
      current_de_results <- modify_DE_results(current_de_results, group1, group2, cluster, "cluster")
      de_results_dataframe <- update_DE_results_df(de_results_dataframe, current_de_results)

      # Generate volcano plots (static PNG and interactive HTML)
      render_volcano_plot(current_de_results, group1, group2, "volcano_plots", "DEG", cluster)
      render_interactive_volcano_plot(current_de_results, group1, group2, "volcano_plots", "DEG", cluster,
                                      fdr_cutoff = opt$volcano_fdr, log2fc_cutoff = opt$volcano_log2fc)
    }
  }

  # Save DE results
  message("\nSaving DE results...")
  result_name <- paste0("pseudobulk_", basename(tools::file_path_sans_ext(opt$input)))
  de_results_dataframe$comparison <- paste0(de_results_dataframe$group1, "_vs_", de_results_dataframe$group2)
  de_results_out <- de_results_dataframe[, c("gene_name", "logFC", "PValue", "FDR", "group1", "group2", "cluster", "comparison", "clust_res")]
  write.table(de_results_out, file = paste0(result_name, "_DE.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  message("Saved: ", result_name, "_DE.tsv")

} else if (opt$fgsea_only) {
  # Load existing DEG results for fgsea rerun
  message("\n", strrep("-", 60))
  message("STEP 3: Loading existing DEG results for fgsea rerun...")
  message(strrep("-", 60))

  deg_files <- list.files(opt$deg_input, pattern = "_DE\\.tsv$", full.names = TRUE)
  message("Loading DEG results from: ", deg_files[1])

  de_results_dataframe <- read.table(deg_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  de_results_dataframe$gene_name <- de_results_dataframe$gene_name
  result_name <- gsub("_DE\\.tsv$", "", basename(deg_files[1]))

  message("Loaded ", nrow(de_results_dataframe), " DEG entries")
  message("Comparisons found: ", paste(unique(de_results_dataframe$comparison), collapse = ", "))

} else {
  message("\n", strrep("-", 60))
  message("STEP 3: SKIPPED (deg module not selected)")
  message(strrep("-", 60))

  # Auto-load existing DEG results if fgsea modules need them
  if (("fgsea_groups" %in% selected_modules || "fgsea_clusters" %in% selected_modules) &&
      nrow(de_results_dataframe) == 0) {
    deg_search_dir <- if (!is.null(opt$deg_input)) opt$deg_input else output_dir
    deg_files <- list.files(deg_search_dir, pattern = "_DE\\.tsv$", full.names = TRUE)
    if (length(deg_files) > 0) {
      message("  Auto-loading existing DEG results from: ", deg_files[1])
      de_results_dataframe <- read.table(deg_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      result_name <- gsub("_DE\\.tsv$", "", basename(deg_files[1]))
      message("  Loaded ", nrow(de_results_dataframe), " DEG entries")
      message("  Comparisons: ", paste(unique(de_results_dataframe$comparison), collapse = ", "))
    }
  }
}

# ==============================================================================
# Step 4: Pathway Enrichment Analysis (conditional on 'fgsea_groups' module)
# ==============================================================================

if ("fgsea_groups" %in% selected_modules) {

  # Initialize checkpoint for fgsea results (saved to output_dir)
  fgsea_checkpoint_file <- file.path(output_dir, "fgsea_checkpoint.rds")
  if (file.exists(fgsea_checkpoint_file)) {
    fgsea_all_results <- readRDS(fgsea_checkpoint_file)
    message("Loaded existing checkpoint: ", fgsea_checkpoint_file)
  } else {
    fgsea_all_results <- list()
  }

  # Check if DEG results are available
  if (nrow(de_results_dataframe) == 0) {
    message("\n", strrep("-", 60))
    message("STEP 4: ERROR - No DEG results available for pathway enrichment")
    message("  Run with 'deg' module first, or use --deg_input to load existing results")
    message(strrep("-", 60))
  } else {

    message("\n", strrep("-", 60))
    message("STEP 4: Pathway enrichment analysis (comparison-level)...")
    message(strrep("-", 60))

    # Create output directories
    if (!dir.exists("fgsea_groups_clusters")) dir.create("fgsea_groups_clusters")
    if (!dir.exists("fgsea_clusters")) dir.create("fgsea_clusters")

# ============================================================================
# IMPORTANT: Gene Symbol vs Entrez ID handling
# - Reactome uses ENTREZ IDs -> we convert gene symbols to Entrez IDs
# - MSigDB uses GENE SYMBOLS -> we keep original gene symbols
# - Custom signatures use GENE SYMBOLS -> same as MSigDB
# ============================================================================

# Prepare dataframes - keep gene symbols for MSigDB and custom signatures
de_df.msigdb <- de_results_dataframe
names(de_df.msigdb)[names(de_df.msigdb) == "gene_name"] <- "feature"

# Reactome pathway analysis (requires Entrez IDs)
message("\nRunning Reactome pathway analysis (converting to Entrez IDs)...")
entrez.db <- org.Hs.eg.db
de_results_entrez <- convert_to_entrez_ids_df(de_results_dataframe, entrez.db)
g.entrezid <- mapIds(entrez.db,
                     keys = de_results_dataframe$gene_name,
                     keytype = "SYMBOL",
                     column = "ENTREZID")

fgsea.set <- reactomePathways(unique(de_results_entrez$feature))
de_results_entrez <- calc_adjusted_pval(de_results_entrez)

# Fully vectorized leadingEdge ID conversion (~3 min for 1.47M entries)
map_genes_parallel <- function(leadingEdge_col, gene_map) {
  n <- length(leadingEdge_col)
  message("  Converting ", n, " leadingEdge entries to gene symbols...")
  message("  Using fully vectorized flatten+lookup+reassemble")

  splits <- strsplit(leadingEdge_col, ",", fixed = TRUE)
  lens <- lengths(splits)
  nonempty <- which(lens > 0)

  result <- rep("", n)

  if (length(nonempty) > 0) {
    flat <- unlist(splits[nonempty], use.names = FALSE)
    mapped <- gene_map[flat]
    mapped[is.na(mapped)] <- "\x01"
    mapped[mapped == ""] <- "\x01"

    grp <- rep.int(seq_along(nonempty), lens[nonempty])
    reassembled <- tapply(mapped, grp, paste, collapse = ",", simplify = TRUE)

    cleaned <- gsub("\x01,|,\x01|\x01", "", unname(reassembled))
    cleaned <- gsub(",+", ",", cleaned)
    cleaned <- gsub("^,|,$", "", cleaned)
    result[nonempty] <- cleaned
  }

  message("  Conversion complete")
  return(result)
}

reactome_file_namings_info <- list(
  folder_name = "fgsea_groups_clusters",
  level = "condition_clusters",
  db_name = "reactome"
)

# Split by comparison and convert to regular R list
# (split on DataFrame returns CompressedSplitDFrameList which is.list() returns FALSE for)
de_results_dataframe_list <- as.list(split(de_results_entrez, f = de_results_entrez$comparison))
# Also convert each element to regular data.frame
de_results_dataframe_list <- lapply(de_results_dataframe_list, as.data.frame)

# If empty, skip Reactome analysis
if (length(de_results_dataframe_list) == 0) {
  message("WARNING: No data for Reactome analysis after Entrez ID conversion - skipping")
} else {
  file.name <- file.path(".", reactome_file_namings_info$folder_name, paste0(reactome_file_namings_info$level, "_", reactome_file_namings_info$db_name))
  de.fgsea <- performGSEAmultilevel_by_groups(de_df_list = de_results_dataframe_list,
                                            fgsea.sets = list(fgsea.set),
                                            file_namings_info = reactome_file_namings_info,
                                            analysis_type = result_name)
  # Checkpoint: Save raw Reactome results (Entrez IDs)
  fgsea_all_results$groups_reactome <- de.fgsea
  saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
  message("Checkpoint saved: groups_reactome (raw)")

  # Convert Entrez IDs to gene symbols in leadingEdge with error handling
  message("Converting Reactome leadingEdge from Entrez IDs to gene symbols...")
  tryCatch({
    gene_map <- setNames(names(g.entrezid), g.entrezid)
    de.fgsea$leadingEdge <- map_genes_parallel(de.fgsea$leadingEdge, gene_map)
  }, error = function(e) {
    message("WARNING: Could not convert leadingEdge to gene symbols: ", e$message)
    message("LeadingEdge will remain as Entrez IDs")
  })
  write.table(de.fgsea, file = paste0(file.name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved Reactome enrichment results: ", paste0(file.name, "_cluster.tsv"))
  # Checkpoint: Update with converted gene symbols
  fgsea_all_results$groups_reactome <- de.fgsea
  saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
  message("Checkpoint updated: groups_reactome (converted)")
}

# MSigDB pathway analysis (uses GENE SYMBOLS)
message("\nRunning MSigDB pathway analysis (uses gene symbols)...")
fgsea.sets <- get_fgsea_sets_msigdbr()
de_df.msigdb <- calc_adjusted_pval(de_df.msigdb)

msigdb_file_namings_info <- list(
  folder_name = "fgsea_groups_clusters",
  level = "condition_clusters",
  db_name = "msigdb"
)

# Split by comparison and convert to regular R list
de_df.msigdb_list <- as.list(split(de_df.msigdb, f = de_df.msigdb$comparison))
de_df.msigdb_list <- lapply(de_df.msigdb_list, as.data.frame)

file.name <- file.path(".", msigdb_file_namings_info$folder_name, paste0(msigdb_file_namings_info$level, "_", msigdb_file_namings_info$db_name))
de.fgsea <- performGSEAmultilevel_by_groups(de_df_list = de_df.msigdb_list,
                                            fgsea.sets = fgsea.sets,
                                            file_namings_info = msigdb_file_namings_info,
                                            analysis_type = result_name)
# Checkpoint: Save MSigDB results (already gene symbols)
fgsea_all_results$groups_msigdb <- de.fgsea
saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
message("Checkpoint saved: groups_msigdb")

write.table(de.fgsea, file = paste0(file.name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
message("Saved MSigDB enrichment results")

# ============================================================================
# CUSTOM SIGNATURE ENRICHMENT (uses GENE SYMBOLS - same as MSigDB)
# Custom signatures are analyzed using the SAME ranked gene list as MSigDB
# ============================================================================
if (!is.null(custom_signatures)) {
  message("\n", strrep("-", 40))
  message("Running CUSTOM SIGNATURE enrichment analysis...")
  message("  Using GENE SYMBOLS (same as MSigDB)")
  message("  Same ranked gene list as MSigDB analysis")
  message(strrep("-", 40))

  custom_file_namings_info <- list(
    folder_name = "fgsea_groups_clusters",
    level = "condition_clusters",
    db_name = "custom"
  )

  # Use de_df.msigdb_list which has gene SYMBOLS (not Entrez IDs)
  file.name <- file.path(".", custom_file_namings_info$folder_name, paste0(custom_file_namings_info$level, "_", custom_file_namings_info$db_name))
  de.fgsea.custom <- performGSEAmultilevel_by_groups(
    de_df_list = de_df.msigdb_list,  # Uses gene symbols
    fgsea.sets = list(custom_signatures),  # Custom signatures with gene symbols
    file_namings_info = custom_file_namings_info,
    analysis_type = result_name
  )
  # Checkpoint: Save custom results
  fgsea_all_results$groups_custom <- de.fgsea.custom
  saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
  message("Checkpoint saved: groups_custom")

  write.table(de.fgsea.custom, file = paste0(file.name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)

  # Also save to root output directory for easy access
  write.table(de.fgsea.custom, file = "fgsea_custom_signatures.tsv", quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved custom signature enrichment results")
  message("  -> fgsea_groups_clusters/condition_clusters_custom_cluster.tsv")
  message("  -> fgsea_custom_signatures.tsv")
}

  }  # End of DEG results check (nrow > 0)
} else {  # fgsea_groups module not selected
  message("\n", strrep("-", 60))
  message("STEP 4: SKIPPED (fgsea_groups module not selected)")
  message(strrep("-", 60))
}

# ==============================================================================
# Step 5: Cluster-level Pathway Analysis (conditional on 'fgsea_clusters' module)
# ==============================================================================

if ("fgsea_clusters" %in% selected_modules && !opt$fgsea_only) {
  message("\n", strrep("-", 60))
  message("STEP 5: Cluster-level pathway analysis...")
  message(strrep("-", 60))

  # Initialize checkpoint if not already done in Step 4
  if (!exists("fgsea_checkpoint_file")) {
    fgsea_checkpoint_file <- file.path(output_dir, "fgsea_checkpoint.rds")
    if (file.exists(fgsea_checkpoint_file)) {
      fgsea_all_results <- readRDS(fgsea_checkpoint_file)
      message("Loaded existing checkpoint: ", fgsea_checkpoint_file)
    } else {
      fgsea_all_results <- list()
    }
  }

# Perform Wilcoxon test on clusters
# Note: presto::wilcoxauc may fail with newer Seurat versions (5.0+)
# due to deprecated 'slot' argument in GetAssayData
x.genes <- tryCatch({
  wilcoxauc(seurat_obj, celltype.obj, seurat_assay = "RNA")
}, error = function(e) {
  message("WARNING: presto::wilcoxauc failed (likely Seurat 5.0+ compatibility issue)")
  message("  Error: ", e$message)
  message("  Falling back to Seurat::FindAllMarkers...")

  # Fallback using Seurat's FindAllMarkers
  # Filter out NA/NaN cell types before setting Idents
  valid_cells <- !is.na(seurat_obj@meta.data[[celltype.obj]]) &
                 seurat_obj@meta.data[[celltype.obj]] != "NaN" &
                 seurat_obj@meta.data[[celltype.obj]] != ""
  seurat_obj_filtered <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
  Idents(seurat_obj_filtered) <- seurat_obj_filtered@meta.data[[celltype.obj]]

  markers <- FindAllMarkers(seurat_obj_filtered, only.pos = FALSE, min.pct = 0.1,
                            logfc.threshold = 0.1, test.use = "wilcox")

  # Convert to presto-like format
  markers$feature <- markers$gene
  markers$group <- markers$cluster
  markers$auc <- (markers$pct.1 - markers$pct.2 + 1) / 2  # Approximate AUC
  markers$pval <- markers$p_val
  markers$padj <- markers$p_val_adj
  markers$logFC <- markers$avg_log2FC

  # Filter out any remaining NaN/NA clusters
  markers <- markers[!is.na(markers$group) & markers$group != "NaN", ]

  markers[, c("feature", "group", "auc", "pval", "padj", "logFC")]
})

if (is.null(x.genes) || nrow(x.genes) == 0) {
  message("WARNING: No cluster markers found. Skipping Step 5.")
} else {
  x.genes.msigdb <- x.genes  # Keep gene symbols for MSigDB

# Map to Entrez IDs for Reactome (only)
entrez.db <- org.Hs.eg.db
g.entrezid <- mapIds(entrez.db,
                     keys = x.genes$feature,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
x.genes$feature <- g.entrezid
x.genes <- x.genes[complete.cases(x.genes), ]

fgsea.set <- reactomePathways(unique(g.entrezid))

# Reactome cluster-level analysis
sample_name <- "./fgsea_clusters/clusters_reactome"
pdf("./fgsea_clusters/clusters_reactome_pathways.pdf", width = 11, height = 5)
df.vec <- do.call(rbind,
                  lapply(seq_along(unique(x.genes$group)),
                         function(n) {
                           performGSEAmultilevel(markers = x.genes, fgsea.set = fgsea.set, metric = "auc",
                                                 cluster = unique(x.genes$group)[n], n.lines = 10, minS = 1, maxS = Inf, sSize = 101, gseaP = 1,
                                                 scoreT = "pos", eps = 1*10^-10, coll.path = FALSE, sample_name = sample_name)
                         }))
dev.off()
# Checkpoint: Save raw Reactome cluster results (Entrez IDs)
fgsea_all_results$clusters_reactome <- df.vec
saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
message("Checkpoint saved: clusters_reactome (raw)")

# Convert Entrez IDs to gene symbols in leadingEdge with error handling
message("Converting cluster-level Reactome leadingEdge to gene symbols...")
tryCatch({
  gene_map <- setNames(names(g.entrezid), g.entrezid)
  df.vec$leadingEdge <- map_genes_parallel(df.vec$leadingEdge, gene_map)
}, error = function(e) {
  message("WARNING: Could not convert leadingEdge to gene symbols: ", e$message)
})
write.table(df.vec, file = paste0(sample_name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
message("Saved cluster-level Reactome results")
# Checkpoint: Update with converted gene symbols
fgsea_all_results$clusters_reactome <- df.vec
saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
message("Checkpoint updated: clusters_reactome (converted)")

# MSigDB C8 cluster-level analysis (uses gene symbols)
# Note: msigdbr 10.0.0+ uses 'collection' instead of deprecated 'category'
msig.df <- msigdbr(species = "Homo sapiens", collection = "C8")
fgsea.set <- msig.df %>% split(x = .$gene_symbol, f = .$gs_name)

pdf("./fgsea_clusters/clusters_msigdbC8_pathways.pdf", width = 11, height = 5)
sample_name <- "./fgsea_clusters/clusters_msigdbC8"
df.vec <- do.call(rbind,
                  lapply(seq_along(unique(x.genes.msigdb$group)),
                         function(n) {
                           performGSEAmultilevel(markers = x.genes.msigdb, fgsea.set = fgsea.set, metric = "auc",
                                                 cluster = unique(x.genes.msigdb$group)[n], n.lines = 10, minS = 1, maxS = Inf, sSize = 101, gseaP = 1,
                                                 scoreT = "pos", eps = 1*10^-10, coll.path = FALSE, sample_name = sample_name)
                         }))
dev.off()
# Checkpoint: Save MSigDB C8 cluster results
fgsea_all_results$clusters_msigdbC8 <- df.vec
saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
message("Checkpoint saved: clusters_msigdbC8")

write.table(df.vec, file = paste0(sample_name, "_clusters.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
message("Saved cluster-level MSigDB C8 results")

# Custom signatures at cluster level (uses gene symbols)
if (!is.null(custom_signatures)) {
  message("\nRunning cluster-level custom signature analysis...")
  pdf("./fgsea_clusters/clusters_custom_pathways.pdf", width = 11, height = 5)
  sample_name <- "./fgsea_clusters/clusters_custom"
  df.vec.custom <- do.call(rbind,
                           lapply(seq_along(unique(x.genes.msigdb$group)),
                                  function(n) {
                                    performGSEAmultilevel(markers = x.genes.msigdb, fgsea.set = custom_signatures, metric = "auc",
                                                          cluster = unique(x.genes.msigdb$group)[n], n.lines = 10, minS = 1, maxS = Inf, sSize = 101, gseaP = 1,
                                                          scoreT = "pos", eps = 1*10^-10, coll.path = FALSE, sample_name = sample_name)
                                  }))
  dev.off()
  # Checkpoint: Save custom cluster results
  fgsea_all_results$clusters_custom <- df.vec.custom
  saveRDS(fgsea_all_results, file = fgsea_checkpoint_file)
  message("Checkpoint saved: clusters_custom")

  write.table(df.vec.custom, file = paste0(sample_name, "_clusters.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved cluster-level custom signature results")
}

  }  # End of if (x.genes) check for Step 5
} else {
  message("\n", strrep("-", 60))
  message("STEP 5: SKIPPED (fgsea_clusters module not selected or fgsea_only mode)")
  message(strrep("-", 60))
}

# ==============================================================================
# Step 6: Cluster Markers (conditional on 'markers' module)
# ==============================================================================

if ("markers" %in% selected_modules && !opt$fgsea_only) {
  message("\n", strrep("-", 60))
  message("STEP 6: Cluster marker identification...")
  message(strrep("-", 60))

  # Filter out NA/NaN cell types and check for 2+ unique groups
  valid_celltypes <- sce[[celltype.obj]][!is.na(sce[[celltype.obj]]) &
                                          sce[[celltype.obj]] != "NaN" &
                                          sce[[celltype.obj]] != ""]
  unique_celltypes <- unique(valid_celltypes)

  if (length(unique_celltypes) < 2) {
    message("WARNING: Need at least 2 unique cell types for marker identification.")
    message("  Found ", length(unique_celltypes), " unique cell type(s). Skipping Step 6.")
  } else {
    # Filter SCE to only include valid cell types
    valid_cells <- !is.na(sce[[celltype.obj]]) &
                   sce[[celltype.obj]] != "NaN" &
                   sce[[celltype.obj]] != ""
    sce_filtered <- sce[, valid_cells]

    # BPCells: use Seurat::FindAllMarkers (supports on-disk matrices natively)
    # scran::findMarkers requires in-memory dgCMatrix
    use_seurat_markers <- requireNamespace("BPCells", quietly = TRUE) &&
                          inherits(counts(sce_filtered), "IterableMatrix")

    if (use_seurat_markers) {
      message("  BPCells detected: using Seurat::FindAllMarkers instead of scran::findMarkers")
      Idents(seurat_obj) <- seurat_obj@meta.data[[celltype.obj]]
      valid_seurat <- !is.na(Idents(seurat_obj)) & Idents(seurat_obj) != "NaN" & Idents(seurat_obj) != ""
      seurat_filt <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_seurat])
      sam <- FindAllMarkers(seurat_filt, only.pos = FALSE, min.pct = 0.1,
                            logfc.threshold = 0.1, test.use = "wilcox")
      combined_med_markers <- data.frame(
        feature = sam$gene,
        cluster = sam$cluster,
        p.adjusted = sam$p_val_adj,
        p.value = sam$p_val,
        summary.logFC = sam$avg_log2FC,
        Top = 0,
        self.average = 0,
        stringsAsFactors = FALSE
      )
    } else {
      if (!is.null(batch.obj) && batch.obj != "NULL") {
        med_markers <- scran::findMarkers(sce_filtered, sce_filtered[[celltype.obj]], block = sce_filtered[[batch.obj]])
      } else {
        med_markers <- scran::findMarkers(sce_filtered, sce_filtered[[celltype.obj]])
      }

      for (i in 1:length(names(med_markers))) {
        table <- as.data.frame(med_markers[[i]][, 1:4])
        table$cluster <- names(med_markers)[i]
        table <- table %>% rownames_to_column((var <- "feature"))
        if (i == 1) {
          old_table <- table
        } else {
          new_table <- rbind(old_table, table)
          old_table <- new_table
        }
      }
      combined_med_markers <- new_table
      rownames(combined_med_markers) <- NULL

      p.val <- as.numeric(combined_med_markers[, 3])
      p.adj <- p.adjust(p.val, method = "bonferroni")
      combined_med_markers$p.adjusted <- p.adj

      combined_med_markers <- combined_med_markers %>% dplyr::select(1, 6, 7, 3, 4, 5, 2)
    }

  write.table(combined_med_markers,
              file = "combined_markers_clusters.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Saved cluster markers")
  }  # End of 2+ unique cell types check
} else {
  message("\n", strrep("-", 60))
  message("STEP 6: SKIPPED (markers module not selected or fgsea_only mode)")
  message(strrep("-", 60))
}

# ==============================================================================
# Step 7: Cell Composition Analysis (Speckle) - conditional on 'speckle' module
# ==============================================================================

if ("speckle" %in% selected_modules && !opt$fgsea_only) {
  message("\n", strrep("-", 60))
  message("STEP 7: Cell composition analysis (Speckle)...")
  message(strrep("-", 60))

  if (!dir.exists("speckle_diffprop")) dir.create("speckle_diffprop")

  # Always use sample as the grouping variable for cell proportions
  # Batch will be used as a covariate in the design matrix if provided
  group_var.obj <- sample.obj

  # Get transformed proportions
  props <- getTransformedProps(sce[[celltype.obj]], sce[[group_var.obj]], transform = "logit")

  # Create metadata dataframe
  metadata_df <- as.data.frame(colData(sce))

  # Include condition, sample, and all covariates/batch
  meta_cols <- c(condition.obj, group_var.obj)
  if (!is.null(batch.obj) && batch.obj != "NULL") meta_cols <- c(meta_cols, batch.obj)
  if (length(covariate.obj_list) > 0) meta_cols <- c(meta_cols, covariate.obj_list)
  meta_cols <- unique(meta_cols)

  subset_df <- metadata_df[, meta_cols, drop = FALSE]
  unique_values_df <- subset_df %>%
    group_by(!!!syms(group_var.obj)) %>%
    summarize_all(list(~ unique(.)[1]))
  unique_values_df <- unique_values_df %>% distinct(!!!syms(group_var.obj), .keep_all = TRUE)

  # Apply make.names to character/factor columns
  unique_values_df[] <- lapply(unique_values_df, function(col) {
    if (is.character(col) || is.factor(col)) {
      return(make.names(col))
    } else {
      return(col)
    }
  })

  # Add speckle.condition as a factor column for model.matrix
  unique_values_df$speckle.condition <- factor(make.names(unique_values_df[[condition.obj]]))

  # Create design matrix - include all covariates
  covariates_in_formula <- c()
  if (!is.null(batch.obj) && batch.obj != "NULL") {
    covariates_in_formula <- c(covariates_in_formula, batch.obj)
  }
  if (length(covariate.obj_list) > 0) {
    covariates_in_formula <- c(covariates_in_formula, covariate.obj_list)
  }
  covariates_in_formula <- unique(covariates_in_formula)

  if (length(covariates_in_formula) == 0) {
    formula_str <- "~ 0 + speckle.condition"
  } else {
    formula_str <- paste("~ 0 + speckle.condition +", paste(covariates_in_formula, collapse = " + "))
  }

  design <- model.matrix(as.formula(formula_str), data = unique_values_df)

  # Check for rank deficiency in speckle design matrix and auto-drop
  speckle_rank <- qr(design)$rank
  if (speckle_rank < ncol(design)) {
    message("  WARNING: Speckle design matrix rank-deficient (rank=", speckle_rank,
            ", ncol=", ncol(design), "). Auto-dropping ", ncol(design) - speckle_rank, " redundant columns.")
    pivot_cols <- qr(design)$pivot[1:speckle_rank]
    dropped_cols <- colnames(design)[setdiff(seq_len(ncol(design)), pivot_cols)]
    message("  Dropped columns: ", paste(dropped_cols, collapse = ", "))
    design <- design[, pivot_cols, drop = FALSE]
  }

  # Debug: print design matrix column names
  message("  Design matrix columns: ", paste(colnames(design), collapse = ", "))

  # Run propeller for each comparison
  for (i in 1:dim(group_combinations)[2]) {
    group1 <- group_combinations[1, i]
    group2 <- group_combinations[2, i]

    message("Running Speckle for: ", group2, " vs ", group1)

    # Construct contrast string using exact column names from design matrix
    col1 <- paste0("speckle.condition", make.names(group1))
    col2 <- paste0("speckle.condition", make.names(group2))
    contrast_str <- paste0(col1, "-", col2)  # group1 vs group2 (group2 is reference)
    message("  Contrast: ", contrast_str)

    # Wrap in tryCatch to handle empty data or missing cell types
    tryCatch({
      mycontr <- makeContrasts(contrasts = contrast_str, levels = design)
      results <- propeller.ttest(props, design, contrasts = mycontr, robust = TRUE, trend = FALSE, sort = TRUE)
      write.table(results, file = glue("speckle_diffprop/comb_condition_{group1}_vs_{group2}.tsv"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    }, error = function(e) {
      message("  WARNING: Speckle failed for ", group2, " vs ", group1, ": ", e$message)
      message("  Skipping this comparison (likely no cells or empty proportion matrix)")
    })
  }
  message("Saved cell composition results")

  # ==========================================================================
  # Generate cell proportion barplot with significance annotations
  # ==========================================================================
  message("Generating cell proportion barplot...")

  # Get unique conditions from group_combinations
  all_conditions <- unique(as.vector(group_combinations))

  # Read all speckle result files and extract proportions
  speckle_files <- list.files("speckle_diffprop", pattern = "^comb_condition_.*\\.tsv$", full.names = TRUE)

  if (length(speckle_files) > 0) {
    # Initialize storage for proportions and significance
    prop_list <- list()
    sig_list <- list()

    for (f in speckle_files) {
      df <- read.table(f, sep = "\t", header = TRUE, row.names = 1)

      # Extract comparison name from filename
      fname <- basename(f)
      # Pattern: comb_condition_Group1_vs_Group2.tsv
      match <- regmatches(fname, regexec("comb_condition_(.+)_vs_(.+)\\.tsv", fname))[[1]]
      if (length(match) == 3) {
        group1 <- match[2]
        group2 <- match[3]

        # Extract proportions
        prop_cols <- grep("^PropMean", colnames(df), value = TRUE)
        for (pc in prop_cols) {
          # Extract condition name from column name (PropMean.speckle.conditionXXX)
          cond_match <- regmatches(pc, regexec("PropMean\\.speckle\\.condition(.+)", pc))[[1]]
          if (length(cond_match) == 2) {
            cond_name <- cond_match[2]
            for (ct in rownames(df)) {
              key <- paste(ct, cond_name, sep = "||")
              prop_list[[key]] <- df[ct, pc]
            }
          }
        }

        # Extract significance (FDR)
        for (ct in rownames(df)) {
          sig_key <- paste(ct, group1, group2, sep = "||")
          sig_list[[sig_key]] <- df[ct, "FDR"]
        }
      }
    }

    # Build proportion dataframe
    prop_entries <- names(prop_list)
    prop_data <- do.call(rbind, lapply(prop_entries, function(x) {
      parts <- strsplit(x, "\\|\\|")[[1]]
      data.frame(CellType = parts[1], Condition = parts[2], Proportion = prop_list[[x]], stringsAsFactors = FALSE)
    }))

    # Remove duplicates (same cell type + condition from different comparison files)
    prop_data <- prop_data %>% distinct(CellType, Condition, .keep_all = TRUE)

    # Get cell types
    cell_types <- unique(prop_data$CellType)
    conditions <- unique(prop_data$Condition)

    # Build significance dataframe
    sig_entries <- names(sig_list)
    sig_data <- do.call(rbind, lapply(sig_entries, function(x) {
      parts <- strsplit(x, "\\|\\|")[[1]]
      data.frame(CellType = parts[1], Group1 = parts[2], Group2 = parts[3], FDR = sig_list[[x]], stringsAsFactors = FALSE)
    }))

    # Function to convert FDR to asterisks
    get_asterisks <- function(fdr) {
      if (is.na(fdr)) return("")
      if (fdr < 0.001) return("***")
      if (fdr < 0.01) return("**")
      if (fdr < 0.05) return("*")
      return("")
    }

    sig_data$Asterisks <- sapply(sig_data$FDR, get_asterisks)
    sig_data$Significant <- sig_data$FDR < 0.05

    # Calculate total proportion per cell type for ordering
    total_prop <- prop_data %>%
      group_by(CellType) %>%
      summarize(Total = sum(Proportion, na.rm = TRUE)) %>%
      arrange(desc(Total))

    # Set factor levels for ordering
    prop_data$CellType <- factor(prop_data$CellType, levels = total_prop$CellType)
    prop_data$Condition <- factor(prop_data$Condition, levels = conditions)

    # Get max proportion per cell type for bracket positioning
    max_props <- prop_data %>%
      group_by(CellType) %>%
      summarize(MaxProp = max(Proportion, na.rm = TRUE))

    # Filter significant comparisons
    sig_filtered <- sig_data %>% filter(Significant)

    # Define colors (up to 5 conditions)
    color_palette <- c("#2ecc71", "#3498db", "#e74c3c", "#9b59b6", "#f39c12")
    names(color_palette) <- conditions[1:min(length(conditions), 5)]

    # Create base plot
    p <- ggplot(prop_data, aes(x = CellType, y = Proportion, fill = Condition)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = color_palette) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Cell Type",
        y = "Mean Proportion",
        title = "Cell Type Proportions Across Conditions",
        subtitle = "* FDR<0.05, ** FDR<0.01, *** FDR<0.001",
        fill = "Condition"
      )

    # Add significance brackets if there are significant comparisons
    if (nrow(sig_filtered) > 0) {
      cell_type_levels <- levels(prop_data$CellType)
      n_conditions <- length(conditions)

      # Prepare annotation data
      annot_list <- list()
      for (ct in cell_type_levels) {
        ct_sig <- sig_filtered %>% filter(CellType == ct)
        if (nrow(ct_sig) == 0) next

        ct_max <- max_props$MaxProp[max_props$CellType == ct]
        ct_idx <- which(cell_type_levels == ct)

        # Bar positions within group
        bar_offset <- 0.8 / n_conditions
        condition_positions <- setNames(
          ct_idx + ((seq_along(conditions) - 1) - (n_conditions - 1) / 2) * bar_offset,
          conditions
        )

        # Sort by span width (narrower first)
        ct_sig$Span <- sapply(1:nrow(ct_sig), function(i) {
          g1_idx <- which(conditions == ct_sig$Group1[i])
          g2_idx <- which(conditions == ct_sig$Group2[i])
          abs(g1_idx - g2_idx)
        })
        ct_sig <- ct_sig %>% arrange(Span)

        base_y <- ct_max * 1.10
        y_step <- ct_max * 0.55

        for (i in 1:nrow(ct_sig)) {
          g1 <- ct_sig$Group1[i]
          g2 <- ct_sig$Group2[i]
          ast <- ct_sig$Asterisks[i]

          x1 <- condition_positions[g1]
          x2 <- condition_positions[g2]
          if (is.na(x1) || is.na(x2)) next

          y_pos <- base_y + (i - 1) * y_step

          annot_list[[length(annot_list) + 1]] <- data.frame(
            xmin = min(x1, x2),
            xmax = max(x1, x2),
            y = y_pos,
            label = ast,
            CellType = ct
          )
        }
      }

      if (length(annot_list) > 0) {
        annot_df <- do.call(rbind, annot_list)
        bracket_height <- max(prop_data$Proportion, na.rm = TRUE) * 0.015

        p <- p +
          # Horizontal line
          geom_segment(data = annot_df,
                       aes(x = xmin, xend = xmax, y = y + bracket_height, yend = y + bracket_height),
                       inherit.aes = FALSE, linewidth = 0.4) +
          # Left vertical
          geom_segment(data = annot_df,
                       aes(x = xmin, xend = xmin, y = y, yend = y + bracket_height),
                       inherit.aes = FALSE, linewidth = 0.4) +
          # Right vertical
          geom_segment(data = annot_df,
                       aes(x = xmax, xend = xmax, y = y, yend = y + bracket_height),
                       inherit.aes = FALSE, linewidth = 0.4) +
          # Asterisks (proximal to bracket)
          geom_text(data = annot_df,
                    aes(x = (xmin + xmax) / 2, y = y + bracket_height * 1.3, label = label),
                    inherit.aes = FALSE, size = 3, fontface = "bold", vjust = 0)

        # Extend y-axis to fit annotations
        y_max <- max(c(prop_data$Proportion, annot_df$y + max(prop_data$Proportion, na.rm = TRUE) * 0.05), na.rm = TRUE)
        p <- p + coord_cartesian(ylim = c(0, y_max * 1.1))
      }
    }

    # Save plot
    ggsave("speckle_diffprop/speckle_barplot_proportions.png", p, width = 16, height = 9, dpi = 150, bg = "white")
    message("Saved: speckle_diffprop/speckle_barplot_proportions.png")
  } else {
    message("Warning: No speckle result files found for plotting")
  }

} else {
  message("\n", strrep("-", 60))
  message("STEP 7: SKIPPED (speckle module not selected or fgsea_only mode)")
  message(strrep("-", 60))
}

# ==============================================================================
# Pipeline Complete
# ==============================================================================

message("\n", strrep("=", 80))
message("PIPELINE COMPLETE")
message(strrep("=", 80))
message("\nModules executed: ", paste(selected_modules, collapse = ", "))
message("\nOutput files generated in: ", output_dir)

if ("deg" %in% selected_modules) {
  message("  - pseudobulk_*_DE.tsv (differential expression)")
  message("  - deg_completion_report.tsv (per-cluster DEG status)")
  message("  - deg_filtering_stats.tsv (filterByExpr breakdown)")
  message("  - deg_summary.txt (human-readable DEG summary)")
  message("  - volcano_plots/*.png (static volcano plots)")
  message("  - volcano_plots/*.html (interactive volcano plots with hover)")
}

if ("fgsea_groups" %in% selected_modules) {
  message("  - fgsea_groups_clusters/*.tsv (comparison-level pathway enrichment)")
  if (!is.null(custom_signatures)) {
    message("  - fgsea_custom_signatures.tsv (custom signature enrichment)")
  }
}

if ("fgsea_clusters" %in% selected_modules) {
  message("  - fgsea_clusters/*.tsv + *.pdf (cluster-level pathway enrichment)")
  if (!is.null(custom_signatures)) {
    message("  - fgsea_clusters/clusters_custom_clusters.tsv (cluster-level custom)")
  }
}

if ("markers" %in% selected_modules) {
  message("  - combined_markers_clusters.tsv (cluster markers)")
}

if ("speckle" %in% selected_modules) {
  message("  - speckle_diffprop/*.tsv (cell composition)")
  message("  - speckle_diffprop/speckle_barplot_proportions.png (cell proportion barplot)")
}

message("\n")
