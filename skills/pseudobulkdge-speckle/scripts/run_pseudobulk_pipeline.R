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
              help = "|log2FC| cutoff for volcano plot significance [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) stop("--input is required")
if (is.null(opt$condition)) stop("--condition is required")
if (is.null(opt$sample)) stop("--sample is required")
if (is.null(opt$celltype)) stop("--celltype is required")
if (is.null(opt$comparisons)) stop("--comparisons is required")

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
                     scoreType = scoreType, eps = eps, nproc = 1)

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
  main.pathways$leadingEdge <- sapply(main.pathways$leadingEdge, function(x) paste(unlist(x), collapse = ","))

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
setwd(output_dir)

# Load custom signatures if provided (uses GENE SYMBOLS)
custom_signatures <- load_custom_signatures(opt$custom_signatures)

# Store column names
condition.obj <- opt$condition
celltype.obj <- opt$celltype
sample.obj <- opt$sample
batch.obj <- opt$batch
covariate.obj <- opt$covariate

# Create safe variable names for formulas
condition <- make.names(condition.obj)
celltype <- make.names(celltype.obj)
sample <- make.names(sample.obj)
batch <- if (!is.null(batch.obj) && batch.obj != "NULL") make.names(batch.obj) else NULL
covariate <- if (!is.null(covariate.obj) && covariate.obj != "NULL") make.names(covariate.obj) else NULL

# Print configuration
message("\n=== Analysis Configuration ===")
message("Input file: ", opt$input)
message("Output directory: ", output_dir)
message("Condition column: ", condition.obj, " -> ", condition)
message("Cell type column: ", celltype.obj, " -> ", celltype)
message("Sample column: ", sample.obj, " -> ", sample)
message("Batch column: ", ifelse(is.null(batch.obj) || batch.obj == "NULL", "None", batch.obj))
message("Covariate column: ", ifelse(is.null(covariate.obj) || covariate.obj == "NULL", "None", covariate.obj))
message("Number of comparisons: ", ncol(group_combinations))
for (i in 1:ncol(group_combinations)) {
  message("  Comparison ", i, ": ", group_combinations[1, i], " vs ", group_combinations[2, i])
}
message("filterByExpr: min.count=", opt$min_count, ", min.total.count=", opt$min_total_count, ", min.prop=", opt$min_prop)
message("Volcano plot cutoffs: FDR < ", opt$volcano_fdr, ", |log2FC| > ", opt$volcano_log2fc)
message("Custom signatures: ", ifelse(is.null(custom_signatures), "None", paste(names(custom_signatures), collapse = ", ")))

# ==============================================================================
# Step 1: Load Data
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 1: Loading input data...")
message(strrep("-", 60))

seurat_obj <- readRDS(opt$input)
message("Loaded Seurat object with ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " genes")

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
message("Converted to SingleCellExperiment")

# ==============================================================================
# Step 2: Pseudobulk Aggregation
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 2: Pseudobulk aggregation...")
message(strrep("-", 60))

# Build aggregation columns
agg_cols.obj <- c(condition.obj, celltype.obj, sample.obj)
if (!is.null(covariate.obj) && covariate.obj != "NULL") agg_cols.obj <- c(agg_cols.obj, covariate.obj)

# Aggregate cell counts
sce_aggregated <- aggregateAcrossCells(
  sce,
  id = colData(sce)[, agg_cols.obj],
  use.assay.type = "counts"
)
message("Aggregated to ", ncol(sce_aggregated), " pseudobulk samples")

# QC check
table_clusters <- as.data.frame.matrix(table(sce_aggregated[[condition.obj]], sce_aggregated[[celltype.obj]]))
write.table(table_clusters, "pseudobulk_aggregation_qc.tsv", quote = FALSE, sep = "\t")
message("Saved aggregation QC table")

# ==============================================================================
# Step 3: Pseudobulk Differential Expression
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 3: Pseudobulk differential expression analysis...")
message(strrep("-", 60))

# Initialize result containers
de_results_dataframe <- data.frame()
de_metrics_list <- list()

# Build design formula
if (!is.null(covariate) && covariate != "NULL") {
  formula <- as.formula(paste("~", paste(c(condition, covariate), collapse = " + ")))
} else {
  formula <- as.formula(paste("~", condition))
}

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
  current_sce[[condition]] <- relevel(factor(current_sce[[condition.obj]]), ref = group1)

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
    coef = paste0(condition, make.names(group2)),
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
de_results_dataframe$comparison <- paste0(de_results_dataframe$group2, "vs", de_results_dataframe$group1)
de_results_out <- de_results_dataframe[, c("gene_name", "logFC", "PValue", "FDR", "group1", "group2", "cluster", "comparison", "clust_res")]
write.table(de_results_out, file = paste0(result_name, "_DE.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved: ", result_name, "_DE.tsv")

# ==============================================================================
# Step 4: Pathway Enrichment Analysis
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 4: Pathway enrichment analysis...")
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

map_genes <- function(numbers_string, gene_map) {
  # Handle NA, NULL, or empty strings
  if (is.null(numbers_string) || is.na(numbers_string) || numbers_string == "" || length(numbers_string) == 0) {
    return("")
  }
  numbers <- unlist(strsplit(as.character(numbers_string), ","))
  numbers <- numbers[!is.na(numbers) & numbers != ""]
  if (length(numbers) == 0) return("")
  genes <- gene_map[numbers]
  genes <- genes[!is.na(genes)]
  paste(genes, collapse = ",")
}

# Vectorized version for better performance on large dataframes
map_genes_vectorized <- function(leadingEdge_col, gene_map) {
  message("  Converting ", length(leadingEdge_col), " leadingEdge entries to gene symbols...")
  # Pre-allocate result vector
  results <- character(length(leadingEdge_col))

  # Process in chunks for memory efficiency
  chunk_size <- 1000
  n_chunks <- ceiling(length(leadingEdge_col) / chunk_size)

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(leadingEdge_col))

    for (j in start_idx:end_idx) {
      val <- leadingEdge_col[j]
      if (is.null(val) || is.na(val) || val == "" || length(val) == 0) {
        results[j] <- ""
      } else {
        numbers <- unlist(strsplit(as.character(val), ","))
        numbers <- numbers[!is.na(numbers) & numbers != ""]
        if (length(numbers) == 0) {
          results[j] <- ""
        } else {
          genes <- gene_map[numbers]
          genes <- genes[!is.na(genes)]
          results[j] <- paste(genes, collapse = ",")
        }
      }
    }

    if (i %% 10 == 0 || i == n_chunks) {
      message("    Processed ", end_idx, "/", length(leadingEdge_col), " entries")
    }
  }

  return(results)
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
  # Convert Entrez IDs to gene symbols in leadingEdge with error handling
  message("Converting Reactome leadingEdge from Entrez IDs to gene symbols...")
  tryCatch({
    gene_map <- setNames(names(g.entrezid), g.entrezid)
    de.fgsea$leadingEdge <- map_genes_vectorized(de.fgsea$leadingEdge, gene_map)
  }, error = function(e) {
    message("WARNING: Could not convert leadingEdge to gene symbols: ", e$message)
    message("LeadingEdge will remain as Entrez IDs")
  })
  write.table(de.fgsea, file = paste0(file.name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved Reactome enrichment results: ", paste0(file.name, "_cluster.tsv"))
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
  write.table(de.fgsea.custom, file = paste0(file.name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)

  # Also save to root output directory for easy access
  write.table(de.fgsea.custom, file = "fgsea_custom_signatures.tsv", quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved custom signature enrichment results")
  message("  -> fgsea_groups_clusters/condition_clusters_custom_cluster.tsv")
  message("  -> fgsea_custom_signatures.tsv")
}

# ==============================================================================
# Step 5: Cluster-level Pathway Analysis
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 5: Cluster-level pathway analysis...")
message(strrep("-", 60))

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
  Idents(seurat_obj) <- seurat_obj@meta.data[[celltype.obj]]
  markers <- FindAllMarkers(seurat_obj, only.pos = FALSE, min.pct = 0.1,
                            logfc.threshold = 0.1, test.use = "wilcox")

  # Convert to presto-like format
  markers$feature <- markers$gene
  markers$group <- markers$cluster
  markers$auc <- (markers$pct.1 - markers$pct.2 + 1) / 2  # Approximate AUC
  markers$pval <- markers$p_val
  markers$padj <- markers$p_val_adj
  markers$logFC <- markers$avg_log2FC

  markers[, c("feature", "group", "auc", "pval", "padj", "logFC")]
})

if (is.null(x.genes) || nrow(x.genes) == 0) {
  message("WARNING: No cluster markers found. Skipping Step 5.")
} else {
  x.genes.msigdb <- x.genes  # Keep gene symbols for MSigDB

# Map to Entrez IDs for Reactome (only)
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
# Convert Entrez IDs to gene symbols in leadingEdge with error handling
message("Converting cluster-level Reactome leadingEdge to gene symbols...")
tryCatch({
  gene_map <- setNames(names(g.entrezid), g.entrezid)
  df.vec$leadingEdge <- map_genes_vectorized(df.vec$leadingEdge, gene_map)
}, error = function(e) {
  message("WARNING: Could not convert leadingEdge to gene symbols: ", e$message)
})
write.table(df.vec, file = paste0(sample_name, "_cluster.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
message("Saved cluster-level Reactome results")

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
  write.table(df.vec.custom, file = paste0(sample_name, "_clusters.tsv"), quote = FALSE, sep = '\t', row.names = FALSE)
  message("Saved cluster-level custom signature results")
}

}  # End of if (x.genes) check for Step 5

# ==============================================================================
# Step 6: Cluster Markers
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 6: Cluster marker identification...")
message(strrep("-", 60))

if (!is.null(batch.obj) && batch.obj != "NULL") {
  med_markers <- scran::findMarkers(sce, sce[[celltype.obj]], block = sce[[batch.obj]])
} else {
  med_markers <- scran::findMarkers(sce, sce[[celltype.obj]])
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

write.table(combined_med_markers,
            file = "combined_markers_clusters.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved cluster markers")

# ==============================================================================
# Step 7: Cell Composition Analysis (Speckle)
# ==============================================================================

message("\n", strrep("-", 60))
message("STEP 7: Cell composition analysis (Speckle)...")
message(strrep("-", 60))

if (!dir.exists("speckle_diffprop")) dir.create("speckle_diffprop")

# Determine grouping variable
if (!is.null(batch.obj) && batch.obj != "NULL") {
  group_var.obj <- batch.obj
} else {
  group_var.obj <- sample.obj
}

# Get transformed proportions
props <- getTransformedProps(sce[[celltype.obj]], sce[[group_var.obj]], transform = "logit")

# Create metadata dataframe
metadata_df <- as.data.frame(colData(sce))

meta_cols <- c(condition.obj, group_var.obj)
if (!is.null(covariate.obj) && covariate.obj != "NULL") meta_cols <- c(meta_cols, covariate.obj)

subset_df <- metadata_df[, meta_cols, drop = FALSE]
unique_values_df <- subset_df %>%
  group_by(!!!syms(group_var.obj)) %>%
  summarize_all(list(~ unique(.)[1]))
unique_values_df <- unique_values_df %>% distinct(!!!syms(group_var.obj), .keep_all = TRUE)

# Apply make.names
unique_values_df[] <- lapply(unique_values_df, function(col) {
  if (is.character(col) || is.factor(col)) {
    return(make.names(col))
  } else {
    return(col)
  }
})

speckle.condition <- as.character(unique_values_df[[condition.obj]])

# Create design matrix
if (is.null(covariate) || covariate == "NULL") {
  formula_str <- "~ 0 + speckle.condition"
} else {
  formula_str <- paste("~ 0 + speckle.condition +", paste(covariate, collapse = " + "))
}

design <- model.matrix(as.formula(formula_str), data = unique_values_df)

# Run propeller for each comparison
for (i in 1:dim(group_combinations)[2]) {
  group1 <- group_combinations[1, i]
  group2 <- group_combinations[2, i]

  message("Running Speckle for: ", group2, " vs ", group1)

  mycontr <- makeContrasts((glue("speckle.condition{make.names(group2)}-speckle.condition{make.names(group1)}")), levels = design)
  results <- propeller.ttest(props, design, contrasts = mycontr, robust = TRUE, trend = FALSE, sort = TRUE)

  write.table(results, file = glue("speckle_diffprop/comb_condition_{group2}_vs_{group1}.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)
}
message("Saved cell composition results")

# ==============================================================================
# Pipeline Complete
# ==============================================================================

message("\n", strrep("=", 80))
message("PIPELINE COMPLETE")
message(strrep("=", 80))
message("\nOutput files generated in: ", output_dir)
message("  - pseudobulk_*_DE.tsv (differential expression)")
message("  - fgsea_groups_clusters/*.tsv (pathway enrichment)")
message("  - fgsea_clusters/*.tsv + *.pdf (cluster-level enrichment)")
message("  - combined_markers_clusters.tsv (cluster markers)")
message("  - speckle_diffprop/*.tsv (cell composition)")
message("  - volcano_plots/*.png (static volcano plots)")
message("  - volcano_plots/*.html (interactive volcano plots with hover)")
if (!is.null(custom_signatures)) {
  message("  - fgsea_custom_signatures.tsv (custom signature enrichment)")
  message("  - fgsea_clusters/clusters_custom_clusters.tsv (cluster-level custom)")
}
message("\n")
