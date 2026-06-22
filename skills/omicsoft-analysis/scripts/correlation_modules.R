# =============================================================================
# Correlation Analysis Modules
# =============================================================================
# Module 1: Gene-Gene Correlation (pooled across all samples, no group faceting)
# Module 1b: GSVA-GSVA Correlation (pairwise between signatures)
# Module 2: Target Gene vs Continuous Metadata Scatter
# Module 3: GSVA Signature vs Continuous Metadata Scatter
# =============================================================================

suppressPackageStartupMessages({
  if (!require(ggplot2, quietly = TRUE)) stop("ggplot2 required")
  if (!require(dplyr, quietly = TRUE)) stop("dplyr required")
  if (!require(plotly, quietly = TRUE)) stop("plotly required")
  if (!require(htmlwidgets, quietly = TRUE)) stop("htmlwidgets required")
})

format_cor_annotation <- function(r_val, p_val) {
  if (is.na(r_val) || is.na(p_val)) return("")
  asterisk <- if (p_val < 0.0001) "****"
              else if (p_val < 0.001) "***"
              else if (p_val < 0.01) "**"
              else if (p_val < 0.05) "*"
              else ""
  paste0(sprintf("%.2f", r_val), asterisk)
}

compute_pairwise_cor <- function(mat, method = "spearman") {
  genes <- rownames(mat)
  n <- length(genes)
  if (n < 2) return(NULL)

  r_mat <- matrix(NA, n, n, dimnames = list(genes, genes))
  p_mat <- matrix(NA, n, n, dimnames = list(genes, genes))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        r_mat[i, j] <- 1
        p_mat[i, j] <- 0
      } else if (j > i) {
        x <- mat[i, ]
        y <- mat[j, ]
        valid <- !is.na(x) & !is.na(y)
        if (sum(valid) >= 5) {
          if (sd(x[valid]) == 0 || sd(y[valid]) == 0) {
            r_mat[i, j] <- NA; p_mat[i, j] <- NA
            r_mat[j, i] <- NA; p_mat[j, i] <- NA
            next
          }
          ct <- tryCatch(
            cor.test(x[valid], y[valid], method = method),
            error = function(e) NULL
          )
          if (!is.null(ct)) {
            r_mat[i, j] <- ct$estimate; p_mat[i, j] <- ct$p.value
            r_mat[j, i] <- ct$estimate; p_mat[j, i] <- ct$p.value
          }
        }
      }
    }
  }
  list(r = r_mat, p = p_mat)
}

save_cor_plot <- function(p, output_path_base, width = 12, height = 8) {
  png_path <- paste0(output_path_base, ".png")
  html_path <- paste0(output_path_base, ".html")

  ggsave(png_path, p, width = width, height = height, dpi = 300)
  message(paste("  Saved:", png_path))

  p_int <- tryCatch({
    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(hoverlabel = list(bgcolor = "white"))
  }, error = function(e) {
    plotly::ggplotly(p)
  })
  htmlwidgets::saveWidget(p_int, html_path, selfcontained = TRUE)
  message(paste("  Saved:", html_path))
}

# =============================================================================
# MODULE 1: Gene-Gene Correlation (pooled across all samples)
# Direct port of target_query_plot.R lines 1548-1637
# =============================================================================
run_gene_gene_correlation <- function(expr_mat, metadata, target_genes, output_dir,
                                       study_name = "", cor_method = "spearman") {
  message(paste("\n--- Module 1: Gene-Gene Correlation ---", study_name))

  if (all(target_genes %in% rownames(expr_mat))) {
    gene_mat <- expr_mat[intersect(target_genes, rownames(expr_mat)), , drop = FALSE]
  } else if (all(target_genes %in% colnames(expr_mat))) {
    gene_mat <- t(expr_mat[, intersect(target_genes, colnames(expr_mat)), drop = FALSE])
  } else {
    found <- target_genes[target_genes %in% rownames(expr_mat) | target_genes %in% colnames(expr_mat)]
    if (length(found) < 2) { message("  WARNING: < 2 target genes found. Skipping."); return(invisible(NULL)) }
    if (sum(target_genes %in% rownames(expr_mat)) >= sum(target_genes %in% colnames(expr_mat))) {
      gene_mat <- expr_mat[intersect(target_genes, rownames(expr_mat)), , drop = FALSE]
    } else {
      gene_mat <- t(expr_mat[, intersect(target_genes, colnames(expr_mat)), drop = FALSE])
    }
  }

  found_genes <- rownames(gene_mat)
  if (length(found_genes) < 2) { message("  WARNING: < 2 target genes. Skipping."); return(invisible(NULL)) }
  message(paste("  Genes:", length(found_genes), "of", length(target_genes), "found"))

  common_samples <- intersect(colnames(gene_mat), rownames(metadata))
  if (length(common_samples) < 5) { message("  WARNING: < 5 common samples. Skipping."); return(invisible(NULL)) }
  gene_mat <- gene_mat[, common_samples, drop = FALSE]

  cor_result <- compute_pairwise_cor(gene_mat, method = cor_method)
  if (is.null(cor_result)) { message("  WARNING: Correlation computation failed. Skipping."); return(invisible(NULL)) }

  genes <- rownames(cor_result$r)
  stats_rows <- list()
  for (i in seq_along(genes)) {
    for (j in seq_along(genes)) {
      if (j <= i) next
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        gene1 = genes[i], gene2 = genes[j],
        spearman_r = cor_result$r[i, j], pvalue = cor_result$p[i, j],
        n = length(common_samples), stringsAsFactors = FALSE
      )
    }
  }
  if (length(stats_rows) == 0) { message("  WARNING: No correlation results. Skipping."); return(invisible(NULL)) }
  stats_df <- do.call(rbind, stats_rows)
  stats_df$annotation <- mapply(format_cor_annotation, stats_df$spearman_r, stats_df$pvalue)

  # Build full symmetric matrix for heatmap
  stats_mirror <- stats_df
  stats_mirror$gene1_tmp <- stats_mirror$gene2
  stats_mirror$gene2 <- stats_mirror$gene1
  stats_mirror$gene1 <- stats_mirror$gene1_tmp
  stats_mirror$gene1_tmp <- NULL
  diag_df <- data.frame(
    gene1 = found_genes, gene2 = found_genes,
    spearman_r = 1, pvalue = 0, n = NA_integer_, annotation = "1.00",
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(stats_df, stats_mirror, diag_df)
  plot_df$gene1 <- factor(plot_df$gene1, levels = found_genes)
  plot_df$gene2 <- factor(plot_df$gene2, levels = rev(found_genes))
  plot_df$hover_text <- paste0("Genes: ", plot_df$gene1, " vs ", plot_df$gene2,
    "\nrho = ", sprintf("%.3f", plot_df$spearman_r),
    "\np = ", format(plot_df$pvalue, digits = 3, scientific = TRUE))

  # Single heatmap (no faceting)
  p_heat <- ggplot(plot_df, aes(x = gene1, y = gene2, fill = spearman_r, text = hover_text)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    geom_text(aes(label = annotation), size = 2.8, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
    labs(title = paste0("Gene-Gene Correlation — ", study_name),
         subtitle = paste0("Spearman r, pooled (n=", length(common_samples),
                           "; * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001)"),
         x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9),
          legend.position = "right")

  cor_dir <- file.path(output_dir, "correlations", "gene_gene")
  dir.create(cor_dir, recursive = TRUE, showWarnings = FALSE)
  n_genes <- length(found_genes)
  plot_width <- max(8, n_genes * 1.8)
  plot_height <- max(6, n_genes * 1.8)
  save_cor_plot(p_heat, file.path(cor_dir, paste0(study_name, "_gene_gene_heatmap")),
                width = plot_width, height = plot_height)

  csv_path <- file.path(cor_dir, paste0(study_name, "_gene_gene_stats.csv"))
  write.csv(stats_df[, c("gene1", "gene2", "spearman_r", "pvalue", "n")], csv_path, row.names = FALSE)
  message(paste("  Saved:", csv_path, "(", nrow(stats_df), "rows)"))
  invisible(list(heatmap = p_heat, stats_df = stats_df))
}

# =============================================================================
# MODULE 1b: GSVA-GSVA Correlation (pooled across all samples)
# Direct port of target_query_plot.R lines 1639-1754
# =============================================================================
run_gsva_gsva_correlation <- function(gsva_scores, output_dir, study_name = "",
                                       cor_method = "spearman") {
  if (is.null(gsva_scores) || nrow(gsva_scores) < 2) {
    message("  GSVA-GSVA correlation: < 2 signatures, skipping.")
    return(invisible(NULL))
  }
  message(paste("\n--- Module 1b: GSVA-GSVA Correlation ---", study_name))
  n_sigs <- nrow(gsva_scores)
  n_samples <- ncol(gsva_scores)
  if (n_samples < 5) { message("  WARNING: < 5 samples. Skipping."); return(invisible(NULL)) }
  message(paste("  Signatures:", n_sigs, "| Samples:", n_samples))

  cor_result <- compute_pairwise_cor(gsva_scores, method = cor_method)
  if (is.null(cor_result)) { message("  WARNING: Correlation computation failed. Skipping."); return(invisible(NULL)) }

  sigs <- rownames(cor_result$r)
  stats_rows <- list()
  for (i in seq_along(sigs)) {
    for (j in seq_along(sigs)) {
      if (j <= i) next
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        sig1 = sigs[i], sig2 = sigs[j],
        spearman_r = cor_result$r[i, j], pvalue = cor_result$p[i, j],
        n = n_samples, stringsAsFactors = FALSE
      )
    }
  }
  if (length(stats_rows) == 0) { message("  WARNING: No results. Skipping."); return(invisible(NULL)) }
  stats_df <- do.call(rbind, stats_rows)
  stats_df$annotation <- mapply(format_cor_annotation, stats_df$spearman_r, stats_df$pvalue)

  stats_mirror <- stats_df
  stats_mirror$sig1_tmp <- stats_mirror$sig2
  stats_mirror$sig2 <- stats_mirror$sig1
  stats_mirror$sig1 <- stats_mirror$sig1_tmp
  stats_mirror$sig1_tmp <- NULL
  diag_df <- data.frame(
    sig1 = sigs, sig2 = sigs,
    spearman_r = 1, pvalue = 0, n = NA_integer_, annotation = "1.00",
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(stats_df, stats_mirror, diag_df)
  plot_df$sig1 <- factor(plot_df$sig1, levels = sigs)
  plot_df$sig2 <- factor(plot_df$sig2, levels = rev(sigs))
  plot_df$hover_text <- paste0("Signatures: ", plot_df$sig1, " vs ", plot_df$sig2,
    "\nrho = ", sprintf("%.3f", plot_df$spearman_r),
    "\np = ", format(plot_df$pvalue, digits = 3, scientific = TRUE))

  p_heat <- ggplot(plot_df, aes(x = sig1, y = sig2, fill = spearman_r, text = hover_text)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    geom_text(aes(label = annotation), size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
    labs(title = paste0("GSVA Signature Correlation — ", study_name),
         subtitle = paste0("Spearman r, pooled (n=", n_samples, ")"), x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")

  cor_dir <- file.path(output_dir, "correlations", "gsva_gsva")
  dir.create(cor_dir, recursive = TRUE, showWarnings = FALSE)
  plot_width <- max(8, n_sigs * 2)
  plot_height <- max(6, n_sigs * 2)
  save_cor_plot(p_heat, file.path(cor_dir, paste0(study_name, "_gsva_gsva_heatmap")),
                width = plot_width, height = plot_height)

  # Pairwise scatter plots
  scatter_data <- list()
  for (idx in seq_len(nrow(stats_df))) {
    s1 <- stats_df$sig1[idx]
    s2 <- stats_df$sig2[idx]
    pair_df <- data.frame(
      x = gsva_scores[s1, ], y = gsva_scores[s2, ],
      pair = paste0(s1, " vs ", s2), stringsAsFactors = FALSE
    )
    pair_df$hover_text <- paste0("x (", s1, "): ", round(pair_df$x, 3),
                                 "\ny (", s2, "): ", round(pair_df$y, 3))
    pair_df$annotation <- paste0("rho = ", sprintf("%.2f", stats_df$spearman_r[idx]),
                                 ", p = ", format(stats_df$pvalue[idx], digits = 2, scientific = TRUE))
    scatter_data[[idx]] <- pair_df
  }
  scatter_df <- do.call(rbind, scatter_data)
  scatter_df$pair <- factor(scatter_df$pair, levels = unique(scatter_df$pair))
  n_pairs <- nrow(stats_df)
  n_cols_scatter <- min(n_pairs, 3)
  n_rows_scatter <- ceiling(n_pairs / n_cols_scatter)
  annot_df <- scatter_df[!duplicated(scatter_df$pair), c("pair", "annotation")]

  p_scatter <- ggplot(scatter_df, aes(x = x, y = y, text = hover_text)) +
    geom_point(alpha = 0.6, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, level = 0.95,
                color = "darkred", linewidth = 0.8, fill = "pink", alpha = 0.2) +
    geom_text(data = annot_df, aes(x = -Inf, y = Inf, label = annotation),
              hjust = -0.05, vjust = 1.3, size = 3, color = "black",
              inherit.aes = FALSE) +
    facet_wrap(~ pair, scales = "free", ncol = n_cols_scatter) +
    labs(title = paste0("GSVA Signature Scatter — ", study_name),
         subtitle = paste0("Spearman r with linear fit + 95% CI (n=", n_samples, ")"),
         x = "GSVA Score (Signature 1)", y = "GSVA Score (Signature 2)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(size = 9, face = "bold"))

  scatter_width <- max(10, n_cols_scatter * 4)
  scatter_height <- max(6, n_rows_scatter * 3.5)
  save_cor_plot(p_scatter, file.path(cor_dir, paste0(study_name, "_gsva_gsva_scatter")),
                width = scatter_width, height = scatter_height)

  csv_path <- file.path(cor_dir, paste0(study_name, "_gsva_gsva_stats.csv"))
  write.csv(stats_df[, c("sig1", "sig2", "spearman_r", "pvalue", "n")], csv_path, row.names = FALSE)
  message(paste("  Saved:", csv_path, "(", nrow(stats_df), "rows)"))
  invisible(list(heatmap = p_heat, scatter = p_scatter, stats_df = stats_df))
}

# =============================================================================
# MODULE 2: Target Gene vs Continuous Metadata Scatter
# Direct port of target_query_plot.R lines 1757-1821
# =============================================================================
run_target_vs_continuous <- function(expr_mat, metadata, target_genes,
                                      continuous_cols, output_dir, study_name = "") {
  if (is.null(continuous_cols) || length(continuous_cols) == 0) return(invisible(NULL))
  message(paste("\n--- Module 2: Target vs Continuous ---", study_name))

  if (sum(target_genes %in% rownames(expr_mat)) >= sum(target_genes %in% colnames(expr_mat))) {
    gene_mat <- expr_mat[intersect(target_genes, rownames(expr_mat)), , drop = FALSE]
  } else {
    gene_mat <- t(expr_mat[, intersect(target_genes, colnames(expr_mat)), drop = FALSE])
  }

  found_genes <- rownames(gene_mat)
  if (length(found_genes) == 0) { message("  WARNING: No target genes found. Skipping."); return(invisible(NULL)) }

  common_samples <- intersect(colnames(gene_mat), rownames(metadata))
  if (length(common_samples) < 5) { message("  WARNING: < 5 common samples. Skipping."); return(invisible(NULL)) }
  gene_mat <- gene_mat[, common_samples, drop = FALSE]
  meta_sub <- metadata[common_samples, , drop = FALSE]

  valid_cols <- continuous_cols[continuous_cols %in% colnames(meta_sub)]
  if (length(valid_cols) == 0) { message("  WARNING: No valid continuous columns. Skipping."); return(invisible(NULL)) }

  cor_dir <- file.path(output_dir, "correlations", "target_vs_clinical")
  dir.create(cor_dir, recursive = TRUE, showWarnings = FALSE)

  all_stats <- list()

  for (var_name in valid_cols) {
    var_vals <- suppressWarnings(as.numeric(meta_sub[[var_name]]))
    names(var_vals) <- rownames(meta_sub)

    non_na_vals <- var_vals[!is.na(var_vals)]
    if (length(non_na_vals) < 3 || sd(non_na_vals) == 0) {
      message(paste("  WARNING:", var_name, "has zero variance or <3 non-NA values. Skipping."))
      next
    }

    scatter_rows <- list()
    annot_rows <- list()

    for (gene in found_genes) {
      expr_vals <- gene_mat[gene, ]
      valid <- !is.na(expr_vals) & !is.na(var_vals)
      n_valid <- sum(valid)

      if (n_valid >= 5) {
        spearman <- cor.test(expr_vals[valid], var_vals[valid], method = "spearman")
        pearson <- cor.test(expr_vals[valid], var_vals[valid], method = "pearson")
        all_stats[[length(all_stats) + 1]] <- data.frame(
          gene = gene, continuous_var = var_name,
          spearman_rho = as.numeric(spearman$estimate), spearman_pval = spearman$p.value,
          pearson_r = as.numeric(pearson$estimate), pearson_pval = pearson$p.value,
          n = n_valid, stringsAsFactors = FALSE
        )
        for (idx in which(valid)) {
          scatter_rows[[length(scatter_rows) + 1]] <- data.frame(
            Gene = gene, x_val = var_vals[idx], y_val = expr_vals[idx], stringsAsFactors = FALSE
          )
        }
        annot_rows[[length(annot_rows) + 1]] <- data.frame(
          Gene = gene,
          label = paste0("rho = ", sprintf("%.2f", spearman$estimate), ", p = ",
                         format(spearman$p.value, digits = 2, scientific = TRUE), "\nn = ", n_valid),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(scatter_rows) == 0) next
    scatter_df <- do.call(rbind, scatter_rows)
    scatter_df$Gene <- factor(scatter_df$Gene, levels = found_genes)
    annot_df <- do.call(rbind, annot_rows)
    annot_df$Gene <- factor(annot_df$Gene, levels = found_genes)
    annot_pos <- scatter_df %>% group_by(Gene) %>%
      summarize(x_pos = min(x_val, na.rm = TRUE) + diff(range(x_val, na.rm = TRUE)) * 0.05,
                y_pos = max(y_val, na.rm = TRUE) - diff(range(y_val, na.rm = TRUE)) * 0.05, .groups = "drop")
    annot_df <- annot_df %>% left_join(annot_pos, by = "Gene")
    scatter_df$hover_text <- paste0("Gene: ", scatter_df$Gene, "\n", var_name, ": ",
      sprintf("%.2f", scatter_df$x_val), "\nExpression: ", sprintf("%.3f", scatter_df$y_val))
    p <- ggplot(scatter_df, aes(x = x_val, y = y_val, text = hover_text)) +
      geom_point(alpha = 0.6, size = 2, color = "#2166AC") +
      geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "#B2182B", fill = "#FDDBC7", linewidth = 0.8) +
      geom_text(data = annot_df, aes(x = x_pos, y = y_pos, label = label), hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
      facet_wrap(~ Gene, scales = "free_y") +
      labs(title = paste0("Target Expression vs ", var_name, " — ", study_name),
           subtitle = "Spearman rho shown; shaded = 95% CI (linear fit)", x = var_name, y = "Normalized Expression") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9),
            strip.text = element_text(size = 11, face = "bold"))
    n_panels <- length(found_genes); n_cols_p <- min(n_panels, 4); n_rows_p <- ceiling(n_panels / n_cols_p)
    save_cor_plot(p, file.path(cor_dir, paste0(study_name, "_target_vs_", var_name)),
                  width = max(10, n_cols_p * 4), height = max(6, n_rows_p * 3.5))
  }
  if (length(all_stats) > 0) {
    stats_combined <- do.call(rbind, all_stats)
    csv_path <- file.path(cor_dir, paste0(study_name, "_target_vs_clinical_stats.csv"))
    write.csv(stats_combined, csv_path, row.names = FALSE)
    message(paste("  Saved:", csv_path, "(", nrow(stats_combined), "rows)"))
    invisible(list(stats_df = stats_combined))
  } else { invisible(NULL) }
}

# =============================================================================
# MODULE 3: GSVA Signature vs Continuous Metadata Scatter
# Direct port of target_query_plot.R lines 1825-1888
# =============================================================================
run_gsva_vs_continuous <- function(gsva_scores, metadata, continuous_cols,
                                    output_dir, study_name = "") {
  if (is.null(continuous_cols) || length(continuous_cols) == 0) return(invisible(NULL))
  if (is.null(gsva_scores) || nrow(gsva_scores) == 0) return(invisible(NULL))
  message(paste("\n--- Module 3: GSVA vs Continuous ---", study_name))

  common_samples <- intersect(colnames(gsva_scores), rownames(metadata))
  if (length(common_samples) < 5) { message("  WARNING: < 5 common samples. Skipping."); return(invisible(NULL)) }
  gsva_sub <- gsva_scores[, common_samples, drop = FALSE]
  meta_sub <- metadata[common_samples, , drop = FALSE]

  sig_names <- rownames(gsva_sub)
  if (length(sig_names) == 0) { message("  WARNING: No signatures. Skipping."); return(invisible(NULL)) }

  valid_cols <- continuous_cols[continuous_cols %in% colnames(meta_sub)]
  if (length(valid_cols) == 0) { message("  WARNING: No valid continuous columns. Skipping."); return(invisible(NULL)) }

  cor_dir <- file.path(output_dir, "correlations", "gsva_vs_clinical")
  dir.create(cor_dir, recursive = TRUE, showWarnings = FALSE)

  all_stats <- list()

  for (var_name in valid_cols) {
    var_vals <- suppressWarnings(as.numeric(meta_sub[[var_name]]))
    names(var_vals) <- rownames(meta_sub)

    non_na_vals <- var_vals[!is.na(var_vals)]
    if (length(non_na_vals) < 3 || sd(non_na_vals) == 0) next

    scatter_rows <- list()
    annot_rows <- list()

    for (sig in sig_names) {
      score_vals <- gsva_sub[sig, ]
      valid <- !is.na(score_vals) & !is.na(var_vals)
      n_valid <- sum(valid)
      if (n_valid >= 5) {
        spearman <- cor.test(score_vals[valid], var_vals[valid], method = "spearman")
        pearson <- cor.test(score_vals[valid], var_vals[valid], method = "pearson")
        all_stats[[length(all_stats) + 1]] <- data.frame(
          signature = sig, continuous_var = var_name,
          spearman_rho = as.numeric(spearman$estimate), spearman_pval = spearman$p.value,
          pearson_r = as.numeric(pearson$estimate), pearson_pval = pearson$p.value,
          n = n_valid, stringsAsFactors = FALSE
        )
        for (idx in which(valid)) {
          scatter_rows[[length(scatter_rows) + 1]] <- data.frame(
            Signature = sig, x_val = var_vals[idx], y_val = score_vals[idx], stringsAsFactors = FALSE
          )
        }
        annot_rows[[length(annot_rows) + 1]] <- data.frame(
          Signature = sig,
          label = paste0("rho = ", sprintf("%.2f", spearman$estimate), ", p = ",
                         format(spearman$p.value, digits = 2, scientific = TRUE), "\nn = ", n_valid),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(scatter_rows) == 0) next
    scatter_df <- do.call(rbind, scatter_rows)
    scatter_df$Signature <- factor(scatter_df$Signature, levels = sig_names)
    annot_df <- do.call(rbind, annot_rows)
    annot_df$Signature <- factor(annot_df$Signature, levels = sig_names)
    annot_pos <- scatter_df %>% group_by(Signature) %>%
      summarize(x_pos = min(x_val, na.rm = TRUE) + diff(range(x_val, na.rm = TRUE)) * 0.05,
                y_pos = max(y_val, na.rm = TRUE) - diff(range(y_val, na.rm = TRUE)) * 0.05, .groups = "drop")
    annot_df <- annot_df %>% left_join(annot_pos, by = "Signature")
    scatter_df$hover_text <- paste0("Signature: ", scatter_df$Signature, "\n", var_name, ": ",
      sprintf("%.2f", scatter_df$x_val), "\nGSVA Score: ", sprintf("%.3f", scatter_df$y_val))
    p <- ggplot(scatter_df, aes(x = x_val, y = y_val, text = hover_text)) +
      geom_point(alpha = 0.6, size = 2, color = "#2166AC") +
      geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "#B2182B", fill = "#FDDBC7", linewidth = 0.8) +
      geom_text(data = annot_df, aes(x = x_pos, y = y_pos, label = label), hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
      facet_wrap(~ Signature, scales = "free_y") +
      labs(title = paste0("GSVA Score vs ", var_name, " — ", study_name),
           subtitle = "Spearman rho shown; shaded = 95% CI (linear fit)", x = var_name, y = "GSVA Score") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9),
            strip.text = element_text(size = 11, face = "bold"))
    n_panels <- length(sig_names); n_cols_p <- min(n_panels, 3); n_rows_p <- ceiling(n_panels / n_cols_p)
    save_cor_plot(p, file.path(cor_dir, paste0(study_name, "_gsva_vs_", var_name)),
                  width = max(10, n_cols_p * 4.5), height = max(5, n_rows_p * 3.5))
  }
  if (length(all_stats) > 0) {
    stats_combined <- do.call(rbind, all_stats)
    csv_path <- file.path(cor_dir, paste0(study_name, "_gsva_vs_clinical_stats.csv"))
    write.csv(stats_combined, csv_path, row.names = FALSE)
    message(paste("  Saved:", csv_path, "(", nrow(stats_combined), "rows)"))
    invisible(list(stats_df = stats_combined))
  } else { invisible(NULL) }
}
