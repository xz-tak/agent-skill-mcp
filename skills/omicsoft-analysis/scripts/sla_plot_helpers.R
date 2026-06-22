# =============================================================================
# Sample-Level Analysis: Plot Helper Functions
# =============================================================================
# Direct port of target_query_plot.R plotting functions.
# Layout: genes on x-axis, groups dodged within each gene.
# Dodge width = 0.75 consistently across all layers.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# =============================================================================
# Utility: Asterisk significance labels
# =============================================================================
get_asterisk <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.0001) return("****")
  if (pval < 0.001) return("***")
  if (pval < 0.01) return("**")
  if (pval < 0.05) return("*")
  return("")
}

# =============================================================================
# Utility: Sanitize path component
# =============================================================================
sanitize_path <- function(x) {
  gsub("[^A-Za-z0-9_]", "_", x)
}

# =============================================================================
# Save plot as PNG + interactive HTML with DEG hover enrichment
# Direct port of target_query_plot.R lines 1261-1487
# =============================================================================
save_plot_pair <- function(p, output_dir, base_name, deg_df = NULL,
                           group_col = NULL, annotation_map = NULL,
                           control_group = NULL, sig_threshold_label = NULL,
                           gsva_stats_df = NULL) {
  ggsave(file.path(output_dir, paste0(base_name, ".png")), p, width = 12, height = 8, dpi = 300)
  message(paste("Saved:", paste0(base_name, ".png")))

  # Fallback ggplotly (ref lines 1267-1280)
  p_int <- tryCatch(
    plotly::ggplotly(p),
    error = function(e) {
      message(paste("  NOTE: plotly conversion failed:", conditionMessage(e), "- using basic ggplotly"))
      tryCatch(
        plotly::ggplotly(p + guides(colour = "none")),
        error = function(e2) {
          message(paste("  NOTE: fallback ggplotly also failed:", conditionMessage(e2), "- skipping HTML"))
          return(NULL)
        }
      )
    }
  )
  if (is.null(p_int)) return(invisible(NULL))

  # Subtitle reconstruction (ref lines 1282-1291)
  plot_subtitle <- p$labels$subtitle
  if (!is.null(plot_subtitle) && nzchar(plot_subtitle)) {
    plot_title <- p$labels$title
    combined_title <- paste0(
      "<b>", plot_title, "</b>",
      "<br><span style='font-size:11px;color:grey'>",
      gsub("\\*", "<span style='font-size:18px;font-weight:bold'>*</span>",
           gsub("\n", "<br>", plot_subtitle)), "</span>"
    )
    p_int <- p_int %>% plotly::layout(title = list(text = combined_title, x = 0.02, xanchor = "left"))
  }

  # Legend fix (ref lines 1293-1310): hide text traces, re-enable marker legendgroup
  sig_legendgroups <- c()
  for (i in seq_along(p_int$x$data)) {
    tr <- p_int$x$data[[i]]
    if (!is.null(tr$mode) && grepl("text", tr$mode)) {
      p_int$x$data[[i]]$showlegend <- FALSE
      if (!is.null(tr$legendgroup)) sig_legendgroups <- c(sig_legendgroups, tr$legendgroup)
    }
  }
  if (length(sig_legendgroups) > 0) {
    for (i in seq_along(p_int$x$data)) {
      tr <- p_int$x$data[[i]]
      if (!is.null(tr$mode) && grepl("markers", tr$mode) &&
          !is.null(tr$legendgroup) && tr$legendgroup %in% sig_legendgroups) {
        p_int$x$data[[i]]$showlegend <- TRUE
      }
    }
  }

  # DEG hover enrichment (ref lines 1311-1389)
  if (!is.null(deg_df) && nrow(deg_df) > 0 && !is.null(group_col)) {
    pval_col <- if ("pval" %in% colnames(deg_df)) "pval" else
                if ("pvalue" %in% colnames(deg_df)) "pvalue" else NULL
    has_contrast <- "comparison_contrast" %in% colnames(deg_df)

    deg_base <- deg_df %>%
      select(gene, comparison_id, log2fc, padj,
             any_of(c(pval_col, "comparison_contrast"))) %>%
      distinct()

    lookup_rows <- list()
    if (!is.null(annotation_map)) {
      for (i in seq_len(nrow(deg_base))) {
        comp_id <- deg_base$comparison_id[i]
        for (pattern in names(annotation_map)) {
          if (grepl(pattern, comp_id, fixed = TRUE)) {
            matched <- unlist(annotation_map[[pattern]])
            if (is.character(matched)) {
              for (mg in matched) {
                if (!is.null(control_group) && mg == control_group) next
                lookup_rows[[length(lookup_rows) + 1]] <- data.frame(
                  gene = deg_base$gene[i], mapped_group = mg,
                  log2fc = deg_base$log2fc[i], padj = deg_base$padj[i],
                  pval = if (!is.null(pval_col) && pval_col %in% colnames(deg_base)) deg_base[[pval_col]][i] else NA_real_,
                  contrast = if (has_contrast) deg_base$comparison_contrast[i] else NA_character_,
                  stringsAsFactors = FALSE
                )
              }
            }
            break
          }
        }
      }
    } else {
      for (i in seq_len(nrow(deg_base))) {
        lookup_rows[[length(lookup_rows) + 1]] <- data.frame(
          gene = deg_base$gene[i], mapped_group = "Case",
          log2fc = deg_base$log2fc[i], padj = deg_base$padj[i],
          pval = if (!is.null(pval_col) && pval_col %in% colnames(deg_base)) deg_base[[pval_col]][i] else NA_real_,
          contrast = if (has_contrast) deg_base$comparison_contrast[i] else NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(lookup_rows) > 0) {
      deg_lookup <- do.call(rbind, lookup_rows)

      for (i in seq_along(p_int$x$data)) {
        tr <- p_int$x$data[[i]]
        if (!is.null(tr$mode) && grepl("markers", tr$mode) && !is.null(tr$text)) {
          new_text <- sapply(tr$text, function(txt) {
            gene_match <- trimws(regmatches(txt, regexpr("(?<=Gene: )[^\n<]+", txt, perl = TRUE)))
            grp_match <- trimws(regmatches(txt, regexpr(paste0("(?<=", group_col, ": )[^\n<]+"), txt, perl = TRUE)))
            if (length(grp_match) == 0)
              grp_match <- trimws(regmatches(txt, regexpr("(?<=fill: )[^\n<]+", txt, perl = TRUE)))
            if (length(gene_match) > 0 && length(grp_match) > 0) {
              hit <- deg_lookup[which(deg_lookup$gene == gene_match & deg_lookup$mapped_group == grp_match), ]
              if (nrow(hit) > 0) {
                extra <- ""
                if (!is.na(hit$contrast[1]))
                  extra <- paste0(extra, "<br>Comparison: ", hit$contrast[1])
                extra <- paste0(extra, "<br>log2FC: ", round(hit$log2fc[1], 3))
                if (!is.na(hit$pval[1]))
                  extra <- paste0(extra, "<br>pvalue: ", signif(hit$pval[1], 3))
                extra <- paste0(extra, "<br>padj: ", signif(hit$padj[1], 3))
                if (!is.null(sig_threshold_label))
                  extra <- paste0(extra, "<br>Cutoff: ", gsub("<", "&lt;", sig_threshold_label))
                return(paste0(txt, extra))
              }
            }
            txt
          }, USE.NAMES = FALSE)
          p_int$x$data[[i]]$text <- new_text
        }
      }
    }
  }

  # GSVA stats hover (ref lines 1390-1427)
  if (!is.null(gsva_stats_df) && nrow(gsva_stats_df) > 0) {
    for (i in seq_along(p_int$x$data)) {
      tr <- p_int$x$data[[i]]
      if (!is.null(tr$mode) && grepl("markers", tr$mode) && !is.null(tr$text)) {
        new_text <- sapply(tr$text, function(txt) {
          sig_match <- trimws(regmatches(txt, regexpr("(?<=Signature: )[^\n<]+", txt, perl = TRUE)))
          grp_match <- trimws(regmatches(txt, regexpr("(?<=grp_val: |fill: )[^\n<]+", txt, perl = TRUE)))
          if (length(sig_match) == 0)
            sig_match <- trimws(regmatches(txt, regexpr("(?<=x: )[^\n<]+", txt, perl = TRUE)))
          if (length(sig_match) > 0 && length(grp_match) > 0) {
            hit <- gsva_stats_df[which(gsva_stats_df$signature == sig_match &
                                        gsva_stats_df$group == grp_match), ]
            if (nrow(hit) == 0 && "treatment" %in% colnames(gsva_stats_df)) {
              hit <- gsva_stats_df[which(gsva_stats_df$signature == sig_match), ]
            }
            if (nrow(hit) > 0) {
              pval_col <- if ("P.Value" %in% colnames(hit)) "P.Value"
                          else if ("pvalue" %in% colnames(hit)) "pvalue" else NULL
              padj_col <- if ("adj.P.Val" %in% colnames(hit)) "adj.P.Val"
                          else if ("padj" %in% colnames(hit)) "padj" else NULL
              extra <- ""
              if ("logFC" %in% colnames(hit))
                extra <- paste0(extra, "<br>log2FC: ", round(hit$logFC[1], 3))
              if (!is.null(pval_col))
                extra <- paste0(extra, "<br>pvalue: ", signif(hit[[pval_col]][1], 3))
              if (!is.null(padj_col))
                extra <- paste0(extra, "<br>padj: ", signif(hit[[padj_col]][1], 3))
              if (!is.null(padj_col))
                extra <- paste0(extra, "<br>Cutoff: padj&lt;0.05")
              return(paste0(txt, extra))
            }
          }
          txt
        }, USE.NAMES = FALSE)
        p_int$x$data[[i]]$text <- new_text
      }
    }
  }

  # Numeric x-axis conversion (ref lines 1428-1483)
  box_indices <- c()
  group_names_ordered <- c()

  for (i in seq_along(p_int$x$data)) {
    tr <- p_int$x$data[[i]]
    if (!is.null(tr$type) && tr$type == "box") {
      box_indices <- c(box_indices, i)
      if (!is.null(tr$name) && nchar(tr$name) > 0 && !(tr$name %in% group_names_ordered))
        group_names_ordered <- c(group_names_ordered, tr$name)
    }
  }

  if (length(box_indices) > 0) {
    categories <- p_int$x$layout$xaxis$categoryarray
    if (is.null(categories)) {
      all_cat_vals <- c()
      for (i in box_indices) all_cat_vals <- c(all_cat_vals, unique(p_int$x$data[[i]]$x))
      categories <- sort(unique(all_cat_vals))
    }

    n_groups <- length(group_names_ordered)
    dodge_width <- 0.75

    group_offsets <- setNames(
      dodge_width * (seq_len(n_groups) - (n_groups + 1) / 2) / n_groups,
      group_names_ordered
    )

    for (i in box_indices) {
      tr <- p_int$x$data[[i]]
      grp_name <- tr$name
      offset <- if (!is.null(grp_name) && grp_name %in% names(group_offsets))
                  group_offsets[[grp_name]] else 0

      p_int$x$data[[i]]$x <- tr$x + offset
      p_int$x$data[[i]]$width <- dodge_width / n_groups * 0.85
    }

    axis_config <- list(
      tickvals = seq_along(categories),
      ticktext = categories,
      tickmode = "array"
    )
    layout_names <- names(p_int$x$layout)
    xaxis_names <- layout_names[grepl("^xaxis", layout_names)]
    if (length(xaxis_names) == 0) xaxis_names <- "xaxis"
    for (ax_name in xaxis_names) {
      if (is.null(p_int$x$layout[[ax_name]])) p_int$x$layout[[ax_name]] <- list()
      p_int$x$layout[[ax_name]]$tickvals <- axis_config$tickvals
      p_int$x$layout[[ax_name]]$ticktext <- axis_config$ticktext
      p_int$x$layout[[ax_name]]$tickmode <- axis_config$tickmode
    }
  }

  # Enlarge Comparison legend markers in HTML (asterisk, ~3x default)
  for (i in seq_along(p_int$x$data)) {
    tr <- p_int$x$data[[i]]
    if (!is.null(tr$mode) && grepl("markers", tr$mode) &&
        !is.null(tr$marker$opacity) && tr$marker$opacity == 0 &&
        isTRUE(tr$showlegend)) {
      p_int$x$data[[i]]$marker$size <- 18
      p_int$x$data[[i]]$marker$symbol <- "asterisk"
      p_int$x$data[[i]]$marker$opacity <- 1
      p_int$x$data[[i]]$marker$line <- list(width = 3, color = p_int$x$data[[i]]$marker$color)
    }
  }

  htmlwidgets::saveWidget(p_int, file.path(output_dir, paste0(base_name, ".html")), selfcontained = TRUE)
  message(paste("Saved:", paste0(base_name, ".html")))
}

# =============================================================================
# Expression Boxplot Builder
# Direct port of target_query_plot.R lines 385-519
# =============================================================================
build_expr_plot <- function(expr_long, group_col, group_levels, group_colors,
                            control_group, sig_annot, plot_title, plot_subtitle,
                            gene_cols_in_targets, facet_col = NULL,
                            colored_sig_annot = NULL, n_comparisons = NULL,
                            all_comparisons = NULL) {
  p <- ggplot(expr_long, aes(x = Gene, y = Expression, fill = .data[[group_col]])) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
               alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = group_colors) +
    labs(title = plot_title, subtitle = plot_subtitle,
         x = "Gene", y = "Normalized Expression", fill = "Group") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), scales = "free_y")
  }

  if (!is.null(colored_sig_annot) || (!is.null(all_comparisons) && length(all_comparisons) > 1)) {
    n_comps <- if (!is.null(n_comparisons)) n_comparisons
               else if (!is.null(all_comparisons)) length(all_comparisons)
               else if (!is.null(colored_sig_annot) && nrow(colored_sig_annot) > 0) length(unique(colored_sig_annot$comp_label))
               else 0

    if (n_comps == 1 && !is.null(colored_sig_annot) && nrow(colored_sig_annot) > 0) {
      p <- p + geom_text(
        data = colored_sig_annot,
        aes(x = x_pos, y = y_pos, label = label),
        colour = "black", hjust = 0.5, vjust = 0, size = 8, fontface = "bold",
        inherit.aes = FALSE, show.legend = FALSE
      ) +
      coord_cartesian(ylim = c(NA, max(max(colored_sig_annot$y_pos),
                                        max(expr_long$Expression, na.rm = TRUE)) * 1.1), clip = "off")
    } else if (n_comps > 1) {
      if (!is.null(all_comparisons)) {
        comp_colors <- setNames(
          sapply(all_comparisons, function(c) c$color),
          sapply(all_comparisons, function(c) c$label)
        )
        comp_colors <- comp_colors[order(names(comp_colors))]
      } else if (!is.null(colored_sig_annot) && nrow(colored_sig_annot) > 0) {
        comp_colors <- setNames(colored_sig_annot$comp_color, colored_sig_annot$comp_label)
        comp_colors <- comp_colors[!duplicated(names(comp_colors))]
      } else {
        comp_colors <- character(0)
      }

      if (!is.null(colored_sig_annot) && nrow(colored_sig_annot) > 0) {
        p <- p + geom_text(
          data = colored_sig_annot,
          aes(x = x_pos, y = y_pos, label = label, colour = comp_label),
          hjust = 0.5, vjust = 0, size = 8, fontface = "bold",
          inherit.aes = FALSE, show.legend = FALSE
        ) +
        coord_cartesian(ylim = c(NA, max(max(colored_sig_annot$y_pos),
                                          max(expr_long$Expression, na.rm = TRUE)) * 1.1), clip = "off")
      }

      if (length(comp_colors) > 0) {
        legend_df <- data.frame(
          x = rep(gene_cols_in_targets[1], length(comp_colors)),
          y = rep(-Inf, length(comp_colors)),
          comp_label = names(comp_colors),
          stringsAsFactors = FALSE
        )
        p <- p + geom_point(
          data = legend_df,
          aes(x = x, y = y, colour = comp_label),
          size = 0, alpha = 0, inherit.aes = FALSE, show.legend = TRUE
        ) +
        scale_colour_manual(values = comp_colors, name = "Comparison",
                            breaks = sort(names(comp_colors)),
                            guide = guide_legend(override.aes = list(size = 10, alpha = 1)))
      }
    }
  } else if (!is.null(sig_annot) && nrow(sig_annot) > 0) {
    # Complete grid expansion for correct position_dodge alignment
    complete_grid <- expand.grid(
      Gene = gene_cols_in_targets,
      group = group_levels,
      stringsAsFactors = FALSE
    )
    if (!is.null(facet_col) && facet_col %in% colnames(sig_annot)) {
      facet_vals <- unique(sig_annot[[facet_col]])
      complete_grid <- expand.grid(
        Gene = gene_cols_in_targets,
        group = group_levels,
        facet_var = facet_vals,
        stringsAsFactors = FALSE
      )
      colnames(complete_grid)[colnames(complete_grid) == "facet_var"] <- facet_col
      complete_grid <- merge(complete_grid, sig_annot, by = c("Gene", "group", facet_col), all.x = TRUE)
    } else {
      complete_grid <- merge(complete_grid, sig_annot[, c("Gene", "group", "y_pos", "label")],
                             by = c("Gene", "group"), all.x = TRUE)
    }
    complete_grid$label[is.na(complete_grid$label)] <- ""
    # Fill y_pos for empty cells from gene max (per-facet if applicable)
    if (!is.null(facet_col) && facet_col %in% colnames(expr_long) && facet_col %in% colnames(complete_grid)) {
      gene_max_fill <- expr_long %>%
        group_by(Gene, .data[[facet_col]]) %>%
        summarize(max_expr = max(Expression, na.rm = TRUE), .groups = "drop")
      complete_grid <- merge(complete_grid, gene_max_fill, by = c("Gene", facet_col), all.x = TRUE)
    } else {
      gene_max_fill <- expr_long %>%
        group_by(Gene) %>%
        summarize(max_expr = max(Expression, na.rm = TRUE), .groups = "drop")
      complete_grid <- merge(complete_grid, gene_max_fill, by = "Gene", all.x = TRUE)
    }
    complete_grid$y_pos[is.na(complete_grid$y_pos)] <- complete_grid$max_expr[is.na(complete_grid$y_pos)] * 1.005
    complete_grid$Gene <- factor(complete_grid$Gene, levels = gene_cols_in_targets)
    complete_grid$group <- factor(complete_grid$group, levels = group_levels)

    p <- p + geom_text(
      data = complete_grid,
      aes(x = Gene, y = y_pos, label = label, group = group),
      position = position_dodge(width = 0.75),
      vjust = 0, size = 6, fontface = "bold",
      inherit.aes = FALSE, show.legend = FALSE
    )
  }
  return(p)
}

# =============================================================================
# GSVA Boxplot Builder
# Direct port of target_query_plot.R lines 522-627
# =============================================================================
build_gsva_plot <- function(gsva_long, gsva_sig_annot, group_levels, group_colors,
                            control_group, gsva_title, gsva_subtitle,
                            facet_col = NULL, all_comparisons = NULL) {
  p <- ggplot(gsva_long, aes(x = Signature, y = GSVA_Score, fill = grp_val)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
               alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = group_colors) +
    labs(title = gsva_title, subtitle = gsva_subtitle,
         x = "Signature", y = "GSVA Score", fill = "Group") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), scales = "free_y")
  }

  if (!is.null(gsva_sig_annot) && nrow(gsva_sig_annot) > 0) {
    sig_names <- levels(gsva_long$Signature)
    has_colored <- "comp_label" %in% colnames(gsva_sig_annot) && "comp_color" %in% colnames(gsva_sig_annot)

    if (has_colored && !is.null(all_comparisons) && length(all_comparisons) > 1) {
      colored_annot <- gsva_sig_annot[gsva_sig_annot$label != "", ]
      if (nrow(colored_annot) > 0) {
        p <- p + geom_text(
          data = colored_annot,
          aes(x = Signature, y = y_pos, label = label, group = grp_val, colour = comp_label),
          position = position_dodge(width = 0.75),
          vjust = 0, size = 6, fontface = "bold",
          inherit.aes = FALSE, show.legend = FALSE
        )
      }
    } else {
      # Standard black annotations with complete grid for position_dodge
      if (!is.null(facet_col) && facet_col %in% colnames(gsva_sig_annot)) {
        facet_vals <- unique(gsva_sig_annot[[facet_col]])
        complete_grid <- expand.grid(
          Signature = sig_names, grp_val = group_levels, facet_var = facet_vals,
          stringsAsFactors = FALSE
        )
        colnames(complete_grid)[colnames(complete_grid) == "facet_var"] <- facet_col
        complete_grid <- merge(complete_grid, gsva_sig_annot,
                               by = c("Signature", "grp_val", facet_col), all.x = TRUE)
      } else {
        complete_grid <- expand.grid(
          Signature = sig_names, grp_val = group_levels, stringsAsFactors = FALSE
        )
        merge_cols <- intersect(c("Signature", "grp_val", "y_pos", "label"), colnames(gsva_sig_annot))
        complete_grid <- merge(complete_grid, gsva_sig_annot[, merge_cols, drop = FALSE],
                               by = c("Signature", "grp_val"), all.x = TRUE)
      }
      complete_grid$label[is.na(complete_grid$label)] <- ""
      sig_max_fill <- gsva_long %>%
        group_by(Signature) %>%
        summarize(max_score = max(GSVA_Score, na.rm = TRUE), .groups = "drop")
      complete_grid <- merge(complete_grid, sig_max_fill, by = "Signature", all.x = TRUE)
      complete_grid$y_pos[is.na(complete_grid$y_pos)] <- complete_grid$max_score[is.na(complete_grid$y_pos)] * 1.005
      complete_grid$Signature <- factor(complete_grid$Signature, levels = sig_names)
      complete_grid$grp_val <- factor(complete_grid$grp_val, levels = group_levels)

      p <- p + geom_text(
        data = complete_grid,
        aes(x = Signature, y = y_pos, label = label, group = grp_val),
        position = position_dodge(width = 0.75),
        vjust = 0, size = 6, fontface = "bold",
        inherit.aes = FALSE, show.legend = FALSE
      )
    }
  }

  # Comparison legend (always visible when multiple comparisons)
  if (!is.null(all_comparisons) && length(all_comparisons) > 1) {
    comp_colors <- setNames(
      sapply(all_comparisons, function(c) c$color),
      sapply(all_comparisons, function(c) c$label)
    )
    comp_colors <- comp_colors[order(names(comp_colors))]
    sig_names <- levels(gsva_long$Signature)
    legend_df <- data.frame(
      x = rep(sig_names[1], length(comp_colors)),
      y = rep(-Inf, length(comp_colors)),
      comp_label = names(comp_colors),
      stringsAsFactors = FALSE
    )
    p <- p + geom_point(
      data = legend_df,
      aes(x = x, y = y, colour = comp_label),
      size = 0, alpha = 0, inherit.aes = FALSE, show.legend = TRUE
    ) +
    scale_colour_manual(values = comp_colors, name = "Comparison",
                        breaks = sort(names(comp_colors)),
                        guide = guide_legend(override.aes = list(size = 10, alpha = 1)))
  }

  return(p)
}

# =============================================================================
# Significance Annotation Builders
# =============================================================================

# Standard annotation (Engitix_FFPE, SPARC, external)
# Direct port of target_query_plot.R lines 630-711
build_sig_annotations <- function(deg_for_annot, expr_long, gene_cols_in_targets,
                                  group_levels, control_group, annotation_map,
                                  grouping_mode, sig_threshold = 0.05) {
  sig_annot <- data.frame(Gene = character(), group = character(),
                          y_pos = numeric(), label = character(), stringsAsFactors = FALSE)

  if (is.null(deg_for_annot) || nrow(deg_for_annot) == 0 || !"padj" %in% colnames(deg_for_annot)) {
    return(sig_annot)
  }

  gene_max <- expr_long %>%
    group_by(Gene) %>%
    summarize(max_expr = max(Expression, na.rm = TRUE), .groups = "drop")

  for (i in seq_len(nrow(deg_for_annot))) {
    gene <- deg_for_annot$gene[i]
    padj_val <- deg_for_annot$padj[i]
    if (is.na(padj_val) || padj_val >= sig_threshold) next
    if (!gene %in% gene_cols_in_targets) next

    asterisk <- get_asterisk(padj_val)
    if (asterisk == "") next

    matched_groups <- NULL

    if (!is.null(annotation_map)) {
      comp_id <- deg_for_annot$comparison_id[i]
      for (pattern in names(annotation_map)) {
        if (grepl(paste0("(^|[._])", pattern, "$"), comp_id)) {
          target_groups <- unlist(annotation_map[[pattern]])
          if (is.character(target_groups)) {
            matched_groups <- target_groups
          }
          break
        }
      }
    } else if (grouping_mode == "comparison") {
      matched_groups <- "Case"
    } else {
      contrast <- deg_for_annot$comparison_contrast[i]
      case_group_raw <- NA_character_
      if (grepl(" vs ", contrast)) {
        case_group_raw <- trimws(strsplit(contrast, " vs ")[[1]][1])
      } else if (grepl("_vs_", contrast)) {
        case_group_raw <- trimws(gsub("_", " ", strsplit(contrast, "_vs_")[[1]][1]))
      }
      if (!is.na(case_group_raw) && case_group_raw != "") {
        for (gl in group_levels) {
          if (is.na(gl) || gl == "") next
          if (tolower(case_group_raw) == tolower(gl) ||
              grepl(tolower(case_group_raw), tolower(gl), fixed = TRUE) ||
              grepl(tolower(gl), tolower(case_group_raw), fixed = TRUE)) {
            matched_groups <- gl
            break
          }
        }
      }
    }

    if (is.null(matched_groups)) next

    gene_max_val <- gene_max$max_expr[gene_max$Gene == gene]
    if (length(gene_max_val) == 0) next

    for (mg in matched_groups) {
      if (!is.null(control_group) && mg == control_group) next
      if (!mg %in% group_levels) next
      sig_annot <- rbind(sig_annot, data.frame(
        Gene = gene, group = mg, y_pos = gene_max_val * 1.005,
        label = asterisk, stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(sig_annot) > 0) {
    sig_annot <- sig_annot[!duplicated(sig_annot[, c("Gene", "group")]), ]
    sig_annot$Gene <- factor(sig_annot$Gene, levels = gene_cols_in_targets)
    sig_annot$group <- factor(sig_annot$group, levels = group_levels)
  }

  return(sig_annot)
}

# Varsity-specific significance annotations (dual thresholds, per-facet y_pos)
# Direct port of target_query_plot.R lines 714-784
build_varsity_sig_annotations <- function(deg_df, expr_long, gene_cols_in_targets,
                                          group_levels, annotation_map, varsity_rules,
                                          facet_col = "varsity_treatment") {
  sig_annot <- data.frame(Gene = character(), group = character(),
                          y_pos = numeric(), label = character(),
                          varsity_treatment = character(), stringsAsFactors = FALSE)

  if (is.null(deg_df) || nrow(deg_df) == 0) return(sig_annot)

  gene_max <- expr_long %>%
    group_by(Gene, .data[[facet_col]]) %>%
    summarize(max_expr = max(Expression, na.rm = TRUE), .groups = "drop")

  for (i in seq_len(nrow(deg_df))) {
    gene <- deg_df$gene[i]
    if (!gene %in% gene_cols_in_targets) next

    comp_id <- deg_df$comparison_id[i]

    if (grepl("Adalimumab", comp_id)) {
      sig_col <- varsity_rules$ada_sig_col
      threshold <- varsity_rules$ada_threshold
      treatment <- "Adalimumab"
    } else if (grepl("Vedolizumab", comp_id)) {
      sig_col <- varsity_rules$vedo_sig_col
      threshold <- varsity_rules$vedo_threshold
      treatment <- "Vedolizumab"
    } else {
      next
    }

    sig_val <- if (sig_col %in% colnames(deg_df)) deg_df[[sig_col]][i] else deg_df$padj[i]
    if (is.na(sig_val) || sig_val >= threshold) next

    asterisk <- get_asterisk(sig_val)
    if (asterisk == "") next

    matched_group <- NULL
    if (!is.null(annotation_map)) {
      for (pattern in names(annotation_map)) {
        if (grepl(pattern, comp_id, fixed = TRUE)) {
          matched_group <- annotation_map[[pattern]]
          break
        }
      }
    }
    if (is.null(matched_group)) next

    for (mg in matched_group) {
      if (!mg %in% group_levels) next
      gene_max_val <- gene_max$max_expr[gene_max$Gene == gene &
                                         gene_max[[facet_col]] == treatment]
      if (length(gene_max_val) == 0) next
      sig_annot <- rbind(sig_annot, data.frame(
        Gene = gene, group = mg, y_pos = gene_max_val[1] * 1.005,
        label = asterisk, varsity_treatment = treatment, stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(sig_annot) > 0) {
    sig_annot <- sig_annot[!duplicated(sig_annot[, c("Gene", "group", "varsity_treatment")]), ]
    sig_annot$Gene <- factor(sig_annot$Gene, levels = gene_cols_in_targets)
    sig_annot$group <- factor(sig_annot$group, levels = group_levels)
    colnames(sig_annot)[colnames(sig_annot) == "varsity_treatment"] <- facet_col
  }

  return(sig_annot)
}

# Yokohama-specific colored significance annotations (midpoint of numerator stages)
# Direct port of target_query_plot.R lines 787-850
build_yokohama_sig_annotations <- function(deg_df, expr_long, gene_cols_in_targets,
                                           group_levels, comparisons) {
  sig_annot <- data.frame(Gene = character(), x_pos = numeric(), y_pos = numeric(),
                          label = character(), comp_label = character(),
                          comp_color = character(), stringsAsFactors = FALSE)

  if (is.null(deg_df) || nrow(deg_df) == 0 || length(comparisons) == 0) return(sig_annot)

  padj_threshold <- 0.05

  gene_max <- expr_long %>%
    group_by(Gene) %>%
    summarize(max_expr = max(Expression, na.rm = TRUE), .groups = "drop")

  n_groups <- length(group_levels)
  dodge_width <- 0.75
  group_offsets <- setNames(
    dodge_width * (seq_len(n_groups) - (n_groups + 1) / 2) / n_groups,
    group_levels
  )

  annotation_count <- list()

  for (comp in comparisons) {
    comp_id <- comp$comparison_id
    comp_label <- comp$label
    comp_color <- comp$color
    numerator_stages <- unlist(comp$numerator_stages)

    comp_deg <- deg_df[grepl(paste0("(^|[._])", comp_id, "$"), deg_df$comparison_id), ]
    if (nrow(comp_deg) == 0) next

    for (gene in gene_cols_in_targets) {
      gene_deg <- comp_deg[comp_deg$gene == gene, ]
      if (nrow(gene_deg) == 0) next

      padj_val <- gene_deg$padj[1]
      if (is.na(padj_val) || padj_val >= padj_threshold) next

      asterisk <- get_asterisk(padj_val)
      if (asterisk == "") next

      gene_idx <- which(gene_cols_in_targets == gene)
      valid_stages <- intersect(numerator_stages, group_levels)
      if (length(valid_stages) == 0) next

      stage_offsets <- group_offsets[valid_stages]
      x_pos <- gene_idx + mean(stage_offsets)

      gene_max_val <- gene_max$max_expr[gene_max$Gene == gene]
      if (length(gene_max_val) == 0) next

      count_key <- gene
      existing_count <- if (is.null(annotation_count[[count_key]])) 0 else annotation_count[[count_key]]
      y_pos <- gene_max_val * 1.005 + existing_count * gene_max_val * 0.04
      annotation_count[[count_key]] <- existing_count + 1

      sig_annot <- rbind(sig_annot, data.frame(
        Gene = gene, x_pos = x_pos, y_pos = y_pos, label = asterisk,
        comp_label = comp_label, comp_color = comp_color, stringsAsFactors = FALSE
      ))
    }
  }

  return(sig_annot)
}

# =============================================================================
# Enrich expr_long with DEG stats for plotly hover
# Direct port of target_query_plot.R lines 995-1080
# =============================================================================
enrich_hover_text <- function(expr_long, deg_df, group_col, annotation_map = NULL,
                              control_group = NULL, sig_threshold_label = NULL) {
  expr_long$hover_text <- paste0(
    "Gene: ", expr_long$Gene,
    "\nGroup: ", expr_long[[group_col]],
    "\nExpression: ", round(expr_long$Expression, 3)
  )
  if (is.null(deg_df) || nrow(deg_df) == 0) return(expr_long)
  if (!all(c("gene", "log2fc", "padj") %in% colnames(deg_df))) return(expr_long)

  pval_col <- if ("pval" %in% colnames(deg_df)) "pval" else "pvalue"
  has_contrast <- "comparison_contrast" %in% colnames(deg_df)

  deg_summary <- deg_df %>%
    select(gene, comparison_id, log2fc, padj, any_of(c(pval_col, "comparison_contrast"))) %>%
    distinct()

  hover_rows <- list()
  if (!is.null(annotation_map)) {
    for (i in seq_len(nrow(deg_summary))) {
      comp_id <- deg_summary$comparison_id[i]
      for (pattern in names(annotation_map)) {
        if (grepl(pattern, comp_id, fixed = TRUE)) {
          matched <- unlist(annotation_map[[pattern]])
          if (is.character(matched) && length(matched) >= 1) {
            for (mg in matched) {
              if (!is.null(control_group) && mg == control_group) next
              hover_rows[[length(hover_rows) + 1]] <- data.frame(
                gene = deg_summary$gene[i],
                mapped_group = mg,
                log2fc = deg_summary$log2fc[i],
                padj = deg_summary$padj[i],
                pval = if (pval_col %in% colnames(deg_summary)) deg_summary[[pval_col]][i] else NA_real_,
                contrast = if (has_contrast) deg_summary$comparison_contrast[i] else NA_character_,
                stringsAsFactors = FALSE
              )
            }
          }
          break
        }
      }
    }
  } else {
    for (i in seq_len(nrow(deg_summary))) {
      hover_rows[[length(hover_rows) + 1]] <- data.frame(
        gene = deg_summary$gene[i],
        mapped_group = "Case",
        log2fc = deg_summary$log2fc[i],
        padj = deg_summary$padj[i],
        pval = if (pval_col %in% colnames(deg_summary)) deg_summary[[pval_col]][i] else NA_real_,
        contrast = if (has_contrast) deg_summary$comparison_contrast[i] else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(hover_rows) == 0) return(expr_long)
  hover_df <- do.call(rbind, hover_rows)

  for (i in seq_len(nrow(hover_df))) {
    gene_val <- hover_df$gene[i]
    grp_val <- hover_df$mapped_group[i]
    if (is.na(grp_val) || is.na(gene_val)) next
    mask <- expr_long$Gene == gene_val & expr_long[[group_col]] == grp_val
    mask[is.na(mask)] <- FALSE
    if (!any(mask)) next
    contrast_str <- if (!is.na(hover_df$contrast[i])) {
      paste0("\nComparison: ", hover_df$contrast[i])
    } else ""
    pval_str <- if (!is.na(hover_df$pval[i])) {
      paste0("\npvalue: ", signif(hover_df$pval[i], 3))
    } else ""
    threshold_str <- if (!is.null(sig_threshold_label)) {
      paste0("\nCutoff: ", sig_threshold_label)
    } else ""
    expr_long$hover_text[mask] <- paste0(
      expr_long$hover_text[mask],
      contrast_str,
      "\nlog2FC: ", round(hover_df$log2fc[i], 3),
      pval_str,
      "\nFDR: ", signif(hover_df$padj[i], 3),
      threshold_str
    )
  }
  return(expr_long)
}

# =============================================================================
# GSVA Analysis Runner
# Direct port of target_query_plot.R lines 1083-1258
# =============================================================================
run_gsva_analysis <- function(expr_df_subset, group_col, group_levels, group_colors,
                              control_group, gene_cols_for_gsva, custom_signatures,
                              comparisons = list()) {
  if (length(custom_signatures) == 0 || length(gene_cols_for_gsva) < 5) {
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  if (!requireNamespace("GSVA", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    message("  WARNING: GSVA or limma not available. Skipping.")
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  num_cols <- gene_cols_for_gsva[sapply(expr_df_subset[, gene_cols_for_gsva, drop = FALSE], is.numeric)]
  gsva_mat <- as.matrix(t(expr_df_subset[, num_cols, drop = FALSE]))
  native_mask <- rowSums(abs(gsva_mat)) > 0
  gsva_mat <- gsva_mat[native_mask, , drop = FALSE]
  gsva_mat <- gsva_mat[!duplicated(rownames(gsva_mat)), , drop = FALSE]
  message(paste("  GSVA matrix:", nrow(gsva_mat), "genes x", ncol(gsva_mat), "samples"))

  sigs_filtered <- lapply(custom_signatures, function(genes) {
    intersect(genes, rownames(gsva_mat))
  })
  sigs_filtered <- sigs_filtered[sapply(sigs_filtered, length) >= 3]

  if (length(sigs_filtered) == 0) {
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  for (sig_name in names(sigs_filtered)) {
    message(paste("  Signature", sig_name, ":",
                  length(sigs_filtered[[sig_name]]), "of",
                  length(custom_signatures[[sig_name]]), "genes in data;",
                  nrow(gsva_mat) - length(sigs_filtered[[sig_name]]), "background genes"))
  }

  gsva_params <- GSVA::gsvaParam(exprData = gsva_mat, geneSets = sigs_filtered, kcdf = "Gaussian")
  gsva_scores <- GSVA::gsva(gsva_params, verbose = FALSE)

  # Limma-based significance testing
  sample_groups <- factor(expr_df_subset[[group_col]], levels = group_levels)
  names(sample_groups) <- rownames(expr_df_subset)

  gsva_stats_rows <- list()

  if (length(comparisons) > 0) {
    for (comp in comparisons) {
      num_stages <- unlist(comp$numerator_stages)
      comp_label <- comp$label
      sample_binary <- ifelse(sample_groups %in% num_stages, "num", "denom")
      sample_binary <- factor(sample_binary, levels = c("denom", "num"))
      if (min(table(sample_binary)) < 2) next

      design <- model.matrix(~ 0 + sample_binary)
      colnames(design) <- c("denom", "num")
      fit <- limma::lmFit(gsva_scores, design)
      contrasts_mat <- limma::makeContrasts(contrasts = "num-denom", levels = design)
      fit2 <- limma::contrasts.fit(fit, contrasts_mat)
      fit2 <- tryCatch(limma::eBayes(fit2), error = function(e) NULL)
      if (is.null(fit2)) next
      results <- limma::topTable(fit2, number = Inf, sort.by = "none")
      if (nrow(results) == nrow(gsva_scores)) {
        rownames(results) <- rownames(gsva_scores)
      }

      for (sig_name in rownames(results)) {
        gsva_stats_rows[[length(gsva_stats_rows) + 1]] <- data.frame(
          signature = sig_name, group = comp_label,
          vs_control = paste(setdiff(group_levels, num_stages), collapse = ","),
          n_group = sum(sample_binary == "num"),
          n_control = sum(sample_binary == "denom"),
          logFC = results[sig_name, "logFC"],
          P.Value = results[sig_name, "P.Value"],
          adj.P.Val = results[sig_name, "adj.P.Val"],
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    design <- model.matrix(~ 0 + sample_groups)
    colnames(design) <- levels(sample_groups)
    fit <- limma::lmFit(gsva_scores, design)

    for (grp in setdiff(group_levels, control_group)) {
      if (sum(sample_groups == grp) < 2) next
      contrast_str <- paste0("`", grp, "` - `", control_group, "`")
      contrasts_mat <- tryCatch(
        limma::makeContrasts(contrasts = contrast_str, levels = design),
        error = function(e) NULL
      )
      if (is.null(contrasts_mat)) next

      fit2 <- limma::contrasts.fit(fit, contrasts_mat)
      fit2 <- tryCatch(limma::eBayes(fit2), error = function(e) NULL)
      if (is.null(fit2)) next
      results <- limma::topTable(fit2, number = Inf, sort.by = "none")
      if (nrow(results) == nrow(gsva_scores)) {
        rownames(results) <- rownames(gsva_scores)
      }

      for (sig_name in rownames(results)) {
        gsva_stats_rows[[length(gsva_stats_rows) + 1]] <- data.frame(
          signature = sig_name, group = grp, vs_control = control_group,
          n_group = sum(sample_groups == grp),
          n_control = sum(sample_groups == control_group),
          logFC = results[sig_name, "logFC"],
          P.Value = results[sig_name, "P.Value"],
          adj.P.Val = results[sig_name, "adj.P.Val"],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  gsva_stats_df <- if (length(gsva_stats_rows) > 0) do.call(rbind, gsva_stats_rows) else {
    message("  NOTE: No GSVA limma contrasts could be computed (insufficient samples per group)")
    NULL
  }

  gsva_long <- gsva_scores %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Signature") %>%
    pivot_longer(cols = -Signature, names_to = "Sample", values_to = "GSVA_Score")

  sample_groups_df <- data.frame(Sample = rownames(expr_df_subset),
                              grp_val = as.character(expr_df_subset[[group_col]]),
                              stringsAsFactors = FALSE)
  gsva_long <- gsva_long %>% left_join(sample_groups_df, by = "Sample")
  gsva_long$Signature <- factor(gsva_long$Signature, levels = rownames(gsva_scores))
  gsva_long$grp_val <- factor(gsva_long$grp_val, levels = group_levels)

  # GSVA significance annotations
  gsva_sig_annot <- data.frame(Signature = character(), grp_val = character(),
                               y_pos = numeric(), label = character(),
                               comp_label = character(), comp_color = character(),
                               stringsAsFactors = FALSE)
  if (!is.null(gsva_stats_df) && nrow(gsva_stats_df) > 0) {
    sig_max <- gsva_long %>%
      group_by(Signature) %>%
      summarize(max_score = max(GSVA_Score, na.rm = TRUE), .groups = "drop")

    comp_to_group <- list()
    comp_to_color <- list()
    if (length(comparisons) > 0) {
      for (comp in comparisons) {
        first_num <- intersect(unlist(comp$numerator_stages), group_levels)[1]
        if (!is.na(first_num)) comp_to_group[[comp$label]] <- first_num
        comp_to_color[[comp$label]] <- comp$color
      }
    }

    for (i in seq_len(nrow(gsva_stats_df))) {
      padj <- gsva_stats_df$adj.P.Val[i]
      asterisk <- get_asterisk(padj)
      if (asterisk == "") next
      sig_name <- gsva_stats_df$signature[i]
      grp <- gsva_stats_df$group[i]
      annot_grp <- if (!is.null(comp_to_group[[grp]])) comp_to_group[[grp]] else grp
      if (!annot_grp %in% group_levels) next
      max_val <- sig_max$max_score[sig_max$Signature == sig_name]
      if (length(max_val) == 0) next
      annot_color <- if (!is.null(comp_to_color[[grp]])) comp_to_color[[grp]] else "black"
      annot_comp_label <- if (!is.null(comp_to_color[[grp]])) grp else ""
      gsva_sig_annot <- rbind(gsva_sig_annot, data.frame(
        Signature = sig_name, grp_val = annot_grp,
        y_pos = max_val * 1.005, label = asterisk,
        comp_label = annot_comp_label, comp_color = annot_color,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(gsva_sig_annot) > 0) {
    gsva_sig_annot$Signature <- factor(gsva_sig_annot$Signature, levels = rownames(gsva_scores))
    gsva_sig_annot$grp_val <- factor(gsva_sig_annot$grp_val, levels = group_levels)
  }

  return(list(gsva_scores = gsva_scores, limma_results = gsva_stats_df,
              p_gsva_long = gsva_long, gsva_sig_annot = gsva_sig_annot,
              gsva_stats_df = gsva_stats_df))
}

# =============================================================================
# Varsity GSVA: limma-based per treatment x week analysis
# Direct port of target_query_plot.R lines 853-992
# =============================================================================
run_varsity_gsva <- function(expr_df_subset, group_col, group_levels, group_colors,
                             custom_signatures, gene_cols_for_gsva, varsity_rules,
                             facet_col = "varsity_treatment") {
  if (length(custom_signatures) == 0 || length(gene_cols_for_gsva) < 5) {
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  if (!requireNamespace("GSVA", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  num_cols <- gene_cols_for_gsva[sapply(expr_df_subset[, gene_cols_for_gsva, drop = FALSE], is.numeric)]
  gsva_mat <- as.matrix(t(expr_df_subset[, num_cols, drop = FALSE]))
  native_mask <- rowSums(abs(gsva_mat)) > 0
  gsva_mat <- gsva_mat[native_mask, , drop = FALSE]
  gsva_mat <- gsva_mat[!duplicated(rownames(gsva_mat)), , drop = FALSE]
  message(paste("  GSVA matrix (Varsity):", nrow(gsva_mat), "genes x", ncol(gsva_mat), "samples"))

  sigs_filtered <- lapply(custom_signatures, function(genes) {
    intersect(genes, rownames(gsva_mat))
  })
  sigs_filtered <- sigs_filtered[sapply(sigs_filtered, length) >= 3]
  if (length(sigs_filtered) == 0) {
    return(list(p_gsva_long = NULL, gsva_sig_annot = NULL, gsva_stats_df = NULL))
  }

  gsva_params <- GSVA::gsvaParam(exprData = gsva_mat, geneSets = sigs_filtered, kcdf = "Gaussian")
  gsva_scores <- GSVA::gsva(gsva_params, verbose = FALSE)

  # Build long format with treatment info
  gsva_long <- gsva_scores %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Signature") %>%
    pivot_longer(cols = -Signature, names_to = "Sample", values_to = "GSVA_Score")

  sample_info <- data.frame(
    Sample = rownames(expr_df_subset),
    grp_val = as.character(expr_df_subset[[group_col]]),
    treatment = as.character(expr_df_subset[[facet_col]]),
    stringsAsFactors = FALSE
  )
  gsva_long <- gsva_long %>% left_join(sample_info, by = "Sample")
  gsva_long$Signature <- factor(gsva_long$Signature, levels = rownames(gsva_scores))
  gsva_long$grp_val <- factor(gsva_long$grp_val, levels = group_levels)

  # Limma per treatment x per week: Yes.Yes vs No.No within each week
  gsva_stats_rows <- list()
  treatments <- unique(sample_info$treatment[!is.na(sample_info$treatment)])
  weeks <- c("Wk0", "Wk14", "Wk52")

  for (trt in treatments) {
    for (wk in weeks) {
      wk_pattern <- paste0("^", wk, "_")
      wk_trt_samples <- sample_info$Sample[
        sample_info$treatment == trt & grepl(wk_pattern, sample_info$grp_val)
      ]
      if (length(wk_trt_samples) < 4) next

      wk_scores <- gsva_scores[, colnames(gsva_scores) %in% wk_trt_samples, drop = FALSE]
      wk_groups <- sample_info$grp_val[sample_info$Sample %in% colnames(wk_scores)]

      resp_label <- paste0(wk, "_Yes.Yes")
      nr_label <- paste0(wk, "_No.No")
      if (!resp_label %in% wk_groups || !nr_label %in% wk_groups) next

      binary_factor <- factor(ifelse(wk_groups == resp_label, "Resp", "NR"),
                              levels = c("NR", "Resp"))
      if (length(levels(binary_factor)) < 2 || min(table(binary_factor)) < 2) next

      design <- model.matrix(~ 0 + binary_factor)
      colnames(design) <- c("NR", "Resp")

      fit <- limma::lmFit(wk_scores, design)
      contrasts_mat <- limma::makeContrasts(contrasts = "NR-Resp", levels = design)
      fit2 <- limma::contrasts.fit(fit, contrasts_mat)
      fit2 <- tryCatch(limma::eBayes(fit2), error = function(e) NULL)
      if (is.null(fit2)) next

      results <- limma::topTable(fit2, number = Inf, sort.by = "none")
      if (nrow(results) == nrow(gsva_scores)) {
        rownames(results) <- rownames(gsva_scores)
      }

      for (sig_name in rownames(results)) {
        pval <- results[sig_name, "P.Value"]
        padj_val <- results[sig_name, "adj.P.Val"]
        gsva_stats_rows[[length(gsva_stats_rows) + 1]] <- data.frame(
          signature = sig_name, treatment = trt, timepoint = wk,
          P.Value = pval, adj.P.Val = padj_val,
          pvalue = pval, padj = padj_val,
          logFC = results[sig_name, "logFC"], sig_val = padj_val,
          n_signatures = nrow(results),
          threshold = 0.05,
          significant = padj_val < 0.05,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  gsva_stats_df <- if (length(gsva_stats_rows) > 0) do.call(rbind, gsva_stats_rows) else NULL

  # Build significance annotations per-facet
  gsva_sig_annot <- data.frame(Signature = character(), grp_val = character(),
                               y_pos = numeric(), label = character(),
                               varsity_treatment = character(), stringsAsFactors = FALSE)

  if (!is.null(gsva_stats_df) && nrow(gsva_stats_df) > 0) {
    sig_max <- gsva_long %>%
      group_by(Signature, treatment) %>%
      summarize(max_score = max(GSVA_Score, na.rm = TRUE), .groups = "drop")

    sig_rows <- gsva_stats_df[gsva_stats_df$significant == TRUE, ]
    for (i in seq_len(nrow(sig_rows))) {
      sig_name <- sig_rows$signature[i]
      trt <- sig_rows$treatment[i]
      wk <- sig_rows$timepoint[i]
      pval <- sig_rows$sig_val[i]
      asterisk <- get_asterisk(pval)
      if (asterisk == "") next

      max_val <- sig_max$max_score[sig_max$Signature == sig_name & sig_max$treatment == trt]
      if (length(max_val) == 0) next

      nr_grp <- paste0(wk, "_No.No")
      if (!nr_grp %in% group_levels) next
      gsva_sig_annot <- rbind(gsva_sig_annot, data.frame(
        Signature = sig_name, grp_val = nr_grp, y_pos = max_val[1] * 1.005,
        label = asterisk, varsity_treatment = trt, stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(gsva_sig_annot) > 0) {
    gsva_sig_annot$Signature <- factor(gsva_sig_annot$Signature, levels = rownames(gsva_scores))
    gsva_sig_annot$grp_val <- factor(gsva_sig_annot$grp_val, levels = group_levels)
    colnames(gsva_sig_annot)[colnames(gsva_sig_annot) == "varsity_treatment"] <- "treatment"
    gsva_sig_annot <- gsva_sig_annot[!duplicated(gsva_sig_annot[, c("Signature", "grp_val", "treatment")]), ]
  }

  return(list(p_gsva_long = gsva_long, gsva_sig_annot = gsva_sig_annot,
              gsva_stats_df = gsva_stats_df))
}
