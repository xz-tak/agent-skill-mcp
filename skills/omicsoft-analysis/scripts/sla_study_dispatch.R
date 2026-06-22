# =============================================================================
# Sample-Level Analysis: Study Dispatch Logic
# =============================================================================
# Per-study dispatch for internal and external studies.
# Direct port of target_query_plot.R study loop logic.
# =============================================================================

# =============================================================================
# Group Derivation Functions
# =============================================================================

derive_varsity_groups <- function(df) {
  message("  Deriving Varsity week_response groups...")
  if (!"metadata_Visit.Type" %in% colnames(df)) {
    message("    WARNING: metadata_Visit.Type not found. Skipping derivation.")
    df$week_response <- NA_character_
    return(df)
  }
  df$varsity_week <- dplyr::case_when(
    df$metadata_Visit.Type == "Visit 1" ~ "Wk0",
    df$metadata_Visit.Type == "Visit 9" ~ "Wk14",
    df$metadata_Visit.Type == "Visit 28" ~ "Wk52"
  )
  df$varsity_treatment <- dplyr::case_when(
    grepl("Adalimumab", df$metadata_Planned.Treatment.for.Period.01) ~ "Adalimumab",
    grepl("Vedolizumab", df$metadata_Planned.Treatment.for.Period.01) ~ "Vedolizumab"
  )
  df$varsity_response <- dplyr::case_when(
    df$metadata_Clinical.Response.at.Week.14 == "Yes" &
      df$metadata_Clinical.Remission.at.Week.52 == "Yes" ~ "Yes.Yes",
    df$metadata_Clinical.Response.at.Week.14 == "No" &
      df$metadata_Clinical.Remission.at.Week.52 == "No" ~ "No.No"
  )
  df$week_response <- paste0(df$varsity_week, "_", df$varsity_response)
  valid_levels <- c("Wk0_No.No", "Wk0_Yes.Yes", "Wk14_No.No",
                    "Wk14_Yes.Yes", "Wk52_No.No", "Wk52_Yes.Yes")
  df$week_response[!df$week_response %in% valid_levels] <- NA
  n_valid <- sum(!is.na(df$week_response))
  message(paste("    Valid samples:", n_valid, "/", nrow(df)))
  df
}

derive_yokohama_rna_fibrosis <- function(df) {
  df$yokohama_rna_fibrosis <- paste0("F", as.integer(df$metadata_stage))
  valid_levels <- c("F0", "F1", "F2", "F3", "F4")
  df$yokohama_rna_fibrosis[!df$yokohama_rna_fibrosis %in% valid_levels] <- NA
  df
}

derive_yokohama_rna_nash <- function(df) {
  nas_score <- suppressWarnings(as.numeric(df$metadata_nas_score))
  stage <- suppressWarnings(as.numeric(df$metadata_stage))
  df$yokohama_rna_nash <- NA_character_
  df$yokohama_rna_nash[nas_score >= 4 & stage >= 2] <- "at_risk"
  df$yokohama_rna_nash[!is.na(nas_score) & !(nas_score >= 4 & stage >= 2)] <- "control"
  df
}

derive_yokohama_rna_diagnosis <- function(df) {
  df$yokohama_rna_diagnosis <- NA_character_
  valid_diag <- df$metadata_diagnosis %in% c("NAFL", "MASH")
  df$yokohama_rna_diagnosis[valid_diag] <- df$metadata_diagnosis[valid_diag]
  df
}

derive_yokohama_prot_fibrosis <- function(df) {
  df$yokohama_prot_fibrosis <- as.character(df$metadata_Fibrosis)
  valid_levels <- c("Healthy", "F0", "F1", "F2", "F3", "F4")
  df$yokohama_prot_fibrosis[!df$yokohama_prot_fibrosis %in% valid_levels] <- NA
  df
}

derive_yokohama_prot_nash <- function(df) {
  nas_val <- suppressWarnings(as.numeric(df$metadata_NAS))
  fibrosis <- as.character(df$metadata_Fibrosis)
  df$yokohama_prot_nash <- NA_character_
  at_risk_mask <- nas_val >= 4 & fibrosis %in% c("F2", "F3", "F4")
  control_mask <- !is.na(nas_val) & !at_risk_mask
  df$yokohama_prot_nash[at_risk_mask] <- "at_risk"
  df$yokohama_prot_nash[control_mask] <- "control"
  df
}

derive_yokohama_prot_diagnosis <- function(df) {
  df$yokohama_prot_diagnosis <- NA_character_
  valid_diag <- df$metadata_diagnosis %in% c("NAFL", "MASH")
  df$yokohama_prot_diagnosis[valid_diag] <- df$metadata_diagnosis[valid_diag]
  df
}

derive_sparc_tissue <- function(df) {
  if ("metadata_CHARACTERISTICS_BIO_MATERIAL" %in% colnames(df)) {
    df$sparc_tissue <- ifelse(
      grepl("ileum|Ileum", df$metadata_CHARACTERISTICS_BIO_MATERIAL, ignore.case = TRUE),
      "ileum", "nonileum"
    )
  } else {
    df$sparc_tissue <- df$meta_tissue
  }
  df
}

subset_sparc_disease <- function(df) {
  diagnosis <- as.character(df$metadata_DIAGNOSIS)
  valid_mask <- !is.na(diagnosis) & diagnosis != "" &
                diagnosis != "-1" & !grepl("Unclassified", diagnosis, ignore.case = TRUE)
  df_valid <- df[valid_mask, , drop = FALSE]
  disease_label_map <- c(
    "Crohn's Disease" = "CD", "crohn's disease (CD)" = "CD",
    "Ulcerative Colitis" = "UC", "ulcerative colitis (UC)" = "UC"
  )
  diseases <- unique(df_valid$metadata_DIAGNOSIS)
  diseases <- diseases[!diseases %in% c("", "-1", "NA", "\\N", "IBD Unclassified")]
  subsets <- list()
  for (d in diseases) {
    short <- if (d %in% names(disease_label_map)) disease_label_map[[d]] else d
    mask <- df_valid$metadata_DIAGNOSIS == d
    if (sum(mask) > 0) subsets[[short]] <- df_valid[mask, , drop = FALSE]
  }
  subsets
}

# =============================================================================
# Internal Study Handler
# Direct port of target_query_plot.R lines 2026-2601
# =============================================================================
run_internal_study <- function(study_name, study_df, expr_mat, deg_df,
                               config, target_genes, signatures,
                               gene_cols, output_dir) {

  message(paste("\n========================================"))
  message(paste("Processing Internal Study:", study_name))
  message(paste("========================================"))

  study_output <- file.path(output_dir, sanitize_path(study_name))
  dir.create(study_output, recursive = TRUE, showWarnings = FALSE)

  message(paste("Study:", study_name, "| Samples:", nrow(study_df)))

  # Derive Varsity groups if needed
  if (isTRUE(config$derive_groups) && study_name == "Varsity") {
    study_df <- derive_varsity_groups(study_df)
    study_df <- study_df[study_df$week_response %in% unlist(config$group_levels), ]
    message(paste("  Varsity: derived week_response, remaining samples:", nrow(study_df)))
  }

  # Determine grouping configs (multi-grouping for Yokohama)
  if (!is.null(config$groupings)) {
    grouping_list <- config$groupings
  } else {
    grouping_list <- list(default = list(
      group_col = config$group_col,
      group_levels = config$group_levels,
      control = config$control,
      group_colors = config$group_colors,
      annotation_map = config$annotation_map,
      facet_col = config$facet_col,
      disease_col = config$disease_col
    ))
  }

  all_gsva_stats <- list()
  gene_cols_in_targets <- intersect(target_genes, colnames(study_df))
  if (length(gene_cols_in_targets) == 0) {
    gene_cols_in_targets <- intersect(target_genes, rownames(expr_mat))
  }

  for (grp_name in names(grouping_list)) {
    grp_cfg <- grouping_list[[grp_name]]
    annotation_map <- grp_cfg$annotation_map
    group_col <- grp_cfg$group_col
    group_levels <- unlist(grp_cfg$group_levels)
    control_group <- grp_cfg$control
    group_colors <- unlist(grp_cfg$group_colors)
    facet_col <- grp_cfg$facet_col
    disease_col <- grp_cfg$disease_col

    grp_suffix <- if (grp_name == "default") "" else paste0("_", grp_name)

    # Apply group derivation
    view_df <- study_df
    if (isTRUE(grp_cfg$derive_fibrosis_from_stage)) {
      view_df <- derive_yokohama_rna_fibrosis(view_df)
    } else if (isTRUE(grp_cfg$derive_nash_nas4_stage2)) {
      view_df <- derive_yokohama_rna_nash(view_df)
    } else if (isTRUE(grp_cfg$derive_diagnosis_from_metadata)) {
      view_df <- derive_yokohama_rna_diagnosis(view_df)
    } else if (isTRUE(grp_cfg$derive_fibrosis_from_metadata)) {
      view_df <- derive_yokohama_prot_fibrosis(view_df)
    } else if (isTRUE(grp_cfg$derive_nash_nas4_fib2)) {
      view_df <- derive_yokohama_prot_nash(view_df)
    } else if (isTRUE(grp_cfg$derive_diagnosis_from_prot_metadata)) {
      view_df <- derive_yokohama_prot_diagnosis(view_df)
    }

    if (is.null(group_col) || !group_col %in% colnames(view_df)) {
      message(paste("    SKIP grouping", grp_name, "- group_col not found:", group_col))
      next
    }

    grp_expr <- view_df[view_df[[group_col]] %in% group_levels, ]
    if (nrow(grp_expr) == 0) next
    grp_expr[[group_col]] <- factor(grp_expr[[group_col]], levels = group_levels)
    message(paste("  Grouping:", grp_name, "|", paste(group_levels, collapse = " | "),
                  "| Samples:", nrow(grp_expr)))

    # Derive tissue facet for SPARC
    if (!is.null(facet_col) && facet_col == "sparc_tissue") {
      grp_expr <- derive_sparc_tissue(grp_expr)
    }

    # Determine disease subsets for SPARC
    disease_subsets <- list(list(label = NULL, filter = NULL))
    if (!is.null(disease_col) && disease_col %in% colnames(grp_expr)) {
      subsets <- subset_sparc_disease(grp_expr)
      if (length(subsets) > 1) {
        disease_subsets <- lapply(names(subsets), function(nm) {
          list(label = nm, filter = nm, data = subsets[[nm]])
        })
      }
    }

    if (nrow(grp_expr) < 4) {
      message(paste("  SKIP:", study_name, grp_name, "- too few samples after grouping"))
      next
    }

    for (ds in disease_subsets) {
      ds_expr <- grp_expr
      file_suffix <- grp_suffix
      title_suffix <- if (grp_name == "default") "" else paste0(" [", grp_name, "]")

      if (!is.null(ds$label)) {
        if (!is.null(ds$data)) {
          ds_expr <- ds$data
        } else {
          ds_expr <- grp_expr[grp_expr[[disease_col]] == ds$filter, ]
        }
        file_suffix <- paste0(grp_suffix, "_", ds$label)
        title_suffix <- paste0(title_suffix, " (", ds$label, ")")
        message(paste("    Disease subset:", ds$label, "| Samples:", nrow(ds_expr)))
        if (nrow(ds_expr) < 4) {
          message(paste("    SKIP:", ds$label, "- too few samples"))
          next
        }
      }

      # Identify gene columns available in this subset
      found_targets <- intersect(target_genes, colnames(ds_expr))
      if (length(found_targets) == 0) {
        message("    WARNING: No target genes in data. Skipping.")
        next
      }
      gene_cols_for_gsva <- intersect(gene_cols, colnames(ds_expr))

      # Build expression long format
      keep_cols <- c(found_targets, group_col)
      if (!is.null(facet_col) && facet_col %in% colnames(ds_expr)) {
        keep_cols <- c(keep_cols, facet_col)
      }
      expr_long <- ds_expr[, keep_cols, drop = FALSE] %>%
        tibble::rownames_to_column("Sample") %>%
        pivot_longer(cols = all_of(found_targets), names_to = "Gene", values_to = "Expression")
      expr_long$Gene <- factor(expr_long$Gene, levels = found_targets)

      # Filter DEG
      ds_deg <- deg_df
      if (!is.null(ds$label) && !is.null(ds_deg) && nrow(ds_deg) > 0 &&
          "comparison_id" %in% colnames(ds_deg)) {
        ds_pattern <- ds$label
        ds_deg_filtered <- ds_deg[grepl(ds_pattern, ds_deg$comparison_id, ignore.case = TRUE), ]
        if (nrow(ds_deg_filtered) > 0) ds_deg <- ds_deg_filtered
      }

      # SPARC tissue-split DEG
      if (!is.null(facet_col) && facet_col == "sparc_tissue" && !is.null(ds_deg) && nrow(ds_deg) > 0) {
        tissue_vals <- unique(ds_expr[[facet_col]])
        tissue_vals <- tissue_vals[!is.na(tissue_vals)]
        has_tissue_in_comp <- any(grepl("ileum|nonileum", ds_deg$comparison_id, ignore.case = TRUE))
        if (has_tissue_in_comp && length(tissue_vals) > 0) {
          ds_deg_by_tissue <- list()
          for (tv in tissue_vals) {
            tissue_pattern <- if (tv == "ileum") "_ileum_" else "_nonileum_"
            tv_deg <- ds_deg[grepl(tissue_pattern, ds_deg$comparison_id, fixed = TRUE), ]
            if (nrow(tv_deg) > 0) ds_deg_by_tissue[[tv]] <- tv_deg
          }
          attr(ds_deg, "tissue_split") <- ds_deg_by_tissue
        }
      }

      # Significance annotations (three paths)
      varsity_rules <- config$varsity_sig_rules
      colored_sig_annot <- NULL

      sig_threshold_val <- if (!is.null(config$sig_threshold)) config$sig_threshold else 0.05

      if (!is.null(varsity_rules)) {
        sig_annot <- build_varsity_sig_annotations(
          ds_deg, expr_long, found_targets, group_levels,
          annotation_map, varsity_rules, facet_col = facet_col
        )
      } else if (!is.null(grp_cfg$comparisons) && length(grp_cfg$comparisons) > 0) {
        sig_annot <- data.frame(Gene = character(), group = character(),
                                y_pos = numeric(), label = character(), stringsAsFactors = FALSE)
        colored_sig_annot <- build_yokohama_sig_annotations(
          ds_deg, expr_long, found_targets, group_levels, grp_cfg$comparisons
        )
      } else {
        tissue_split <- attr(ds_deg, "tissue_split")
        if (!is.null(tissue_split) && length(tissue_split) > 0 &&
            !is.null(facet_col) && facet_col %in% colnames(expr_long)) {
          sig_annot_list <- list()
          for (tv in names(tissue_split)) {
            tv_expr <- expr_long[expr_long[[facet_col]] == tv, ]
            if (nrow(tv_expr) == 0) next
            tv_annot <- build_sig_annotations(tissue_split[[tv]], tv_expr, found_targets,
                                               group_levels, control_group, annotation_map,
                                               "internal_metadata", sig_threshold = sig_threshold_val)
            if (nrow(tv_annot) > 0) {
              tv_annot[[facet_col]] <- tv
              sig_annot_list[[length(sig_annot_list) + 1]] <- tv_annot
            }
          }
          sig_annot <- if (length(sig_annot_list) > 0) do.call(rbind, sig_annot_list)
                       else data.frame(Gene = character(), group = character(),
                                       y_pos = numeric(), label = character(), stringsAsFactors = FALSE)
        } else {
          sig_annot <- build_sig_annotations(ds_deg, expr_long, found_targets,
                                             group_levels, control_group, annotation_map,
                                             "internal_metadata", sig_threshold = sig_threshold_val)
        }
      }

      # Study-specific threshold label
      sig_threshold_label <- if (!is.null(varsity_rules)) {
        "Ada: pval<0.01 | Vedo: padj<0.05"
      } else {
        sig_col_name <- if (!is.null(config$sig_col)) config$sig_col else "padj"
        paste0(sig_col_name, "<", sig_threshold_val)
      }

      # Enrich hover text
      expr_long <- enrich_hover_text(expr_long, ds_deg, group_col, annotation_map,
                                     control_group = control_group,
                                     sig_threshold_label = sig_threshold_label)

      plot_title <- paste0("Target Gene Expression — ", study_name, title_suffix)
      ctrl_label <- if (is.null(control_group)) "baseline" else control_group
      plot_subtitle <- if (!is.null(varsity_rules)) {
        paste0("Ada: * pval<", varsity_rules$ada_threshold,
               " | Vedo: * padj<", varsity_rules$vedo_threshold, " (R vs NR, * on NR)")
      } else if (!is.null(grp_cfg$comparisons) && length(grp_cfg$comparisons) > 0) {
        paste0("Colored * = ", sig_threshold_label, " (per comparison)")
      } else {
        paste0("* ", sig_threshold_label, " vs ", ctrl_label)
      }

      p_expr <- build_expr_plot(expr_long, group_col, group_levels, group_colors,
                                control_group, sig_annot, plot_title, plot_subtitle,
                                found_targets, facet_col = facet_col,
                                colored_sig_annot = colored_sig_annot,
                                n_comparisons = length(grp_cfg$comparisons),
                                all_comparisons = grp_cfg$comparisons)
      save_plot_pair(p_expr, study_output, paste0("target_expression", file_suffix),
                     deg_df = ds_deg, group_col = group_col,
                     annotation_map = annotation_map, control_group = control_group,
                     sig_threshold_label = sig_threshold_label)

      # GSVA routing (three paths)
      gsva_result <- NULL
      if (study_name == "Varsity" && !is.null(varsity_rules)) {
        gsva_result <- run_varsity_gsva(ds_expr, group_col, group_levels, group_colors,
                                        signatures, gene_cols_for_gsva, varsity_rules,
                                        facet_col = facet_col)
        if (!is.null(gsva_result$p_gsva_long)) {
          gsva_title <- paste0("GSVA Signature Scores — ", study_name, title_suffix)
          gsva_subtitle <- "* adj.P.Val < 0.05 (limma, per week NR vs R)"
          p_gsva <- build_gsva_plot(gsva_result$p_gsva_long, gsva_result$gsva_sig_annot,
                                    group_levels, group_colors, control_group,
                                    gsva_title, gsva_subtitle, facet_col = "treatment")
          save_plot_pair(p_gsva, study_output, paste0("signature_gsva", file_suffix),
                         gsva_stats_df = gsva_result$gsva_stats_df)
        }
      } else if (!is.null(control_group) && !is.null(facet_col) && facet_col == "sparc_tissue" &&
                 facet_col %in% colnames(ds_expr)) {
        # SPARC: compute GSVA on ALL samples, then run limma per-tissue
        tissue_vals <- unique(ds_expr[[facet_col]])
        tissue_vals <- tissue_vals[!is.na(tissue_vals)]
        num_cols <- gene_cols_for_gsva[sapply(ds_expr[, gene_cols_for_gsva, drop = FALSE], is.numeric)]
        gsva_mat <- as.matrix(t(ds_expr[, num_cols, drop = FALSE]))
        native_mask <- rowSums(abs(gsva_mat)) > 0
        gsva_mat <- gsva_mat[native_mask, , drop = FALSE]
        gsva_mat <- gsva_mat[!duplicated(rownames(gsva_mat)), , drop = FALSE]
        sigs_filtered <- lapply(signatures, function(genes) intersect(genes, rownames(gsva_mat)))
        sigs_filtered <- sigs_filtered[sapply(sigs_filtered, length) >= 3]

        if (length(sigs_filtered) > 0 && ncol(gsva_mat) >= 5 &&
            requireNamespace("GSVA", quietly = TRUE) && requireNamespace("limma", quietly = TRUE)) {
          message(paste("  GSVA matrix (SPARC):", nrow(gsva_mat), "genes x", ncol(gsva_mat), "samples"))
          gsva_params <- GSVA::gsvaParam(exprData = gsva_mat, geneSets = sigs_filtered, kcdf = "Gaussian")
          gsva_scores_mat <- GSVA::gsva(gsva_params, verbose = FALSE)

          sample_info <- data.frame(
            Sample = rownames(ds_expr),
            grp_val = as.character(ds_expr[[group_col]]),
            tissue = as.character(ds_expr[[facet_col]]),
            stringsAsFactors = FALSE
          )
          gsva_long_full <- gsva_scores_mat %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Signature") %>%
            pivot_longer(cols = -Signature, names_to = "Sample", values_to = "GSVA_Score") %>%
            left_join(sample_info, by = "Sample")
          gsva_long_full$Signature <- factor(gsva_long_full$Signature, levels = rownames(gsva_scores_mat))
          gsva_long_full$grp_val <- factor(gsva_long_full$grp_val, levels = group_levels)
          colnames(gsva_long_full)[colnames(gsva_long_full) == "tissue"] <- facet_col

          # Per-tissue limma
          gsva_stats_parts <- list()
          gsva_annot_parts <- list()
          for (tv in tissue_vals) {
            tv_samples <- sample_info$Sample[sample_info$tissue == tv]
            tv_groups <- sample_info$grp_val[sample_info$tissue == tv]
            names(tv_groups) <- tv_samples
            tv_groups <- tv_groups[!is.na(tv_groups)]
            if (length(unique(tv_groups)) < 2 || sum(tv_groups == control_group) < 2) next
            valid_samples <- intersect(names(tv_groups), colnames(gsva_scores_mat))
            if (length(valid_samples) < 2) next
            tv_groups <- tv_groups[valid_samples]
            panel_gsva <- gsva_scores_mat[, valid_samples, drop = FALSE]
            present_levels <- intersect(group_levels, unique(tv_groups))
            grp_factor <- factor(tv_groups, levels = present_levels)
            design <- model.matrix(~ 0 + grp_factor)
            safe_levels <- make.names(present_levels)
            colnames(design) <- safe_levels
            fit <- limma::lmFit(panel_gsva, design)
            test_groups <- setdiff(present_levels, control_group)
            for (grp in test_groups) {
              safe_grp <- make.names(grp)
              safe_ctrl <- make.names(control_group)
              contrast_str <- paste0(safe_grp, " - ", safe_ctrl)
              contrasts_mat <- tryCatch(limma::makeContrasts(contrasts = contrast_str, levels = design), error = function(e) NULL)
              if (is.null(contrasts_mat)) next
              fit2 <- limma::contrasts.fit(fit, contrasts_mat)
              fit2 <- tryCatch(limma::eBayes(fit2), error = function(e) NULL)
              if (is.null(fit2)) next
              results <- limma::topTable(fit2, number = Inf, sort.by = "none")
              if (nrow(results) == nrow(gsva_scores_mat)) rownames(results) <- rownames(gsva_scores_mat)
              for (sig_name in rownames(results)) {
                padj_val <- results[sig_name, "adj.P.Val"]
                gsva_stats_parts[[length(gsva_stats_parts) + 1]] <- data.frame(
                  signature = sig_name, group = grp, vs_control = control_group,
                  tissue = tv, logFC = results[sig_name, "logFC"],
                  P.Value = results[sig_name, "P.Value"], adj.P.Val = padj_val,
                  n_group = sum(tv_groups == grp), n_control = sum(tv_groups == control_group),
                  stringsAsFactors = FALSE
                )
                asterisk <- get_asterisk(padj_val)
                if (asterisk != "") {
                  sig_max <- max(gsva_long_full$GSVA_Score[gsva_long_full$Signature == sig_name &
                                                            gsva_long_full[[facet_col]] == tv], na.rm = TRUE)
                  gsva_annot_parts[[length(gsva_annot_parts) + 1]] <- data.frame(
                    Signature = sig_name, grp_val = grp,
                    y_pos = sig_max * 1.005, label = asterisk,
                    sparc_tissue = tv, stringsAsFactors = FALSE
                  )
                }
              }
            }
          }

          combined_gsva_annot <- if (length(gsva_annot_parts) > 0) {
            annot_df <- do.call(rbind, gsva_annot_parts)
            annot_df$Signature <- factor(annot_df$Signature, levels = rownames(gsva_scores_mat))
            annot_df$grp_val <- factor(annot_df$grp_val, levels = group_levels)
            annot_df
          } else {
            data.frame(Signature = character(), grp_val = character(),
                       y_pos = numeric(), label = character(),
                       sparc_tissue = character(), stringsAsFactors = FALSE)
          }
          combined_gsva_stats <- if (length(gsva_stats_parts) > 0) do.call(rbind, gsva_stats_parts) else NULL

          gsva_title <- paste0("GSVA Signature Scores — ", study_name, title_suffix)
          gsva_subtitle <- paste0("* adj.P.Val < 0.05 vs ", control_group, " (limma, per tissue)")
          p_gsva <- build_gsva_plot(gsva_long_full, combined_gsva_annot,
                                    group_levels, group_colors, control_group,
                                    gsva_title, gsva_subtitle, facet_col = facet_col)
          save_plot_pair(p_gsva, study_output, paste0("signature_gsva", file_suffix),
                         gsva_stats_df = combined_gsva_stats)
          gsva_result <- list(p_gsva_long = gsva_long_full, gsva_sig_annot = combined_gsva_annot,
                              gsva_stats_df = combined_gsva_stats, gsva_scores = gsva_scores_mat)
        }
      } else if (!is.null(control_group)) {
        # Standard limma GSVA (Yokohama, Engitix)
        gsva_comparisons <- if (!is.null(grp_cfg$comparisons)) grp_cfg$comparisons else list()
        gsva_result <- run_gsva_analysis(ds_expr, group_col, group_levels, group_colors,
                                         control_group, gene_cols_for_gsva, signatures,
                                         comparisons = gsva_comparisons)
        if (!is.null(gsva_result$p_gsva_long)) {
          gsva_title <- paste0("GSVA Signature Scores — ", study_name, title_suffix)
          gsva_subtitle <- paste0("* adj.P.Val < 0.05 vs ", control_group, " (limma)")
          p_gsva <- build_gsva_plot(gsva_result$p_gsva_long, gsva_result$gsva_sig_annot,
                                    group_levels, group_colors, control_group,
                                    gsva_title, gsva_subtitle,
                                    all_comparisons = grp_cfg$comparisons)
          save_plot_pair(p_gsva, study_output, paste0("signature_gsva", file_suffix),
                         gsva_stats_df = gsva_result$gsva_stats_df)
        }
      }

      if (!is.null(gsva_result$gsva_stats_df)) {
        gsva_result$gsva_stats_df$view <- grp_name
        all_gsva_stats[[length(all_gsva_stats) + 1]] <- gsva_result$gsva_stats_df
      }
    }
  }

  # Stats CSVs
  if (!is.null(deg_df) && nrow(deg_df) > 0) {
    write.csv(deg_df, file.path(study_output, "target_comparison_stats.csv"), row.names = FALSE)
    message(paste("  Saved: target_comparison_stats.csv (", nrow(deg_df), "rows)"))
  }
  if (length(all_gsva_stats) > 0) {
    combined_gsva_stats <- do.call(rbind, all_gsva_stats)
    write.csv(combined_gsva_stats, file.path(study_output, "gsva_comparison_stats.csv"), row.names = FALSE)
    message(paste("  Saved: gsva_comparison_stats.csv"))
  }

  # Correlation modules (pooled across all study samples)
  message("\n  --- Running Correlation Modules ---")
  found_targets_all <- intersect(target_genes, colnames(study_df))
  if (nrow(study_df) >= 5 && length(found_targets_all) >= 2) {
    run_gene_gene_correlation(
      expr_mat = t(study_df[, found_targets_all, drop = FALSE]),
      metadata = study_df,
      target_genes = found_targets_all,
      output_dir = study_output,
      study_name = study_name
    )
  }

  # GSVA-GSVA correlation
  if (!exists("gsva_result")) gsva_result <- NULL
  if (!is.null(gsva_result) && !is.null(gsva_result$p_gsva_long) &&
      length(unique(gsva_result$p_gsva_long$Signature)) >= 2) {
    gsva_scores_for_corr <- gsva_result$p_gsva_long %>%
      dplyr::select(Signature, Sample, GSVA_Score) %>%
      pivot_wider(names_from = Sample, values_from = GSVA_Score) %>%
      tibble::column_to_rownames("Signature") %>%
      as.matrix()
    run_gsva_gsva_correlation(
      gsva_scores = gsva_scores_for_corr,
      output_dir = study_output,
      study_name = study_name
    )
  }

  # Target vs continuous
  continuous_cols <- unlist(config$continuous_cols)
  if (!is.null(continuous_cols) && length(continuous_cols) > 0 && length(found_targets_all) >= 1) {
    run_target_vs_continuous(
      expr_mat = t(study_df[, found_targets_all, drop = FALSE]),
      metadata = study_df,
      target_genes = found_targets_all,
      continuous_cols = continuous_cols,
      output_dir = study_output,
      study_name = study_name
    )
    if (!is.null(gsva_result) && !is.null(gsva_result$gsva_scores)) {
      run_gsva_vs_continuous(
        gsva_scores = gsva_result$gsva_scores,
        metadata = study_df,
        continuous_cols = continuous_cols,
        output_dir = study_output,
        study_name = study_name
      )
    }
  }
}

# =============================================================================
# External Study Handler
# Direct port of target_query_plot.R lines 2606-2719
# Flat output: files at study_output root, not subdirectories
# =============================================================================
run_external_study <- function(study_name, study_df, expr_mat, deg_df,
                               target_genes, signatures, gene_cols, output_dir) {

  message(paste("\n========================================"))
  message(paste("Processing External Study:", study_name))
  message(paste("========================================"))

  study_output <- file.path(output_dir, sanitize_path(study_name))
  dir.create(study_output, recursive = TRUE, showWarnings = FALSE)

  group_levels <- c("Control", "Case")
  control_group <- "Control"
  group_colors <- c("Control" = "lightblue", "Case" = "#F8766D")
  ext_padj_threshold <- 0.05

  if (!"meta_comparison_group" %in% colnames(study_df)) {
    message(paste("  SKIP:", study_name, "- no meta_comparison_group for comparison-based grouping"))
    return(invisible(NULL))
  }

  if (is.null(deg_df) || nrow(deg_df) == 0) {
    message(paste("  SKIP:", study_name, "- no DEG data for comparisons"))
    return(invisible(NULL))
  }

  comparisons <- unique(deg_df$comparison_id)
  message(paste("  Comparisons:", length(comparisons)))

  all_gsva_stats <- list()
  found_targets <- intersect(target_genes, colnames(study_df))
  if (length(found_targets) == 0) {
    message(paste("  SKIP:", study_name, "- no target genes in data"))
    return(invisible(NULL))
  }
  gene_cols_for_gsva <- intersect(gene_cols, colnames(study_df))

  for (comp_id in comparisons) {
    comp_contrast <- deg_df$comparison_contrast[deg_df$comparison_id == comp_id][1]
    comp_label <- gsub("[^A-Za-z0-9_]", "_", comp_contrast)
    comp_label <- gsub("_+", "_", comp_label)
    comp_label <- gsub("^_|_$", "", comp_label)
    if (nchar(comp_label) > 80) comp_label <- substr(comp_label, 1, 80)
    message(paste("    Processing comparison:", comp_contrast))

    group_col <- "query_group"
    comp_expr <- study_df
    comp_expr[[group_col]] <- assign_comparison_roles(
      comp_expr$meta_comparison_group, comp_id
    )
    comp_expr <- comp_expr[!is.na(comp_expr[[group_col]]), ]
    comp_expr[[group_col]] <- factor(comp_expr[[group_col]], levels = group_levels)

    if (nrow(comp_expr) < 4) {
      message(paste("    SKIP comparison:", comp_contrast, "- too few samples"))
      next
    }
    message(paste("    Samples: Case=", sum(comp_expr[[group_col]] == "Case"),
                  "Control=", sum(comp_expr[[group_col]] == "Control")))

    # Expression long format
    expr_long <- comp_expr[, c(found_targets, group_col), drop = FALSE] %>%
      tibble::rownames_to_column("Sample") %>%
      pivot_longer(cols = all_of(found_targets), names_to = "Gene", values_to = "Expression")
    expr_long$Gene <- factor(expr_long$Gene, levels = found_targets)

    # Sig annotations for this comparison
    comp_deg <- deg_df[deg_df$comparison_id == comp_id, ]
    sig_annot <- build_sig_annotations(comp_deg, expr_long, found_targets,
                                       group_levels, control_group, NULL, "comparison",
                                       sig_threshold = ext_padj_threshold)

    # Enrich hover text
    expr_long <- enrich_hover_text(expr_long, comp_deg, group_col, NULL,
                                   control_group = control_group,
                                   sig_threshold_label = paste0("padj<", ext_padj_threshold))

    plot_title <- paste("Target Gene Expression —", study_name)
    plot_subtitle <- paste0(comp_contrast, "\n* padj < ", ext_padj_threshold, " vs Control")

    p_expr <- build_expr_plot(expr_long, group_col, group_levels, group_colors,
                              control_group, sig_annot, plot_title, plot_subtitle,
                              found_targets)
    # Flat output naming at study_output root
    save_plot_pair(p_expr, study_output, paste0("target_expression_", comp_label),
                   deg_df = comp_deg, group_col = group_col, annotation_map = NULL,
                   control_group = control_group,
                   sig_threshold_label = paste0("padj<", ext_padj_threshold))

    # GSVA per comparison
    gsva_result <- run_gsva_analysis(comp_expr, group_col, group_levels, group_colors,
                                     control_group, gene_cols_for_gsva, signatures)
    if (!is.null(gsva_result$p_gsva_long)) {
      gsva_title <- paste("GSVA Signature Scores —", study_name)
      gsva_subtitle <- paste0(comp_contrast, "\n* adj.P.Val < 0.05 vs Control (limma)")
      p_gsva <- build_gsva_plot(gsva_result$p_gsva_long, gsva_result$gsva_sig_annot,
                                group_levels, group_colors, control_group,
                                gsva_title, gsva_subtitle)
      save_plot_pair(p_gsva, study_output, paste0("signature_gsva_", comp_label),
                     gsva_stats_df = gsva_result$gsva_stats_df)
      if (!is.null(gsva_result$gsva_stats_df)) {
        all_gsva_stats[[length(all_gsva_stats) + 1]] <- gsva_result$gsva_stats_df
      }
    }
  }

  # Stats CSVs at study level
  if (length(all_gsva_stats) > 0) {
    combined_gsva_stats <- do.call(rbind, all_gsva_stats)
    write.csv(combined_gsva_stats, file.path(study_output, "gsva_comparison_stats.csv"), row.names = FALSE)
    message(paste("  Saved: gsva_comparison_stats.csv"))
  }
  write.csv(deg_df, file.path(study_output, "target_comparison_stats.csv"), row.names = FALSE)
  message(paste("  Saved: target_comparison_stats.csv (", nrow(deg_df), "rows, all comparisons)"))

  # Gene-Gene Correlation ONCE pooled per study (not per-comparison)
  message("\n  --- Running Correlation Modules ---")
  if (nrow(study_df) >= 5 && length(found_targets) >= 2) {
    run_gene_gene_correlation(
      expr_mat = t(study_df[, found_targets, drop = FALSE]),
      metadata = study_df,
      target_genes = found_targets,
      output_dir = study_output,
      study_name = study_name
    )
  }
}

# =============================================================================
# Comparison Role Assignment for External Studies
# =============================================================================
assign_comparison_roles <- function(comparison_group_col, comparison_id) {
  comp_base <- comparison_id
  comp_base_escaped <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", comp_base)
  roles <- rep(NA_character_, length(comparison_group_col))
  for (i in seq_along(comparison_group_col)) {
    cg <- comparison_group_col[i]
    if (is.na(cg) || cg == "" || cg == "\\N") next
    if (grepl(paste0(comp_base_escaped, "@case"), cg)) {
      roles[i] <- "Case"
    } else if (grepl(paste0(comp_base_escaped, "@control"), cg)) {
      roles[i] <- "Control"
    }
  }
  return(roles)
}

# =============================================================================
# Utility: sanitize_path (if not already sourced)
# =============================================================================
if (!exists("sanitize_path")) {
  sanitize_path <- function(x) {
    gsub("[^A-Za-z0-9_]", "_", x)
  }
}
