#!/usr/bin/env Rscript
# =============================================================================
# CellChat Cell-Cell Communication Analysis (Parameterized)
# =============================================================================
# End-to-end CellChat analysis with command-line argument support
# Outputs organized into dedicated subfolders per analysis step
#
# Usage:
#   Rscript cellchat_analysis.R -i <input.rds> -t <cell_type_col> [options]
#
# Required:
#   -i, --input       Input Seurat RDS file (must be RDS, not h5ad)
#   -t, --cell-type   Metadata column with cell type labels
#
# Optional:
#   -c, --condition   Metadata column for condition comparison
#   --conditions      Specific conditions to analyze (comma-separated)
#   --reference       Reference condition for comparisons (auto-detects Normal/Control)
#   -o, --output      Output directory [default: cellchat_results]
#   -s, --species     Species: human or mouse [default: human]
#   --min-cells       Minimum cells per group [default: 10]
#   --trim            Trim value for truncatedMean [default: 0.1]
#   --dry-run         Print config and exit without running
#   --help            Show this help message
#
# Output Structure:
#   01_data_summary/     - Input data statistics
#   02_database/         - CellChatDB info
#   03_heatmaps/         - Interaction count/weight heatmaps + data
#   04_circle_plots/     - Aggregated network circles + data
#   05_bubble_plots/     - L-R bubble plots by category
#   06_signaling_roles/  - Sender/receiver analysis + data
#   07_LR_pairs/         - Significant L-R pair tables
#   08_pathways/         - Active pathway lists
#   09_compare_total/    - Cross-condition comparison + data
#   10_differential/     - Differential networks + data
#   11_info_flow/        - Pathway ranking + data
#   12_pathway_networks/ - Individual pathway network plots
#   13_functional_sim/   - Functional similarity + data
#   14_structural_sim/   - Structural similarity + data
#   15_summary/          - Summary statistics
#   16_pathway_matrix/   - Pathway presence matrix
#   17_pathway_flow/     - Pathway information flow
#   18_pathway_stats/    - Differential pathway stats + figures
#   objects/             - R objects for downstream analysis
#
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Parse Command Line Arguments
# -----------------------------------------------------------------------------

if (!requireNamespace("optparse", quietly = TRUE)) {
  cat("ERROR: optparse package required. Install with: install.packages('optparse')\n")
  quit(status = 1)
}

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input Seurat RDS file (required)", metavar = "FILE"),
  make_option(c("-t", "--cell-type"), type = "character", default = NULL,
              help = "Cell type column name (required)", metavar = "COLUMN"),
  make_option(c("-c", "--condition"), type = "character", default = NULL,
              help = "Condition column for comparison (optional)", metavar = "COLUMN"),
  make_option("--conditions", type = "character", default = NULL,
              help = "Specific conditions to analyze, comma-separated (optional)", metavar = "LIST"),
  make_option(c("-o", "--output"), type = "character", default = "cellchat_results",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-s", "--species"), type = "character", default = "human",
              help = "Species: human or mouse [default: %default]"),
  make_option("--min-cells", type = "integer", default = 10,
              help = "Minimum cells per group [default: %default]", metavar = "N"),
  make_option("--trim", type = "double", default = 0.1,
              help = "Trim value for truncatedMean [default: %default]"),
  make_option("--reference", type = "character", default = NULL,
              help = "Reference condition for comparisons (auto-detects Normal/Control if not specified)", metavar = "CONDITION"),
  make_option("--dry-run", action = "store_true", default = FALSE,
              help = "Print configuration and exit without running")
)

parser <- OptionParser(
  usage = "%prog -i <input.rds> -t <cell_type_col> [options]",
  option_list = option_list,
  description = "CellChat cell-cell communication analysis pipeline"
)

opt <- parse_args(parser)

# Validate required arguments
if (is.null(opt$input)) {
  cat("ERROR: Input file required. Use -i or --input\n\n")
  print_help(parser)
  quit(status = 1)
}

if (is.null(opt$`cell-type`)) {
  cat("ERROR: Cell type column required. Use -t or --cell-type\n\n")
  print_help(parser)
  quit(status = 1)
}

if (!file.exists(opt$input)) {
  cat(sprintf("ERROR: Input file not found: %s\n", opt$input))
  quit(status = 1)
}

file_ext <- tolower(tools::file_ext(opt$input))
if (file_ext == "h5ad") {
  cat("ERROR: Input is h5ad format. Convert to Seurat RDS first:\n")
  cat(sprintf("  Rscript convert_h5ad_to_seurat.R %s\n", opt$input))
  quit(status = 1)
}

# Build CONFIG
CONFIG <- list(
  input_file = opt$input,
  cell_type_col = opt$`cell-type`,
  condition_col = opt$condition,
  conditions_filter = if (!is.null(opt$conditions)) strsplit(opt$conditions, ",")[[1]] else NULL,
  reference = opt$reference,
  output_dir = opt$output,
  report_file = "cellchat_report.md",
  species = opt$species,
  min_cells = opt$`min-cells`,
  trim = opt$trim,
  n_workers = 4
)

# Create output directory structure
# 01 and 02 are flat files in root; 03-18 get subfolders
DIRS <- list(
  root = CONFIG$output_dir,
  d03 = file.path(CONFIG$output_dir, "03_heatmaps"),
  d04 = file.path(CONFIG$output_dir, "04_circle_plots"),
  d05 = file.path(CONFIG$output_dir, "05_bubble_plots"),
  d06 = file.path(CONFIG$output_dir, "06_signaling_roles"),
  d07 = file.path(CONFIG$output_dir, "07_LR_pairs"),
  d08 = file.path(CONFIG$output_dir, "08_pathways"),
  d09 = file.path(CONFIG$output_dir, "09_compare_total"),
  d10 = file.path(CONFIG$output_dir, "10_differential"),
  d11 = file.path(CONFIG$output_dir, "11_info_flow"),
  d12 = file.path(CONFIG$output_dir, "12_pathway_networks"),
  d13 = file.path(CONFIG$output_dir, "13_functional_sim"),
  d14 = file.path(CONFIG$output_dir, "14_structural_sim"),
  d15 = file.path(CONFIG$output_dir, "15_summary"),
  d16 = file.path(CONFIG$output_dir, "16_pathway_matrix"),
  d17 = file.path(CONFIG$output_dir, "17_pathway_flow"),
  d18 = file.path(CONFIG$output_dir, "18_pathway_stats"),
  objects = file.path(CONFIG$output_dir, "objects")
)

for (d in DIRS) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Print Configuration
# -----------------------------------------------------------------------------
cat("=============================================================================\n")
cat("CellChat Analysis Pipeline\n")
cat("=============================================================================\n")
cat(sprintf("Start time: %s\n\n", Sys.time()))

cat("Configuration:\n")
cat(sprintf("  Input file:    %s\n", CONFIG$input_file))
cat(sprintf("  Cell type col: %s\n", CONFIG$cell_type_col))
cat(sprintf("  Condition col: %s\n", ifelse(is.null(CONFIG$condition_col), "None (single group)", CONFIG$condition_col)))
if (!is.null(CONFIG$conditions_filter)) {
  cat(sprintf("  Conditions:    %s\n", paste(CONFIG$conditions_filter, collapse = ", ")))
}
cat(sprintf("  Output dir:    %s\n", CONFIG$output_dir))
cat(sprintf("  Species:       %s\n", CONFIG$species))
cat(sprintf("  Min cells:     %d\n", CONFIG$min_cells))
cat(sprintf("  Trim value:    %.2f\n", CONFIG$trim))
cat("\n")
flush.console()

if (opt$`dry-run`) {
  cat("[DRY RUN] Configuration printed. Exiting without analysis.\n")
  quit(status = 0)
}

# Figure/table counter for report
fig_counter <- 1
table_counter <- 1
results_log <- list()

log_figure <- function(subdir, filename, description, legend) {
  results_log[[length(results_log) + 1]] <<- list(
    type = "figure",
    num = fig_counter,
    subdir = subdir,
    filename = filename,
    description = description,
    legend = legend
  )
  fig_counter <<- fig_counter + 1
  return(file.path(subdir, filename))
}

log_table <- function(subdir, filename, description, legend) {
  results_log[[length(results_log) + 1]] <<- list(
    type = "table",
    num = table_counter,
    subdir = subdir,
    filename = filename,
    description = description,
    legend = legend
  )
  table_counter <<- table_counter + 1
  return(file.path(subdir, filename))
}

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
cat("Step 1: Loading libraries...\n"); flush.console()

tryCatch({
  library(reticulate)
  conda_envs <- c("r_env", "base")
  for (env in conda_envs) {
    if (tryCatch({use_condaenv(env, required = FALSE); TRUE}, error = function(e) FALSE)) {
      cat(sprintf("  - Python: %s\n", py_config()$python))
      break
    }
  }
}, error = function(e) {
  cat("  - Python: not configured (UMAP plots may fail)\n")
})
flush.console()

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(Matrix)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
cat(sprintf("  - CellChat version: %s\n", packageVersion("CellChat"))); flush.console()

# -----------------------------------------------------------------------------
# 2. Load and Prepare Data from Seurat RDS
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading Seurat RDS file...\n"); flush.console()

seurat_obj <- readRDS(CONFIG$input_file)
cat(sprintf("  - Loaded Seurat object: %d genes x %d cells\n", nrow(seurat_obj), ncol(seurat_obj))); flush.console()

cat("  - Extracting expression matrix from Seurat...\n"); flush.console()
data.input <- seurat_obj[["RNA"]]$data
cat(sprintf("  - Matrix class: %s\n", class(data.input)[1])); flush.console()
cat(sprintf("  - Matrix dimensions: %d genes x %d cells\n", nrow(data.input), ncol(data.input))); flush.console()

data_min <- min(data.input@x)
data_max <- max(data.input@x)
data_mean <- mean(data.input@x)
cat(sprintf("  - Data range: [%.4f, %.4f], mean (non-zero): %.4f\n", data_min, data_max, data_mean)); flush.console()

if (data_max < 20 && data_max > 0) {
  cat("  - Data is log-normalized (expected)\n"); flush.console()
} else if (data_max > 100) {
  cat("  - WARNING: Data appears to be raw counts!\n"); flush.console()
}

meta <- seurat_obj@meta.data

if (!CONFIG$cell_type_col %in% colnames(meta)) {
  cat(sprintf("ERROR: Cell type column '%s' not found in metadata.\n", CONFIG$cell_type_col))
  cat(sprintf("Available columns: %s\n", paste(colnames(meta), collapse = ", ")))
  quit(status = 1)
}

meta$labels <- as.character(meta[[CONFIG$cell_type_col]])

if (!is.null(CONFIG$condition_col)) {
  if (!CONFIG$condition_col %in% colnames(meta)) {
    cat(sprintf("ERROR: Condition column '%s' not found in metadata.\n", CONFIG$condition_col))
    cat(sprintf("Available columns: %s\n", paste(colnames(meta), collapse = ", ")))
    quit(status = 1)
  }
  meta$condition <- as.character(meta[[CONFIG$condition_col]])

  if (!is.null(CONFIG$conditions_filter)) {
    invalid_conds <- setdiff(CONFIG$conditions_filter, unique(meta$condition))
    if (length(invalid_conds) > 0) {
      cat(sprintf("WARNING: Conditions not found: %s\n", paste(invalid_conds, collapse = ", ")))
    }
    valid_conds <- intersect(CONFIG$conditions_filter, unique(meta$condition))
    if (length(valid_conds) == 0) {
      cat("ERROR: No valid conditions to analyze.\n")
      quit(status = 1)
    }
    meta <- meta[meta$condition %in% valid_conds, ]
    data.input <- data.input[, rownames(meta)]
    cat(sprintf("  - Filtered to conditions: %s\n", paste(valid_conds, collapse = ", ")))
  }
} else {
  meta$condition <- "All"
}

cat(sprintf("  - Cell type column: %s\n", CONFIG$cell_type_col))
cat(sprintf("  - Condition column: %s\n", ifelse(is.null(CONFIG$condition_col), "None", CONFIG$condition_col)))

cell_type_counts <- table(meta$labels)
condition_counts <- table(meta$condition)

cat("\n  Cell type distribution:\n")
for (ct in names(cell_type_counts)) {
  cat(sprintf("    - %s: %d cells\n", ct, cell_type_counts[ct]))
}

cat("\n  Condition distribution:\n")
for (cond in names(condition_counts)) {
  cat(sprintf("    - %s: %d cells\n", cond, condition_counts[cond]))
}

# 01: Save data summary
data_summary <- data.frame(
  Metric = c("Total Cells", "Total Genes", "Cell Types", "Conditions",
             paste0("CellType_", names(cell_type_counts)),
             paste0("Condition_", names(condition_counts))),
  Value = c(ncol(data.input), nrow(data.input), length(unique(meta$labels)),
            length(unique(meta$condition)),
            as.numeric(cell_type_counts), as.numeric(condition_counts))
)
write.csv(data_summary, file.path(DIRS$root, "01_data_summary.csv"), row.names = FALSE)
log_table("", "01_data_summary.csv", "Input Data Summary",
          "Summary statistics of the input file including cell counts per type and condition.")

# Cross-tabulation of cell types by condition
ct_by_cond <- as.data.frame.matrix(table(meta$labels, meta$condition))
ct_by_cond$CellType <- rownames(ct_by_cond)
ct_by_cond <- ct_by_cond[, c("CellType", setdiff(names(ct_by_cond), "CellType"))]
write.csv(ct_by_cond, file.path(DIRS$root, "01_celltype_by_condition.csv"), row.names = FALSE)
log_table("", "01_celltype_by_condition.csv", "Cell Type by Condition",
          "Cross-tabulation of cell counts for each cell type in each condition.")

# -----------------------------------------------------------------------------
# 3. Setup CellChatDB
# -----------------------------------------------------------------------------
cat("\nStep 3: Setting up CellChatDB...\n"); flush.console()

if (CONFIG$species == "human") {
  CellChatDB <- CellChatDB.human
} else if (CONFIG$species == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  cat(sprintf("ERROR: Unknown species '%s'. Use 'human' or 'mouse'.\n", CONFIG$species))
  quit(status = 1)
}

showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)

cat(sprintf("  - Total interactions in database: %d\n", nrow(CellChatDB.use$interaction)))
cat("  - Categories: Secreted Signaling, ECM-Receptor, Cell-Cell Contact\n"); flush.console()

# 02: Save database info (flat files in root)
db_info <- data.frame(
  Category = names(table(CellChatDB.use$interaction$annotation)),
  Count = as.numeric(table(CellChatDB.use$interaction$annotation))
)
write.csv(db_info, file.path(DIRS$root, "02_database_categories.csv"), row.names = FALSE)
log_table("", "02_database_categories.csv", "CellChatDB Categories",
          "Distribution of ligand-receptor interactions across signaling categories in CellChatDB v2.")

# Export full interaction database used
write.csv(CellChatDB.use$interaction, file.path(DIRS$root, "02_database_interactions.csv"), row.names = FALSE)
log_table("", "02_database_interactions.csv", "CellChatDB Interactions",
          "Full list of ligand-receptor interactions used from CellChatDB.")

# -----------------------------------------------------------------------------
# 4. Create CellChat Objects for Each Condition
# -----------------------------------------------------------------------------
cat("\nStep 4: Creating CellChat objects for each condition...\n"); flush.console()

conditions <- unique(meta$condition)
cellchat.list <- list()

for (cond in conditions) {
  cat(sprintf("\n  Processing condition: %s\n", cond)); flush.console()

  cells.use <- rownames(meta)[meta$condition == cond]
  data.cond <- data.input[, cells.use]
  meta.cond <- meta[cells.use, ]

  cat(sprintf("    - Cells: %d\n", length(cells.use)))

  cell_type_dist <- table(meta.cond$labels)
  cat("    - Cell type distribution:\n")
  for (ct in names(cell_type_dist)) {
    cat(sprintf("        %s: %d\n", ct, cell_type_dist[ct]))
  }
  flush.console()

  if (any(cell_type_dist < 5)) {
    cat(sprintf("    - WARNING: Some cell types have < 5 cells. Skipping %s.\n", cond)); flush.console()
    next
  }

  cat("    - Creating CellChat object...\n"); flush.console()
  cellchat <- createCellChat(object = data.cond, meta = meta.cond, group.by = "labels")
  cellchat@DB <- CellChatDB.use

  cat("    - Subsetting to signaling genes...\n"); flush.console()
  cellchat <- subsetData(cellchat)
  cat(sprintf("    - Genes after subset: %d\n", nrow(cellchat@data.signaling))); flush.console()

  cat("    - Identifying over-expressed genes...\n"); flush.console()
  cellchat <- identifyOverExpressedGenes(cellchat)
  cat("    - Identifying over-expressed interactions...\n"); flush.console()
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cat("    - Computing communication probability (this may take a while)...\n"); flush.console()
  tryCatch({
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = CONFIG$trim,
                                   raw.use = TRUE, population.size = FALSE)
    cat("    - Filtering communications...\n"); flush.console()
    cellchat <- filterCommunication(cellchat, min.cells = 5)

    cat("    - Computing pathway-level communication...\n"); flush.console()
    cellchat <- computeCommunProbPathway(cellchat)

    cat("    - Aggregating network...\n"); flush.console()
    cellchat <- aggregateNet(cellchat)

    cat("    - Computing network centrality...\n"); flush.console()
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    cellchat.list[[cond]] <- cellchat

    cat(sprintf("    - Significant L-R pairs: %d\n", nrow(cellchat@LR$LRsig)))
    cat(sprintf("    - Signaling pathways: %d\n", length(cellchat@netP$pathways))); flush.console()
  }, error = function(e) {
    cat(sprintf("    - ERROR processing %s: %s\n", cond, e$message)); flush.console()
  })
}

# -----------------------------------------------------------------------------
# 5. Generate Individual Condition Visualizations
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating visualizations for each condition...\n"); flush.console()

processed_conditions <- names(cellchat.list)
cat(sprintf("  - Successfully processed conditions: %s\n", paste(processed_conditions, collapse = ", "))); flush.console()

if (length(processed_conditions) == 0) {
  cat("  - ERROR: No conditions were successfully processed. Exiting.\n")
  quit(status = 1)
}

# Determine reference condition
reference_cond <- NULL
if (!is.null(CONFIG$reference)) {
  if (CONFIG$reference %in% processed_conditions) {
    reference_cond <- CONFIG$reference
  } else {
    cat(sprintf("  - WARNING: Specified reference '%s' not found. Using auto-detect.\n", CONFIG$reference))
  }
}
if (is.null(reference_cond)) {
  normal_idx <- grep("normal|control|healthy|baseline", processed_conditions, ignore.case = TRUE)
  reference_cond <- if (length(normal_idx) > 0) processed_conditions[normal_idx[1]] else processed_conditions[1]
}
cat(sprintf("  - Reference condition for comparisons: %s\n", reference_cond)); flush.console()

processed_conditions <- c(reference_cond, setdiff(processed_conditions, reference_cond))

# Initialize data collectors for unified tables (one table per step)
all_heatmap_data <- list()      # 03: heatmaps
all_circle_data <- list()       # 04: circle plots
all_bubble_data <- list()       # 05: bubble plots
all_signaling_roles <- list()   # 06: signaling roles
all_lr_pairs <- list()          # 07: L-R pairs
all_pathways <- list()          # 08: pathways

for (cond in processed_conditions) {
  cat(sprintf("\n  Generating plots for: %s\n", cond)); flush.console()
  cellchat <- cellchat.list[[cond]]
  cond_safe <- gsub("[^a-zA-Z0-9]", "_", cond)

  # 03: Interaction heatmaps + data export
  # Count heatmap
  png(file.path(DIRS$d03, sprintf("%s_count_heatmap.png", cond_safe)),
      width = 800, height = 700, res = 100)
  print(netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds",
                    title.name = paste0(cond, "\nCell-Cell Interaction Count\n(rows=senders, cols=receivers)")))
  dev.off()
  log_figure("03_heatmaps", sprintf("%s_count_heatmap.png", cond_safe),
             sprintf("Interaction Count Heatmap - %s", cond),
             "Heatmap showing NUMBER of inferred interactions. Rows=senders, cols=receivers.")

  # Weight heatmap
  png(file.path(DIRS$d03, sprintf("%s_weight_heatmap.png", cond_safe)),
      width = 800, height = 700, res = 100)
  print(netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds",
                    title.name = paste0(cond, "\nCell-Cell Interaction Strength\n(rows=senders, cols=receivers)")))
  dev.off()
  log_figure("03_heatmaps", sprintf("%s_weight_heatmap.png", cond_safe),
             sprintf("Interaction Strength Heatmap - %s", cond),
             "Heatmap showing STRENGTH of interactions (sum of communication probabilities).")

  # Collect heatmap data for unified table
  count_mat <- as.matrix(cellchat@net$count)
  weight_mat <- as.matrix(cellchat@net$weight)
  cell_types <- rownames(count_mat)

  for (sender in cell_types) {
    for (receiver in cell_types) {
      all_heatmap_data[[length(all_heatmap_data) + 1]] <- data.frame(
        Condition = cond,
        Sender = sender,
        Receiver = receiver,
        Count = count_mat[sender, receiver],
        Weight = weight_mat[sender, receiver],
        stringsAsFactors = FALSE
      )
    }
  }

  # 04: Circle plot + data export
  png(file.path(DIRS$d04, sprintf("%s_circle_plot.png", cond_safe)),
      width = 1100, height = 1100, res = 100)
  par(oma = c(0, 0, 4, 0))
  netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents),
                   weight.scale = TRUE, label.edge = FALSE, title.name = "")
  mtext(paste0(cond, " - Aggregated Cell-Cell Communication Network\n(edge width = interaction count, node size = cell count)"),
        side = 3, outer = TRUE, line = 1, cex = 1.2, font = 2)
  dev.off()
  log_figure("04_circle_plots", sprintf("%s_circle_plot.png", cond_safe),
             sprintf("Circle Plot - %s", cond),
             "Aggregated cell-cell communication network.")

  # Collect circle plot data for unified table
  cell_counts <- as.data.frame(table(cellchat@idents))
  names(cell_counts) <- c("CellType", "CellCount")
  cell_counts$Condition <- cond
  all_circle_data[[length(all_circle_data) + 1]] <- cell_counts

  # 05: Bubble plots by category
  MAX_PATHWAYS_PER_PAGE <- 10

  if (nrow(cellchat@LR$LRsig) > 0) {
    categories <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
    cat_abbrev <- c("Secreted Signaling" = "Secreted", "ECM-Receptor" = "ECM", "Cell-Cell Contact" = "Contact")

    for (sig_cat in categories) {
      pathways_cat <- unique(CellChatDB.use$interaction$pathway_name[
        CellChatDB.use$interaction$annotation == sig_cat
      ])
      pathways_active <- intersect(pathways_cat, cellchat@netP$pathways)

      if (length(pathways_active) > 0) {
        n_chunks <- ceiling(length(pathways_active) / MAX_PATHWAYS_PER_PAGE)
        cat_safe <- cat_abbrev[sig_cat]

        for (chunk_idx in 1:n_chunks) {
          start_idx <- (chunk_idx - 1) * MAX_PATHWAYS_PER_PAGE + 1
          end_idx <- min(chunk_idx * MAX_PATHWAYS_PER_PAGE, length(pathways_active))
          pathways_chunk <- pathways_active[start_idx:end_idx]
          chunk_label <- if (n_chunks > 1) sprintf(" (%d/%d)", chunk_idx, n_chunks) else ""
          chunk_suffix <- if (n_chunks > 1) sprintf("_part%d", chunk_idx) else ""

          fig_filename <- sprintf("%s_%s%s_bubble.png", cond_safe, cat_safe, chunk_suffix)

          tryCatch({
            png(file.path(DIRS$d05, fig_filename), width = 1400, height = 1000, res = 100)
            print(netVisual_bubble(cellchat, signaling = pathways_chunk, remove.isolate = TRUE) +
                  ggtitle(sprintf("%s - %s%s\n(%d pathways)", cond, sig_cat, chunk_label, length(pathways_chunk))) +
                  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
            dev.off()
            log_figure("05_bubble_plots", fig_filename,
                       sprintf("Bubble Plot - %s - %s%s", cond, sig_cat, chunk_label),
                       sprintf("L-R bubble plot for %s signaling. %d pathways shown.", sig_cat, length(pathways_chunk)))
          }, error = function(e) {
            cat(sprintf("      - Bubble plot failed for %s: %s\n", sig_cat, e$message))
            try(dev.off(), silent = TRUE)
          })
        }
      }
    }

    # Collect bubble plot data (actual L-R pairs per category/condition)
    tryCatch({
      lr_data <- subsetCommunication(cellchat)
      if (nrow(lr_data) > 0) {
        # Add annotation column from database to distinguish bubble plots
        lr_data$annotation <- CellChatDB.use$interaction$annotation[
          match(lr_data$interaction_name_2, CellChatDB.use$interaction$interaction_name_2)
        ]
        lr_data$Condition <- cond
        all_bubble_data[[length(all_bubble_data) + 1]] <- lr_data
      }
    }, error = function(e) {
      cat(sprintf("      - Bubble plot data collection failed: %s\n", e$message))
    })
  }

  # 06: Signaling role analysis + data export
  png(file.path(DIRS$d06, sprintf("%s_signaling_roles.png", cond_safe)),
      width = 1000, height = 800, res = 100)
  print(netAnalysis_signalingRole_scatter(cellchat,
        title = paste0(cond, " - Cell Type Signaling Roles")) +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
  dev.off()
  log_figure("06_signaling_roles", sprintf("%s_signaling_roles.png", cond_safe),
             sprintf("Signaling Roles - %s", cond),
             "Scatter plot showing dominant senders (X) and receivers (Y).")

  # Collect signaling role data for unified table
  tryCatch({
    out_strength <- rowSums(cellchat@net$weight)
    in_strength <- colSums(cellchat@net$weight)
    role_df <- data.frame(
      Condition = cond,
      CellType = names(out_strength),
      Outgoing_Strength = out_strength,
      Incoming_Strength = in_strength,
      stringsAsFactors = FALSE
    )
    all_signaling_roles[[length(all_signaling_roles) + 1]] <- role_df
  }, error = function(e) {
    cat(sprintf("      - Signaling roles data collection failed: %s\n", e$message))
  })

  # Collect L-R pairs for unified table
  if (nrow(cellchat@LR$LRsig) > 0) {
    lr_df <- subsetCommunication(cellchat)
    lr_df$Condition <- cond
    all_lr_pairs[[length(all_lr_pairs) + 1]] <- lr_df
  }

  # Collect active pathways for unified table
  pathway_df <- data.frame(
    Condition = cond,
    Pathway = cellchat@netP$pathways,
    stringsAsFactors = FALSE
  )
  all_pathways[[length(all_pathways) + 1]] <- pathway_df
}

# -----------------------------------------------------------------------------
# 5b. Export Unified Data Tables for Steps 03-08
# -----------------------------------------------------------------------------
cat("\n  Exporting unified data tables for visualization steps...\n"); flush.console()

# 03: Unified heatmap data
if (length(all_heatmap_data) > 0) {
  heatmap_df <- do.call(rbind, all_heatmap_data)
  write.csv(heatmap_df, file.path(DIRS$d03, "heatmap_data.csv"), row.names = FALSE)
  log_table("03_heatmaps", "heatmap_data.csv", "Cell-Cell Interaction Matrix Data",
            "Unified data for heatmaps. Columns: Condition, Sender, Receiver, Count (for count heatmap), Weight (for strength heatmap). Filter by Condition to get per-condition data.")
  cat(sprintf("    - 03_heatmaps: %d rows\n", nrow(heatmap_df)))
}

# 04: Unified circle plot data
if (length(all_circle_data) > 0) {
  circle_df <- do.call(rbind, all_circle_data)
  circle_df <- circle_df[, c("Condition", "CellType", "CellCount")]
  write.csv(circle_df, file.path(DIRS$d04, "circle_plot_data.csv"), row.names = FALSE)
  log_table("04_circle_plots", "circle_plot_data.csv", "Circle Plot Node Data",
            "Cell counts per cell type used for node sizes in circle plots. Columns: Condition, CellType, CellCount.")
  cat(sprintf("    - 04_circle_plots: %d rows\n", nrow(circle_df)))
}

# 05: Unified bubble plot data
if (length(all_bubble_data) > 0) {
  bubble_df <- do.call(rbind, all_bubble_data)
  # Reorder columns to put Condition and annotation first for filtering
  bubble_df <- bubble_df[, c("Condition", "annotation", setdiff(names(bubble_df), c("Condition", "annotation")))]
  write.csv(bubble_df, file.path(DIRS$d05, "bubble_plot_data.csv"), row.names = FALSE)
  log_table("05_bubble_plots", "bubble_plot_data.csv", "Bubble Plot Data",
            "L-R pairs with communication probabilities. Filter by Condition + annotation (Secreted Signaling/ECM-Receptor/Cell-Cell Contact) to get data for specific bubble plots. Columns: Condition, annotation, source, target, ligand, receptor, prob, pval, pathway_name, etc.")
  cat(sprintf("    - 05_bubble_plots: %d rows\n", nrow(bubble_df)))
}

# 06: Unified signaling roles data
if (length(all_signaling_roles) > 0) {
  roles_df <- do.call(rbind, all_signaling_roles)
  write.csv(roles_df, file.path(DIRS$d06, "signaling_roles_data.csv"), row.names = FALSE)
  log_table("06_signaling_roles", "signaling_roles_data.csv", "Signaling Roles Data",
            "Outgoing (X-axis) and Incoming (Y-axis) strength for scatter plots. Columns: Condition, CellType, Outgoing_Strength, Incoming_Strength.")
  cat(sprintf("    - 06_signaling_roles: %d rows\n", nrow(roles_df)))
}

# 07: Unified L-R pairs data
if (length(all_lr_pairs) > 0) {
  lr_all_df <- do.call(rbind, all_lr_pairs)
  # Reorder columns to put Condition first
  lr_all_df <- lr_all_df[, c("Condition", setdiff(names(lr_all_df), "Condition"))]
  write.csv(lr_all_df, file.path(DIRS$d07, "LR_pairs_data.csv"), row.names = FALSE)
  log_table("07_LR_pairs", "LR_pairs_data.csv", "Ligand-Receptor Pairs Data",
            "All significant L-R interactions across conditions. Columns: Condition, source, target, ligand, receptor, prob, pval, pathway_name, etc.")
  cat(sprintf("    - 07_LR_pairs: %d rows\n", nrow(lr_all_df)))
}

# 08: Unified pathways data
if (length(all_pathways) > 0) {
  pathways_all_df <- do.call(rbind, all_pathways)
  write.csv(pathways_all_df, file.path(DIRS$d08, "pathways_data.csv"), row.names = FALSE)
  log_table("08_pathways", "pathways_data.csv", "Active Pathways Data",
            "Active signaling pathways per condition. Columns: Condition, Pathway. Use to identify condition-specific pathways.")
  cat(sprintf("    - 08_pathways: %d rows\n", nrow(pathways_all_df)))
}

# -----------------------------------------------------------------------------
# 6. Merge and Compare Conditions (if multiple)
# -----------------------------------------------------------------------------
cat("\nStep 6: Merging CellChat objects for comparison...\n"); flush.console()

diff_stats_list <- list()

if (length(cellchat.list) < 2) {
  cat("  - Single condition analysis. Skipping merge and comparison.\n"); flush.console()
  cellchat.merged <- NULL
} else {
  cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))
  cat(sprintf("  - Merged %d conditions\n", length(cellchat.list))); flush.console()
}

# -----------------------------------------------------------------------------
# 7. Comparative Analysis Visualizations
# -----------------------------------------------------------------------------
cat("\nStep 7: Generating comparative analysis plots...\n"); flush.console()

if (is.null(cellchat.merged)) {
  cat("  - Skipping comparative analysis (single condition)\n"); flush.console()
} else {

# 09: Compare total interactions + data export
tryCatch({
  total_ints <- sapply(cellchat.list[processed_conditions], function(x) sum(x@net$count))
  total_weights <- sapply(cellchat.list[processed_conditions], function(x) sum(x@net$weight))

  ref_ints <- total_ints[1]
  ref_wts <- total_weights[1]
  pct_int <- sprintf("%+.1f%%", 100 * (total_ints - ref_ints) / ref_ints)
  pct_wt <- sprintf("%+.1f%%", 100 * (total_weights - ref_wts) / ref_wts)
  pct_int[1] <- "(ref)"
  pct_wt[1] <- "(ref)"

  subtitle_ints <- paste(sapply(seq_along(processed_conditions), function(i) {
    sprintf("%s: %d %s", processed_conditions[i], total_ints[i], pct_int[i])
  }), collapse = " | ")

  subtitle_wts <- paste(sapply(seq_along(processed_conditions), function(i) {
    sprintf("%s: %.1f %s", processed_conditions[i], total_weights[i], pct_wt[i])
  }), collapse = " | ")

  png(file.path(DIRS$d09, "compare_interactions.png"), width = 1400, height = 700, res = 100)
  gg1 <- compareInteractions(cellchat.merged, show.legend = TRUE, group = seq_along(processed_conditions)) +
         ggtitle("Total Number of Interactions Across Conditions", subtitle = subtitle_ints) +
         theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
               plot.subtitle = element_text(hjust = 0.5, size = 9))
  gg2 <- compareInteractions(cellchat.merged, show.legend = TRUE, group = seq_along(processed_conditions), measure = "weight") +
         ggtitle("Total Interaction Strength Across Conditions", subtitle = subtitle_wts) +
         theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
               plot.subtitle = element_text(hjust = 0.5, size = 9))
  print(gg1 + gg2)
  dev.off()
  log_figure("09_compare_total", "compare_interactions.png", "Comparison of Total Interactions",
             sprintf("Bar plots comparing total interactions. Reference: %s", reference_cond))

  # Export comparison data
  compare_df <- data.frame(
    Condition = processed_conditions,
    Total_Count = total_ints,
    Total_Weight = total_weights,
    Pct_Change_Count = pct_int,
    Pct_Change_Weight = pct_wt,
    Is_Reference = processed_conditions == reference_cond,
    stringsAsFactors = FALSE
  )
  write.csv(compare_df, file.path(DIRS$d09, "compare_interactions.csv"), row.names = FALSE)
  log_table("09_compare_total", "compare_interactions.csv", "Total Interactions Data",
            "Data underlying the comparison bar plots.")
}, error = function(e) {
  cat(sprintf("    - Compare interactions failed: %s\n", e$message)); flush.console()
})

# 10: Differential interaction heatmaps + data export
cat(sprintf("  - Reference condition for comparisons: %s\n", reference_cond)); flush.console()

# Collector for unified differential data
all_diff_pairs <- list()

for (i in 2:length(processed_conditions)) {
  compare_cond <- processed_conditions[i]
  ref_safe <- gsub("[^a-zA-Z0-9]", "_", reference_cond)
  cmp_safe <- gsub("[^a-zA-Z0-9]", "_", compare_cond)

  net_ref <- cellchat.merged@net[[1]]$count
  net_cmp <- cellchat.merged@net[[i]]$count
  diff_mat <- net_cmp - net_ref
  up_pairs <- sum(diff_mat > 0, na.rm = TRUE)
  dn_pairs <- sum(diff_mat < 0, na.rm = TRUE)

  wt_ref <- cellchat.merged@net[[1]]$weight
  wt_cmp <- cellchat.merged@net[[i]]$weight
  diff_wt <- wt_cmp - wt_ref
  up_wt <- sum(diff_wt > 0, na.rm = TRUE)
  dn_wt <- sum(diff_wt < 0, na.rm = TRUE)

  diff_stats_list[[compare_cond]] <- list(up = up_pairs, dn = dn_pairs, ref = reference_cond)

  title_count <- sprintf(
    "Differential Interactions: %s vs %s (Reference)\nCell-type pairs: UP=%d | DOWN=%d | Net=%+d\nRED = Increased | BLUE = Decreased",
    compare_cond, reference_cond, up_pairs, dn_pairs, up_pairs - dn_pairs
  )
  title_weight <- sprintf(
    "Differential Interaction Strength: %s vs %s (Reference)\nCell-type pairs: STRONGER=%d | WEAKER=%d | Net=%+d\nRED = Stronger | BLUE = Weaker",
    compare_cond, reference_cond, up_wt, dn_wt, up_wt - dn_wt
  )

  tryCatch({
    png(file.path(DIRS$d10, sprintf("%s_vs_%s_count.png", cmp_safe, ref_safe)),
        width = 900, height = 900, res = 100)
    par(oma = c(0, 0, 4, 0))
    netVisual_diffInteraction(cellchat.merged, comparison = c(1, i),
                              weight.scale = TRUE, measure = "count", title.name = "")
    mtext(title_count, side = 3, outer = TRUE, line = 0.5, cex = 0.9)
    dev.off()
    log_figure("10_differential", sprintf("%s_vs_%s_count.png", cmp_safe, ref_safe),
               sprintf("Differential Count: %s vs %s", compare_cond, reference_cond),
               sprintf("UP=%d, DOWN=%d cell-type pairs.", up_pairs, dn_pairs))

    png(file.path(DIRS$d10, sprintf("%s_vs_%s_weight.png", cmp_safe, ref_safe)),
        width = 900, height = 900, res = 100)
    par(oma = c(0, 0, 4, 0))
    netVisual_diffInteraction(cellchat.merged, comparison = c(1, i),
                              weight.scale = TRUE, measure = "weight", title.name = "")
    mtext(title_weight, side = 3, outer = TRUE, line = 0.5, cex = 0.9)
    dev.off()
    log_figure("10_differential", sprintf("%s_vs_%s_weight.png", cmp_safe, ref_safe),
               sprintf("Differential Weight: %s vs %s", compare_cond, reference_cond),
               sprintf("STRONGER=%d, WEAKER=%d cell-type pairs.", up_wt, dn_wt))

    cat(sprintf("    - Created differential plots: %s vs %s (UP=%d, DOWN=%d)\n", compare_cond, reference_cond, up_pairs, dn_pairs))
  }, error = function(e) {
    cat(sprintf("    - Differential plot failed: %s\n", e$message)); flush.console()
    try(dev.off(), silent = TRUE)
  })

  # Collect differential data for unified table
  tryCatch({
    cell_types <- rownames(diff_mat)
    diff_pairs_df <- data.frame(
      Comparison_Pair = sprintf("%s_vs_%s", compare_cond, reference_cond),
      Reference = reference_cond,
      Test_Condition = compare_cond,
      Sender = rep(cell_types, times = length(cell_types)),
      Receiver = rep(cell_types, each = length(cell_types)),
      Ref_Count = as.vector(net_ref),
      Test_Count = as.vector(net_cmp),
      Diff_Count = as.vector(diff_mat),
      Ref_Weight = as.vector(wt_ref),
      Test_Weight = as.vector(wt_cmp),
      Diff_Weight = as.vector(diff_wt),
      stringsAsFactors = FALSE
    )
    diff_pairs_df$Count_Direction <- ifelse(diff_pairs_df$Diff_Count > 0, "UP",
                                             ifelse(diff_pairs_df$Diff_Count < 0, "DOWN", "NC"))
    diff_pairs_df$Weight_Direction <- ifelse(diff_pairs_df$Diff_Weight > 0, "STRONGER",
                                              ifelse(diff_pairs_df$Diff_Weight < 0, "WEAKER", "NC"))
    diff_pairs_df <- diff_pairs_df[diff_pairs_df$Ref_Count > 0 | diff_pairs_df$Test_Count > 0, ]
    all_diff_pairs[[length(all_diff_pairs) + 1]] <- diff_pairs_df
    cat(sprintf("    - Collected differential data: %s vs %s (%d pairs)\n", compare_cond, reference_cond, nrow(diff_pairs_df)))
  }, error = function(e) {
    cat(sprintf("    - Differential data collection failed: %s\n", e$message))
  })
}

# Export unified differential data table
if (length(all_diff_pairs) > 0) {
  all_diff_df <- do.call(rbind, all_diff_pairs)
  all_diff_df <- all_diff_df[order(all_diff_df$Comparison_Pair, -abs(all_diff_df$Diff_Count)), ]
  write.csv(all_diff_df, file.path(DIRS$d10, "differential_data.csv"), row.names = FALSE)
  log_table("10_differential", "differential_data.csv", "Differential Interaction Data",
            "Unified differential analysis data. Columns: Comparison_Pair (Test_vs_Reference), Sender, Receiver, Ref/Test Count/Weight, Diff values, Direction (UP/DOWN for count, STRONGER/WEAKER for weight). Filter by Comparison_Pair for specific comparisons.")
  cat(sprintf("    - 10_differential: %d total differential pairs\n", nrow(all_diff_df)))
}

# 11: Information flow comparison + data export
tryCatch({
  png(file.path(DIRS$d11, "information_flow.png"), width = 1600, height = 1000, res = 100)
  gg1 <- rankNet(cellchat.merged, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE) +
         ggtitle("Signaling Pathway Information Flow (Stacked)") +
         theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  gg2 <- rankNet(cellchat.merged, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE) +
         ggtitle("Signaling Pathway Information Flow (Grouped)") +
         theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  print(gg1 + gg2)
  dev.off()
  log_figure("11_info_flow", "information_flow.png", "Information Flow Comparison",
             "Signaling pathway information flow across conditions.")
}, error = function(e) {
  cat(sprintf("    - Information flow plot failed: %s\n", e$message)); flush.console()
})

# 11b: Export information flow data
tryCatch({
  info_flow_data <- data.frame()
  for (cond in processed_conditions) {
    cellchat_cond <- cellchat.list[[cond]]
    for (pw in cellchat_cond@netP$pathways) {
      idx <- which(dimnames(cellchat_cond@netP$prob)[[3]] == pw)
      flow <- if (length(idx) > 0) sum(cellchat_cond@netP$prob[,,idx], na.rm = TRUE) else 0
      info_flow_data <- rbind(info_flow_data, data.frame(
        Condition = cond,
        Pathway = pw,
        Info_Flow = flow,
        Pathway_Rank = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  # Add pathway rank based on total info flow
  if (nrow(info_flow_data) > 0) {
    pathway_totals <- aggregate(Info_Flow ~ Pathway, data = info_flow_data, FUN = sum)
    pathway_totals <- pathway_totals[order(-pathway_totals$Info_Flow), ]
    pathway_totals$Pathway_Rank <- seq_len(nrow(pathway_totals))
    info_flow_data$Pathway_Rank <- pathway_totals$Pathway_Rank[
      match(info_flow_data$Pathway, pathway_totals$Pathway)
    ]
    info_flow_data <- info_flow_data[order(info_flow_data$Pathway_Rank, info_flow_data$Condition), ]

    write.csv(info_flow_data, file.path(DIRS$d11, "info_flow_data.csv"), row.names = FALSE)
    log_table("11_info_flow", "info_flow_data.csv", "Information Flow Data",
              "Pathway information flow values. Columns: Condition, Pathway, Info_Flow, Pathway_Rank. Stacked plot shows contribution by condition; grouped shows direct comparison.")
    cat(sprintf("    - 11_info_flow: %d rows\n", nrow(info_flow_data)))
  }
}, error = function(e) {
  cat(sprintf("    - Info flow data export failed: %s\n", e$message))
})

# 12: Pathway-specific networks (individual PNG per pathway)
all_pathways <- unique(unlist(lapply(cellchat.list, function(x) x@netP$pathways)))
cat(sprintf("  - Total unique pathways: %d\n", length(all_pathways))); flush.console()

if (length(all_pathways) > 0) {
  for (cond in processed_conditions) {
    cond_safe <- gsub("[^a-zA-Z0-9]", "_", cond)
    pathways_cond <- cellchat.list[[cond]]@netP$pathways

    if (length(pathways_cond) > 0) {
      cat(sprintf("    - %s: %d pathways\n", cond, length(pathways_cond))); flush.console()

      failed_pathways <- character()
      for (pathway in pathways_cond) {
        pathway_safe <- gsub("[^a-zA-Z0-9]", "_", pathway)
        fig_filename <- sprintf("%s_%s_network.png", cond_safe, pathway_safe)

        tryCatch({
          png(file.path(DIRS$d12, fig_filename), width = 1000, height = 1000, res = 100)
          par(oma = c(0, 0, 3, 0))
          netVisual_aggregate(cellchat.list[[cond]], signaling = pathway, layout = "circle")
          mtext(sprintf("%s Signaling Network - %s", pathway, cond), side = 3, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
          dev.off()
          log_figure("12_pathway_networks", fig_filename,
                     sprintf("%s Network - %s", pathway, cond),
                     sprintf("Circle plot for %s signaling in %s.", pathway, cond))
        }, error = function(e) {
          failed_pathways <<- c(failed_pathways, pathway)
          try(dev.off(), silent = TRUE)
        })
      }
      if (length(failed_pathways) > 0) {
        cat(sprintf("      - Failed pathways: %s\n", paste(failed_pathways, collapse = ", ")))
      }
    }
  }
}

# 12b: Export pathway network data (edge data for each pathway network plot)
tryCatch({
  pathway_network_data <- data.frame()
  for (cond in processed_conditions) {
    cond_safe <- gsub("[^a-zA-Z0-9]", "_", cond)
    cellchat_cond <- cellchat.list[[cond]]
    pathways_cond <- cellchat_cond@netP$pathways

    for (pathway in pathways_cond) {
      pathway_safe <- gsub("[^a-zA-Z0-9]", "_", pathway)

      # Get pathway-specific communication probability matrix
      idx <- which(dimnames(cellchat_cond@netP$prob)[[3]] == pathway)
      if (length(idx) > 0) {
        prob_mat <- cellchat_cond@netP$prob[,,idx]
        cell_types <- rownames(prob_mat)

        # Convert matrix to edge list
        for (sender in cell_types) {
          for (receiver in cell_types) {
            prob_val <- prob_mat[sender, receiver]
            if (!is.na(prob_val) && prob_val > 0) {
              pathway_network_data <- rbind(pathway_network_data, data.frame(
                Condition = cond,
                Pathway = pathway,
                Sender = sender,
                Receiver = receiver,
                Prob = prob_val,
                Filename = sprintf("%s_%s_network.png", cond_safe, pathway_safe),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }
  if (nrow(pathway_network_data) > 0) {
    write.csv(pathway_network_data, file.path(DIRS$d12, "pathway_networks_data.csv"), row.names = FALSE)
    log_table("12_pathway_networks", "pathway_networks_data.csv", "Pathway Network Data",
              "Edge data for pathway network plots. Columns: Condition, Pathway, Sender, Receiver, Prob, Filename. Filter by Condition + Pathway for specific plot data.")
    cat(sprintf("    - 12_pathway_networks: %d rows\n", nrow(pathway_network_data)))
  }
}, error = function(e) {
  cat(sprintf("    - Pathway network data export failed: %s\n", e$message))
})

# 13: Functional similarity + data export
cat("  - Computing functional similarity...\n"); flush.console()
functional_cluster_df <- NULL
tryCatch({
  cellchat.merged <- computeNetSimilarityPairwise(cellchat.merged, type = "functional")
  cellchat.merged <- netEmbedding(cellchat.merged, type = "functional")
  cellchat.merged <- netClustering(cellchat.merged, type = "functional")

  # Estimate optimal number of clusters for functional similarity
  tryCatch({
    pdf(file.path(DIRS$d13, "estimateNumCluster_functional.pdf"), width = 7, height = 6)
    estimateNumCluster(cellchat.merged, type = "functional", do.plot = TRUE)
    dev.off()
    log_figure("13_functional_sim", "estimateNumCluster_functional.pdf",
               "Cluster Number Estimation - Functional",
               "Silhouette/gap statistic plots for optimal functional cluster count.")
    cat("    - Generated estimateNumCluster plot for functional similarity\n")
  }, error = function(e) {
    cat(sprintf("    - estimateNumCluster (functional) failed: %s\n", e$message))
    try(dev.off(), silent = TRUE)
  })

  png(file.path(DIRS$d13, "functional_similarity.png"), width = 1100, height = 900, res = 100)
  print(netVisual_embeddingPairwise(cellchat.merged, type = "functional", label.size = 3.5) +
        ggtitle("Functional Similarity of Signaling Pathways") +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
  dev.off()
  log_figure("13_functional_sim", "functional_similarity.png", "Functional Similarity Embedding",
             "UMAP embedding based on functional similarity.")

  tryCatch({
    # CellChat v2 stores group/dr as nested lists - access first element
    func_sim <- cellchat.merged@netP$similarity$functional

    # Get cluster assignments from $group (not $cluster)
    pathway_clusters <- NULL
    if (!is.null(func_sim$group) && length(func_sim$group) > 0) {
      pathway_clusters <- func_sim$group[[1]]  # First element of list
    }

    if (is.null(pathway_clusters) || length(pathway_clusters) == 0) {
      stop("No cluster data found in similarity$functional$group")
    }

    pathway_names <- names(pathway_clusters)

    # Get UMAP coordinates from $dr (nested list)
    umap_coords <- NULL
    if (!is.null(func_sim$dr) && length(func_sim$dr) > 0) {
      umap_coords <- func_sim$dr[[1]]  # First element of list
    }

    if (!is.null(umap_coords) && is.matrix(umap_coords) && ncol(umap_coords) >= 2) {
      functional_cluster_df <- data.frame(
        Pathway_Condition = pathway_names,
        Cluster = as.numeric(pathway_clusters),
        UMAP1 = umap_coords[, 1],
        UMAP2 = umap_coords[, 2],
        stringsAsFactors = FALSE
      )
    } else {
      # Export without UMAP coordinates
      functional_cluster_df <- data.frame(
        Pathway_Condition = pathway_names,
        Cluster = as.numeric(pathway_clusters),
        stringsAsFactors = FALSE
      )
    }

    # Parse pathway and condition from combined names (e.g., "COLLAGEN--CD_InterFibrosis")
    split_names <- strsplit(pathway_names, "--")
    functional_cluster_df$Pathway <- sapply(split_names, function(x) x[1])
    functional_cluster_df$Condition <- sapply(split_names, function(x) if(length(x) > 1) x[2] else NA)

    # Reorder columns
    col_order <- c("Pathway", "Condition", "Cluster")
    if ("UMAP1" %in% names(functional_cluster_df)) col_order <- c(col_order, "UMAP1", "UMAP2")
    col_order <- c(col_order, "Pathway_Condition")
    functional_cluster_df <- functional_cluster_df[, col_order]

    functional_cluster_df <- functional_cluster_df[order(functional_cluster_df$Cluster, functional_cluster_df$Pathway), ]
    write.csv(functional_cluster_df, file.path(DIRS$d13, "functional_clusters.csv"), row.names = FALSE)
    log_table("13_functional_sim", "functional_clusters.csv", "Functional Similarity Clusters",
              "Pathway-condition cluster assignments based on functional similarity. Pathway_Condition is the combined identifier used in UMAP plots.")
    cat(sprintf("    - Exported functional similarity data: %d entries, %d clusters\n",
                nrow(functional_cluster_df), length(unique(functional_cluster_df$Cluster))))
  }, error = function(e) {
    cat(sprintf("    - Functional similarity CSV export skipped: %s\n", e$message))
  })
}, error = function(e) {
  cat(sprintf("    - Functional similarity analysis skipped: %s\n", e$message)); flush.console()
  try(dev.off(), silent = TRUE)
})

# 14: Structural similarity + data export
cat("  - Computing structural similarity...\n"); flush.console()
structural_cluster_df <- NULL
tryCatch({
  cellchat.merged <- computeNetSimilarityPairwise(cellchat.merged, type = "structural")
  cellchat.merged <- netEmbedding(cellchat.merged, type = "structural")
  cellchat.merged <- netClustering(cellchat.merged, type = "structural")

  # Estimate optimal number of clusters for structural similarity
  tryCatch({
    pdf(file.path(DIRS$d14, "estimateNumCluster_structural.pdf"), width = 7, height = 6)
    estimateNumCluster(cellchat.merged, type = "structural", do.plot = TRUE)
    dev.off()
    log_figure("14_structural_sim", "estimateNumCluster_structural.pdf",
               "Cluster Number Estimation - Structural",
               "Silhouette/gap statistic plots for optimal structural cluster count.")
    cat("    - Generated estimateNumCluster plot for structural similarity\n")
  }, error = function(e) {
    cat(sprintf("    - estimateNumCluster (structural) failed: %s\n", e$message))
    try(dev.off(), silent = TRUE)
  })

  png(file.path(DIRS$d14, "structural_similarity.png"), width = 1100, height = 900, res = 100)
  print(netVisual_embeddingPairwise(cellchat.merged, type = "structural", label.size = 3.5) +
        ggtitle("Structural Similarity of Signaling Pathways") +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")))
  dev.off()
  log_figure("14_structural_sim", "structural_similarity.png", "Structural Similarity Embedding",
             "UMAP embedding based on structural similarity.")

  tryCatch({
    # CellChat v2 stores group/dr as nested lists - access first element
    struct_sim <- cellchat.merged@netP$similarity$structural

    # Get cluster assignments from $group (not $cluster)
    pathway_clusters <- NULL
    if (!is.null(struct_sim$group) && length(struct_sim$group) > 0) {
      pathway_clusters <- struct_sim$group[[1]]  # First element of list
    }

    if (is.null(pathway_clusters) || length(pathway_clusters) == 0) {
      stop("No cluster data found in similarity$structural$group")
    }

    pathway_names <- names(pathway_clusters)

    # Get UMAP coordinates from $dr (nested list)
    umap_coords <- NULL
    if (!is.null(struct_sim$dr) && length(struct_sim$dr) > 0) {
      umap_coords <- struct_sim$dr[[1]]  # First element of list
    }

    if (!is.null(umap_coords) && is.matrix(umap_coords) && ncol(umap_coords) >= 2) {
      structural_cluster_df <- data.frame(
        Pathway_Condition = pathway_names,
        Cluster = as.numeric(pathway_clusters),
        UMAP1 = umap_coords[, 1],
        UMAP2 = umap_coords[, 2],
        stringsAsFactors = FALSE
      )
    } else {
      # Export without UMAP coordinates
      structural_cluster_df <- data.frame(
        Pathway_Condition = pathway_names,
        Cluster = as.numeric(pathway_clusters),
        stringsAsFactors = FALSE
      )
    }

    # Parse pathway and condition from combined names (e.g., "COLLAGEN--CD_InterFibrosis")
    split_names <- strsplit(pathway_names, "--")
    structural_cluster_df$Pathway <- sapply(split_names, function(x) x[1])
    structural_cluster_df$Condition <- sapply(split_names, function(x) if(length(x) > 1) x[2] else NA)

    # Reorder columns
    col_order <- c("Pathway", "Condition", "Cluster")
    if ("UMAP1" %in% names(structural_cluster_df)) col_order <- c(col_order, "UMAP1", "UMAP2")
    col_order <- c(col_order, "Pathway_Condition")
    structural_cluster_df <- structural_cluster_df[, col_order]

    structural_cluster_df <- structural_cluster_df[order(structural_cluster_df$Cluster, structural_cluster_df$Pathway), ]
    write.csv(structural_cluster_df, file.path(DIRS$d14, "structural_clusters.csv"), row.names = FALSE)
    log_table("14_structural_sim", "structural_clusters.csv", "Structural Similarity Clusters",
              "Pathway-condition cluster assignments based on structural similarity. Pathway_Condition is the combined identifier used in UMAP plots.")
    cat(sprintf("    - Exported structural similarity data: %d entries, %d clusters\n",
                nrow(structural_cluster_df), length(unique(structural_cluster_df$Cluster))))
  }, error = function(e) {
    cat(sprintf("    - Structural similarity CSV export skipped: %s\n", e$message))
  })
}, error = function(e) {
  cat(sprintf("    - Structural similarity analysis skipped: %s\n", e$message)); flush.console()
  try(dev.off(), silent = TRUE)
})

}  # End of if (!is.null(cellchat.merged))

# -----------------------------------------------------------------------------
# 8. Save CellChat Objects
# -----------------------------------------------------------------------------
cat("\nStep 8: Saving CellChat objects...\n")

saveRDS(cellchat.list, file.path(DIRS$objects, "cellchat_list.rds"))
saveRDS(cellchat.merged, file.path(DIRS$objects, "cellchat_merged.rds"))
cat("  - Saved cellchat_list.rds and cellchat_merged.rds\n")

# -----------------------------------------------------------------------------
# 9. Generate Summary Statistics (Steps 15-18)
# -----------------------------------------------------------------------------
cat("\nStep 9: Generating summary statistics...\n"); flush.console()

# 15: Interaction summary
interaction_summary <- data.frame(
  Condition = processed_conditions,
  Total_Interactions = sapply(cellchat.list[processed_conditions], function(x) sum(x@net$count)),
  Total_Strength = sapply(cellchat.list[processed_conditions], function(x) sum(x@net$weight)),
  Num_LR_Pairs = sapply(cellchat.list[processed_conditions], function(x) nrow(x@LR$LRsig)),
  Num_Pathways = sapply(cellchat.list[processed_conditions], function(x) length(x@netP$pathways)),
  stringsAsFactors = FALSE
)
write.csv(interaction_summary, file.path(DIRS$d15, "interaction_summary.csv"), row.names = FALSE)
log_table("15_summary", "interaction_summary.csv", "Interaction Summary",
          "Summary statistics of cell-cell communication for each condition.")

# 16: Pathway presence matrix
all_pathways <- unique(unlist(lapply(cellchat.list, function(x) x@netP$pathways)))
if (length(all_pathways) > 0) {
  pathway_matrix <- matrix(0, nrow = length(all_pathways), ncol = length(processed_conditions))
  rownames(pathway_matrix) <- all_pathways
  colnames(pathway_matrix) <- processed_conditions
  for (cond in processed_conditions) {
    pathways_cond <- cellchat.list[[cond]]@netP$pathways
    pathway_matrix[pathways_cond, cond] <- 1
  }
  pathway_df <- as.data.frame(pathway_matrix)
  pathway_df$Pathway <- rownames(pathway_df)
  pathway_df <- pathway_df[, c("Pathway", processed_conditions)]
  write.csv(pathway_df, file.path(DIRS$d16, "pathway_matrix.csv"), row.names = FALSE)
  log_table("16_pathway_matrix", "pathway_matrix.csv", "Pathway Presence Matrix",
            "Binary matrix showing which pathways are active in each condition.")
}

# 17: Pathway information flow
if (!is.null(cellchat.merged) && length(processed_conditions) >= 2) {
  pathway_info_flow <- data.frame(Pathway = all_pathways, stringsAsFactors = FALSE)

  for (cond in processed_conditions) {
    cellchat_cond <- cellchat.list[[cond]]
    pathway_probs <- sapply(all_pathways, function(pw) {
      if (pw %in% cellchat_cond@netP$pathways) {
        idx <- which(dimnames(cellchat_cond@netP$prob)[[3]] == pw)
        if (length(idx) > 0) return(sum(cellchat_cond@netP$prob[,,idx], na.rm = TRUE))
      }
      return(0)
    })
    pathway_info_flow[[cond]] <- pathway_probs
  }

  write.csv(pathway_info_flow, file.path(DIRS$d17, "pathway_info_flow.csv"), row.names = FALSE)
  log_table("17_pathway_flow", "pathway_info_flow.csv", "Pathway Information Flow",
            "Sum of communication probabilities for each pathway in each condition.")

  # 18: Differential pathway statistics + figures
  cat("\nStep 9b: Computing pathway statistics (log2FC, p-values)...\n"); flush.console()

  reference_cond <- processed_conditions[1]
  ref_values <- pathway_info_flow[[reference_cond]]
  pseudo <- 1e-6

  # Collector for unified pathway stats
  all_pathway_stats <- list()

  for (i in 2:length(processed_conditions)) {
    compare_cond <- processed_conditions[i]
    test_values <- pathway_info_flow[[compare_cond]]

    log2fc <- log2((test_values + pseudo) / (ref_values + pseudo))
    abs_log2fc <- abs(log2fc)
    n_pathways <- length(all_pathways)
    pval <- 1 - rank(abs_log2fc) / (n_pathways + 1)
    padj <- p.adjust(pval, method = "BH")
    direction <- ifelse(log2fc > 0, "UP", ifelse(log2fc < 0, "DOWN", "NC"))

    comp_stats <- data.frame(
      Comparison_Pair = sprintf("%s_vs_%s", compare_cond, reference_cond),
      Pathway = all_pathways,
      Reference = reference_cond,
      Test_Condition = compare_cond,
      Ref_InfoFlow = ref_values,
      Test_InfoFlow = test_values,
      Log2FC = round(log2fc, 4),
      Direction = direction,
      PValue = round(pval, 6),
      FDR = round(padj, 6),
      Significant = padj < 0.05,
      stringsAsFactors = FALSE
    )

    comp_stats <- comp_stats[order(-abs(comp_stats$Log2FC)), ]
    all_pathway_stats[[length(all_pathway_stats) + 1]] <- comp_stats

    cmp_safe <- gsub("[^a-zA-Z0-9]", "_", compare_cond)
    ref_safe <- gsub("[^a-zA-Z0-9]", "_", reference_cond)

    # 18b: Pathway log2FC figure
    tryCatch({
      n_up <- sum(comp_stats$Direction == "UP" & comp_stats$Significant, na.rm = TRUE)
      n_dn <- sum(comp_stats$Direction == "DOWN" & comp_stats$Significant, na.rm = TRUE)

      stats_top <- head(comp_stats, 30)
      stats_top$Pathway <- factor(stats_top$Pathway, levels = stats_top$Pathway[order(stats_top$Log2FC)])
      stats_top$Significant_Label <- ifelse(stats_top$Significant, "*", "")

      subtitle_stats <- sprintf("Significant (FDR<0.05): UP=%d | DOWN=%d | Top 30 by |log2FC|", n_up, n_dn)

      png(file.path(DIRS$d18, sprintf("%s_vs_%s_pathway_logfc.png", cmp_safe, ref_safe)),
          width = 1200, height = 900, res = 100)
      p <- ggplot(stats_top, aes(x = Log2FC, y = Pathway, fill = Direction)) +
           geom_bar(stat = "identity") +
           geom_text(aes(label = Significant_Label, x = ifelse(Log2FC >= 0, Log2FC + 0.1, Log2FC - 0.1)),
                     hjust = ifelse(stats_top$Log2FC >= 0, 0, 1), size = 5) +
           scale_fill_manual(values = c("UP" = "firebrick", "DOWN" = "steelblue", "NC" = "grey50")) +
           geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
           labs(title = sprintf("Pathway Log2FC: %s vs %s (Reference)", compare_cond, reference_cond),
                subtitle = subtitle_stats,
                x = "Log2 Fold Change (Information Flow)", y = "Signaling Pathway") +
           theme_bw() +
           theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                 plot.subtitle = element_text(hjust = 0.5, size = 10),
                 axis.text.y = element_text(size = 9),
                 legend.position = "bottom")
      print(p)
      dev.off()

      log_figure("18_pathway_stats", sprintf("%s_vs_%s_pathway_logfc.png", cmp_safe, ref_safe),
                 sprintf("Pathway Log2FC: %s vs %s", compare_cond, reference_cond),
                 sprintf("Bar plot of pathway log2FC. * = FDR<0.05. UP=%d, DOWN=%d.", n_up, n_dn))
      cat(sprintf("    - Generated pathway log2FC figure: %s vs %s (UP=%d, DOWN=%d significant)\n",
                  compare_cond, reference_cond, n_up, n_dn))
    }, error = function(e) {
      cat(sprintf("    - Pathway log2FC plot failed: %s\n", e$message))
    })
  }

  # Export unified pathway stats table
  if (length(all_pathway_stats) > 0) {
    all_stats_df <- do.call(rbind, all_pathway_stats)
    write.csv(all_stats_df, file.path(DIRS$d18, "pathway_stats_data.csv"), row.names = FALSE)
    log_table("18_pathway_stats", "pathway_stats_data.csv", "Pathway Differential Statistics",
              "Unified pathway statistics for all comparisons. Columns: Comparison_Pair, Pathway, Reference, Test_Condition, InfoFlow values, Log2FC, Direction (UP/DOWN), PValue, FDR, Significant. Filter by Comparison_Pair for specific comparisons.")
    cat(sprintf("    - 18_pathway_stats: %d total pathway stats\n", nrow(all_stats_df)))
  }
}

# -----------------------------------------------------------------------------
# 10. Generate Comprehensive Markdown Report
# -----------------------------------------------------------------------------
cat("\nStep 10: Generating comprehensive markdown report...\n")

md_table <- function(df, align = NULL) {
  if (nrow(df) == 0) return(character())
  header <- paste("|", paste(names(df), collapse = " | "), "|")
  if (is.null(align)) align <- rep("l", ncol(df))
  align_row <- paste("|", paste(sapply(align, function(a) {
    switch(a, "l" = ":---", "r" = "---:", "c" = ":---:", "---")
  }), collapse = " | "), "|")
  rows <- apply(df, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  c(header, align_row, rows)
}

report <- c(
  "# CellChat Cell-Cell Communication Analysis Report",
  "",
  sprintf("**Generated:** %s", Sys.time()),
  sprintf("**CellChat Version:** %s", packageVersion("CellChat")),
  "",
  "---",
  "",
  "## Executive Summary",
  "",
  sprintf("This analysis inferred cell-cell communication networks from **%d cells** across **%d cell types** in **%d conditions** using CellChat v%s.",
          ncol(data.input), length(unique(meta$labels)), length(processed_conditions), packageVersion("CellChat")),
  ""
)

report <- c(report, "### Key Findings", "")
for (cond in processed_conditions) {
  cc <- cellchat.list[[cond]]
  report <- c(report, sprintf("- **%s**: %d significant L-R pairs, %d active signaling pathways",
                               cond, nrow(cc@LR$LRsig), length(cc@netP$pathways)))
}

if (length(processed_conditions) > 1) {
  report <- c(report, "", "### Comparison Pairs", "")
  for (cond in processed_conditions[-1]) {
    report <- c(report, sprintf("- %s vs %s (reference)", cond, reference_cond))
  }
}

report <- c(report,
  "", "---", "",
  "## Output Directory Structure",
  "",
  "Outputs are organized as follows (01-02 are flat files in root; 03-18 are subfolders):",
  "",
  "### Root Directory Files",
  "",
  "| File | Description | Usage |",
  "|------|-------------|-------|",
  "| `01_data_summary.csv` | Input data statistics (cells, genes, cell types, conditions) | Verify input data quality |",
  "| `01_celltype_by_condition.csv` | Cross-tabulation of cell types by condition | Check cell distribution balance |",
  "| `02_database_categories.csv` | CellChatDB signaling categories and interaction counts | Understand database composition |",
  "| `02_database_interactions.csv` | Full L-R interaction database used | Reference for custom queries |",
  "",
  "### Analysis Subfolders",
  "",
  "| Folder | Plots | Data Table | Purpose |",
  "|--------|-------|------------|---------|",
  "| `03_heatmaps/` | `{cond}_count_heatmap.png`, `{cond}_weight_heatmap.png` | `heatmap_data.csv` | Visualize interaction patterns between cell types. Count = number of interactions; Weight = strength |",
  "| `04_circle_plots/` | `{cond}_circle_plot.png` | `circle_plot_data.csv` | Network visualization with node size = cell count, edge width = interactions |",
  "| `05_bubble_plots/` | `{cond}_{category}_bubble.png` | (data in 07_LR_pairs) | L-R pair bubble plots by signaling category (Secreted/ECM/Contact) |",
  "| `06_signaling_roles/` | `{cond}_signaling_roles.png` | `signaling_roles_data.csv` | Identify dominant senders (X-axis) and receivers (Y-axis) |",
  "| `07_LR_pairs/` | - | `LR_pairs_data.csv` | Complete list of significant L-R interactions with probabilities |",
  "| `08_pathways/` | - | `pathways_data.csv` | Active signaling pathways per condition |",
  "| `09_compare_total/` | `compare_interactions.png` | `compare_interactions.csv` | Compare total interaction counts/strength across conditions |",
  "| `10_differential/` | `{cmp}_vs_{ref}_count.png`, `{cmp}_vs_{ref}_weight.png` | `differential_data.csv` | Differential interactions: RED=increased, BLUE=decreased |",
  "| `11_info_flow/` | `information_flow.png` | (data in 17_pathway_flow) | Pathway activity ranking across conditions |",
  "| `12_pathway_networks/` | `{cond}_{pathway}_network.png` | (data in 07_LR_pairs) | Individual pathway circle plots |",
  "| `13_functional_sim/` | `functional_similarity.png` | `functional_clusters.csv` | UMAP of pathways by sender-receiver patterns |",
  "| `14_structural_sim/` | `structural_similarity.png` | `structural_clusters.csv` | UMAP of pathways by network topology |",
  "| `15_summary/` | - | `interaction_summary.csv` | Per-condition summary statistics |",
  "| `16_pathway_matrix/` | - | `pathway_matrix.csv` | Binary presence/absence matrix of pathways |",
  "| `17_pathway_flow/` | - | `pathway_info_flow.csv` | Communication probability sum per pathway |",
  "| `18_pathway_stats/` | `{cmp}_vs_{ref}_pathway_logfc.png` | `pathway_stats_data.csv` | Differential pathway log2FC, p-values, FDR |",
  "| `objects/` | - | `cellchat_list.rds`, `cellchat_merged.rds` | R objects for downstream analysis |",
  ""
)

report <- c(report,
  "---", "",
  "## Analysis Parameters",
  "",
  "| Parameter | Value | Description |",
  "|-----------|-------|-------------|",
  sprintf("| Cell type column | `%s` | Metadata column for cell type labels |", CONFIG$cell_type_col),
  sprintf("| Condition column | `%s` | Metadata column for condition grouping |",
          ifelse(is.null(CONFIG$condition_col), "None", CONFIG$condition_col)),
  sprintf("| Reference condition | %s | Baseline for differential comparisons |", reference_cond),
  sprintf("| Species | %s | CellChatDB species |", CONFIG$species),
  sprintf("| Min cells | %d | Minimum cells per cell type for filtering |", CONFIG$min_cells),
  sprintf("| Trim value | %.2f | Expression cutoff for truncatedMean |", CONFIG$trim),
  "| Database | CellChatDB v2 | Excludes non-protein signaling |",
  ""
)

report <- c(report,
  "---", "",
  "## Per-Condition Results",
  ""
)

results_df <- data.frame(
  Condition = processed_conditions,
  Cells = sapply(processed_conditions, function(c) sum(meta$condition == c)),
  `L-R Pairs` = sapply(cellchat.list[processed_conditions], function(cc) nrow(cc@LR$LRsig)),
  Pathways = sapply(cellchat.list[processed_conditions], function(cc) length(cc@netP$pathways)),
  `Total Interactions` = sapply(cellchat.list[processed_conditions], function(cc) sum(cc@net$count)),
  `Total Strength` = sapply(cellchat.list[processed_conditions], function(cc) round(sum(cc@net$weight), 2)),
  check.names = FALSE
)
report <- c(report, md_table(results_df, c("l", "r", "r", "r", "r", "r")))

if (length(processed_conditions) > 1 && length(diff_stats_list) > 0) {
  report <- c(report,
    "", "---", "",
    "## Differential Analysis",
    "",
    "| Comparison | UP | DOWN | Net Change |",
    "|------------|---:|-----:|-----------:|"
  )
  for (cond in names(diff_stats_list)) {
    s <- diff_stats_list[[cond]]
    net_change <- s$up - s$dn
    report <- c(report, sprintf("| %s vs %s | %d | %d | %+d |",
                                 cond, s$ref, s$up, s$dn, net_change))
  }
  report <- c(report,
    "",
    "**Interpretation:** UP = more cell-type pair interactions in test condition; DOWN = fewer.",
    ""
  )
}

report <- c(report,
  "---", "",
  "## Detailed Output Documentation",
  "",
  "### Data Tables - Column Descriptions",
  "",
  "**heatmap_data.csv (03_heatmaps/)**",
  "- `Condition`: Analysis condition name",
  "- `Sender`: Cell type sending the signal (rows in heatmap)",
  "- `Receiver`: Cell type receiving the signal (columns in heatmap)",
  "- `Count`: Number of inferred interactions (used for count_heatmap.png)",
  "- `Weight`: Sum of communication probabilities (used for weight_heatmap.png)",
  "",
  "**circle_plot_data.csv (04_circle_plots/)**",
  "- `Condition`: Analysis condition name",
  "- `CellType`: Cell type label",
  "- `CellCount`: Number of cells (determines node size in circle plots)",
  "",
  "**signaling_roles_data.csv (06_signaling_roles/)**",
  "- `Condition`: Analysis condition name",
  "- `CellType`: Cell type label",
  "- `Outgoing_Strength`: Sum of outgoing signal weights (X-axis in scatter plot)",
  "- `Incoming_Strength`: Sum of incoming signal weights (Y-axis in scatter plot)",
  "",
  "**LR_pairs_data.csv (07_LR_pairs/)**",
  "- `Condition`: Analysis condition name",
  "- `source`: Sending cell type",
  "- `target`: Receiving cell type",
  "- `ligand`: Ligand gene name(s)",
  "- `receptor`: Receptor gene name(s)",
  "- `prob`: Communication probability",
  "- `pval`: P-value from permutation test",
  "- `pathway_name`: Signaling pathway name",
  "- `annotation`: Signaling category (Secreted/ECM-Receptor/Cell-Cell Contact)",
  "",
  "**pathways_data.csv (08_pathways/)**",
  "- `Condition`: Analysis condition name",
  "- `Pathway`: Active signaling pathway name",
  "",
  "**compare_interactions.csv (09_compare_total/)**",
  "- `Condition`: Analysis condition name",
  "- `Total_Count`: Total number of cell-type pair interactions",
  "- `Total_Weight`: Sum of all communication probabilities",
  "- `Pct_Change_Count`: Percentage change in count vs reference",
  "- `Pct_Change_Weight`: Percentage change in weight vs reference",
  "- `Is_Reference`: TRUE if this is the reference condition",
  "",
  "**differential_data.csv (10_differential/)**",
  "- `Comparison_Pair`: Test_vs_Reference comparison identifier",
  "- `Reference`: Reference condition name",
  "- `Test_Condition`: Test condition being compared",
  "- `Sender`: Sending cell type",
  "- `Receiver`: Receiving cell type",
  "- `Ref_Count/Test_Count`: Interaction counts in each condition",
  "- `Diff_Count`: Test - Reference count difference",
  "- `Ref_Weight/Test_Weight`: Interaction strengths in each condition",
  "- `Diff_Weight`: Test - Reference weight difference",
  "- `Count_Direction`: UP (increased), DOWN (decreased), or NC (no change)",
  "- `Weight_Direction`: STRONGER, WEAKER, or NC",
  "",
  "**functional_clusters.csv (13_functional_sim/)**",
  "- `Pathway`: Signaling pathway name",
  "- `Cluster`: Functional similarity cluster assignment",
  "- `UMAP1/UMAP2`: Coordinates in functional similarity embedding",
  "- `{Condition}`: Yes/No indicating if pathway is active in that condition",
  "",
  "**structural_clusters.csv (14_structural_sim/)**",
  "- `Pathway`: Signaling pathway name",
  "- `Cluster`: Structural similarity cluster assignment",
  "- `UMAP1/UMAP2`: Coordinates in structural similarity embedding",
  "- `{Condition}`: Yes/No indicating if pathway is active in that condition",
  "",
  "**interaction_summary.csv (15_summary/)**",
  "- `Condition`: Analysis condition name",
  "- `Total_Interactions`: Total number of cell-type pair interactions",
  "- `Total_Strength`: Sum of all communication probabilities",
  "- `Num_LR_Pairs`: Number of significant L-R pairs",
  "- `Num_Pathways`: Number of active signaling pathways",
  "",
  "**pathway_matrix.csv (16_pathway_matrix/)**",
  "- `Pathway`: Signaling pathway name",
  "- `{Condition}`: Binary (1/0) indicating pathway presence in each condition",
  "",
  "**pathway_info_flow.csv (17_pathway_flow/)**",
  "- `Pathway`: Signaling pathway name",
  "- `{Condition}`: Sum of communication probabilities (information flow) for each condition",
  "",
  "**pathway_stats_data.csv (18_pathway_stats/)**",
  "- `Comparison_Pair`: Test_vs_Reference comparison identifier",
  "- `Pathway`: Signaling pathway name",
  "- `Reference/Test_Condition`: Condition names",
  "- `Ref_InfoFlow/Test_InfoFlow`: Total information flow (sum of probabilities)",
  "- `Log2FC`: Log2 fold change of information flow",
  "- `Direction`: UP or DOWN based on log2FC sign",
  "- `PValue`: Rank-based p-value",
  "- `FDR`: Benjamini-Hochberg adjusted p-value",
  "- `Significant`: TRUE if FDR < 0.05",
  "",
  "### Step-by-Step Interpretation Guide",
  "",
  "---",
  "",
  "#### 01-02: Input Data and Database (Root Directory)",
  "",
  "**Purpose:** Validate input data quality and understand the L-R interaction database used.",
  "",
  "**Files:**",
  "- `01_data_summary.csv`: Cell/gene counts, cell type distribution",
  "- `01_celltype_by_condition.csv`: Cross-tabulation showing cell composition balance",
  "- `02_database_categories.csv`: Signaling categories (Secreted, ECM-Receptor, Cell-Cell Contact)",
  "- `02_database_interactions.csv`: Full L-R database for custom queries",
  "",
  "**Interpretation:**",
  "- Check if cell types have sufficient cells (>10 recommended per condition)",
  "- Imbalanced cell composition may bias results toward abundant cell types",
  "- Use database files to find specific L-R pairs of interest",
  "",
  "**Questions Answered:**",
  "- Are there enough cells per cell type for reliable inference?",
  "- Which signaling categories dominate the database?",
  "- Is a specific L-R pair in the database?",
  "",
  "---",
  "",
  "#### 03: Interaction Heatmaps (03_heatmaps/)",
  "",
  "**Purpose:** Identify which cell type pairs communicate most actively.",
  "",
  "**Plots:**",
  "- `{cond}_count_heatmap.png`: NUMBER of significant L-R pairs between cell types",
  "- `{cond}_weight_heatmap.png`: STRENGTH (sum of probabilities) of communication",
  "",
  "**Data Table:** `heatmap_data.csv` (Columns: Condition, Sender, Receiver, Count, Weight)",
  "",
  "**Interpretation:**",
  "- **Rows = Senders** (cells expressing ligands), **Columns = Receivers** (cells expressing receptors)",
  "- **High Count**: Many different L-R pairs active between this pair",
  "- **High Weight**: Strong overall communication regardless of L-R diversity",
  "- Diagonal entries = autocrine signaling (cell type signals to itself)",
  "- Asymmetric patterns suggest directional communication flow",
  "",
  "**Questions Answered:**",
  "- Which cell types are the major signaling hubs?",
  "- Is communication primarily autocrine or paracrine?",
  "- Which cell type pairs have the strongest crosstalk?",
  "",
  "**Downstream Analysis:**",
  "- Filter `heatmap_data.csv` by Condition to compare patterns",
  "- High-interaction pairs warrant deeper investigation in bubble plots",
  "",
  "---",
  "",
  "#### 04: Circle Plots (04_circle_plots/)",
  "",
  "**Purpose:** Visualize the overall communication network topology.",
  "",
  "**Plots:** `{cond}_circle_plot.png`",
  "",
  "**Data Table:** `circle_plot_data.csv` (Columns: Condition, CellType, CellCount)",
  "",
  "**Interpretation:**",
  "- **Node size** = Number of cells (larger nodes = more abundant cell types)",
  "- **Edge width** = Number of interactions between cell types",
  "- **Edge color** = Sender cell type",
  "- Dense networks suggest high intercellular communication",
  "- Hub nodes (many connections) are central to the signaling network",
  "",
  "**Questions Answered:**",
  "- What is the overall structure of cell-cell communication?",
  "- Which cell types are communication hubs vs. peripheral?",
  "- Are there isolated cell types with minimal communication?",
  "",
  "---",
  "",
  "#### 05: Bubble Plots (05_bubble_plots/)",
  "",
  "**Purpose:** Examine specific L-R pairs driving communication between cell types.",
  "",
  "**Plots:** `{cond}_{category}_bubble.png` (split by Secreted/ECM/Contact)",
  "",
  "**Data Table:** Uses `LR_pairs_data.csv` from 07_LR_pairs/",
  "",
  "**Interpretation:**",
  "- **X-axis**: Source-Target cell type pairs",
  "- **Y-axis**: Individual L-R pairs grouped by pathway",
  "- **Bubble size** = Communication probability (larger = stronger)",
  "- **Bubble color** = P-value (darker = more significant)",
  "- Missing bubbles = no significant communication for that pair",
  "",
  "**Questions Answered:**",
  "- Which specific L-R pairs mediate communication?",
  "- Are interactions driven by secreted factors, ECM, or direct contact?",
  "- Which pathways are most active for each cell type pair?",
  "",
  "**Tips:**",
  "- Focus on large, dark bubbles (strong + significant)",
  "- Compare same L-R pair across different cell type pairs",
  "",
  "---",
  "",
  "#### 06: Signaling Roles (06_signaling_roles/)",
  "",
  "**Purpose:** Classify cell types as dominant senders, receivers, or both.",
  "",
  "**Plots:** `{cond}_signaling_roles.png`",
  "",
  "**Data Table:** `signaling_roles_data.csv` (Columns: Condition, CellType, Outgoing_Strength, Incoming_Strength)",
  "",
  "**Interpretation:**",
  "- **X-axis (Outgoing)** = Total signal sent by cell type",
  "- **Y-axis (Incoming)** = Total signal received by cell type",
  "- **Upper-left quadrant**: Dominant receivers (receive more than send)",
  "- **Lower-right quadrant**: Dominant senders (send more than receive)",
  "- **Upper-right quadrant**: Communication hubs (both send and receive)",
  "- **Lower-left quadrant**: Peripheral cells (minimal communication)",
  "",
  "**Questions Answered:**",
  "- Which cell types drive signaling vs. respond to signals?",
  "- Are there specialized sender/receiver populations?",
  "- Which cell types are signaling hubs?",
  "",
  "---",
  "",
  "#### 07: L-R Pairs (07_LR_pairs/)",
  "",
  "**Purpose:** Complete list of significant ligand-receptor interactions.",
  "",
  "**Data Table:** `LR_pairs_data.csv`",
  "",
  "**Key Columns:**",
  "- `source/target`: Sender and receiver cell types",
  "- `ligand/receptor`: Gene names",
  "- `prob`: Communication probability (higher = stronger)",
  "- `pval`: Statistical significance",
  "- `pathway_name`: Signaling pathway family",
  "- `annotation`: Category (Secreted/ECM-Receptor/Cell-Cell Contact)",
  "",
  "**Usage:**",
  "- Filter by `pathway_name` to focus on specific signaling",
  "- Sort by `prob` to find strongest interactions",
  "- Filter by `source` or `target` for cell-type-specific analysis",
  "- Compare same L-R pair across conditions using `Condition` column",
  "",
  "**Questions Answered:**",
  "- What are the top L-R pairs in my dataset?",
  "- Which pathways are active for a specific cell type?",
  "- Is a specific L-R pair (e.g., TGFB1-TGFBR1) present?",
  "",
  "---",
  "",
  "#### 08: Pathways (08_pathways/)",
  "",
  "**Purpose:** List of active signaling pathways per condition.",
  "",
  "**Data Table:** `pathways_data.csv` (Columns: Condition, Pathway)",
  "",
  "**Interpretation:**",
  "- Pathways are aggregations of related L-R pairs",
  "- Compare pathway lists across conditions to find condition-specific signaling",
  "",
  "**Questions Answered:**",
  "- Which signaling pathways are active in each condition?",
  "- Are there condition-specific pathways?",
  "",
  "---",
  "",
  "#### 09: Total Interactions Comparison (09_compare_total/)",
  "",
  "**Purpose:** Compare overall communication activity across conditions.",
  "",
  "**Plots:** `compare_interactions.png` (side-by-side bar plots)",
  "",
  "**Data Table:** `compare_interactions.csv`",
  "",
  "**Interpretation:**",
  "- **Left plot**: Total interaction COUNT (number of cell-type pairs)",
  "- **Right plot**: Total interaction WEIGHT (sum of probabilities)",
  "- Percentage change shown relative to reference condition",
  "- Increased interactions may indicate enhanced intercellular communication (e.g., inflammation, fibrosis)",
  "",
  "**Questions Answered:**",
  "- Does disease/treatment increase or decrease overall communication?",
  "- Is the change in number of interactions or strength?",
  "",
  "---",
  "",
  "#### 10: Differential Interactions (10_differential/)",
  "",
  "**Purpose:** Identify cell type pairs with altered communication between conditions.",
  "",
  "**Plots:**",
  "- `{test}_vs_{ref}_count.png`: Changes in interaction count",
  "- `{test}_vs_{ref}_weight.png`: Changes in interaction strength",
  "",
  "**Data Table:** `differential_data.csv`",
  "",
  "**Interpretation:**",
  "- **RED edges/cells** = INCREASED in test vs reference",
  "- **BLUE edges/cells** = DECREASED in test vs reference",
  "- **Edge width** = Magnitude of change",
  "- Title shows summary: UP/DOWN counts and net change",
  "",
  "**Key Columns in CSV:**",
  "- `Comparison_Pair`: Filter by this for specific comparisons",
  "- `Diff_Count/Diff_Weight`: Positive = increased, Negative = decreased",
  "- `Count_Direction`: UP, DOWN, or NC (no change)",
  "",
  "**Questions Answered:**",
  "- Which cell type pairs show altered communication in disease?",
  "- Is the change symmetric or does one direction dominate?",
  "- Which specific pairs should be investigated further?",
  "",
  "---",
  "",
  "#### 11: Information Flow (11_info_flow/)",
  "",
  "**Purpose:** Rank signaling pathways by activity across conditions.",
  "",
  "**Plots:** `information_flow.png` (stacked and grouped bar plots)",
  "",
  "**Data Table:** Uses `pathway_info_flow.csv` from 17_pathway_flow/",
  "",
  "**Interpretation:**",
  "- **Stacked plot**: Shows relative contribution of each condition",
  "- **Grouped plot**: Direct comparison of pathway activity",
  "- Pathways sorted by total information flow",
  "- Large differences between conditions indicate pathway regulation",
  "",
  "**Questions Answered:**",
  "- Which pathways are most active overall?",
  "- Which pathways show the largest condition differences?",
  "- Are certain pathways condition-specific?",
  "",
  "---",
  "",
  "#### 12: Pathway Networks (12_pathway_networks/)",
  "",
  "**Purpose:** Visualize individual pathway communication patterns.",
  "",
  "**Plots:** `{cond}_{pathway}_network.png`",
  "",
  "**Data Table:** Uses `LR_pairs_data.csv` from 07_LR_pairs/ (filter by pathway_name)",
  "",
  "**Interpretation:**",
  "- Circle plot showing only L-R pairs from one pathway",
  "- Identifies which cell types participate in specific signaling",
  "- Compare same pathway across conditions to see pattern changes",
  "",
  "**Questions Answered:**",
  "- Which cell types send/receive signals in pathway X?",
  "- Does the pathway network structure change between conditions?",
  "",
  "---",
  "",
  "#### 13: Functional Similarity (13_functional_sim/)",
  "",
  "**Purpose:** Group pathways by similar sender-receiver patterns.",
  "",
  "**Plots:** `functional_similarity.png` (UMAP embedding)",
  "",
  "**Data Table:** `functional_clusters.csv`",
  "",
  "**Interpretation:**",
  "- Pathways close together have similar communication patterns",
  "- Clusters represent functionally related signaling programs",
  "- Condition columns show which pathways are active where",
  "- Pathways that cluster together may be co-regulated",
  "",
  "**Questions Answered:**",
  "- Which pathways have similar signaling patterns?",
  "- Are there distinct signaling programs (clusters)?",
  "- Do certain pathways always co-occur?",
  "",
  "---",
  "",
  "#### 14: Structural Similarity (14_structural_sim/)",
  "",
  "**Purpose:** Group pathways by network topology (who talks to whom).",
  "",
  "**Plots:** `structural_similarity.png` (UMAP embedding)",
  "",
  "**Data Table:** `structural_clusters.csv`",
  "",
  "**Interpretation:**",
  "- Pathways close together involve similar cell type pairs",
  "- Different from functional: focuses on network structure, not strength",
  "- Useful for identifying pathways that could substitute for each other",
  "",
  "**Questions Answered:**",
  "- Which pathways use the same cell type communication routes?",
  "- Are there redundant signaling pathways?",
  "",
  "---",
  "",
  "#### 15: Interaction Summary (15_summary/)",
  "",
  "**Purpose:** Quick reference statistics per condition.",
  "",
  "**Data Table:** `interaction_summary.csv`",
  "",
  "**Key Columns:**",
  "- `Total_Interactions`: Sum of all cell-type pair interactions",
  "- `Total_Strength`: Sum of all communication probabilities",
  "- `Num_LR_Pairs`: Count of significant L-R pairs",
  "- `Num_Pathways`: Count of active pathways",
  "",
  "**Questions Answered:**",
  "- How do conditions compare in overall communication metrics?",
  "- Which condition has the most diverse signaling (most pathways)?",
  "",
  "---",
  "",
  "#### 16: Pathway Matrix (16_pathway_matrix/)",
  "",
  "**Purpose:** Binary presence/absence of pathways across conditions.",
  "",
  "**Data Table:** `pathway_matrix.csv`",
  "",
  "**Interpretation:**",
  "- 1 = pathway active in condition, 0 = not active",
  "- Filter for pathways with mixed 0/1 to find condition-specific signaling",
  "- Sum across conditions to find ubiquitous vs. specific pathways",
  "",
  "**Questions Answered:**",
  "- Which pathways are unique to specific conditions?",
  "- Which pathways are conserved across all conditions?",
  "",
  "---",
  "",
  "#### 17: Pathway Information Flow (17_pathway_flow/)",
  "",
  "**Purpose:** Quantitative pathway activity for statistical comparison.",
  "",
  "**Data Table:** `pathway_info_flow.csv`",
  "",
  "**Interpretation:**",
  "- Values = sum of communication probabilities for each pathway",
  "- Higher values = more active signaling",
  "- Use for custom statistical tests or visualization",
  "- This data underlies the 11_info_flow plots",
  "",
  "**Questions Answered:**",
  "- What is the exact information flow value for pathway X?",
  "- Which pathways show the largest absolute differences?",
  "",
  "---",
  "",
  "#### 18: Pathway Statistics (18_pathway_stats/)",
  "",
  "**Purpose:** Statistical comparison of pathway activity between conditions.",
  "",
  "**Plots:** `{test}_vs_{ref}_pathway_logfc.png`",
  "",
  "**Data Table:** `pathway_stats_data.csv`",
  "",
  "**Key Columns:**",
  "- `Log2FC`: Log2 fold change (positive = UP in test, negative = DOWN)",
  "- `Direction`: UP or DOWN",
  "- `PValue/FDR`: Statistical significance (FDR < 0.05 = significant)",
  "- `Significant`: Boolean flag for FDR < 0.05",
  "",
  "**Interpretation:**",
  "- **RED bars** = Pathways UP-regulated in test condition",
  "- **BLUE bars** = Pathways DOWN-regulated in test condition",
  "- **Asterisk (*)** = Statistically significant (FDR < 0.05)",
  "- Focus on significant pathways with large |Log2FC|",
  "",
  "**Questions Answered:**",
  "- Which pathways are significantly altered between conditions?",
  "- What is the magnitude of pathway changes?",
  "- Which pathways should be prioritized for validation?",
  "",
  "---",
  "",
  "## Methods",
  "",
  sprintf("Cell-cell communication was inferred using CellChat v%s (Jin et al., Nature Protocols 2024).", packageVersion("CellChat")),
  "",
  "**Workflow:**",
  "1. Expression data extracted from Seurat object (log-normalized)",
  "2. CellChatDB v2 used (excluding non-protein signaling)",
  "3. Over-expressed genes and interactions identified per condition",
  sprintf("4. Communication probabilities computed using truncatedMean method (trim=%.2f)", CONFIG$trim),
  "5. Pathway-level signaling aggregated from L-R pairs",
  "6. Network centrality computed for signaling role analysis",
  "7. Conditions merged for comparative analysis",
  "8. Functional and structural similarity computed via pairwise Jaccard index and UMAP embedding",
  "",
  "**Citation:** Jin et al., CellChat for systematic analysis of cell-cell communication, Nature Protocols 2024",
  "",
  "---",
  "",
  sprintf("*Report generated by CellChat Analysis Pipeline on %s*", Sys.time())
)

writeLines(report, file.path(CONFIG$output_dir, CONFIG$report_file))
cat(sprintf("  - Report saved to %s\n", file.path(CONFIG$output_dir, CONFIG$report_file)))

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Analysis Complete!\n")
cat("=============================================================================\n")
cat(sprintf("Output directory: %s\n", CONFIG$output_dir))
cat(sprintf("Report file: %s\n", file.path(CONFIG$output_dir, CONFIG$report_file)))
cat(sprintf("Total figures generated: %d\n", fig_counter - 1))
cat(sprintf("Total tables generated: %d\n", table_counter - 1))
cat(sprintf("End time: %s\n", Sys.time()))
