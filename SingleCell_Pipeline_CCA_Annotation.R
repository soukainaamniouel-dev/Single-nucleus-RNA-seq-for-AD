##########################################################
# 0) PACKAGES
############################################################
setwd("H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq/")

# Core
library(Seurat)
library(tidyverse)   # dplyr, ggplot2, tibble, etc.
library(patchwork)
library(cowplot)
library(magrittr)
library(data.table)

# Normalization / models
library(glmGamPoi)
library(future)

#Doublets
library(DoubletFinder)
library(scDblFinder)   # loaded in original; not used here (DoubletFinder is primary)

# Annotation / references
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scCustomize)

#DE / abundance (guarded for replicates; descriptive path used for n=1)
library(DESeq2)
# library(speckle) 

# Viz / misc
library(viridis)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(lisi)
library(qs)

# Optional downstream (uncomment if used): CellChat, NMF, plot1cell, metap, multtest

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB; raise if you have RAM

set.seed(1234)  # reproducibility

############################################################
# RESTART CCA CONSERVED-MARKER ANALYSIS FROM A SAVED QS FILE
#
# Uses the saved CCA cluster assignments, but tests markers
# on the joined, log-normalized RNA assay.
#
# A "strong conserved" gene must satisfy ALL three samples:
#   - adjusted p-value <= 0.05
#   - |avg_log2FC| > 1
#   - |pct.1 - pct.2| >= 0.20
#   - the direction must be consistent across all samples
#
# Positive conserved markers:
#   avg_log2FC > 1 and pct.1 - pct.2 >= 0.20 in all samples
#
# Negative conserved markers:
#   avg_log2FC < -1 and pct.1 - pct.2 <= -0.20 in all samples
############################################################

rm(list = ls())
gc()

############################################################
# 1) PACKAGES AND SETTINGS
############################################################

library(Seurat)
library(qs)
library(dplyr)

if (!requireNamespace("metap", quietly = TRUE)) {
  stop(
    "The 'metap' package is required. Install it, restart R, ",
    "and rerun this script.\n",
    "Suggested commands:\n",
    "  install.packages('BiocManager')\n",
    "  BiocManager::install('multtest')\n",
    "  install.packages('metap')"
  )
}

set.seed(1234)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 8 * 1024^3)

# Samples that must all support a marker.
required_samples <- c("Ctrl", "Fibrils", "FibJ8")

# Strict final thresholds.
min_abs_log2fc       <- 1.00
max_adjusted_p       <- 0.05
min_abs_pct_diff     <- 0.25  # 0.20 = 20 percentage points
min_cells_per_group  <- 10L
min_detection        <- 0.10
n_top_markers        <- 25L


############################################################
# 2) LOCATE AND LOAD THE SAVED CCA OBJECT
############################################################

project_dir <- paste0(
  "H:/Documents/Qi_Projects/AD_scRNAseq/",
  "AD_Qi_scRNAseq"
)

# The first existing file in this list will be loaded.
cca_candidates <- c(
  file.path(
    project_dir,
    "output_new_process_scRNAseq",
    "integrated_CCA_with_cluster_markers.qs"
  ),
  file.path(
    project_dir,
    "output_new_process_scRNAseq",
    "integrated_CCA.qs"
  ),
  file.path(
    project_dir,
    "output_new_process_scRNAseq",
    "v11_output",
    "integrated_CCA_with_cluster_markers.qs"
  ),
  file.path(
    project_dir,
    "output_new_process_scRNAseq",
    "v11_output",
    "integrated_CCA.qs"
  )
)

existing_cca_files <- cca_candidates[file.exists(cca_candidates)]

if (length(existing_cca_files) == 0L) {
  stop(
    "No saved CCA QS file was found. Checked:\n",
    paste(cca_candidates, collapse = "\n")
  )
}

cca_file <- existing_cca_files[[1]]

message("Loading CCA object from: ", cca_file)
cca <- qs::qread(cca_file)

if (!inherits(cca, "Seurat")) {
  stop("The loaded object is not a Seurat object.")
}

message(
  "Loaded ",
  ncol(cca),
  " cells and ",
  nrow(cca),
  " features."
)

# Write all new outputs into a separate folder so that old files
# are not overwritten.
results_dir <- file.path(
  dirname(cca_file),
  "cca_conserved_markers_restart"
)

by_cluster_dir <- file.path(
  results_dir,
  "positive_markers_by_cluster"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  by_cluster_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3) RESTORE CCA CLUSTER IDENTITIES
############################################################

metadata_columns <- colnames(cca[[]])

if ("cca_cluster" %in% metadata_columns) {
  cluster_source <- "cca_cluster"
} else if ("seurat_clusters" %in% metadata_columns) {
  cluster_source <- "seurat_clusters"
} else {
  stop(
    "Neither 'cca_cluster' nor 'seurat_clusters' is present ",
    "in the saved CCA object's metadata."
  )
}

cluster_values <- as.character(
  cca[[]][[cluster_source]]
)

if (anyNA(cluster_values) || any(cluster_values == "")) {
  stop("The saved CCA object has missing cluster labels.")
}

cluster_levels <- unique(cluster_values)

# Keep numeric clusters in natural numeric order.
if (all(grepl("^[0-9]+$", cluster_levels))) {
  cluster_levels <- as.character(
    sort(as.integer(cluster_levels))
  )
} else {
  cluster_levels <- sort(cluster_levels)
}

cca$cca_cluster <- factor(
  cluster_values,
  levels = cluster_levels
)

Idents(cca) <- "cca_cluster"

message(
  "Using CCA clusters: ",
  paste(levels(Idents(cca)), collapse = ", ")
)

############################################################
# 4) VERIFY THE THREE SAMPLES
############################################################

if (!"sample" %in% metadata_columns) {
  stop("Metadata column 'sample' is missing from the CCA object.")
}

sample_values <- as.character(
  cca[[]][["sample"]]
)

if (anyNA(sample_values) || any(sample_values == "")) {
  stop("The sample metadata contains missing or empty values.")
}

samples_present <- sort(unique(sample_values))

if (!setequal(samples_present, required_samples)) {
  stop(
    "Expected exactly these samples: ",
    paste(required_samples, collapse = ", "),
    ". Found: ",
    paste(samples_present, collapse = ", ")
  )
}

cca$sample <- factor(
  sample_values,
  levels = required_samples
)

print(table(cca$sample))
print(table(cca$cca_cluster, cca$sample))

############################################################
# 5) PREPARE RNA FOR DIFFERENTIAL EXPRESSION
############################################################

# Obtain assay names directly from the Seurat object.
# This avoids a possible function-name conflict involving Assays().
assay_names <- names(cca@assays)

message(
  "Assays in the saved CCA object: ",
  paste(assay_names, collapse = ", ")
)

if (!("RNA" %in% assay_names)) {
  stop(
    "The saved CCA object does not contain an RNA assay. ",
    "Assays found: ",
    paste(assay_names, collapse = ", ")
  )
}

SeuratObject::DefaultAssay(cca) <- "RNA"

# Seurat v5 assays can contain multiple expression layers.
# Join them before differential-expression testing.
if (inherits(cca[["RNA"]], "Assay5")) {
  
  cca <- SeuratObject::JoinLayers(
    object = cca,
    assay = "RNA"
  )
  
  rna_layers <- SeuratObject::Layers(
    cca[["RNA"]]
  )
  
  message(
    "RNA layers after JoinLayers: ",
    paste(rna_layers, collapse = ", ")
  )
  
  if (!("counts" %in% rna_layers)) {
    stop(
      "The joined RNA assay does not contain a counts layer. ",
      "Layers found: ",
      paste(rna_layers, collapse = ", ")
    )
  }
  
  # Create normalized RNA data only when it is absent.
  if (!("data" %in% rna_layers)) {
    
    message(
      "RNA data layer is absent; running NormalizeData()."
    )
    
    cca <- Seurat::NormalizeData(
      object = cca,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
    
  } else {
    
    message(
      "Existing normalized RNA data layer will be used."
    )
  }
  
} else {
  
  # Compatibility path for an older Seurat RNA Assay object.
  message(
    "RNA is an older Assay object; running NormalizeData()."
  )
  
  cca <- Seurat::NormalizeData(
    object = cca,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

# Restore the settings used for marker analysis.
SeuratObject::DefaultAssay(cca) <- "RNA"
SeuratObject::Idents(cca) <- "cca_cluster"

message(
  "Default assay: ",
  SeuratObject::DefaultAssay(cca)
)

message(
  "Number of CCA clusters: ",
  length(levels(SeuratObject::Idents(cca)))
)

############################################################
# 6) IDENTIFY CLUSTERS REPRESENTED IN EVERY SAMPLE
############################################################

cluster_sample_n <- table(
  cluster = factor(
    as.character(cca$cca_cluster),
    levels = cluster_levels
  ),
  sample = factor(
    as.character(cca$sample),
    levels = required_samples
  )
)

sample_totals <- colSums(cluster_sample_n)

eligible_clusters <- rownames(cluster_sample_n)[
  apply(
    cluster_sample_n,
    MARGIN = 1,
    FUN = function(cluster_counts) {
      enough_in_cluster <- all(
        cluster_counts >= min_cells_per_group
      )
      
      enough_outside_cluster <- all(
        sample_totals - cluster_counts >= min_cells_per_group
      )
      
      enough_in_cluster && enough_outside_cluster
    }
  )
]

if (length(eligible_clusters) == 0L) {
  stop(
    "No cluster has at least ",
    min_cells_per_group,
    " cells inside and outside the cluster in every sample."
  )
}

message(
  "Clusters eligible in all three samples: ",
  paste(eligible_clusters, collapse = ", ")
)

write.csv(
  as.data.frame.matrix(cluster_sample_n),
  file.path(
    results_dir,
    "cca_cluster_cell_counts_by_sample.csv"
  ),
  row.names = TRUE
)

############################################################
# 7) RUN FindConservedMarkers() FOR EACH CCA CLUSTER
############################################################

conserved_marker_list_raw <- setNames(
  vector("list", length(eligible_clusters)),
  eligible_clusters
)

for (cluster_id in eligible_clusters) {
  
  message(
    "Finding conserved markers for CCA cluster ",
    cluster_id
  )
  
  result <- tryCatch(
    FindConservedMarkers(
      object = cca,
      ident.1 = cluster_id,
      ident.2 = NULL,
      grouping.var = "sample",
      assay = "RNA",
      slot = "data",
      test.use = "wilcox",
      only.pos = FALSE,
      
      # These are efficient candidate filters. The code below
      # applies the strict final filters again to every sample.
      logfc.threshold = min_abs_log2fc,
      min.pct = min_detection,
      min.diff.pct = min_abs_pct_diff,
      
      min.cells.group = min_cells_per_group,
      base = 2,
      random.seed = 1234,
      meta.method = metap::minimump,
      verbose = FALSE
    ),
    error = function(e) {
      warning(
        "Cluster ",
        cluster_id,
        " failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  if (is.null(result) || nrow(result) == 0L) {
    warning(
      "No conserved candidates were returned for cluster ",
      cluster_id,
      ".",
      call. = FALSE
    )
    
    next
  }
  
  result <- as.data.frame(
    result,
    check.names = FALSE
  )
  
  result$gene <- rownames(result)
  result$cluster <- as.character(cluster_id)
  rownames(result) <- NULL
  
  conserved_marker_list_raw[[cluster_id]] <- result
}

conserved_marker_list_raw <- Filter(
  Negate(is.null),
  conserved_marker_list_raw
)

if (length(conserved_marker_list_raw) == 0L) {
  stop("No cluster returned a conserved-marker table.")
}

############################################################
# 8) STRICTLY FILTER EVERY SAMPLE-SPECIFIC RESULT
############################################################

add_strict_filters <- function(marker_table) {
  
  if (is.null(marker_table) || nrow(marker_table) == 0L) {
    return(NULL)
  }
  
  # Seurat v5 normally uses avg_log2FC. Older versions may use
  # avg_logFC; support either naming convention.
  fc_cols_v5 <- paste0(
    required_samples,
    "_avg_log2FC"
  )
  
  fc_cols_old <- paste0(
    required_samples,
    "_avg_logFC"
  )
  
  if (all(fc_cols_v5 %in% colnames(marker_table))) {
    fc_cols <- fc_cols_v5
  } else if (all(fc_cols_old %in% colnames(marker_table))) {
    fc_cols <- fc_cols_old
  } else {
    stop(
      "Could not find all sample-specific log-fold-change columns.\n",
      "Available columns:\n",
      paste(colnames(marker_table), collapse = ", ")
    )
  }
  
  p_adj_cols <- paste0(
    required_samples,
    "_p_val_adj"
  )
  
  pct1_cols <- paste0(
    required_samples,
    "_pct.1"
  )
  
  pct2_cols <- paste0(
    required_samples,
    "_pct.2"
  )
  
  needed_cols <- c(
    fc_cols,
    p_adj_cols,
    pct1_cols,
    pct2_cols
  )
  
  missing_cols <- setdiff(
    needed_cols,
    colnames(marker_table)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(marker_table), collapse = ", ")
    )
  }
  
  to_numeric_matrix <- function(column_names) {
    values <- lapply(
      marker_table[, column_names, drop = FALSE],
      function(x) as.numeric(as.character(x))
    )
    
    output <- as.matrix(
      as.data.frame(
        values,
        check.names = FALSE
      )
    )
    
    colnames(output) <- column_names
    output
  }
  
  fc_matrix <- to_numeric_matrix(fc_cols)
  p_adj_matrix <- to_numeric_matrix(p_adj_cols)
  pct1_matrix <- to_numeric_matrix(pct1_cols)
  pct2_matrix <- to_numeric_matrix(pct2_cols)
  
  delta_pct_matrix <- pct1_matrix - pct2_matrix
  
  # Save sample-specific pct differences into the output table.
  for (i in seq_along(required_samples)) {
    marker_table[[
      paste0(
        required_samples[[i]],
        "_delta_pct"
      )
    ]] <- delta_pct_matrix[, i]
  }
  
  all_rows_pass <- function(value_matrix, condition_function) {
    apply(
      value_matrix,
      MARGIN = 1,
      FUN = function(values) {
        all(
          !is.na(values) &
            is.finite(values) &
            condition_function(values)
        )
      }
    )
  }
  
  # Statistical significance in Ctrl, Fibrils, and FibJ8.
  marker_table$pass_p_adj_all_samples <- all_rows_pass(
    p_adj_matrix,
    function(x) x <= max_adjusted_p
  )
  
  # Positive markers: higher expression and higher detection
  # in the cluster in all three samples.
  marker_table$pass_positive_all_samples <- (
    all_rows_pass(
      fc_matrix,
      function(x) x > min_abs_log2fc
    ) &
      all_rows_pass(
        delta_pct_matrix,
        function(x) x >= min_abs_pct_diff
      )
  )
  
  # Negative markers: lower expression and lower detection
  # in the cluster in all three samples.
  marker_table$pass_negative_all_samples <- (
    all_rows_pass(
      fc_matrix,
      function(x) x < -min_abs_log2fc
    ) &
      all_rows_pass(
        delta_pct_matrix,
        function(x) x <= -min_abs_pct_diff
      )
  )
  
  marker_table$direction <- ifelse(
    marker_table$pass_positive_all_samples,
    "Higher_in_cluster",
    ifelse(
      marker_table$pass_negative_all_samples,
      "Lower_in_cluster",
      NA_character_
    )
  )
  
  marker_table$pass_strong_conserved <- (
    marker_table$pass_p_adj_all_samples &
      (
        marker_table$pass_positive_all_samples |
          marker_table$pass_negative_all_samples
      )
  )
  
  # Worst-case values across the three samples. These are useful
  # for ranking because a conserved marker is only as strong as
  # its weakest sample.
  marker_table$minimum_abs_log2FC <- apply(
    abs(fc_matrix),
    MARGIN = 1,
    FUN = min
  )
  
  marker_table$minimum_abs_pct_difference <- apply(
    abs(delta_pct_matrix),
    MARGIN = 1,
    FUN = min
  )
  
  marker_table$maximum_sample_p_adj <- apply(
    p_adj_matrix,
    MARGIN = 1,
    FUN = max
  )
  
  marker_table
}

conserved_marker_list_annotated <- lapply(
  conserved_marker_list_raw,
  add_strict_filters
)

cca_conserved_all_candidates <- dplyr::bind_rows(
  conserved_marker_list_annotated
)

cca_conserved_strong <- cca_conserved_all_candidates[
  cca_conserved_all_candidates$pass_strong_conserved %in% TRUE,
  ,
  drop = FALSE
]

cca_conserved_positive <- cca_conserved_all_candidates[
  cca_conserved_all_candidates$pass_p_adj_all_samples %in% TRUE &
    cca_conserved_all_candidates$pass_positive_all_samples %in% TRUE,
  ,
  drop = FALSE
]

cca_conserved_negative <- cca_conserved_all_candidates[
  cca_conserved_all_candidates$pass_p_adj_all_samples %in% TRUE &
    cca_conserved_all_candidates$pass_negative_all_samples %in% TRUE,
  ,
  drop = FALSE
]

############################################################
# 9) SORT RESULTS
############################################################

sort_marker_table <- function(marker_table) {
  
  if (is.null(marker_table) || nrow(marker_table) == 0L) {
    return(marker_table)
  }
  
  cluster_rank <- match(
    as.character(marker_table$cluster),
    cluster_levels
  )
  
  marker_table[
    order(
      cluster_rank,
      -marker_table$minimum_abs_log2FC,
      -marker_table$minimum_abs_pct_difference,
      marker_table$maximum_sample_p_adj
    ),
    ,
    drop = FALSE
  ]
}

cca_conserved_all_candidates <- sort_marker_table(
  cca_conserved_all_candidates
)

cca_conserved_strong <- sort_marker_table(
  cca_conserved_strong
)

cca_conserved_positive <- sort_marker_table(
  cca_conserved_positive
)

cca_conserved_negative <- sort_marker_table(
  cca_conserved_negative
)

###########################################################
# 10) TOP POSITIVE MARKERS FOR CLUSTER ANNOTATION
############################################################

if (nrow(cca_conserved_positive) > 0L) {
  positive_split <- split(
    cca_conserved_positive,
    cca_conserved_positive$cluster
  )
  
  top_positive_list <- lapply(
    positive_split,
    function(marker_table) {
      head(marker_table, n_top_markers)
    }
  )
  
  cca_top_positive_conserved <- dplyr::bind_rows(
    top_positive_list
  )
} else {
  positive_split <- list()
  cca_top_positive_conserved <- cca_conserved_positive
}

############################################################
# 11) SUMMARY COUNTS
############################################################

count_for_cluster <- function(marker_table, cluster_id, pass_column = NULL) {
  
  if (is.null(marker_table) || nrow(marker_table) == 0L) {
    return(0L)
  }
  
  cluster_rows <- as.character(marker_table$cluster) == cluster_id
  
  if (is.null(pass_column)) {
    return(as.integer(sum(cluster_rows)))
  }
  
  as.integer(
    sum(
      cluster_rows &
        marker_table[[pass_column]] %in% TRUE
    )
  )
}

marker_summary <- data.frame(
  cluster = eligible_clusters,
  n_candidates = vapply(
    eligible_clusters,
    function(cluster_id) {
      count_for_cluster(
        cca_conserved_all_candidates,
        cluster_id
      )
    },
    integer(1)
  ),
  n_strong_positive = vapply(
    eligible_clusters,
    function(cluster_id) {
      count_for_cluster(
        cca_conserved_all_candidates,
        cluster_id,
        "pass_positive_all_samples"
      )
    },
    integer(1)
  ),
  n_strong_negative = vapply(
    eligible_clusters,
    function(cluster_id) {
      count_for_cluster(
        cca_conserved_all_candidates,
        cluster_id,
        "pass_negative_all_samples"
      )
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)

# The significance condition is common to the final positive and
# negative tables, so report final retained counts explicitly.
marker_summary$n_final_positive <- vapply(
  eligible_clusters,
  function(cluster_id) {
    as.integer(
      sum(
        as.character(cca_conserved_positive$cluster) == cluster_id
      )
    )
  },
  integer(1)
)

marker_summary$n_final_negative <- vapply(
  eligible_clusters,
  function(cluster_id) {
    as.integer(
      sum(
        as.character(cca_conserved_negative$cluster) == cluster_id
      )
    )
  },
  integer(1)
)

for (sample_id in required_samples) {
  marker_summary[[paste0("n_cells_", sample_id)]] <- as.integer(
    cluster_sample_n[
      eligible_clusters,
      sample_id
    ]
  )
}

print(marker_summary)

############################################################
# 12) SAVE TABLES AND PER-CLUSTER LISTS
############################################################

write.csv(
  cca_conserved_all_candidates,
  file.path(
    results_dir,
    "cca_conserved_all_candidates.csv"
  ),
  row.names = FALSE
)

write.csv(
  cca_conserved_strong,
  file.path(
    results_dir,
    paste0(
      "cca_conserved_strong_both_directions_",
      "absLog2FC_gt1_",
      "absPctDiff_ge20_",
      "padj_le05_all_samples.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cca_conserved_positive,
  file.path(
    results_dir,
    paste0(
      "cca_conserved_positive_",
      "log2FC_gt1_",
      "pctDiff_ge20_",
      "padj_le05_all_samples.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cca_conserved_negative,
  file.path(
    results_dir,
    paste0(
      "cca_conserved_negative_",
      "log2FC_lt_minus1_",
      "pctDiff_le_minus20_",
      "padj_le05_all_samples.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  cca_top_positive_conserved,
  file.path(
    results_dir,
    "cca_top10_positive_conserved_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

write.csv(
  marker_summary,
  file.path(
    results_dir,
    "cca_conserved_marker_summary_by_cluster.csv"
  ),
  row.names = FALSE
)

saveRDS(
  conserved_marker_list_annotated,
  file.path(
    results_dir,
    "cca_conserved_candidates_by_cluster.rds"
  )
)

saveRDS(
  split(
    cca_conserved_positive,
    cca_conserved_positive$cluster
  ),
  file.path(
    results_dir,
    "cca_conserved_positive_by_cluster.rds"
  )
)

saveRDS(
  split(
    cca_conserved_negative,
    cca_conserved_negative$cluster
  ),
  file.path(
    results_dir,
    "cca_conserved_negative_by_cluster.rds"
  )
)

# Also write one positive-marker CSV for each cluster.
if (length(positive_split) > 0L) {
  for (cluster_id in names(positive_split)) {
    write.csv(
      positive_split[[cluster_id]],
      file.path(
        by_cluster_dir,
        paste0(
          "cca_cluster_",
          cluster_id,
          "_positive_conserved_markers.csv"
        )
      ),
      row.names = FALSE
    )
  }
}

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "sessionInfo.txt"
  )
)

############################################################
# 13) FINAL CONSOLE REPORT
############################################################

message("Conserved-marker analysis complete.")
message(
  "Strong positive markers retained: ",
  nrow(cca_conserved_positive)
)
message(
  "Strong negative markers retained: ",
  nrow(cca_conserved_negative)
)
message(
  "Results saved to: ",
  normalizePath(
    results_dir,
    winslash = "/",
    mustWork = FALSE
  )
)

# Main table for cluster annotation:
#   cca_conserved_positive
#
# Top 10 positive markers per cluster:
#   cca_top_positive_conserved
#
# Negative/depleted genes:
#   cca_conserved_negative

###################################################
# ignore cluster 8 , 9 , 12, 10, 14


#####################################################
# Finalizing markers
####################################################

neuron_markers <- list(
  
  cl0 = c("STK32B","CDH13","MCTP1","DOK5","NEO1"), #STK32B⁺/CDH13⁺ deep-layer intratelencephalic-like excitatory neurons #STK32B⁺/CDH13⁺ L6 IT-like excitatory neurons
  
  cl1 = c("DPY19L1","PRSS12", "PLXNA2", "SEMA3C","PTPN4","FRMD4B", "ST3GAL6"), #DPY19L1⁺/SEMA3C⁺ migrating cortical glutamatergic neurons
  
  cl2 = c("PDGFRB","PDE3A","TNC","FBN1","COL11A1","PTN","PTPRM","IL33", "GLI3","LIPG","SLCO1C1", "INTU","DYNC2H1"), #PDGFRB⁺/PDE3A⁺ perivascular mesenchymal cells
  
  cl3 = c("LRRC4C", "RORB","GRM7","KCNH5","NTNG1","CDH18","TENM1","LUZP2"), #RORB⁺ mid-layer excitatory neurons #RORB⁺ mid-layer IT-like excitatory neurons #RORB⁺/KCNH5⁺ mid-layer excitatory neurons
  
  cl4 = c("GAD1", "PDZRN3", "NRXN3","GAD2","DLX6-AS1", "ST18", "DCLK2"), #DLX6-AS1⁺/ST18⁺ developing inhibitory neurons
  
  cl5 = c("COPG2IT1", "PAX6","MEST", "ST18", "DLX6-AS1","GAD2", "ARX", "NRXN3","GAD1","SOX2-OT"), #PAX6⁺/ARX⁺ local immature GABAergic neurons
  
  cl6 = c("OLIG1", "FERMT1","SMOC1", "PDGFRA", "PLP1", "SOX10", "CLDN11", "MBP",
          "DLL3", "CSPG4", "SLC24A3", "CEROX1", "MAP3K1", "POLR2F"), #OLIG1⁺ early OPC, #OLIG1⁺ early oligodendroglial progenitors
  
  cl7 =  c("UNC5D", "SOX11","DCC", "NRP1",  "LRP8","EPHA3",  "DOK6","KCNQ3","CNR1","CDH4"), #UNC5D⁺ immature upper-layer excitatory neurons
  
  cl11 = c("HS3ST4","TRPM3", "TLE4", "FOXP2","ZFPM2", "SOX5","ADAMTSL1", "SEMA3E"), #FOXP2⁺/TLE4⁺ L6 corticothalamic-like neurons
  
  cl13 = c("AQP4", "SLC4A4", "GFAP", "MEGF10", "SPARCL1", "CD44", "GLIS3", "RFX4", "ID4", "EDNRB", "LIFR",
           "LAMA1", "ITGA6", "NID1"), #AQP4⁺/SLC4A4⁺ developmental astrocytes #endfoot-like astrocytes
  
  cl15 = c("SGCZ", "PEX5L","RALYL", "CD36","NWD2","ST6GALNAC3","RELN",
           "HCN1","RYR2", "KCNIP4","CACNA2D3", "CADPS2","GRIK2","GRM8","GRM5") #RELN⁺ L5 deep-layer excitatory neurons

)      

p_dot_neuron <- DotPlot(
  cca,
  features = unique(unlist(neuron_markers))
) + RotatedAxis()

p_dot_neuron

#####################################################################################################################
############################################################
# CREATE TWO CURATED CCA OBJECTS FROM A SAVED QS FILE
#
# Object 1:
#   remove clusters 8, 9, 10, 12, and 14
#
# Object 2:
#   start from Object 1, then merge clusters 4 and 5
#
# The selected source QS file is never overwritten.
############################################################

library(Seurat)
library(qs)
library(ggplot2)

############################################################
# 1) SELECT AND LOAD THE UNTOUCHED CCA OBJECT
############################################################

# In the file chooser, select either:
#   integrated_CCA.qs
# or
#   integrated_CCA_with_cluster_markers.qs
source_qs <- file.choose()

cca_source <- qs::qread(source_qs)

if (!inherits(cca_source, "Seurat")) {
  stop("The selected QS file does not contain a Seurat object.")
}

message("Loaded: ", normalizePath(source_qs, winslash = "/", mustWork = TRUE))
message("Cells in source object: ", ncol(cca_source))
############################################################
# 2) SETTINGS
############################################################

clusters_to_remove <- c("8", "9", "10", "12", "14")
clusters_to_merge  <- c("4", "5")
merged_cluster_name <- "4_5"

# Save both new objects in a separate folder next to the source QS.
output_dir <- file.path(
  dirname(source_qs),
  "curated_cluster_objects"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_dir <- file.path(output_dir, "plots")
dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3) IDENTIFY THE ORIGINAL CCA CLUSTER COLUMN
############################################################

metadata_columns <- colnames(cca_source[[]])

# Prefer seurat_clusters because FindClusters() created this column.
# This avoids accidentally reusing a previously edited cca_cluster column.
if ("seurat_clusters" %in% metadata_columns) {
  cluster_source_column <- "seurat_clusters"
} else if ("cca_cluster_original" %in% metadata_columns) {
  cluster_source_column <- "cca_cluster_original"
} else if ("cca_cluster" %in% metadata_columns) {
  cluster_source_column <- "cca_cluster"
} else {
  stop(
    "No cluster column was found. Expected one of: ",
    "seurat_clusters, cca_cluster_original, or cca_cluster."
  )
}

source_cluster_raw <- cca_source[[]][[cluster_source_column]]
source_cluster_chr <- as.character(source_cluster_raw)

if (anyNA(source_cluster_chr) || any(source_cluster_chr == "")) {
  stop("The source cluster column contains missing or empty labels.")
}

# Preserve the original cluster order when possible.
if (is.factor(source_cluster_raw)) {
  source_cluster_levels <- levels(droplevels(source_cluster_raw))
} else {
  source_cluster_levels <- unique(source_cluster_chr)
  
  # Use numeric ordering when all labels are integers.
  if (all(grepl("^-?[0-9]+$", source_cluster_levels))) {
    source_cluster_levels <- as.character(
      sort(as.integer(source_cluster_levels))
    )
  }
}

source_cluster_levels <- source_cluster_levels[
  source_cluster_levels %in% source_cluster_chr
]

available_clusters <- source_cluster_levels

message("Original cluster column: ", cluster_source_column)
message("Original clusters: ", paste(available_clusters, collapse = ", "))

required_clusters <- unique(c(clusters_to_remove, clusters_to_merge))
missing_clusters <- setdiff(required_clusters, available_clusters)

if (length(missing_clusters) > 0L) {
  stop(
    "The following requested clusters were not found: ",
    paste(missing_clusters, collapse = ", "),
    ". Available clusters: ",
    paste(available_clusters, collapse = ", ")
  )
}

############################################################
# 4) PRESERVE ORIGINAL LABELS AND SAVE SOURCE COUNTS
############################################################

cca_source$cca_cluster_original <- factor(
  source_cluster_chr,
  levels = source_cluster_levels
)

source_cluster_counts <- as.data.frame(
  table(
    cluster = cca_source$cca_cluster_original,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(source_cluster_counts)[2] <- "n_cells"
source_cluster_counts$version <- "source"

write.csv(
  source_cluster_counts,
  file.path(output_dir, "source_cluster_cell_counts.csv"),
  row.names = FALSE
)

############################################################
# 5) OBJECT 1: REMOVE CLUSTERS 8, 9, 10, 12, AND 14
############################################################

remove_cell_flag <- source_cluster_chr %in% clusters_to_remove
removed_cells <- colnames(cca_source)[remove_cell_flag]
kept_cells <- colnames(cca_source)[!remove_cell_flag]

if (length(removed_cells) == 0L) {
  stop("No cells matched the requested clusters to remove.")
}

removed_cell_manifest <- data.frame(
  cell = removed_cells,
  cluster = source_cluster_chr[remove_cell_flag],
  sample = if ("sample" %in% metadata_columns) {
    as.character(cca_source$sample[remove_cell_flag])
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)

write.csv(
  removed_cell_manifest,
  file.path(output_dir, "removed_cells_clusters_8_9_10_12_14.csv"),
  row.names = FALSE
)

cca_removed <- subset(
  x = cca_source,
  cells = kept_cells,
  droplevels.meta.data = TRUE
)

remaining_levels <- source_cluster_levels[
  !source_cluster_levels %in% clusters_to_remove
]

cca_removed$cca_cluster_original <- factor(
  as.character(cca_removed$cca_cluster_original),
  levels = remaining_levels
)

# Working curated identity. For Object 1, labels remain unchanged.
cca_removed$cca_cluster <- factor(
  as.character(cca_removed$cca_cluster_original),
  levels = remaining_levels
)

cca_removed$cca_cluster_curated <- cca_removed$cca_cluster
SeuratObject::Idents(cca_removed) <- "cca_cluster"

# Store curation history inside the object.
cca_removed@misc$cluster_curation <- list(
  source_qs = normalizePath(source_qs, winslash = "/", mustWork = TRUE),
  source_cluster_column = cluster_source_column,
  clusters_removed = clusters_to_remove,
  clusters_merged = character(0),
  merged_cluster_name = NA_character_,
  created = as.character(Sys.time())
)

if (any(as.character(cca_removed$cca_cluster_original) %in% clusters_to_remove)) {
  stop("At least one requested cluster remains in Object 1.")
}

expected_n_after_removal <- ncol(cca_source) - length(removed_cells)

if (ncol(cca_removed) != expected_n_after_removal) {
  stop("Object 1 has an unexpected number of cells after subsetting.")
}

removed_file <- file.path(
  output_dir,
  "integrated_CCA_removed_clusters_8_9_10_12_14.qs"
)

qs::qsave(
  cca_removed,
  removed_file
)

message(
  "Saved Object 1 (removed clusters only): ",
  normalizePath(removed_file, winslash = "/", mustWork = TRUE)
)

# The untouched source remains on disk. Remove its in-memory copy
# before creating Object 2 to reduce RAM use.
rm(cca_source)
gc()

############################################################
# 6) OBJECT 2: COPY OBJECT 1 AND MERGE CLUSTERS 4 AND 5
############################################################

cca_removed_merge45 <- cca_removed

# Preserve the labels immediately before merging.
cca_removed_merge45$cca_cluster_premerge <- factor(
  as.character(cca_removed_merge45$cca_cluster),
  levels = levels(cca_removed_merge45$cca_cluster)
)

premerge_labels <- as.character(
  cca_removed_merge45$cca_cluster_premerge
)

n_cluster4_before <- sum(premerge_labels == "4")
n_cluster5_before <- sum(premerge_labels == "5")

merged_labels <- premerge_labels
merged_labels[merged_labels %in% clusters_to_merge] <- merged_cluster_name

old_levels <- levels(cca_removed_merge45$cca_cluster_premerge)
first_merge_position <- min(match(clusters_to_merge, old_levels))
new_levels <- old_levels[!old_levels %in% clusters_to_merge]
new_levels <- append(
  new_levels,
  merged_cluster_name,
  after = first_merge_position - 1L
)

cca_removed_merge45$cca_cluster <- factor(
  merged_labels,
  levels = new_levels
)

cca_removed_merge45$cca_cluster_curated <- cca_removed_merge45$cca_cluster
SeuratObject::Idents(cca_removed_merge45) <- "cca_cluster"

cca_removed_merge45@misc$cluster_curation <- list(
  source_qs = normalizePath(source_qs, winslash = "/", mustWork = TRUE),
  source_cluster_column = cluster_source_column,
  clusters_removed = clusters_to_remove,
  clusters_merged = clusters_to_merge,
  merged_cluster_name = merged_cluster_name,
  created = as.character(Sys.time())
)

if (any(as.character(cca_removed_merge45$cca_cluster) %in% clusters_to_merge)) {
  stop("Clusters 4 or 5 remain as separate curated identities in Object 2.")
}

n_merged_after <- sum(
  as.character(cca_removed_merge45$cca_cluster) == merged_cluster_name
)

if (n_merged_after != n_cluster4_before + n_cluster5_before) {
  stop("The merged 4_5 cluster has an unexpected number of cells.")
}

if (!identical(colnames(cca_removed), colnames(cca_removed_merge45))) {
  stop("Object 1 and Object 2 do not contain the same retained cells.")
}

merged_file <- file.path(
  output_dir,
  "integrated_CCA_removed_8_9_10_12_14_merged_4_5.qs"
)

qs::qsave(
  cca_removed_merge45,
  merged_file
)

message(
  "Saved Object 2 (removed clusters plus merged 4 and 5): ",
  normalizePath(merged_file, winslash = "/", mustWork = TRUE)
)

############################################################
# 7) VERIFY AND SAVE CELL-COUNT TABLES
############################################################

make_cluster_count_table <- function(object, version_name) {
  tab <- as.data.frame(
    table(
      cluster = object$cca_cluster,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(tab)[2] <- "n_cells"
  tab$version <- version_name
  tab
}

removed_counts <- make_cluster_count_table(
  cca_removed,
  "removed_only"
)

merged_counts <- make_cluster_count_table(
  cca_removed_merge45,
  "removed_plus_merged_4_5"
)

cluster_count_comparison <- rbind(
  source_cluster_counts[, c("version", "cluster", "n_cells")],
  removed_counts[, c("version", "cluster", "n_cells")],
  merged_counts[, c("version", "cluster", "n_cells")]
)

write.csv(
  cluster_count_comparison,
  file.path(output_dir, "cluster_count_comparison.csv"),
  row.names = FALSE
)

message("\nObject 1 cluster counts:")
print(table(cca_removed$cca_cluster))

message("\nObject 2 cluster counts:")
print(table(cca_removed_merge45$cca_cluster))

message(
  "\nCells removed: ",
  length(removed_cells),
  "; cells retained in each curated object: ",
  ncol(cca_removed)
)

message(
  "Cluster 4 cells before merge: ", n_cluster4_before,
  "; cluster 5 cells before merge: ", n_cluster5_before,
  "; merged 4_5 cells: ", n_merged_after
)

############################################################
# 8) SAMPLE-BY-CLUSTER TABLES
############################################################

if ("sample" %in% colnames(cca_removed[[]])) {
  
  removed_by_sample <- as.data.frame(
    table(
      sample = cca_removed$sample,
      cluster = cca_removed$cca_cluster
    ),
    stringsAsFactors = FALSE
  )
  
  merged_by_sample <- as.data.frame(
    table(
      sample = cca_removed_merge45$sample,
      cluster = cca_removed_merge45$cca_cluster
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    removed_by_sample,
    file.path(output_dir, "removed_only_cluster_counts_by_sample.csv"),
    row.names = FALSE
  )
  
  write.csv(
    merged_by_sample,
    file.path(output_dir, "removed_plus_merged_4_5_counts_by_sample.csv"),
    row.names = FALSE
  )
}

############################################################
# 9) UMAP PLOTS
############################################################

if ("umap" %in% names(cca_removed@reductions)) {
  
  p_removed <- Seurat::DimPlot(
    object = cca_removed,
    reduction = "umap",
    group.by = "cca_cluster",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("CCA: clusters 8, 9, 10, 12, and 14 removed") +
    Seurat::NoLegend()
  
  p_merged <- Seurat::DimPlot(
    object = cca_removed_merge45,
    reduction = "umap",
    group.by = "cca_cluster",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("CCA: removed clusters; clusters 4 and 5 merged") +
    Seurat::NoLegend()
  
  ggsave(
    filename = file.path(plot_dir, "cca_removed_clusters_only.pdf"),
    plot = p_removed,
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = file.path(plot_dir, "cca_removed_clusters_merged_4_5.pdf"),
    plot = p_merged,
    width = 8,
    height = 6
  )
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_compare <- patchwork::wrap_plots(
      p_removed,
      p_merged,
      ncol = 2
    )
    
    ggsave(
      filename = file.path(plot_dir, "cca_curated_objects_side_by_side.pdf"),
      plot = p_compare,
      width = 16,
      height = 6
    )
  }
  
  print(p_removed)
  print(p_merged)
  
} else {
  warning("No UMAP reduction was found; UMAP plots were skipped.")
}

############################################################
# 10) FINAL FILE CHECK
############################################################

if (!file.exists(removed_file) || !file.exists(merged_file)) {
  stop("One or both curated QS files were not written successfully.")
}

message("\nDone. The untouched source QS was not overwritten.")
message("Object 1: ", normalizePath(removed_file, winslash = "/"))
message("Object 2: ", normalizePath(merged_file, winslash = "/"))
message("Output folder: ", normalizePath(output_dir, winslash = "/"))

############################################################
# OBJECTS AVAILABLE IN THE CURRENT R SESSION
############################################################
# cca_removed
#   - clusters 8, 9, 10, 12, and 14 removed
#   - clusters 4 and 5 remain separate
#
# cca_removed_merge45
#   - the same cells as cca_removed
#   - clusters 4 and 5 relabeled as 4_5
#
# The original cluster labels are retained in:
#   cca_cluster_original
#
# In the merged object, the premerge labels are also retained in:
#   cca_cluster_premerge
############################################################
############################################################
# MANUAL CELL-TYPE ANNOTATION FOR TWO CURATED CCA OBJECTS
#
# Object 1: clusters 8, 9, 10, 12, 14 removed; clusters 4 and 5 merged
# Object 2: clusters 8, 9, 10, 12, 14 removed; clusters 4 and 5 separate
#
# This script adds metadata only. It does not rerun CCA, PCA,
# neighbors, clustering, or UMAP.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
})

############################################################
# 1) INPUT DIRECTORY
############################################################

# Edit this directory if your curated QS files are elsewhere.
input_dir <- paste0(
  "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq/",
  "output_new_process_scRNAseq/curated_cluster_objects"
)

merged_file <- file.path(
  input_dir,
  "integrated_CCA_removed_8_9_10_12_14_merged_4_5.qs"
)

separate_file <- file.path(
  input_dir,
  "integrated_CCA_removed_clusters_8_9_10_12_14.qs"
)

# Fall back to interactive selection when a path is not found.
if (!file.exists(merged_file)) {
  message(
    "Merged 4_5 file was not found at the expected path. ",
    "Select integrated_CCA_removed_8_9_10_12_14_merged_4_5.qs"
  )
  merged_file <- file.choose()
}

if (!file.exists(separate_file)) {
  message(
    "Separate 4/5 file was not found at the expected path. ",
    "Select integrated_CCA_removed_clusters_8_9_10_12_14.qs"
  )
  separate_file <- file.choose()
}

############################################################
# 2) CONCISE AND MARKER-RICH LABELS: MERGED 4_5 OBJECT
############################################################

merged_cell_type <- c(
  "0"   = "L6 IT ExN",
  "1"   = "Migrating cortical ExN",
  "2"   = "Perivascular mesenchymal",
  "3"   = "Mid-layer IT ExN",
  "4_5" = "GABAergic neurons",
  "6"   = "Early OPC",
  "7"   = "Immature upper-layer ExN",
  "11"  = "L6 CT-like ExN",
  "13"  = "Developmental astrocytes",
  "15"  = "L5 deep-layer ExN"
)

merged_cell_type_markers <- c(
  "0"   = "STK32B+/CDH13+ L6 IT ExN",
  "1"   = "DPY19L1+/SEMA3C+ migrating cortical ExN",
  "2"   = "PDGFRB+/PDE3A+ perivascular mesenchymal cells",
  "3"   = "RORB+ mid-layer IT-like ExN",
  "4_5" = "GAD1+/DLX6-AS1+ GABAergic neurons",
  "6"   = "OLIG1+ early OPC",
  "7"   = "UNC5D+ immature upper-layer ExN",
  "11"  = "FOXP2+/TLE4+ L6 corticothalamic-like ExN",
  "13"  = "AQP4+/SLC4A4+ developmental astrocytes",
  "15"  = "RELN+ L5 deep-layer ExN"
)

############################################################
# 3) CONCISE AND MARKER-RICH LABELS: SEPARATE 4/5 OBJECT
############################################################

separate_cell_type <- c(
  "0"  = "L6 IT ExN",
  "1"  = "Migrating cortical ExN",
  "2"  = "Perivascular mesenchymal",
  "3"  = "Mid-layer IT ExN",
  "4"  = "Developing GABAergic neurons",
  "5"  = "Immature GABAergic neurons",
  "6"  = "Early OPC",
  "7"  = "Immature upper-layer ExN",
  "11" = "L6 CT-like ExN",
  "13" = "Developmental astrocytes",
  "15" = "L5 deep-layer ExN"
)

separate_cell_type_markers <- c(
  "0"  = "STK32B+/CDH13+ L6 IT ExN",
  "1"  = "DPY19L1+/SEMA3C+ migrating cortical ExN",
  "2"  = "PDGFRB+/PDE3A+ perivascular mesenchymal cells",
  "3"  = "RORB+ mid-layer IT-like ExN",
  "4"  = "DLX6-AS1+/ST18+ developing GABAergic neurons",
  "5"  = "PAX6+/ARX+ local immature GABAergic neurons",
  "6"  = "OLIG1+ early OPC",
  "7"  = "UNC5D+ immature upper-layer ExN",
  "11" = "FOXP2+/TLE4+ L6 corticothalamic-like ExN",
  "13" = "AQP4+/SLC4A4+ developmental astrocytes",
  "15" = "RELN+ L5 deep-layer ExN"
)

############################################################
# 4) HELPER: ANNOTATE, CHECK, PLOT, AND SAVE ONE OBJECT
############################################################

annotate_and_save <- function(
    input_file,
    cell_type_map,
    marker_label_map,
    output_stub,
    plot_title
) {
  
  message("\nLoading: ", input_file)
  obj <- qs::qread(input_file)
  
  if (!inherits(obj, "Seurat")) {
    stop("The selected QS file does not contain a Seurat object: ", input_file)
  }
  
  metadata_columns <- colnames(obj@meta.data)
  
  cluster_column <- if ("cca_cluster" %in% metadata_columns) {
    "cca_cluster"
  } else if ("seurat_clusters" %in% metadata_columns) {
    "seurat_clusters"
  } else {
    stop(
      "Neither 'cca_cluster' nor 'seurat_clusters' was found in: ",
      input_file
    )
  }
  
  cluster_id <- as.character(obj@meta.data[[cluster_column]])
  observed_clusters <- unique(cluster_id)
  
  missing_annotations <- setdiff(
    observed_clusters,
    names(cell_type_map)
  )
  
  if (length(missing_annotations) > 0L) {
    stop(
      "No annotation was supplied for cluster(s): ",
      paste(missing_annotations, collapse = ", "),
      ". Clusters present in the object: ",
      paste(sort(observed_clusters), collapse = ", ")
    )
  }
  
  missing_marker_labels <- setdiff(
    observed_clusters,
    names(marker_label_map)
  )
  
  if (length(missing_marker_labels) > 0L) {
    stop(
      "No marker-rich label was supplied for cluster(s): ",
      paste(missing_marker_labels, collapse = ", ")
    )
  }
  
  cluster_order <- names(cell_type_map)[
    names(cell_type_map) %in% observed_clusters
  ]
  
  # Preserve the cluster ID used for annotation.
  obj$cca_cluster_id <- factor(
    cluster_id,
    levels = cluster_order
  )
  
  # Concise label for UMAPs and downstream grouping.
  obj$cell_type <- unname(
    cell_type_map[cluster_id]
  )
  
  # Detailed label retaining the defining marker genes.
  obj$cell_type_markers <- unname(
    marker_label_map[cluster_id]
  )
  
  # Cluster number plus concise label is useful for review.
  obj$cluster_cell_type <- paste0(
    "C",
    cluster_id,
    ": ",
    obj$cell_type
  )
  
  concise_order <- unname(
    cell_type_map[cluster_order]
  )
  
  marker_order <- unname(
    marker_label_map[cluster_order]
  )
  
  cluster_cell_type_order <- paste0(
    "C",
    cluster_order,
    ": ",
    concise_order
  )
  
  obj$cell_type <- factor(
    obj$cell_type,
    levels = concise_order
  )
  
  obj$cell_type_markers <- factor(
    obj$cell_type_markers,
    levels = marker_order
  )
  
  obj$cluster_cell_type <- factor(
    obj$cluster_cell_type,
    levels = cluster_cell_type_order
  )
  
  obj$annotation_method <- "manual_marker_based"
  obj$annotation_version <- "v1"
  
  # Set concise cell types as the active identity without changing
  # the original cluster metadata columns.
  SeuratObject::Idents(obj) <- "cell_type"
  
  if (anyNA(obj$cell_type)) {
    stop("At least one cell received an NA cell-type label.")
  }
  
  message("Annotation check:")
  print(
    table(
      cluster = obj$cca_cluster_id,
      cell_type = obj$cell_type
    )
  )
  
  output_dir <- file.path(
    dirname(input_file),
    "cell_type_annotations"
  )
  
  plot_dir <- file.path(
    output_dir,
    "plots"
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Annotation table and cell counts.
  annotation_table <- data.frame(
    cluster = cluster_order,
    cell_type = concise_order,
    cell_type_markers = marker_order,
    n_cells = as.integer(
      table(
        factor(
          cluster_id,
          levels = cluster_order
        )
      )
    ),
    stringsAsFactors = FALSE
  )
  
  if ("sample" %in% colnames(obj@meta.data)) {
    
    sample_order <- unique(
      as.character(obj$sample)
    )
    
    cluster_by_sample <- table(
      cluster = factor(
        cluster_id,
        levels = cluster_order
      ),
      sample = factor(
        as.character(obj$sample),
        levels = sample_order
      )
    )
    
    sample_count_table <- as.data.frame.matrix(
      cluster_by_sample
    )
    
    sample_count_table$cluster <- rownames(
      sample_count_table
    )
    
    rownames(sample_count_table) <- NULL
    
    sample_count_table <- sample_count_table[
      match(
        cluster_order,
        sample_count_table$cluster
      ),
      ,
      drop = FALSE
    ]
    
    annotation_table <- merge(
      annotation_table,
      sample_count_table,
      by = "cluster",
      all.x = TRUE,
      sort = FALSE
    )
    
    annotation_table <- annotation_table[
      match(
        cluster_order,
        annotation_table$cluster
      ),
      ,
      drop = FALSE
    ]
  }
  
  annotation_csv <- file.path(
    output_dir,
    paste0(output_stub, "_cell_type_annotation_table.csv")
  )
  
  write.csv(
    annotation_table,
    annotation_csv,
    row.names = FALSE
  )
  
  # Save cell-level metadata for traceability.
  cell_metadata <- data.frame(
    cell = colnames(obj),
    cluster = as.character(obj$cca_cluster_id),
    cell_type = as.character(obj$cell_type),
    cell_type_markers = as.character(obj$cell_type_markers),
    stringsAsFactors = FALSE
  )
  
  if ("sample" %in% colnames(obj@meta.data)) {
    cell_metadata$sample <- as.character(obj$sample)
  }
  
  write.csv(
    cell_metadata,
    file.path(
      output_dir,
      paste0(output_stub, "_cell_level_annotations.csv")
    ),
    row.names = FALSE
  )
  
  # UMAP plots use the existing saved UMAP coordinates.
  if ("umap" %in% names(obj@reductions)) {
    
    p_cell_type <- Seurat::DimPlot(
      object = obj,
      reduction = "umap",
      group.by = "cell_type",
      label = TRUE,
      repel = TRUE,
      label.size = 3,
      pt.size = 0.25
    ) +
      ggplot2::ggtitle(plot_title) +
      Seurat::NoLegend()
    
    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(output_stub, "_umap_cell_types.pdf")
      ),
      plot = p_cell_type,
      width = 10,
      height = 7
    )
    
    p_cluster_cell_type <- Seurat::DimPlot(
      object = obj,
      reduction = "umap",
      group.by = "cluster_cell_type",
      label = TRUE,
      repel = TRUE,
      label.size = 3,
      pt.size = 0.25
    ) +
      ggplot2::ggtitle(
        paste0(plot_title, " (cluster IDs shown)")
      ) +
      Seurat::NoLegend()
    
    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(output_stub, "_umap_cluster_cell_types.pdf")
      ),
      plot = p_cluster_cell_type,
      width = 11,
      height = 7
    )
    
    if ("sample" %in% colnames(obj@meta.data)) {
      
      p_by_sample <- Seurat::DimPlot(
        object = obj,
        reduction = "umap",
        group.by = "cell_type",
        split.by = "sample",
        ncol = 3,
        pt.size = 0.20
      ) +
        ggplot2::ggtitle(
          paste0(plot_title, " by sample")
        )
      
      ggplot2::ggsave(
        filename = file.path(
          plot_dir,
          paste0(output_stub, "_umap_cell_types_by_sample.pdf")
        ),
        plot = p_by_sample,
        width = 18,
        height = 6
      )
    }
    
  } else {
    warning(
      "No UMAP reduction was found in ",
      basename(input_file),
      "; annotation was saved without UMAP plots."
    )
  }
  
  output_qs <- file.path(
    output_dir,
    paste0(output_stub, "_celltypes_annotated.qs")
  )
  
  qs::qsave(
    obj,
    output_qs
  )
  
  message("Saved annotated object: ", output_qs)
  message("Saved annotation table: ", annotation_csv)
  
  list(
    annotated_qs = output_qs,
    annotation_csv = annotation_csv,
    output_dir = output_dir
  )
}

############################################################
# 5) ANNOTATE AND SAVE THE MERGED 4_5 OBJECT
############################################################

merged_outputs <- annotate_and_save(
  input_file = merged_file,
  cell_type_map = merged_cell_type,
  marker_label_map = merged_cell_type_markers,
  output_stub = "integrated_CCA_removed_8_9_10_12_14_merged_4_5",
  plot_title = "CCA cell types: clusters 4 and 5 merged"
)

gc()

############################################################
# 6) ANNOTATE AND SAVE THE SEPARATE 4/5 OBJECT
############################################################

separate_outputs <- annotate_and_save(
  input_file = separate_file,
  cell_type_map = separate_cell_type,
  marker_label_map = separate_cell_type_markers,
  output_stub = "integrated_CCA_removed_clusters_8_9_10_12_14",
  plot_title = "CCA cell types: clusters 4 and 5 separate"
)

gc()

############################################################
# 7) OUTPUT LOCATIONS
############################################################

message("\nMerged-object outputs:")
print(merged_outputs)

message("\nSeparate-object outputs:")
print(separate_outputs)

message("\nAnnotation complete.")

########################################################################################

neuron_markers <- list(
  
  L6_IT_ExN = c("STK32B","CDH13","MCTP1","DOK5","NEO1"), #STK32B⁺/CDH13⁺ deep-layer intratelencephalic-like excitatory neurons #STK32B⁺/CDH13⁺ L6 IT-like excitatory neurons
  
  Migrating_Cortical_ExN = c("DPY19L1","PRSS12", "PLXNA2", "SEMA3C","PTPN4","FRMD4B", "ST3GAL6"), #DPY19L1⁺/SEMA3C⁺ migrating cortical glutamatergic neurons
  
  Perivascular_Mesenchymal = c("PDGFRB","PDE3A","TNC","FBN1","COL11A1","PTN","PTPRM","IL33", "GLI3","LIPG","SLCO1C1", "INTU","DYNC2H1"), #PDGFRB⁺/PDE3A⁺ perivascular mesenchymal cells
  
  Mid_layer_IT_ExN = c("LRRC4C", "RORB","GRM7","KCNH5","NTNG1","CDH18","TENM1","LUZP2"), #RORB⁺ mid-layer excitatory neurons #RORB⁺ mid-layer IT-like excitatory neurons #RORB⁺/KCNH5⁺ mid-layer excitatory neurons
  
  Developing_GABAergic_neurons = c("GAD1", "PDZRN3", "NRXN3","GAD2","DLX6-AS1", "ST18", "DCLK2"), #DLX6-AS1⁺/ST18⁺ developing inhibitory neurons
  
  Immature_GABAergic_neurons = c("COPG2IT1", "PAX6","MEST", "ST18", "DLX6-AS1","GAD2", "ARX", "NRXN3","GAD1","SOX2-OT"), #PAX6⁺/ARX⁺ local immature GABAergic neurons
  
  Early_OPC = c("OLIG1", "FERMT1","SMOC1", "PDGFRA", "PLP1", "SOX10", "CLDN11", "MBP",
          "DLL3", "CSPG4", "SLC24A3", "CEROX1", "MAP3K1", "POLR2F"), #OLIG1⁺ early OPC, #OLIG1⁺ early oligodendroglial progenitors
  
  Immature_upper_layer_ExN =  c("UNC5D", "SOX11","DCC", "NRP1",  "LRP8","EPHA3",  "DOK6","KCNQ3","CNR1","CDH4"), #UNC5D⁺ immature upper-layer excitatory neurons
  
  L6_CT_like_ExN = c("HS3ST4","TRPM3", "TLE4", "FOXP2","ZFPM2", "SOX5","ADAMTSL1", "SEMA3E"), #FOXP2⁺/TLE4⁺ L6 corticothalamic-like neurons
  
  Developmental_astrocytes = c("AQP4", "SLC4A4", "GFAP", "MEGF10", "SPARCL1", "CD44", "GLIS3", "RFX4", "ID4", "EDNRB", "LIFR",
           "LAMA1", "ITGA6", "NID1"), #AQP4⁺/SLC4A4⁺ developmental astrocytes #endfoot-like astrocytes
  
  L5_deep_layer_ExN = c("SGCZ", "PEX5L","RALYL", "CD36","NWD2","ST6GALNAC3","RELN",
           "HCN1","RYR2", "KCNIP4","CACNA2D3", "CADPS2","GRIK2","GRM8","GRM5") #RELN⁺ L5 deep-layer excitatory neurons
  
)      

####################################################################################

