############################################################
setwd("H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2/")

############################################################
# DESCRIPTIVE FIBRILS EFFECT, FIBJ8 REVERSAL,
# AND CROSS-STUDY AD CONCORDANCE ANALYSIS
#
# Dataset:
#   obj2_celltypes_2026_v3.rds
#
# Metadata:
#   cell type column: celltype_final
#   condition labels: Control, Fibrils, FibJ8
#
# Contrasts retained throughout this script:
#   1) Fibrils vs Control
#   2) FibJ8 vs Fibrils
#
# The FibJ8 vs Control contrast is intentionally NOT calculated,
# saved, ranked, summarized, or plotted.
#
# IMPORTANT STATISTICAL LIMITATION
# --------------------------------
# This script does not treat individual cells as biological replicates.
# It does not calculate p-values for the user's Control/Fibrils/FibJ8
# data when there is one biological library per condition.
# Raw RNA counts are summed within each library and selected cell
# population, then library-size-normalized log2-CPM effect sizes are shown.
# With n = 1 library per condition, replicate-aware DE, mixed models,
# and biological-variance estimates cannot be fitted.
#
# Leng and Morabito spreadsheets are used as independent published
# reference summaries. Cross-study results are descriptive concordance
# analyses, not proof of treatment efficacy and not a substitute for
# biological replication.
############################################################

options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 8 * 1024^3)
set.seed(1234)

############################################################
# 0) PACKAGES
############################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "readxl",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readxl)
  library(patchwork)
  library(scales)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_qs <- requireNamespace("qs", quietly = TRUE)

############################################################
# 1) USER SETTINGS
############################################################

# Change this one line only if the Windows drive or parent directory differs.
project_dir <- "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2"

# The script uses these files directly when they exist. If a path is not
# found, a file chooser opens so the analysis can still proceed.
object_file <- file.path(project_dir, "obj2_celltypes_2026_v3.rds")
leng_file <- file.path(project_dir, "Leng_EX_Braak6_vs_0.xlsx")
morabito_file <- file.path(project_dir, "Morabito_EX_AD_vs_Control.xlsx")

# Appropriate default for comparison with excitatory-neuron external studies.
# Options: "all_excitatory", "single_cell_type", "all_cells"
analysis_mode <- "all_excitatory"

# Used only when analysis_mode == "single_cell_type".
# The supplied UMAP uses the broad excitatory label "ExN".
single_cell_type <- "ExN"

# Exact cell-type metadata column requested for this dataset.
annotation_col <- "celltype_final"

# The first existing condition column is used.
condition_candidates <- c(
  "condition",
  "Condition",
  "sample",
  "Sample",
  "group",
  "treatment",
  "orig.ident"
)

# Exact normalized condition order used throughout the script.
condition_levels <- c(
  "Control",
  "Fibrils",
  "FibJ8"
)

#################################################################

############################################################
# FIGURE 5a-d: DATA-DERIVED EXPRESSION ANALYSIS WITH HGNC GENE FILTERS
#
# All expression differences, thresholds, quadrant assignments, rankings,
# and plotted values are calculated only from the input CSV/TSV/XLSX/XLS
# table. No values are copied from the reference figure. Official HGNC
# annotation is used only to classify protein-coding genes and to exclude
# LINC genes, pseudogenes, and chromosome-Y genes.
#
# Definitions calculated directly from the sheet:
#   delta_F = log2CPM_Fibrils - log2CPM_Control
#   delta_J = log2CPM_FibJ8   - log2CPM_Fibrils
#
# Analysis set used in panels a-d:
#   genes in the selected gene scope with abs(delta_F) >= 0.20 and
#   abs(delta_J) >= 0.20.
#
# GENE-SCOPE RULE:
#   gene_scope_mode = "protein_coding" retains HGNC loci classified as
#   protein-coding genes. The requested ENSG-only, LINC, pseudogene, and
#   chromosome-Y exclusions are then applied before the effect thresholds.
############################################################

options(stringsAsFactors = FALSE)
set.seed(1234)

############################################################
# 1) USER SETTINGS
############################################################

project_dir <- "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2"

# CSV/TSV/XLSX/XLS are supported. Change this to your Excel filename if needed.
input_file <- file.path(
  project_dir,
  "all_tested_genes_fibrils_crossstudy_fibj8_reversal.csv"
)
excel_sheet <- 1

output_dir <- file.path(project_dir, "Figure5_filtered_protein_coding_abcd")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Exact effect threshold requested for the Figure 5 analysis set.
effect_threshold <- 0.10

# Data-driven labels in panel a. Set to 0 to suppress gene labels.
panel_a_labels_per_quadrant <- 10L

# Number of Q2 genes shown in panel d.
top_q2_n <- 20L

# Gene-scope behavior:
#   "protein_coding" = use HGNC locus annotation to retain protein-coding genes.
#   "auto"           = same as protein_coding when HGNC annotation is available.
#   "all_rows"       = use all valid rows after the requested exclusions.
# The default preserves the requested panel-5b label: "Protein-coding genes".
gene_scope_mode <- "protein_coding"

# Requested exclusions. These are applied before the |delta| thresholds and
# therefore affect every panel (5a-5d).
exclude_unannotated_ensembl_ids <- TRUE
exclude_lincRNAs <- TRUE
exclude_pseudogenes <- TRUE
exclude_chrY_genes <- TRUE

# HGNC annotation is downloaded once and then reused from this local cache.
# If your workstation cannot access the internet, download this file manually
# and save it at hgnc_annotation_file before running the script.
hgnc_annotation_file <- file.path(project_dir, "hgnc_complete_set.txt")
hgnc_annotation_url <- paste0(
  "https://storage.googleapis.com/public-download-files/hgnc/",
  "tsv/tsv/hgnc_complete_set.txt"
)
download_hgnc_if_missing <- TRUE

figure_width_in <- 12
figure_height_in <- 9
figure_dpi <- 600

############################################################
# 2) PACKAGES
############################################################

required_packages <- c(
  "dplyr", "tibble", "ggplot2", "readr", "readxl",
  "patchwork", "ggrepel", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

############################################################
# 3) HELPERS
############################################################

read_gene_table <- function(path, sheet = 1) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  
  extension <- tolower(tools::file_ext(path))
  
  if (extension == "csv") {
    return(
      readr::read_csv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("tsv", "txt")) {
    return(
      readr::read_tsv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("xlsx", "xls")) {
    return(
      readxl::read_excel(
        path,
        sheet = sheet,
        na = c("", "NA", "NaN", "nan")
      )
    )
  }
  
  stop("Unsupported input extension: .", extension)
}

load_hgnc_annotation <- function(
    cache_file,
    source_url,
    download_if_missing = TRUE
) {
  if (!file.exists(cache_file)) {
    if (!isTRUE(download_if_missing)) {
      stop(
        "HGNC annotation file not found: ", cache_file,
        ". Set download_hgnc_if_missing <- TRUE or place the HGNC complete ",
        "set at that location."
      )
    }
    
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    message("Downloading HGNC annotation to: ", cache_file)
    
    download_result <- tryCatch(
      utils::download.file(
        url = source_url,
        destfile = cache_file,
        mode = "wb",
        quiet = FALSE
      ),
      error = function(e) e
    )
    
    if (
      inherits(download_result, "error") ||
      !file.exists(cache_file) ||
      is.na(file.info(cache_file)$size) ||
      file.info(cache_file)$size < 1000
    ) {
      if (file.exists(cache_file)) {
        unlink(cache_file)
      }
      
      error_message <- if (inherits(download_result, "error")) {
        conditionMessage(download_result)
      } else {
        "the downloaded file was absent or unexpectedly small"
      }
      
      stop(
        "Could not download the HGNC annotation: ", error_message, "\n",
        "Download it manually from:\n", source_url, "\n",
        "and save it as:\n", cache_file
      )
    }
  }
  
  hgnc_raw <- readr::read_tsv(
    cache_file,
    na = c("", "NA", "NaN", "nan"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  required_hgnc_columns <- c(
    "symbol", "name", "locus_group", "locus_type", "location"
  )
  missing_hgnc_columns <- setdiff(required_hgnc_columns, names(hgnc_raw))
  
  if (length(missing_hgnc_columns) > 0L) {
    stop(
      "The HGNC annotation file is missing required columns: ",
      paste(missing_hgnc_columns, collapse = ", ")
    )
  }
  
  hgnc_raw %>%
    transmute(
      gene_key = toupper(trimws(as.character(symbol))),
      hgnc_symbol = trimws(as.character(symbol)),
      hgnc_name = trimws(as.character(name)),
      hgnc_locus_group = tolower(trimws(as.character(locus_group))),
      hgnc_locus_type = tolower(trimws(as.character(locus_type))),
      hgnc_location = trimws(as.character(location))
    ) %>%
    filter(!is.na(gene_key), nzchar(gene_key)) %>%
    distinct(gene_key, .keep_all = TRUE)
}

blank_panel <- function(title, message) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message,
      size = 4,
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.margin = margin(8, 8, 8, 8)
    )
}

############################################################
# 4) READ, ANNOTATE, AND CALCULATE DESCRIPTIVE DIFFERENCES
############################################################

raw_gene_table <- read_gene_table(input_file, sheet = excel_sheet) %>%
  as_tibble()

required_columns <- c(
  "gene",
  "log2cpm_Control",
  "log2cpm_Fibrils",
  "log2cpm_FibJ8"
)

missing_columns <- setdiff(required_columns, names(raw_gene_table))
if (length(missing_columns) > 0L) {
  stop(
    "The input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

valid_gene_scope_modes <- c("auto", "protein_coding", "all_rows")
if (!gene_scope_mode %in% valid_gene_scope_modes) {
  stop(
    "gene_scope_mode must be one of: ",
    paste(valid_gene_scope_modes, collapse = ", ")
  )
}

hgnc_annotation <- load_hgnc_annotation(
  cache_file = hgnc_annotation_file,
  source_url = hgnc_annotation_url,
  download_if_missing = download_hgnc_if_missing
)

prepared_gene_table <- raw_gene_table %>%
  mutate(
    gene = trimws(as.character(gene)),
    gene_key = toupper(gene),
    log2cpm_Control = suppressWarnings(as.numeric(log2cpm_Control)),
    log2cpm_Fibrils = suppressWarnings(as.numeric(log2cpm_Fibrils)),
    log2cpm_FibJ8 = suppressWarnings(as.numeric(log2cpm_FibJ8))
  ) %>%
  left_join(hgnc_annotation, by = "gene_key") %>%
  mutate(
    delta_F = log2cpm_Fibrils - log2cpm_Control,
    delta_J = log2cpm_FibJ8 - log2cpm_Fibrils,
    valid_gene_identifier = !is.na(gene) & nzchar(gene),
    finite_differences = is.finite(delta_F) & is.finite(delta_J),
    hgnc_mapped = !is.na(hgnc_symbol) & nzchar(hgnc_symbol),
    
    # Exact ENSG-only rows are considered unannotated identifiers here.
    is_unannotated_ensembl_id = grepl(
      "^ENSG[0-9]+([.][0-9]+)?$",
      gene,
      ignore.case = TRUE
    ),
    
    # LINC detection is conservative: approved LINC symbols and HGNC names
    # explicitly described as long intergenic non-protein-coding RNA.
    is_lincRNA = grepl("^LINC", gene, ignore.case = TRUE) |
      (
        !is.na(hgnc_name) &
          grepl(
            "long intergenic non[- ]protein coding RNA",
            hgnc_name,
            ignore.case = TRUE
          )
      ),
    
    # HGNC locus classification avoids unreliable suffix rules such as P/P1,
    # which would incorrectly remove legitimate protein-coding symbols.
    is_pseudogene = (
      !is.na(hgnc_locus_group) & hgnc_locus_group == "pseudogene"
    ) | (
      !is.na(hgnc_locus_type) &
        grepl("pseudogene", hgnc_locus_type, ignore.case = TRUE)
    ),
    
    # A location beginning with Y denotes a Y-chromosome locus. X/Y
    # pseudoautosomal entries are not classified as Y-specific by this rule.
    is_chrY_gene = !is.na(hgnc_location) &
      grepl("^Y($|p|q)", hgnc_location, ignore.case = TRUE),
    
    is_hgnc_protein_coding = (
      !is.na(hgnc_locus_group) &
        hgnc_locus_group == "protein-coding gene"
    ) | (
      !is.na(hgnc_locus_type) &
        hgnc_locus_type == "gene with protein product"
    )
  )

if (gene_scope_mode %in% c("auto", "protein_coding")) {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = is_hgnc_protein_coding)
  gene_scope_label <- "Protein-coding genes"
  gene_scope_source <- "HGNC locus_group/locus_type"
} else {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = TRUE)
  gene_scope_label <- "All genes after requested exclusions"
  gene_scope_source <- "all valid rows in the input sheet"
}

prepared_gene_table <- prepared_gene_table %>%
  mutate(
    excluded_unannotated_ensembl =
      isTRUE(exclude_unannotated_ensembl_ids) &
      is_unannotated_ensembl_id,
    excluded_lincRNA = isTRUE(exclude_lincRNAs) & is_lincRNA,
    excluded_pseudogene = isTRUE(exclude_pseudogenes) & is_pseudogene,
    excluded_chrY = isTRUE(exclude_chrY_genes) & is_chrY_gene,
    passes_requested_exclusions = !(
      excluded_unannotated_ensembl |
        excluded_lincRNA |
        excluded_pseudogene |
        excluded_chrY
    ),
    exclusion_reason = paste0(
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    exclusion_reason = sub("; $", "", exclusion_reason),
    exclusion_reason = na_if(exclusion_reason, ""),
    analysis_exclusion_reason = paste0(
      if_else(!is_in_gene_scope, "not HGNC protein-coding; ", ""),
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    analysis_exclusion_reason = sub(
      "; $", "", analysis_exclusion_reason
    ),
    analysis_exclusion_reason = na_if(analysis_exclusion_reason, "")
  )

valid_gene_mask <- prepared_gene_table$valid_gene_identifier
finite_difference_mask <- prepared_gene_table$finite_differences

if (anyDuplicated(prepared_gene_table$gene[valid_gene_mask]) > 0L) {
  duplicate_genes <- unique(
    prepared_gene_table$gene[
      valid_gene_mask & duplicated(prepared_gene_table$gene)
    ]
  )
  
  stop(
    "The gene column contains duplicate identifiers. Resolve duplicates before ",
    "plotting. First duplicated genes: ",
    paste(utils::head(duplicate_genes, 10L), collapse = ", ")
  )
}

# Optional audit: compare newly calculated differences with any pre-existing
# effect columns, but never use those columns to select or rank genes.
delta_consistency <- tibble(
  comparison = character(),
  max_absolute_difference = numeric()
)

if ("fibrils_vs_control_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_F <- suppressWarnings(
    as.numeric(prepared_gene_table$fibrils_vs_control_log2fc)
  )
  max_difference_F <- suppressWarnings(
    max(abs(prepared_gene_table$delta_F - existing_delta_F), na.rm = TRUE)
  )
  if (!is.finite(max_difference_F)) {
    max_difference_F <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_F versus fibrils_vs_control_log2fc",
      max_absolute_difference = max_difference_F
    )
  )
}

if ("fibj8_vs_fibrils_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_J <- suppressWarnings(
    as.numeric(prepared_gene_table$fibj8_vs_fibrils_log2fc)
  )
  max_difference_J <- suppressWarnings(
    max(abs(prepared_gene_table$delta_J - existing_delta_J), na.rm = TRUE)
  )
  if (!is.finite(max_difference_J)) {
    max_difference_J <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_J versus fibj8_vs_fibrils_log2fc",
      max_absolute_difference = max_difference_J
    )
  )
}

if (nrow(delta_consistency) > 0L) {
  readr::write_csv(
    delta_consistency,
    file.path(output_dir, "Figure5_delta_consistency_audit.csv")
  )
}

############################################################
# 5) DEFINE THE FILTERED ANALYSIS SET FROM THE SHEET
############################################################

analysis_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    valid_gene_identifier,
    finite_differences,
    !is.na(is_in_gene_scope),
    is_in_gene_scope,
    passes_requested_exclusions,
    abs(delta_F) >= effect_threshold,
    abs(delta_J) >= effect_threshold
  ) %>%
  dplyr::mutate(
    quadrant = case_when(
      delta_F > 0 & delta_J > 0 ~ "Q1",
      delta_F < 0 & delta_J > 0 ~ "Q2",
      delta_F < 0 & delta_J < 0 ~ "Q3",
      delta_F > 0 & delta_J < 0 ~ "Q4",
      TRUE ~ NA_character_
    ),
    quadrant = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
    direction_class = if_else(
      quadrant %in% c("Q2", "Q4"),
      "Directional reversal (Q2/Q4)",
      "Same direction (Q1/Q3)"
    ),
    direction_class = factor(
      direction_class,
      levels = c(
        "Directional reversal (Q2/Q4)",
        "Same direction (Q1/Q3)"
      )
    ),
    fibrils_direction = if_else(
      delta_F < 0,
      "Down in Fibrils",
      "Up in Fibrils"
    ),
    fibrils_direction = factor(
      fibrils_direction,
      levels = c("Down in Fibrils", "Up in Fibrils")
    ),
    paired_effect = pmin(abs(delta_F), abs(delta_J)),
    total_effect = abs(delta_F) + abs(delta_J),
    reversal_strength = abs(delta_J)
  ) %>%
  dplyr::filter(!is.na(quadrant)) %>%
  arrange(quadrant, gene)

if (nrow(analysis_gene_table) == 0L) {
  stop(
    "No genes remain after the selected gene scope, requested exclusions, ",
    "and abs(delta_F) >= ", effect_threshold,
    " / abs(delta_J) >= ", effect_threshold, " thresholds."
  )
}

quadrant_text <- c(
  Q1 = "Q1: same direction (+,+)",
  Q2 = "Q2: reversal (-,+)",
  Q3 = "Q3: same direction (-,-)",
  Q4 = "Q4: reversal (+,-)"
)

quadrant_counts <- analysis_gene_table %>%
  dplyr::count(quadrant, name = "gene_count", .drop = FALSE) %>%
  dplyr::mutate(
    quadrant_label = factor(
      quadrant_text[as.character(quadrant)],
      levels = unname(quadrant_text)
    )
  )

base_valid_mask <-
  prepared_gene_table$valid_gene_identifier &
  prepared_gene_table$finite_differences
scope_mask <-
  base_valid_mask &
  !is.na(prepared_gene_table$is_in_gene_scope) &
  prepared_gene_table$is_in_gene_scope
post_exclusion_mask <-
  scope_mask & prepared_gene_table$passes_requested_exclusions
post_delta_F_mask <-
  post_exclusion_mask &
  abs(prepared_gene_table$delta_F) >= effect_threshold

filter_audit <- tibble(
  step = c(
    "Rows in input table",
    "Rows with nonempty identifiers and finite delta_F/delta_J",
    "Rows mapped to an approved HGNC symbol",
    "Unannotated ENSG-only IDs identified",
    "LINC RNAs identified",
    "Pseudogenes identified by HGNC",
    "Chromosome-Y genes identified by HGNC",
    paste0("Rows in selected scope: ", gene_scope_label),
    "Rows after requested gene-class exclusions",
    paste0("Rows after |delta_F| >= ", effect_threshold),
    paste0(
      "Final rows after |delta_F| >= ", effect_threshold,
      " and |delta_J| >= ", effect_threshold
    )
  ),
  n = c(
    nrow(prepared_gene_table),
    sum(base_valid_mask, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$hgnc_mapped, na.rm = TRUE),
    sum(
      base_valid_mask & prepared_gene_table$is_unannotated_ensembl_id,
      na.rm = TRUE
    ),
    sum(base_valid_mask & prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_chrY_gene, na.rm = TRUE),
    sum(scope_mask, na.rm = TRUE),
    sum(post_exclusion_mask, na.rm = TRUE),
    sum(post_delta_F_mask, na.rm = TRUE),
    nrow(analysis_gene_table)
  )
)

exclusion_summary <- tibble(
  exclusion = c(
    "Unannotated ENSG-only ID",
    "LINC RNA",
    "Pseudogene",
    "Chromosome Y"
  ),
  enabled = c(
    exclude_unannotated_ensembl_ids,
    exclude_lincRNAs,
    exclude_pseudogenes,
    exclude_chrY_genes
  ),
  n_identified = c(
    sum(prepared_gene_table$is_unannotated_ensembl_id, na.rm = TRUE),
    sum(prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$is_chrY_gene, na.rm = TRUE)
  ),
  n_flagged_for_exclusion = c(
    sum(prepared_gene_table$excluded_unannotated_ensembl, na.rm = TRUE),
    sum(prepared_gene_table$excluded_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$excluded_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$excluded_chrY, na.rm = TRUE)
  )
)

excluded_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    !is_in_gene_scope |
      !passes_requested_exclusions
  ) %>%
  dplyr::select(
    gene,
    hgnc_symbol,
    hgnc_name,
    hgnc_locus_group,
    hgnc_locus_type,
    hgnc_location,
    is_hgnc_protein_coding,
    is_unannotated_ensembl_id,
    is_lincRNA,
    is_pseudogene,
    is_chrY_gene,
    exclusion_reason,
    analysis_exclusion_reason
  ) %>%
  arrange(gene)

readr::write_csv(
  analysis_gene_table,
  file.path(output_dir, "Figure5_sheet_only_analysis_genes.csv")
)
readr::write_csv(
  quadrant_counts,
  file.path(output_dir, "Figure5_sheet_only_quadrant_counts.csv")
)
readr::write_csv(
  filter_audit,
  file.path(output_dir, "Figure5_sheet_only_filter_audit.csv")
)
readr::write_csv(
  exclusion_summary,
  file.path(output_dir, "Figure5_gene_exclusion_summary.csv")
)
readr::write_csv(
  excluded_gene_table,
  file.path(output_dir, "Figure5_excluded_genes_and_reasons.csv")
)

message("Gene scope: ", gene_scope_label)
message("Gene-scope source: ", gene_scope_source)
message("HGNC annotation cache: ", hgnc_annotation_file)
message("Final analysis genes: ", scales::comma(nrow(analysis_gene_table)))
message(
  "Quadrant counts: ",
  paste0(
    as.character(quadrant_counts$quadrant),
    "=",
    quadrant_counts$gene_count,
    collapse = ", "
  )
)

############################################################
# 6) COMMON STYLE
############################################################

delta_symbol <- "\u0394"
greater_equal_symbol <- "\u2265"
minus_symbol <- "\u2212"

reversal_blue <- "#2C7FB8"
same_direction_orange <- "#E67E39"

category_colors <- c(
  "Directional reversal (Q2/Q4)" = reversal_blue,
  "Same direction (Q1/Q3)" = same_direction_orange
)

quadrant_colors <- c(
  Q1 = "#F8766D",
  Q2 = "#7CAE00",
  Q3 = "#00BFC4",
  Q4 = "#C77CFF"
)

box_colors <- c(
  "Down in Fibrils" = "#70B96B",
  "Up in Fibrils" = "#B784C6"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5),
    plot.caption = element_text(hjust = 0, size = 7.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 8),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 7, 7, 7)
  )

############################################################
# 7) PANEL a: DESCRIPTIVE EFFECT-DIRECTION SCATTER
############################################################

if (panel_a_labels_per_quadrant > 0L) {
  panel_a_labels <- analysis_gene_table %>%
    group_by(quadrant) %>%
    arrange(
      desc(paired_effect),
      desc(total_effect),
      gene,
      .by_group = TRUE
    ) %>%
    slice_head(n = panel_a_labels_per_quadrant) %>%
    ungroup()
} else {
  panel_a_labels <- analysis_gene_table[0, , drop = FALSE]
}

readr::write_csv(
  panel_a_labels,
  file.path(output_dir, "Figure5_panel_a_data_selected_labels.csv")
)

panel_a <- ggplot(
  analysis_gene_table,
  aes(
    x = delta_F,
    y = delta_J,
    color = direction_class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_point(size = 1.15, alpha = 0.58) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.70
  ) +
  ggrepel::geom_text_repel(
    data = panel_a_labels,
    aes(label = gene),
    color = "black",
    size = 2.45,
    box.padding = 0.25,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.3,
    max.overlaps = Inf,
    seed = 1234,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Q2",
    hjust = -0.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Q1",
    hjust = 1.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    label = "Q3",
    hjust = -0.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = "Q4",
    hjust = 1.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Fibrils-associated expression and FibJ8 response",
    subtitle = "Descriptive log2CPM differences per condition (single pooled library each).",
    caption = paste0(
      "Labels: top ", panel_a_labels_per_quadrant,
      " genes per quadrant by min(|", delta_symbol, "F|, |", delta_symbol,
      "J|). Black line: descriptive least-squares trend."
    ),
    x = paste0(
      delta_symbol, "F = log2CPM(Fibrils) ", minus_symbol,
      " log2CPM(Control)"
    ),
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    ),
    color = "Effect direction"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.55, "cm"),
    plot.margin = margin(7, 18, 7, 7)
  )

############################################################
# 8) PANEL b: QUADRANT COUNTS FOR THE EXACT ANALYSIS SET
############################################################

panel_b_subtitle <- paste0(
  gene_scope_label, ", |", delta_symbol, "F| ",
  greater_equal_symbol, " ", sprintf("%.2f", effect_threshold),
  " & |", delta_symbol, "J| ", greater_equal_symbol, " ",
  sprintf("%.2f", effect_threshold)
)

panel_b <- ggplot(
  quadrant_counts,
  aes(x = quadrant_label, y = gene_count, fill = quadrant)
) +
  geom_col(width = 0.67) +
  geom_text(
    aes(label = scales::comma(gene_count)),
    vjust = -0.45,
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = quadrant_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = scales::comma
  ) +
  labs(
    title = "Effect-direction quadrant distribution",
    subtitle = panel_b_subtitle,
    caption = paste0(
      "Total genes in analysis set: ", scales::comma(nrow(analysis_gene_table)),
      ". ENSG-only IDs, LINC RNAs, pseudogenes, and chromosome-Y genes were excluded."
    ),
    x = NULL,
    y = "Gene count"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.1),
    plot.margin = margin(7, 7, 18, 7)
  )

############################################################
# 9) PANEL c: GENE-LEVEL DESCRIPTIVE DISTRIBUTION
############################################################

panel_c <- ggplot(
  analysis_gene_table,
  aes(
    x = fibrils_direction,
    y = delta_J,
    fill = fibrils_direction
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.45,
    alpha = 0.78
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1234),
    size = 0.50,
    alpha = 0.22,
    color = "black"
  ) +
  scale_fill_manual(values = box_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "FibJ8 response stratified by Fibrils-effect direction",
    subtitle = "Gene-level description; no inferential test performed.",
    caption = paste0(
      "Boxes and points summarize genes in the same threshold-qualified set: ",
      gene_scope_label, "."
    ),
    x = NULL,
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    )
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = -18, hjust = 0, vjust = 0.5, size = 7.5),
    plot.margin = margin(7, 12, 18, 7)
  )

############################################################
# 10) PANEL d: Q2 RANKED BY abs(delta_J)
############################################################

q2_gene_table <- analysis_gene_table %>%
  filter(quadrant == "Q2") %>%
  arrange(
    desc(reversal_strength),
    desc(abs(delta_F)),
    gene
  )

readr::write_csv(
  q2_gene_table,
  file.path(output_dir, "Figure5_all_Q2_genes_ranked_by_abs_delta_J.csv")
)

panel_d_data <- q2_gene_table %>%
  slice_head(n = top_q2_n) %>%
  arrange(reversal_strength, gene) %>%
  mutate(gene_plot = factor(gene, levels = gene))

readr::write_csv(
  panel_d_data,
  file.path(output_dir, "Figure5_top_Q2_by_abs_delta_J.csv")
)

panel_d_subtitle <- paste0(
  "Reversal strength = |", delta_symbol, "J| (|FibJ8 ",
  minus_symbol, " Fibrils|); Q2 = down in Fibrils, up with J8"
)

if (nrow(panel_d_data) == 0L) {
  panel_d <- blank_panel(
    "Top Q2 directional-reversal genes",
    paste0("No Q2 genes met the selected gene-scope and effect-threshold criteria: ", gene_scope_label, ".")
  )
} else {
  panel_d <- ggplot(
    panel_d_data,
    aes(x = reversal_strength, y = gene_plot)
  ) +
    geom_col(fill = reversal_blue, width = 0.86) +
    geom_text(
      aes(label = scales::number(reversal_strength, accuracy = 0.01)),
      hjust = -0.12,
      size = 2.5
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.16)),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    labs(
      title = "Top Q2 directional-reversal genes",
      subtitle = panel_d_subtitle,
      caption = paste0(
        "Top ", min(top_q2_n, nrow(q2_gene_table)),
        " of ", scales::comma(nrow(q2_gene_table)),
        " Q2 genes; ranking uses only |", delta_symbol, "J|."
      ),
      x = paste0("Reversal strength = |", delta_symbol, "J|"),
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y = element_text(size = 7.2),
      plot.margin = margin(7, 12, 7, 7)
    )
}

############################################################
# 11) COMBINE AND SAVE
############################################################

figure_5_abcd <- (
  panel_a + panel_b + plot_layout(widths = c(1.45, 1.00))
) / (
  panel_c + panel_d + plot_layout(widths = c(1.20, 1.00))
) +
  plot_layout(heights = c(1.00, 1.00)) +
  plot_annotation(
    title = "Figure 5",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 1, size = 13),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

combined_png <- file.path(output_dir, "Figure5_abcd_sheet_only.png")
combined_pdf <- file.path(output_dir, "Figure5_abcd_sheet_only.pdf")

ggsave(
  combined_png,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  combined_pdf,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

panel_list <- list(a = panel_a, b = panel_b, c = panel_c, d = panel_d)
panel_sizes <- list(
  a = c(7.0, 5.2),
  b = c(5.4, 4.7),
  c = c(6.2, 4.9),
  d = c(5.6, 4.9)
)

for (panel_name in names(panel_list)) {
  ggsave(
    file.path(
      output_dir,
      paste0("Figure5_panel_", panel_name, "_sheet_only.png")
    ),
    panel_list[[panel_name]],
    width = panel_sizes[[panel_name]][1],
    height = panel_sizes[[panel_name]][2],
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )
}

message("Figure and audit files written to:")
message("  ", output_dir)
message("Combined PNG: ", combined_png)
message("Combined PDF: ", combined_pdf)
##############################################################
############################################################
# 1) USER SETTINGS
############################################################

project_dir <- "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2"

# CSV/TSV/XLSX/XLS are supported. Change this to your Excel filename if needed.
input_file <- file.path(
  project_dir,
  "all_tested_genes_fibrils_crossstudy_fibj8_reversal.csv"
)
excel_sheet <- 1

output_dir <- file.path(project_dir, "Figure5_filtered_protein_coding_Q2_Q4_abcd")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Exact effect threshold requested for the Figure 5 analysis set.
effect_threshold <- 0.10

# Data-driven labels in panel a. Set to 0 to suppress gene labels.
panel_a_labels_per_quadrant <- 10L

# Number of genes shown from EACH directional-reversal quadrant in panel d.
# For example, 10 displays up to 10 Q2 genes and up to 10 Q4 genes.
top_reversal_n_per_quadrant <- 20L

# Gene-scope behavior:
#   "protein_coding" = use HGNC locus annotation to retain protein-coding genes.
#   "auto"           = same as protein_coding when HGNC annotation is available.
#   "all_rows"       = use all valid rows after the requested exclusions.
# The default preserves the requested panel-5b label: "Protein-coding genes".
gene_scope_mode <- "protein_coding"

# Requested exclusions. These are applied before the |delta| thresholds and
# therefore affect every panel (5a-5d).
exclude_unannotated_ensembl_ids <- TRUE
exclude_lincRNAs <- TRUE
exclude_pseudogenes <- TRUE
exclude_chrY_genes <- TRUE

# HGNC annotation is downloaded once and then reused from this local cache.
# If your workstation cannot access the internet, download this file manually
# and save it at hgnc_annotation_file before running the script.
hgnc_annotation_file <- file.path(project_dir, "hgnc_complete_set.txt")
hgnc_annotation_url <- paste0(
  "https://storage.googleapis.com/public-download-files/hgnc/",
  "tsv/tsv/hgnc_complete_set.txt"
)
download_hgnc_if_missing <- TRUE

figure_width_in <- 12
figure_height_in <- 9
figure_dpi <- 600

############################################################
# 2) PACKAGES
############################################################

required_packages <- c(
  "dplyr", "tibble", "ggplot2", "readr", "readxl",
  "patchwork", "ggrepel", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

############################################################
# 3) HELPERS
############################################################

read_gene_table <- function(path, sheet = 1) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  
  extension <- tolower(tools::file_ext(path))
  
  if (extension == "csv") {
    return(
      readr::read_csv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("tsv", "txt")) {
    return(
      readr::read_tsv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("xlsx", "xls")) {
    return(
      readxl::read_excel(
        path,
        sheet = sheet,
        na = c("", "NA", "NaN", "nan")
      )
    )
  }
  
  stop("Unsupported input extension: .", extension)
}

load_hgnc_annotation <- function(
    cache_file,
    source_url,
    download_if_missing = TRUE
) {
  if (!file.exists(cache_file)) {
    if (!isTRUE(download_if_missing)) {
      stop(
        "HGNC annotation file not found: ", cache_file,
        ". Set download_hgnc_if_missing <- TRUE or place the HGNC complete ",
        "set at that location."
      )
    }
    
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    message("Downloading HGNC annotation to: ", cache_file)
    
    download_result <- tryCatch(
      utils::download.file(
        url = source_url,
        destfile = cache_file,
        mode = "wb",
        quiet = FALSE
      ),
      error = function(e) e
    )
    
    if (
      inherits(download_result, "error") ||
      !file.exists(cache_file) ||
      is.na(file.info(cache_file)$size) ||
      file.info(cache_file)$size < 1000
    ) {
      if (file.exists(cache_file)) {
        unlink(cache_file)
      }
      
      error_message <- if (inherits(download_result, "error")) {
        conditionMessage(download_result)
      } else {
        "the downloaded file was absent or unexpectedly small"
      }
      
      stop(
        "Could not download the HGNC annotation: ", error_message, "\n",
        "Download it manually from:\n", source_url, "\n",
        "and save it as:\n", cache_file
      )
    }
  }
  
  hgnc_raw <- readr::read_tsv(
    cache_file,
    na = c("", "NA", "NaN", "nan"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  required_hgnc_columns <- c(
    "symbol", "name", "locus_group", "locus_type", "location"
  )
  missing_hgnc_columns <- setdiff(required_hgnc_columns, names(hgnc_raw))
  
  if (length(missing_hgnc_columns) > 0L) {
    stop(
      "The HGNC annotation file is missing required columns: ",
      paste(missing_hgnc_columns, collapse = ", ")
    )
  }
  
  hgnc_raw %>%
    transmute(
      gene_key = toupper(trimws(as.character(symbol))),
      hgnc_symbol = trimws(as.character(symbol)),
      hgnc_name = trimws(as.character(name)),
      hgnc_locus_group = tolower(trimws(as.character(locus_group))),
      hgnc_locus_type = tolower(trimws(as.character(locus_type))),
      hgnc_location = trimws(as.character(location))
    ) %>%
    filter(!is.na(gene_key), nzchar(gene_key)) %>%
    distinct(gene_key, .keep_all = TRUE)
}

blank_panel <- function(title, message) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message,
      size = 4,
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.margin = margin(8, 8, 8, 8)
    )
}

############################################################
# 4) READ, ANNOTATE, AND CALCULATE DESCRIPTIVE DIFFERENCES
############################################################

raw_gene_table <- read_gene_table(input_file, sheet = excel_sheet) %>%
  as_tibble()

required_columns <- c(
  "gene",
  "log2cpm_Control",
  "log2cpm_Fibrils",
  "log2cpm_FibJ8"
)

missing_columns <- setdiff(required_columns, names(raw_gene_table))
if (length(missing_columns) > 0L) {
  stop(
    "The input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

valid_gene_scope_modes <- c("auto", "protein_coding", "all_rows")
if (!gene_scope_mode %in% valid_gene_scope_modes) {
  stop(
    "gene_scope_mode must be one of: ",
    paste(valid_gene_scope_modes, collapse = ", ")
  )
}

hgnc_annotation <- load_hgnc_annotation(
  cache_file = hgnc_annotation_file,
  source_url = hgnc_annotation_url,
  download_if_missing = download_hgnc_if_missing
)

prepared_gene_table <- raw_gene_table %>%
  mutate(
    gene = trimws(as.character(gene)),
    gene_key = toupper(gene),
    log2cpm_Control = suppressWarnings(as.numeric(log2cpm_Control)),
    log2cpm_Fibrils = suppressWarnings(as.numeric(log2cpm_Fibrils)),
    log2cpm_FibJ8 = suppressWarnings(as.numeric(log2cpm_FibJ8))
  ) %>%
  left_join(hgnc_annotation, by = "gene_key") %>%
  mutate(
    delta_F = log2cpm_Fibrils - log2cpm_Control,
    delta_J = log2cpm_FibJ8 - log2cpm_Fibrils,
    valid_gene_identifier = !is.na(gene) & nzchar(gene),
    finite_differences = is.finite(delta_F) & is.finite(delta_J),
    hgnc_mapped = !is.na(hgnc_symbol) & nzchar(hgnc_symbol),
    
    # Exact ENSG-only rows are considered unannotated identifiers here.
    is_unannotated_ensembl_id = grepl(
      "^ENSG[0-9]+([.][0-9]+)?$",
      gene,
      ignore.case = TRUE
    ),
    
    # LINC detection is conservative: approved LINC symbols and HGNC names
    # explicitly described as long intergenic non-protein-coding RNA.
    is_lincRNA = grepl("^LINC", gene, ignore.case = TRUE) |
      (
        !is.na(hgnc_name) &
          grepl(
            "long intergenic non[- ]protein coding RNA",
            hgnc_name,
            ignore.case = TRUE
          )
      ),
    
    # HGNC locus classification avoids unreliable suffix rules such as P/P1,
    # which would incorrectly remove legitimate protein-coding symbols.
    is_pseudogene = (
      !is.na(hgnc_locus_group) & hgnc_locus_group == "pseudogene"
    ) | (
      !is.na(hgnc_locus_type) &
        grepl("pseudogene", hgnc_locus_type, ignore.case = TRUE)
    ),
    
    # A location beginning with Y denotes a Y-chromosome locus. X/Y
    # pseudoautosomal entries are not classified as Y-specific by this rule.
    is_chrY_gene = !is.na(hgnc_location) &
      grepl("^Y($|p|q)", hgnc_location, ignore.case = TRUE),
    
    is_hgnc_protein_coding = (
      !is.na(hgnc_locus_group) &
        hgnc_locus_group == "protein-coding gene"
    ) | (
      !is.na(hgnc_locus_type) &
        hgnc_locus_type == "gene with protein product"
    )
  )

if (gene_scope_mode %in% c("auto", "protein_coding")) {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = is_hgnc_protein_coding)
  gene_scope_label <- "Protein-coding genes"
  gene_scope_source <- "HGNC locus_group/locus_type"
} else {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = TRUE)
  gene_scope_label <- "All genes after requested exclusions"
  gene_scope_source <- "all valid rows in the input sheet"
}

prepared_gene_table <- prepared_gene_table %>%
  mutate(
    excluded_unannotated_ensembl =
      isTRUE(exclude_unannotated_ensembl_ids) &
      is_unannotated_ensembl_id,
    excluded_lincRNA = isTRUE(exclude_lincRNAs) & is_lincRNA,
    excluded_pseudogene = isTRUE(exclude_pseudogenes) & is_pseudogene,
    excluded_chrY = isTRUE(exclude_chrY_genes) & is_chrY_gene,
    passes_requested_exclusions = !(
      excluded_unannotated_ensembl |
        excluded_lincRNA |
        excluded_pseudogene |
        excluded_chrY
    ),
    exclusion_reason = paste0(
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    exclusion_reason = sub("; $", "", exclusion_reason),
    exclusion_reason = na_if(exclusion_reason, ""),
    analysis_exclusion_reason = paste0(
      if_else(!is_in_gene_scope, "not HGNC protein-coding; ", ""),
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    analysis_exclusion_reason = sub(
      "; $", "", analysis_exclusion_reason
    ),
    analysis_exclusion_reason = na_if(analysis_exclusion_reason, "")
  )

valid_gene_mask <- prepared_gene_table$valid_gene_identifier
finite_difference_mask <- prepared_gene_table$finite_differences

if (anyDuplicated(prepared_gene_table$gene[valid_gene_mask]) > 0L) {
  duplicate_genes <- unique(
    prepared_gene_table$gene[
      valid_gene_mask & duplicated(prepared_gene_table$gene)
    ]
  )
  
  stop(
    "The gene column contains duplicate identifiers. Resolve duplicates before ",
    "plotting. First duplicated genes: ",
    paste(utils::head(duplicate_genes, 10L), collapse = ", ")
  )
}

# Optional audit: compare newly calculated differences with any pre-existing
# effect columns, but never use those columns to select or rank genes.
delta_consistency <- tibble(
  comparison = character(),
  max_absolute_difference = numeric()
)

if ("fibrils_vs_control_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_F <- suppressWarnings(
    as.numeric(prepared_gene_table$fibrils_vs_control_log2fc)
  )
  max_difference_F <- suppressWarnings(
    max(abs(prepared_gene_table$delta_F - existing_delta_F), na.rm = TRUE)
  )
  if (!is.finite(max_difference_F)) {
    max_difference_F <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_F versus fibrils_vs_control_log2fc",
      max_absolute_difference = max_difference_F
    )
  )
}

if ("fibj8_vs_fibrils_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_J <- suppressWarnings(
    as.numeric(prepared_gene_table$fibj8_vs_fibrils_log2fc)
  )
  max_difference_J <- suppressWarnings(
    max(abs(prepared_gene_table$delta_J - existing_delta_J), na.rm = TRUE)
  )
  if (!is.finite(max_difference_J)) {
    max_difference_J <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_J versus fibj8_vs_fibrils_log2fc",
      max_absolute_difference = max_difference_J
    )
  )
}

if (nrow(delta_consistency) > 0L) {
  readr::write_csv(
    delta_consistency,
    file.path(output_dir, "Figure5_delta_consistency_audit.csv")
  )
}

############################################################
# 5) DEFINE THE FILTERED ANALYSIS SET FROM THE SHEET
############################################################

analysis_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    valid_gene_identifier,
    finite_differences,
    !is.na(is_in_gene_scope),
    is_in_gene_scope,
    passes_requested_exclusions,
    abs(delta_F) >= effect_threshold,
    abs(delta_J) >= effect_threshold
  ) %>%
  dplyr::mutate(
    quadrant = case_when(
      delta_F > 0 & delta_J > 0 ~ "Q1",
      delta_F < 0 & delta_J > 0 ~ "Q2",
      delta_F < 0 & delta_J < 0 ~ "Q3",
      delta_F > 0 & delta_J < 0 ~ "Q4",
      TRUE ~ NA_character_
    ),
    quadrant = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
    direction_class = if_else(
      quadrant %in% c("Q2", "Q4"),
      "Directional reversal (Q2/Q4)",
      "Same direction (Q1/Q3)"
    ),
    direction_class = factor(
      direction_class,
      levels = c(
        "Directional reversal (Q2/Q4)",
        "Same direction (Q1/Q3)"
      )
    ),
    fibrils_direction = if_else(
      delta_F < 0,
      "Down in Fibrils",
      "Up in Fibrils"
    ),
    fibrils_direction = factor(
      fibrils_direction,
      levels = c("Down in Fibrils", "Up in Fibrils")
    ),
    paired_effect = pmin(abs(delta_F), abs(delta_J)),
    total_effect = abs(delta_F) + abs(delta_J),
    reversal_strength = abs(delta_J)
  ) %>%
  dplyr::filter(!is.na(quadrant)) %>%
  arrange(quadrant, gene)

if (nrow(analysis_gene_table) == 0L) {
  stop(
    "No genes remain after the selected gene scope, requested exclusions, ",
    "and abs(delta_F) >= ", effect_threshold,
    " / abs(delta_J) >= ", effect_threshold, " thresholds."
  )
}

quadrant_text <- c(
  Q1 = "Q1: same direction (+,+)",
  Q2 = "Q2: reversal (-,+)",
  Q3 = "Q3: same direction (-,-)",
  Q4 = "Q4: reversal (+,-)"
)

quadrant_counts <- analysis_gene_table %>%
  dplyr::count(quadrant, name = "gene_count", .drop = FALSE) %>%
  dplyr::mutate(
    quadrant_label = factor(
      quadrant_text[as.character(quadrant)],
      levels = unname(quadrant_text)
    )
  )

base_valid_mask <-
  prepared_gene_table$valid_gene_identifier &
  prepared_gene_table$finite_differences
scope_mask <-
  base_valid_mask &
  !is.na(prepared_gene_table$is_in_gene_scope) &
  prepared_gene_table$is_in_gene_scope
post_exclusion_mask <-
  scope_mask & prepared_gene_table$passes_requested_exclusions
post_delta_F_mask <-
  post_exclusion_mask &
  abs(prepared_gene_table$delta_F) >= effect_threshold

filter_audit <- tibble(
  step = c(
    "Rows in input table",
    "Rows with nonempty identifiers and finite delta_F/delta_J",
    "Rows mapped to an approved HGNC symbol",
    "Unannotated ENSG-only IDs identified",
    "LINC RNAs identified",
    "Pseudogenes identified by HGNC",
    "Chromosome-Y genes identified by HGNC",
    paste0("Rows in selected scope: ", gene_scope_label),
    "Rows after requested gene-class exclusions",
    paste0("Rows after |delta_F| >= ", effect_threshold),
    paste0(
      "Final rows after |delta_F| >= ", effect_threshold,
      " and |delta_J| >= ", effect_threshold
    )
  ),
  n = c(
    nrow(prepared_gene_table),
    sum(base_valid_mask, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$hgnc_mapped, na.rm = TRUE),
    sum(
      base_valid_mask & prepared_gene_table$is_unannotated_ensembl_id,
      na.rm = TRUE
    ),
    sum(base_valid_mask & prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_chrY_gene, na.rm = TRUE),
    sum(scope_mask, na.rm = TRUE),
    sum(post_exclusion_mask, na.rm = TRUE),
    sum(post_delta_F_mask, na.rm = TRUE),
    nrow(analysis_gene_table)
  )
)

exclusion_summary <- tibble(
  exclusion = c(
    "Unannotated ENSG-only ID",
    "LINC RNA",
    "Pseudogene",
    "Chromosome Y"
  ),
  enabled = c(
    exclude_unannotated_ensembl_ids,
    exclude_lincRNAs,
    exclude_pseudogenes,
    exclude_chrY_genes
  ),
  n_identified = c(
    sum(prepared_gene_table$is_unannotated_ensembl_id, na.rm = TRUE),
    sum(prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$is_chrY_gene, na.rm = TRUE)
  ),
  n_flagged_for_exclusion = c(
    sum(prepared_gene_table$excluded_unannotated_ensembl, na.rm = TRUE),
    sum(prepared_gene_table$excluded_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$excluded_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$excluded_chrY, na.rm = TRUE)
  )
)

excluded_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    !is_in_gene_scope |
      !passes_requested_exclusions
  ) %>%
  dplyr::select(
    gene,
    hgnc_symbol,
    hgnc_name,
    hgnc_locus_group,
    hgnc_locus_type,
    hgnc_location,
    is_hgnc_protein_coding,
    is_unannotated_ensembl_id,
    is_lincRNA,
    is_pseudogene,
    is_chrY_gene,
    exclusion_reason,
    analysis_exclusion_reason
  ) %>%
  arrange(gene)

readr::write_csv(
  analysis_gene_table,
  file.path(output_dir, "Figure5_sheet_only_analysis_genes.csv")
)
readr::write_csv(
  quadrant_counts,
  file.path(output_dir, "Figure5_sheet_only_quadrant_counts.csv")
)
readr::write_csv(
  filter_audit,
  file.path(output_dir, "Figure5_sheet_only_filter_audit.csv")
)
readr::write_csv(
  exclusion_summary,
  file.path(output_dir, "Figure5_gene_exclusion_summary.csv")
)
readr::write_csv(
  excluded_gene_table,
  file.path(output_dir, "Figure5_excluded_genes_and_reasons.csv")
)

message("Gene scope: ", gene_scope_label)
message("Gene-scope source: ", gene_scope_source)
message("HGNC annotation cache: ", hgnc_annotation_file)
message("Final analysis genes: ", scales::comma(nrow(analysis_gene_table)))
message(
  "Quadrant counts: ",
  paste0(
    as.character(quadrant_counts$quadrant),
    "=",
    quadrant_counts$gene_count,
    collapse = ", "
  )
)

############################################################
# 6) COMMON STYLE
############################################################

delta_symbol <- "\u0394"
greater_equal_symbol <- "\u2265"
minus_symbol <- "\u2212"

reversal_blue <- "#2C7FB8"
same_direction_orange <- "#E67E39"

category_colors <- c(
  "Directional reversal (Q2/Q4)" = reversal_blue,
  "Same direction (Q1/Q3)" = same_direction_orange
)

quadrant_colors <- c(
  Q1 = "#F8766D",
  Q2 = "#7CAE00",
  Q3 = "#00BFC4",
  Q4 = "#C77CFF"
)

box_colors <- c(
  "Down in Fibrils" = "#70B96B",
  "Up in Fibrils" = "#B784C6"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5),
    plot.caption = element_text(hjust = 0, size = 7.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 8),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 7, 7, 7)
  )

############################################################
# 7) PANEL a: DESCRIPTIVE EFFECT-DIRECTION SCATTER
############################################################

if (panel_a_labels_per_quadrant > 0L) {
  panel_a_labels <- analysis_gene_table %>%
    group_by(quadrant) %>%
    arrange(
      desc(paired_effect),
      desc(total_effect),
      gene,
      .by_group = TRUE
    ) %>%
    slice_head(n = panel_a_labels_per_quadrant) %>%
    ungroup()
} else {
  panel_a_labels <- analysis_gene_table[0, , drop = FALSE]
}

readr::write_csv(
  panel_a_labels,
  file.path(output_dir, "Figure5_panel_a_data_selected_labels.csv")
)

panel_a <- ggplot(
  analysis_gene_table,
  aes(
    x = delta_F,
    y = delta_J,
    color = direction_class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_point(size = 1.15, alpha = 0.58) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.70
  ) +
  ggrepel::geom_text_repel(
    data = panel_a_labels,
    aes(label = gene),
    color = "black",
    size = 2.45,
    box.padding = 0.25,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.3,
    max.overlaps = Inf,
    seed = 1234,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Q2",
    hjust = -0.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Q1",
    hjust = 1.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    label = "Q3",
    hjust = -0.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = "Q4",
    hjust = 1.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Fibrils-associated expression and FibJ8 response",
    subtitle = "Descriptive log2CPM differences per condition (single pooled library each).",
    caption = paste0(
      "Labels: top ", panel_a_labels_per_quadrant,
      " genes per quadrant by min(|", delta_symbol, "F|, |", delta_symbol,
      "J|). Black line: descriptive least-squares trend."
    ),
    x = paste0(
      delta_symbol, "F = log2CPM(Fibrils) ", minus_symbol,
      " log2CPM(Control)"
    ),
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    ),
    color = "Effect direction"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.55, "cm"),
    plot.margin = margin(7, 18, 7, 7)
  )

############################################################
# 8) PANEL b: QUADRANT COUNTS FOR THE EXACT ANALYSIS SET
############################################################

panel_b_subtitle <- paste0(
  gene_scope_label, ", |", delta_symbol, "F| ",
  greater_equal_symbol, " ", sprintf("%.2f", effect_threshold),
  " & |", delta_symbol, "J| ", greater_equal_symbol, " ",
  sprintf("%.2f", effect_threshold)
)

panel_b <- ggplot(
  quadrant_counts,
  aes(x = quadrant_label, y = gene_count, fill = quadrant)
) +
  geom_col(width = 0.67) +
  geom_text(
    aes(label = scales::comma(gene_count)),
    vjust = -0.45,
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = quadrant_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = scales::comma
  ) +
  labs(
    title = "Effect-direction quadrant distribution",
    subtitle = panel_b_subtitle,
    caption = paste0(
      "Total genes in analysis set: ", scales::comma(nrow(analysis_gene_table)),
      ". ENSG-only IDs, LINC RNAs, pseudogenes, and chromosome-Y genes were excluded."
    ),
    x = NULL,
    y = "Gene count"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.1),
    plot.margin = margin(7, 7, 18, 7)
  )

############################################################
# 9) PANEL c: GENE-LEVEL DESCRIPTIVE DISTRIBUTION
############################################################

panel_c <- ggplot(
  analysis_gene_table,
  aes(
    x = fibrils_direction,
    y = delta_J,
    fill = fibrils_direction
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.45,
    alpha = 0.78
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1234),
    size = 0.50,
    alpha = 0.22,
    color = "black"
  ) +
  scale_fill_manual(values = box_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "FibJ8 response stratified by Fibrils-effect direction",
    subtitle = "Gene-level description; no inferential test performed.",
    caption = paste0(
      "Boxes and points summarize genes in the same threshold-qualified set: ",
      gene_scope_label, "."
    ),
    x = NULL,
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    )
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = -18, hjust = 0, vjust = 0.5, size = 7.5),
    plot.margin = margin(7, 12, 18, 7)
  )

############################################################
# 10) PANEL d: Q2 AND Q4 RANKED BY abs(delta_J)
############################################################

# Both directional-reversal quadrants are retained:
#   Q2: delta_F < 0 and delta_J > 0
#       down in Fibrils, then up with J8
#   Q4: delta_F > 0 and delta_J < 0
#       up in Fibrils, then down with J8
#
# Ranking is performed separately within Q2 and Q4 so that both directions
# are represented in panel d whenever both contain qualifying genes.
directional_reversal_gene_table <- analysis_gene_table %>%
  filter(quadrant %in% c("Q2", "Q4")) %>%
  arrange(
    quadrant,
    desc(reversal_strength),
    desc(abs(delta_F)),
    gene
  )

readr::write_csv(
  directional_reversal_gene_table,
  file.path(
    output_dir,
    "Figure5_all_Q2_Q4_genes_ranked_by_abs_delta_J.csv"
  )
)

panel_d_data <- directional_reversal_gene_table %>%
  group_by(quadrant) %>%
  arrange(
    desc(reversal_strength),
    desc(abs(delta_F)),
    gene,
    .by_group = TRUE
  ) %>%
  slice_head(n = top_reversal_n_per_quadrant) %>%
  ungroup() %>%
  arrange(reversal_strength, quadrant, gene) %>%
  mutate(
    quadrant = factor(quadrant, levels = c("Q2", "Q4")),
    gene_plot_label = paste0(gene, " [", quadrant, "]"),
    gene_plot = factor(gene_plot_label, levels = gene_plot_label)
  )

readr::write_csv(
  panel_d_data,
  file.path(
    output_dir,
    "Figure5_top_Q2_Q4_by_abs_delta_J.csv"
  )
)

q2_total <- sum(
  as.character(directional_reversal_gene_table$quadrant) == "Q2",
  na.rm = TRUE
)
q4_total <- sum(
  as.character(directional_reversal_gene_table$quadrant) == "Q4",
  na.rm = TRUE
)
q2_shown <- sum(as.character(panel_d_data$quadrant) == "Q2", na.rm = TRUE)
q4_shown <- sum(as.character(panel_d_data$quadrant) == "Q4", na.rm = TRUE)

panel_d_subtitle <- paste0(
  "Reversal strength = |", delta_symbol, "J| (|FibJ8 ",
  minus_symbol, " Fibrils|)\n",
  "Q2 = down in Fibrils, up with J8; ",
  "Q4 = up in Fibrils, down with J8"
)

if (nrow(panel_d_data) == 0L) {
  panel_d <- blank_panel(
    "Top directional-reversal genes (Q2 and Q4)",
    paste0(
      "No Q2 or Q4 genes met the selected gene-scope and ",
      "effect-threshold criteria: ",
      gene_scope_label,
      "."
    )
  )
} else {
  panel_d <- ggplot(
    panel_d_data,
    aes(
      x = reversal_strength,
      y = gene_plot,
      fill = quadrant
    )
  ) +
    geom_col(width = 0.86) +
    geom_text(
      aes(label = scales::number(reversal_strength, accuracy = 0.01)),
      hjust = -0.12,
      size = 2.35
    ) +
    scale_fill_manual(
      values = quadrant_colors[c("Q2", "Q4")],
      breaks = c("Q2", "Q4"),
      labels = c(
        "Q2: down in Fibrils, up with J8",
        "Q4: up in Fibrils, down with J8"
      ),
      name = "Reversal quadrant",
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.18)),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    labs(
      title = "Top directional-reversal genes (Q2 and Q4)",
      subtitle = panel_d_subtitle,
      caption = paste0(
        "Shown: Q2 ", q2_shown, " of ", q2_total,
        "; Q4 ", q4_shown, " of ", q4_total,
        ". Genes are ranked separately within each quadrant by |",
        delta_symbol, "J|; ties are ordered by |", delta_symbol, "F|."
      ),
      x = paste0("Reversal strength = |", delta_symbol, "J|"),
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.margin = margin(7, 12, 7, 7)
    )
}

############################################################
# 11) COMBINE AND SAVE
############################################################

figure_5_abcd <- (
  panel_a + panel_b + plot_layout(widths = c(1.45, 1.00))
) / (
  panel_c + panel_d + plot_layout(widths = c(1.20, 1.00))
) +
  plot_layout(heights = c(1.00, 1.00)) +
  plot_annotation(
    title = "Figure 5",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 1, size = 13),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

combined_png <- file.path(output_dir, "Figure5_abcd_sheet_only.png")
combined_pdf <- file.path(output_dir, "Figure5_abcd_sheet_only.pdf")

ggsave(
  combined_png,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  combined_pdf,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

panel_list <- list(a = panel_a, b = panel_b, c = panel_c, d = panel_d)
panel_sizes <- list(
  a = c(7.0, 5.2),
  b = c(5.4, 4.7),
  c = c(5.2, 4.9),
  d = c(10, 6)
)

for (panel_name in names(panel_list)) {
  ggsave(
    file.path(
      output_dir,
      paste0("Figure5_panel_", panel_name, "_sheet_only.pdf")
    ),
    panel_list[[panel_name]],
    width = panel_sizes[[panel_name]][1],
    height = panel_sizes[[panel_name]][2],
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )
}

message("Figure and audit files written to:")
message("  ", output_dir)
message("Combined PNG: ", combined_png)
message("Combined PDF: ", combined_pdf)
#######################################################

############################################################
# 12) Q2/Q4 GENE LISTS FOR SHINYGO AND ENRICHMENT SETTINGS
############################################################

############################################################
# INSTALL BIOCONDUCTOR ENRICHMENT PACKAGES
# Run this block once, not every time the analysis is run.
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

# Restore the appropriate CRAN and Bioconductor repositories.
options(repos = BiocManager::repositories())

BiocManager::install(
  c(
    "ReactomePA",
    "clusterProfiler",
    "AnnotationDbi",
    "org.Hs.eg.db"
  ),
  ask = FALSE,
  update = FALSE
)

# This section pools both directional-reversal quadrants for enrichment:
#   Q2: delta_F < 0 and delta_J > 0
#   Q4: delta_F > 0 and delta_J < 0
#
# The pooled test asks which pathways/processes are over-represented among
# all directional-reversal genes, regardless of whether the reversal is Q2
# or Q4. Separate Q2 and Q4 symbol lists are also exported for optional
# direction-specific analysis in ShinyGO.

enrichment_output_dir <- file.path(output_dir, "Q2_Q4_enrichment")
dir.create(enrichment_output_dir, recursive = TRUE, showWarnings = FALSE)

enrichment_fdr_cutoff <- 0.1
min_gene_set_size <- 10L
max_gene_set_size <- 500L

top_terms_kegg <- 15L
top_terms_reactome <- 15L
top_terms_go_bp <- 15L

# KEGG uses the current online KEGG annotation when FALSE. This step requires
# internet access while the script is running.
kegg_use_internal_data <- FALSE

q2_q4_gene_table <- analysis_gene_table %>%
  dplyr::filter(quadrant %in% c("Q2", "Q4")) %>%
  dplyr::mutate(quadrant = as.character(quadrant)) %>%
  dplyr::arrange(quadrant, desc(reversal_strength), desc(abs(delta_F)), gene)

if (nrow(q2_q4_gene_table) == 0L) {
  stop(
    "No Q2 or Q4 genes are available for the enrichment analysis after ",
    "the selected gene-scope, exclusion, and effect-size filters."
  )
}

q2_symbols <- q2_q4_gene_table %>%
  dplyr::filter(quadrant == "Q2") %>%
  pull(gene) %>%
  unique() %>%
  sort()

q4_symbols <- q2_q4_gene_table %>%
  dplyr::filter(quadrant == "Q4") %>%
  pull(gene) %>%
  unique() %>%
  sort()

q2_q4_symbols <- sort(unique(c(q2_symbols, q4_symbols)))

# Header-free, one-symbol-per-line files can be pasted directly into ShinyGO.
writeLines(
  q2_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q2_symbols.txt")
)
writeLines(
  q4_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q4_symbols.txt")
)
writeLines(
  q2_q4_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q2_Q4_combined_symbols.txt")
)

readr::write_csv(
  tibble(gene = q2_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q2_symbols.csv")
)
readr::write_csv(
  tibble(gene = q4_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q4_symbols.csv")
)
readr::write_csv(
  tibble(gene = q2_q4_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q2_Q4_combined_symbols.csv")
)

readr::write_csv(
  q2_q4_gene_table %>%
    dplyr::select(
      gene,
      quadrant,
      delta_F,
      delta_J,
      reversal_strength,
      log2cpm_Control,
      log2cpm_Fibrils,
      log2cpm_FibJ8
    ),
  file.path(enrichment_output_dir, "Q2_Q4_gene_table_with_effects.csv")
)

readr::write_csv(
  tibble(
    gene_set = c("Q2", "Q4", "Q2 + Q4 pooled"),
    n_symbols = c(
      length(q2_symbols),
      length(q4_symbols),
      length(q2_q4_symbols)
    )
  ),
  file.path(enrichment_output_dir, "Q2_Q4_gene_list_counts.csv")
)

# The ShinyGO-ready symbol files above are written before this package check.
# Therefore, even when an enrichment package is missing, the exact Q2/Q4
# lists remain available in the output directory.

required_enrichment_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db"
)

missing_enrichment_packages <- required_enrichment_packages[
  !vapply(
    required_enrichment_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_enrichment_packages) > 0L) {
  stop(
    "Install these Bioconductor packages before running enrichment: ",
    paste(missing_enrichment_packages, collapse = ", "),
    "\n\nRun:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"AnnotationDbi\", \"clusterProfiler\", ",
    "\"ReactomePA\", \"org.Hs.eg.db\"))"
  )
}


############################################################
# 13) MAP SYMBOLS TO ENTREZ IDS AND DEFINE THE TESTED BACKGROUND
############################################################

# The background is not the whole human genome. It is every input-sheet gene
# that had finite delta_F/delta_J values, belonged to the selected HGNC gene
# scope, and passed the requested ENSG/LINC/pseudogene/chrY exclusions. The
# effect-size threshold is intentionally NOT applied to the background.
background_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    valid_gene_identifier,
    finite_differences,
    !is.na(is_in_gene_scope),
    is_in_gene_scope,
    passes_requested_exclusions
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(gene)

background_symbols <- sort(unique(background_gene_table$gene))

orgdb <- org.Hs.eg.db::org.Hs.eg.db

map_symbols_to_entrez <- function(symbols, org_db) {
  symbols <- sort(unique(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  
  if (length(symbols) == 0L) {
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }
  
  suppressMessages(
    AnnotationDbi::select(
      org_db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    )
  ) %>%
    as_tibble() %>%
    dplyr::transmute(
      SYMBOL = trimws(as.character(SYMBOL)),
      ENTREZID = trimws(as.character(ENTREZID))
    ) %>%
    dplyr::filter(
      !is.na(SYMBOL),
      nzchar(SYMBOL),
      !is.na(ENTREZID),
      nzchar(ENTREZID)
    ) %>%
    distinct(SYMBOL, ENTREZID)
}

all_symbol_entrez_map <- map_symbols_to_entrez(
  union(background_symbols, q2_q4_symbols),
  orgdb
)

selected_symbol_entrez_map <- all_symbol_entrez_map %>%
  dplyr::filter(SYMBOL %in% q2_q4_symbols)

background_symbol_entrez_map <- all_symbol_entrez_map %>%
  dplyr::filter(SYMBOL %in% background_symbols)

selected_entrez <- sort(unique(selected_symbol_entrez_map$ENTREZID))
background_entrez <- sort(unique(background_symbol_entrez_map$ENTREZID))

if (length(selected_entrez) == 0L) {
  stop(
    "None of the pooled Q2/Q4 symbols mapped to Entrez IDs in org.Hs.eg.db."
  )
}

if (length(background_entrez) == 0L) {
  stop(
    "None of the eligible tested background symbols mapped to Entrez IDs ",
    "in org.Hs.eg.db."
  )
}

# Ensure the selected set is contained in the supplied universe.
background_entrez <- sort(unique(c(background_entrez, selected_entrez)))

symbol_mapping_audit <- tibble(SYMBOL = union(background_symbols, q2_q4_symbols)) %>%
  left_join(
    q2_q4_gene_table %>%
      distinct(gene, quadrant) %>%
      dplyr::rename(SYMBOL = gene),
    by = "SYMBOL"
  ) %>%
  left_join(all_symbol_entrez_map, by = "SYMBOL") %>%
  dplyr::mutate(
    is_Q2_Q4_selected = SYMBOL %in% q2_q4_symbols,
    is_background_gene = SYMBOL %in% background_symbols,
    mapped_to_entrez = !is.na(ENTREZID) & nzchar(ENTREZID)
  ) %>%
  dplyr::arrange(desc(is_Q2_Q4_selected), quadrant, SYMBOL, ENTREZID)

readr::write_csv(
  symbol_mapping_audit,
  file.path(enrichment_output_dir, "Q2_Q4_symbol_to_Entrez_mapping_audit.csv")
)

unmapped_selected_symbols <- setdiff(
  q2_q4_symbols,
  unique(selected_symbol_entrez_map$SYMBOL)
)
writeLines(
  unmapped_selected_symbols,
  file.path(enrichment_output_dir, "Q2_Q4_unmapped_symbols.txt")
)

mapping_summary <- tibble(
  item = c(
    "Q2 symbols",
    "Q4 symbols",
    "Pooled Q2/Q4 symbols",
    "Pooled Q2/Q4 Entrez IDs",
    "Eligible tested background symbols",
    "Eligible tested background Entrez IDs",
    "Unmapped pooled Q2/Q4 symbols"
  ),
  n = c(
    length(q2_symbols),
    length(q4_symbols),
    length(q2_q4_symbols),
    length(selected_entrez),
    length(background_symbols),
    length(background_entrez),
    length(unmapped_selected_symbols)
  )
)

readr::write_csv(
  mapping_summary,
  file.path(enrichment_output_dir, "Q2_Q4_enrichment_input_summary.csv")
)

############################################################
# 14) KEGG, REACTOME, AND GO:BP OVER-REPRESENTATION ANALYSIS
############################################################

# pvalueCutoff and qvalueCutoff are set to 1 so the full result tables are
# retained. Plotting and significant-result exports then apply the explicit
# BH-adjusted p-value threshold enrichment_fdr_cutoff.

kegg_result <- tryCatch(
  clusterProfiler::enrichKEGG(
    gene = selected_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_entrez,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = 1,
    use_internal_data = kegg_use_internal_data
  ),
  error = function(e) {
    warning("KEGG enrichment failed: ", conditionMessage(e))
    NULL
  }
)

reactome_result <- tryCatch(
  ReactomePA::enrichPathway(
    gene = selected_entrez,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_entrez,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = 1,
    organism = "human",
    readable = TRUE
  ),
  error = function(e) {
    warning("Reactome enrichment failed: ", conditionMessage(e))
    NULL
  }
)

go_bp_result <- tryCatch(
  clusterProfiler::enrichGO(
    gene = selected_entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_entrez,
    qvalueCutoff = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    readable = TRUE,
    pool = FALSE
  ),
  error = function(e) {
    warning("GO:BP enrichment failed: ", conditionMessage(e))
    NULL
  }
)

ratio_to_numeric <- function(x) {
  x <- as.character(x)
  
  vapply(
    strsplit(x, "/", fixed = TRUE),
    FUN.VALUE = numeric(1),
    FUN = function(parts) {
      if (length(parts) != 2L) {
        return(NA_real_)
      }
      
      numerator <- suppressWarnings(as.numeric(parts[1]))
      denominator <- suppressWarnings(as.numeric(parts[2]))
      
      if (
        !is.finite(numerator) ||
        !is.finite(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }
      
      numerator / denominator
    }
  )
}

entrez_to_symbol <- setNames(
  all_symbol_entrez_map$SYMBOL,
  all_symbol_entrez_map$ENTREZID
)

translate_gene_id_field <- function(gene_id_string) {
  if (is.na(gene_id_string) || !nzchar(gene_id_string)) {
    return(NA_character_)
  }
  
  ids <- strsplit(as.character(gene_id_string), "/", fixed = TRUE)[[1]]
  ids <- ids[nzchar(ids)]
  
  if (length(ids) == 0L) {
    return(NA_character_)
  }
  
  # readable=TRUE already gives symbols for GO and Reactome. KEGG usually
  # returns numeric Entrez IDs, which are converted here.
  if (all(grepl("^[0-9]+$", ids))) {
    translated <- unname(entrez_to_symbol[ids])
    translated[is.na(translated) | !nzchar(translated)] <-
      ids[is.na(translated) | !nzchar(translated)]
    ids <- translated
  }
  
  paste(unique(ids), collapse = "/")
}

empty_enrichment_table <- function() {
  tibble(
    database = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    GeneRatio_numeric = numeric(),
    BgRatio_numeric = numeric(),
    FoldEnrichment = numeric(),
    gene_symbols = character()
  )
}

as_enrichment_table <- function(result_object, database_name) {
  if (is.null(result_object)) {
    return(empty_enrichment_table())
  }
  
  result_table <- as.data.frame(result_object) %>%
    as_tibble()
  
  if (nrow(result_table) == 0L) {
    return(empty_enrichment_table())
  }
  
  if (!"Count" %in% names(result_table)) {
    result_table$Count <- NA_integer_
  }
  if (!"GeneRatio" %in% names(result_table)) {
    result_table$GeneRatio <- NA_character_
  }
  if (!"BgRatio" %in% names(result_table)) {
    result_table$BgRatio <- NA_character_
  }
  if (!"geneID" %in% names(result_table)) {
    result_table$geneID <- NA_character_
  }
  
  result_table %>%
    mutate(
      database = database_name,
      Count = suppressWarnings(as.integer(Count)),
      GeneRatio_numeric = ratio_to_numeric(GeneRatio),
      BgRatio_numeric = ratio_to_numeric(BgRatio),
      FoldEnrichment = if_else(
        is.finite(GeneRatio_numeric) &
          is.finite(BgRatio_numeric) &
          BgRatio_numeric > 0,
        GeneRatio_numeric / BgRatio_numeric,
        NA_real_
      ),
      gene_symbols = vapply(
        geneID,
        translate_gene_id_field,
        character(1)
      )
    ) %>%
    relocate(database, .before = 1)
}

kegg_table <- as_enrichment_table(kegg_result, "KEGG")
reactome_table <- as_enrichment_table(reactome_result, "Reactome")
go_bp_table <- as_enrichment_table(go_bp_result, "GO:BP")

write_enrichment_outputs <- function(result_table, file_prefix, fdr_cutoff) {
  all_file <- file.path(
    enrichment_output_dir,
    paste0(file_prefix, "_all_results.csv")
  )
  significant_file <- file.path(
    enrichment_output_dir,
    paste0(file_prefix, "_BH_FDR05_results.csv")
  )
  
  readr::write_csv(result_table, all_file)
  
  significant_table <- result_table %>%
    filter(
      is.finite(p.adjust),
      p.adjust <= fdr_cutoff
    ) %>%
    arrange(p.adjust, desc(Count), Description)
  
  readr::write_csv(significant_table, significant_file)
  invisible(significant_table)
}

kegg_significant <- write_enrichment_outputs(
  kegg_table,
  "KEGG_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)
reactome_significant <- write_enrichment_outputs(
  reactome_table,
  "Reactome_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)
go_bp_significant <- write_enrichment_outputs(
  go_bp_table,
  "GO_BP_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)

############################################################
# 15) FIGURE 5e-f STYLE BAR PLOTS FOR THE POOLED Q2 + Q4 SET
############################################################

wrap_plot_labels <- function(x, width = 34L) {
  vapply(
    as.character(x),
    FUN.VALUE = character(1),
    FUN = function(label) {
      paste(strwrap(label, width = width), collapse = "\n")
    }
  )
}

empty_enrichment_plot <- function(title, subtitle, database_name) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0(
        "No ", database_name,
        " terms met BH FDR ", less_than_or_equal_symbol,
        " ", sprintf("%.2f", enrichment_fdr_cutoff), "."
      ),
      size = 3.7,
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = element_text(
        hjust = 0,
        size = 8,
        face = "italic",
        color = "grey35"
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Define the symbol here because it was not needed in panels a-d.
less_than_or_equal_symbol <- "\u2264"

make_enrichment_barplot <- function(
    result_table,
    title,
    database_name,
    x_variable = c("Count", "FoldEnrichment"),
    x_axis_title,
    top_n = 12L,
    bar_fill = "#2B8CBE"
) {
  x_variable <- match.arg(x_variable)
  
  plot_subtitle <- paste0(
    "ORA on ", scales::comma(length(selected_entrez)),
    " mapped Q2 + Q4 genes; background = ",
    scales::comma(length(background_entrez)),
    " eligible tested genes; BH FDR ",
    less_than_or_equal_symbol, " ",
    sprintf("%.2f", enrichment_fdr_cutoff)
  )
  
  plot_table <- result_table %>%
    filter(
      is.finite(p.adjust),
      p.adjust <= enrichment_fdr_cutoff,
      is.finite(.data[[x_variable]]),
      .data[[x_variable]] > 0
    ) %>%
    arrange(p.adjust, desc(Count), Description) %>%
    slice_head(n = top_n)
  
  if (nrow(plot_table) == 0L) {
    return(
      empty_enrichment_plot(
        title = title,
        subtitle = plot_subtitle,
        database_name = database_name
      )
    )
  }
  
  plot_table <- plot_table %>%
    arrange(.data[[x_variable]], desc(Count), Description) %>%
    mutate(
      plot_key = make.unique(paste0(Description, " [", ID, "]")),
      term_label = wrap_plot_labels(Description, width = 34L),
      plot_key = factor(plot_key, levels = plot_key)
    )
  
  y_label_lookup <- setNames(
    plot_table$term_label,
    as.character(plot_table$plot_key)
  )
  
  ggplot(
    plot_table,
    aes(x = .data[[x_variable]], y = plot_key)
  ) +
    geom_col(width = 0.72, fill = bar_fill) +
    scale_y_discrete(labels = y_label_lookup) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = scales::label_number(accuracy = if (x_variable == "Count") 1 else 0.1)
    ) +
    labs(
      title = title,
      subtitle = plot_subtitle,
      caption = paste0(
        "Top ", min(top_n, nrow(plot_table)),
        " adjusted-significant terms ordered by BH-adjusted p-value; ",
        "ties use overlapping-gene count."
      ),
      x = x_axis_title,
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = element_text(
        hjust = 0,
        size = 7.8,
        face = "italic",
        color = "grey35",
        margin = margin(b = 7)
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 6.8,
        color = "grey35",
        margin = margin(t = 7)
      ),
      axis.title.x = element_text(size = 8.5, margin = margin(t = 8)),
      axis.text.x = element_text(size = 7.6, color = "grey30"),
      axis.text.y = element_text(size = 7.2, color = "grey30"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.35,
        linetype = "dashed"
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

panel_e_kegg <- make_enrichment_barplot(
  result_table = kegg_table,
  title = "KEGG pathways (Q2 + Q4 pooled set)",
  database_name = "KEGG",
  x_variable = "Count",
  x_axis_title = "Overlapping genes",
  top_n = top_terms_kegg,
  bar_fill = "#2B8CBE"
) +
  labs(tag = "e") +
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

panel_f_go_bp <- make_enrichment_barplot(
  result_table = go_bp_table,
  title = "GO biological processes (Q2 + Q4 pooled set)",
  database_name = "GO:BP",
  x_variable = "FoldEnrichment",
  x_axis_title = "Fold enrichment",
  top_n = top_terms_go_bp,
  bar_fill = "#6A3D9A"
) +
  labs(tag = "f") +
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

panel_reactome <- make_enrichment_barplot(
  result_table = reactome_table,
  title = "Reactome pathways (Q2 + Q4 pooled set)",
  database_name = "Reactome",
  x_variable = "Count",
  x_axis_title = "Overlapping genes",
  top_n = top_terms_reactome,
  bar_fill = "#3A9D5D"
)

figure_5_ef_q2_q4 <- panel_e_kegg + panel_f_go_bp +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1))

figure_5_three_database_q2_q4 <-
  make_enrichment_barplot(
    result_table = kegg_table,
    title = "KEGG pathways",
    database_name = "KEGG",
    x_variable = "Count",
    x_axis_title = "Overlapping genes",
    top_n = top_terms_kegg,
    bar_fill = "#2B8CBE"
  ) +
  panel_reactome +
  make_enrichment_barplot(
    result_table = go_bp_table,
    title = "GO biological processes",
    database_name = "GO:BP",
    x_variable = "FoldEnrichment",
    x_axis_title = "Fold enrichment",
    top_n = top_terms_go_bp,
    bar_fill = "#6A3D9A"
  ) +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(
    title = "Q2 + Q4 pooled directional-reversal enrichment",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )
  )

# Two-panel Figure 5e-f analogue: KEGG plus GO:BP.
ggsave(
  file.path(
    enrichment_output_dir,
    "Figure5_ef_Q2_Q4_pooled_KEGG_GO_BP.png"
  ),
  figure_5_ef_q2_q4,
  width = 13.5,
  height = 5.6,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    enrichment_output_dir,
    "Figure5_ef_Q2_Q4_pooled_KEGG_GO_BP.pdf"
  ),
  figure_5_ef_q2_q4,
  width = 13.5,
  height = 5.6,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# Reactome is saved separately because the requested two-panel figure mirrors
# the attached KEGG + GO:BP structure.
ggsave(
  file.path(enrichment_output_dir, "Reactome_Q2_Q4_pooled.png"),
  panel_reactome,
  width = 6.8,
  height = 5.6,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(enrichment_output_dir, "Reactome_Q2_Q4_pooled.pdf"),
  panel_reactome,
  width = 6.8,
  height = 5.6,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# A three-panel version is also exported so KEGG, Reactome, and GO:BP can be
# inspected together without changing the requested two-panel Figure 5e-f.
ggsave(
  file.path(
    enrichment_output_dir,
    "Q2_Q4_pooled_KEGG_Reactome_GO_BP_three_panel.png"
  ),
  figure_5_three_database_q2_q4,
  width = 19.5,
  height = 5.8,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    enrichment_output_dir,
    "Q2_Q4_pooled_KEGG_Reactome_GO_BP_three_panel.pdf"
  ),
  figure_5_three_database_q2_q4,
  width = 19.5,
  height = 5.8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

message("Q2/Q4 ShinyGO lists and enrichment outputs written to:")
message("  ", enrichment_output_dir)
message(
  "Pooled Q2/Q4 symbols: ", length(q2_q4_symbols),
  "; mapped Entrez IDs: ", length(selected_entrez)
)
message(
  "Eligible background symbols: ", length(background_symbols),
  "; mapped Entrez IDs: ", length(background_entrez)
)
message(
  "BH FDR <= ", enrichment_fdr_cutoff,
  " terms: KEGG=", nrow(kegg_significant),
  ", Reactome=", nrow(reactome_significant),
  ", GO:BP=", nrow(go_bp_significant)
)
#############################################################################
############################################################
# FIGURE 5a-d: DATA-DERIVED EXPRESSION ANALYSIS WITH HGNC GENE FILTERS
#
# All expression differences, thresholds, quadrant assignments, rankings,
# and plotted values are calculated only from the input CSV/TSV/XLSX/XLS
# table. No values are copied from the reference figure. Official HGNC
# annotation is used only to classify protein-coding genes and to exclude
# LINC genes, pseudogenes, and chromosome-Y genes.
#
# Definitions calculated directly from the sheet:
#   delta_F = log2CPM_Fibrils - log2CPM_Control
#   delta_J = log2CPM_FibJ8   - log2CPM_Fibrils
#
# Base effect-qualified analysis set:
#   genes in the selected gene scope with abs(delta_F) >= 0.20 and
#   abs(delta_J) >= 0.20.
#
# Expression floor used for panels 5a and 5d:
#   mean_CPM = mean(cpm_Control, cpm_Fibrils, cpm_FibJ8)
#   mean_CPM >= minimum_mean_cpm_panels_a_d
#
# Panel d ranks both directional-reversal quadrants, Q2 and Q4,
# separately by abs(delta_J), after the expression floor.
#
# GENE-SCOPE RULE:
#   gene_scope_mode = "protein_coding" retains HGNC loci classified as
#   protein-coding genes. The requested ENSG-only, LINC, pseudogene, and
#   chromosome-Y exclusions are then applied before the effect thresholds.
############################################################

options(stringsAsFactors = FALSE)
set.seed(1234)

############################################################
# 1) USER SETTINGS
############################################################

project_dir <- "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2"

# CSV/TSV/XLSX/XLS are supported. Change this to your Excel filename if needed.
input_file <- file.path(
  project_dir,
  "all_tested_genes_fibrils_crossstudy_fibj8_reversal.csv"
)
excel_sheet <- 1

output_dir <- file.path(
  project_dir,
  "Figure5_filtered_protein_coding_Q2_Q4_expression_floor"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Exact effect threshold requested for the Figure 5 analysis set.
effect_threshold <- 0.10

# Expression floor for panels 5a and 5d. The arithmetic mean is calculated
# from the three raw CPM columns in the input sheet, not from log2CPM.
# Change 1.00 to the exact CPM threshold stated in your methods, if different.
minimum_mean_cpm_panels_a_d <- 5.00

# Keep the Q2/Q4 ShinyGO lists and enrichment foreground consistent with
# panel 5d. When TRUE, the same mean-CPM floor is also applied to the Q2/Q4
# enrichment set and to its eligible tested-gene background.
apply_expression_floor_to_q2_q4_enrichment <- TRUE

# Data-driven labels in panel a. Set to 0 to suppress gene labels.
panel_a_labels_per_quadrant <- 15L

# Number of genes shown from EACH directional-reversal quadrant in panel d.
# For example, 10 displays up to 10 Q2 genes and up to 10 Q4 genes.
top_reversal_n_per_quadrant <- 15L

# Gene-scope behavior:
#   "protein_coding" = use HGNC locus annotation to retain protein-coding genes.
#   "auto"           = same as protein_coding when HGNC annotation is available.
#   "all_rows"       = use all valid rows after the requested exclusions.
# The default preserves the requested panel-5b label: "Protein-coding genes".
gene_scope_mode <- "protein_coding"

# Requested exclusions. These are applied before the |delta| thresholds and
# therefore affect every panel (5a-5d).
exclude_unannotated_ensembl_ids <- TRUE
exclude_lincRNAs <- TRUE
exclude_pseudogenes <- TRUE
exclude_chrY_genes <- TRUE

# HGNC annotation is downloaded once and then reused from this local cache.
# If your workstation cannot access the internet, download this file manually
# and save it at hgnc_annotation_file before running the script.
hgnc_annotation_file <- file.path(project_dir, "hgnc_complete_set.txt")
hgnc_annotation_url <- paste0(
  "https://storage.googleapis.com/public-download-files/hgnc/",
  "tsv/tsv/hgnc_complete_set.txt"
)
download_hgnc_if_missing <- TRUE

figure_width_in <- 12
figure_height_in <- 9
figure_dpi <- 600

############################################################
# 2) PACKAGES
############################################################

required_packages <- c(
  "dplyr", "tibble", "ggplot2", "readr", "readxl",
  "patchwork", "ggrepel", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

############################################################
# 3) HELPERS
############################################################

read_gene_table <- function(path, sheet = 1) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  
  extension <- tolower(tools::file_ext(path))
  
  if (extension == "csv") {
    return(
      readr::read_csv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("tsv", "txt")) {
    return(
      readr::read_tsv(
        path,
        na = c("", "NA", "NaN", "nan"),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (extension %in% c("xlsx", "xls")) {
    return(
      readxl::read_excel(
        path,
        sheet = sheet,
        na = c("", "NA", "NaN", "nan")
      )
    )
  }
  
  stop("Unsupported input extension: .", extension)
}

load_hgnc_annotation <- function(
    cache_file,
    source_url,
    download_if_missing = TRUE
) {
  if (!file.exists(cache_file)) {
    if (!isTRUE(download_if_missing)) {
      stop(
        "HGNC annotation file not found: ", cache_file,
        ". Set download_hgnc_if_missing <- TRUE or place the HGNC complete ",
        "set at that location."
      )
    }
    
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    message("Downloading HGNC annotation to: ", cache_file)
    
    download_result <- tryCatch(
      utils::download.file(
        url = source_url,
        destfile = cache_file,
        mode = "wb",
        quiet = FALSE
      ),
      error = function(e) e
    )
    
    if (
      inherits(download_result, "error") ||
      !file.exists(cache_file) ||
      is.na(file.info(cache_file)$size) ||
      file.info(cache_file)$size < 1000
    ) {
      if (file.exists(cache_file)) {
        unlink(cache_file)
      }
      
      error_message <- if (inherits(download_result, "error")) {
        conditionMessage(download_result)
      } else {
        "the downloaded file was absent or unexpectedly small"
      }
      
      stop(
        "Could not download the HGNC annotation: ", error_message, "\n",
        "Download it manually from:\n", source_url, "\n",
        "and save it as:\n", cache_file
      )
    }
  }
  
  hgnc_raw <- readr::read_tsv(
    cache_file,
    na = c("", "NA", "NaN", "nan"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  
  required_hgnc_columns <- c(
    "symbol", "name", "locus_group", "locus_type", "location"
  )
  missing_hgnc_columns <- setdiff(required_hgnc_columns, names(hgnc_raw))
  
  if (length(missing_hgnc_columns) > 0L) {
    stop(
      "The HGNC annotation file is missing required columns: ",
      paste(missing_hgnc_columns, collapse = ", ")
    )
  }
  
  hgnc_raw %>%
    transmute(
      gene_key = toupper(trimws(as.character(symbol))),
      hgnc_symbol = trimws(as.character(symbol)),
      hgnc_name = trimws(as.character(name)),
      hgnc_locus_group = tolower(trimws(as.character(locus_group))),
      hgnc_locus_type = tolower(trimws(as.character(locus_type))),
      hgnc_location = trimws(as.character(location))
    ) %>%
    filter(!is.na(gene_key), nzchar(gene_key)) %>%
    distinct(gene_key, .keep_all = TRUE)
}

blank_panel <- function(title, message) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message,
      size = 4,
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.margin = margin(8, 8, 8, 8)
    )
}

############################################################
# 4) READ, ANNOTATE, AND CALCULATE DESCRIPTIVE DIFFERENCES
############################################################

raw_gene_table <- read_gene_table(input_file, sheet = excel_sheet) %>%
  as_tibble()

required_columns <- c(
  "gene",
  "cpm_Control",
  "cpm_Fibrils",
  "cpm_FibJ8",
  "log2cpm_Control",
  "log2cpm_Fibrils",
  "log2cpm_FibJ8"
)

missing_columns <- setdiff(required_columns, names(raw_gene_table))
if (length(missing_columns) > 0L) {
  stop(
    "The input table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

valid_gene_scope_modes <- c("auto", "protein_coding", "all_rows")
if (!gene_scope_mode %in% valid_gene_scope_modes) {
  stop(
    "gene_scope_mode must be one of: ",
    paste(valid_gene_scope_modes, collapse = ", ")
  )
}

if (
  length(minimum_mean_cpm_panels_a_d) != 1L ||
  !is.numeric(minimum_mean_cpm_panels_a_d) ||
  !is.finite(minimum_mean_cpm_panels_a_d) ||
  minimum_mean_cpm_panels_a_d < 0
) {
  stop("minimum_mean_cpm_panels_a_d must be one finite number >= 0.")
}

hgnc_annotation <- load_hgnc_annotation(
  cache_file = hgnc_annotation_file,
  source_url = hgnc_annotation_url,
  download_if_missing = download_hgnc_if_missing
)

prepared_gene_table <- raw_gene_table %>%
  mutate(
    gene = trimws(as.character(gene)),
    gene_key = toupper(gene),
    cpm_Control = suppressWarnings(as.numeric(cpm_Control)),
    cpm_Fibrils = suppressWarnings(as.numeric(cpm_Fibrils)),
    cpm_FibJ8 = suppressWarnings(as.numeric(cpm_FibJ8)),
    log2cpm_Control = suppressWarnings(as.numeric(log2cpm_Control)),
    log2cpm_Fibrils = suppressWarnings(as.numeric(log2cpm_Fibrils)),
    log2cpm_FibJ8 = suppressWarnings(as.numeric(log2cpm_FibJ8))
  ) %>%
  left_join(hgnc_annotation, by = "gene_key") %>%
  mutate(
    delta_F = log2cpm_Fibrils - log2cpm_Control,
    delta_J = log2cpm_FibJ8 - log2cpm_Fibrils,
    mean_cpm_across_conditions = (
      cpm_Control + cpm_Fibrils + cpm_FibJ8
    ) / 3,
    valid_gene_identifier = !is.na(gene) & nzchar(gene),
    finite_differences = is.finite(delta_F) & is.finite(delta_J),
    finite_cpm_values = is.finite(cpm_Control) &
      is.finite(cpm_Fibrils) &
      is.finite(cpm_FibJ8) &
      is.finite(mean_cpm_across_conditions),
    passes_panel_expression_floor = finite_cpm_values &
      mean_cpm_across_conditions >= minimum_mean_cpm_panels_a_d,
    hgnc_mapped = !is.na(hgnc_symbol) & nzchar(hgnc_symbol),
    
    # Exact ENSG-only rows are considered unannotated identifiers here.
    is_unannotated_ensembl_id = grepl(
      "^ENSG[0-9]+([.][0-9]+)?$",
      gene,
      ignore.case = TRUE
    ),
    
    # LINC detection is conservative: approved LINC symbols and HGNC names
    # explicitly described as long intergenic non-protein-coding RNA.
    is_lincRNA = grepl("^LINC", gene, ignore.case = TRUE) |
      (
        !is.na(hgnc_name) &
          grepl(
            "long intergenic non[- ]protein coding RNA",
            hgnc_name,
            ignore.case = TRUE
          )
      ),
    
    # HGNC locus classification avoids unreliable suffix rules such as P/P1,
    # which would incorrectly remove legitimate protein-coding symbols.
    is_pseudogene = (
      !is.na(hgnc_locus_group) & hgnc_locus_group == "pseudogene"
    ) | (
      !is.na(hgnc_locus_type) &
        grepl("pseudogene", hgnc_locus_type, ignore.case = TRUE)
    ),
    
    # A location beginning with Y denotes a Y-chromosome locus. X/Y
    # pseudoautosomal entries are not classified as Y-specific by this rule.
    is_chrY_gene = !is.na(hgnc_location) &
      grepl("^Y($|p|q)", hgnc_location, ignore.case = TRUE),
    
    is_hgnc_protein_coding = (
      !is.na(hgnc_locus_group) &
        hgnc_locus_group == "protein-coding gene"
    ) | (
      !is.na(hgnc_locus_type) &
        hgnc_locus_type == "gene with protein product"
    )
  )

if (gene_scope_mode %in% c("auto", "protein_coding")) {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = is_hgnc_protein_coding)
  gene_scope_label <- "Protein-coding genes"
  gene_scope_source <- "HGNC locus_group/locus_type"
} else {
  prepared_gene_table <- prepared_gene_table %>%
    mutate(is_in_gene_scope = TRUE)
  gene_scope_label <- "All genes after requested exclusions"
  gene_scope_source <- "all valid rows in the input sheet"
}

prepared_gene_table <- prepared_gene_table %>%
  mutate(
    excluded_unannotated_ensembl =
      isTRUE(exclude_unannotated_ensembl_ids) &
      is_unannotated_ensembl_id,
    excluded_lincRNA = isTRUE(exclude_lincRNAs) & is_lincRNA,
    excluded_pseudogene = isTRUE(exclude_pseudogenes) & is_pseudogene,
    excluded_chrY = isTRUE(exclude_chrY_genes) & is_chrY_gene,
    passes_requested_exclusions = !(
      excluded_unannotated_ensembl |
        excluded_lincRNA |
        excluded_pseudogene |
        excluded_chrY
    ),
    exclusion_reason = paste0(
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    exclusion_reason = sub("; $", "", exclusion_reason),
    exclusion_reason = na_if(exclusion_reason, ""),
    analysis_exclusion_reason = paste0(
      if_else(!is_in_gene_scope, "not HGNC protein-coding; ", ""),
      if_else(
        excluded_unannotated_ensembl,
        "unannotated ENSG ID; ",
        ""
      ),
      if_else(excluded_lincRNA, "LINC RNA; ", ""),
      if_else(excluded_pseudogene, "pseudogene; ", ""),
      if_else(excluded_chrY, "chromosome Y; ", "")
    ),
    analysis_exclusion_reason = sub(
      "; $", "", analysis_exclusion_reason
    ),
    analysis_exclusion_reason = na_if(analysis_exclusion_reason, "")
  )

valid_gene_mask <- prepared_gene_table$valid_gene_identifier
finite_difference_mask <- prepared_gene_table$finite_differences

if (anyDuplicated(prepared_gene_table$gene[valid_gene_mask]) > 0L) {
  duplicate_genes <- unique(
    prepared_gene_table$gene[
      valid_gene_mask & duplicated(prepared_gene_table$gene)
    ]
  )
  
  stop(
    "The gene column contains duplicate identifiers. Resolve duplicates before ",
    "plotting. First duplicated genes: ",
    paste(utils::head(duplicate_genes, 10L), collapse = ", ")
  )
}

# Optional audit: compare newly calculated differences with any pre-existing
# effect columns, but never use those columns to select or rank genes.
delta_consistency <- tibble(
  comparison = character(),
  max_absolute_difference = numeric()
)

if ("fibrils_vs_control_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_F <- suppressWarnings(
    as.numeric(prepared_gene_table$fibrils_vs_control_log2fc)
  )
  max_difference_F <- suppressWarnings(
    max(abs(prepared_gene_table$delta_F - existing_delta_F), na.rm = TRUE)
  )
  if (!is.finite(max_difference_F)) {
    max_difference_F <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_F versus fibrils_vs_control_log2fc",
      max_absolute_difference = max_difference_F
    )
  )
}

if ("fibj8_vs_fibrils_log2fc" %in% names(prepared_gene_table)) {
  existing_delta_J <- suppressWarnings(
    as.numeric(prepared_gene_table$fibj8_vs_fibrils_log2fc)
  )
  max_difference_J <- suppressWarnings(
    max(abs(prepared_gene_table$delta_J - existing_delta_J), na.rm = TRUE)
  )
  if (!is.finite(max_difference_J)) {
    max_difference_J <- NA_real_
  }
  
  delta_consistency <- bind_rows(
    delta_consistency,
    tibble(
      comparison = "delta_J versus fibj8_vs_fibrils_log2fc",
      max_absolute_difference = max_difference_J
    )
  )
}

if (nrow(delta_consistency) > 0L) {
  readr::write_csv(
    delta_consistency,
    file.path(output_dir, "Figure5_delta_consistency_audit.csv")
  )
}

############################################################
# 5) DEFINE THE FILTERED ANALYSIS SET FROM THE SHEET
############################################################

analysis_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    valid_gene_identifier,
    finite_differences,
    !is.na(is_in_gene_scope),
    is_in_gene_scope,
    passes_requested_exclusions,
    abs(delta_F) >= effect_threshold,
    abs(delta_J) >= effect_threshold
  ) %>%
  dplyr::mutate(
    quadrant = case_when(
      delta_F > 0 & delta_J > 0 ~ "Q1",
      delta_F < 0 & delta_J > 0 ~ "Q2",
      delta_F < 0 & delta_J < 0 ~ "Q3",
      delta_F > 0 & delta_J < 0 ~ "Q4",
      TRUE ~ NA_character_
    ),
    quadrant = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
    direction_class = if_else(
      quadrant %in% c("Q2", "Q4"),
      "Directional reversal (Q2/Q4)",
      "Same direction (Q1/Q3)"
    ),
    direction_class = factor(
      direction_class,
      levels = c(
        "Directional reversal (Q2/Q4)",
        "Same direction (Q1/Q3)"
      )
    ),
    fibrils_direction = if_else(
      delta_F < 0,
      "Down in Fibrils",
      "Up in Fibrils"
    ),
    fibrils_direction = factor(
      fibrils_direction,
      levels = c("Down in Fibrils", "Up in Fibrils")
    ),
    paired_effect = pmin(abs(delta_F), abs(delta_J)),
    total_effect = abs(delta_F) + abs(delta_J),
    reversal_strength = abs(delta_J)
  ) %>%
  dplyr::filter(!is.na(quadrant)) %>%
  dplyr::arrange(quadrant, gene)

if (nrow(analysis_gene_table) == 0L) {
  stop(
    "No genes remain after the selected gene scope, requested exclusions, ",
    "and abs(delta_F) >= ", effect_threshold,
    " / abs(delta_J) >= ", effect_threshold, " thresholds."
  )
}

# Panels 5a, 5b, and 5d use the additional expression floor. Panel 5c
# intentionally retains the effect-qualified set without the additional CPM
# floor, so its gene-level distribution remains as previously specified.
high_expression_analysis_gene_table <- analysis_gene_table %>%
  dplyr::filter(passes_panel_expression_floor) %>%
  dplyr::arrange(quadrant, gene)

if (nrow(high_expression_analysis_gene_table) == 0L) {
  stop(
    "No genes remain for panels 5a, 5b, and 5d after requiring mean CPM >= ",
    minimum_mean_cpm_panels_a_d,
    " across Control, Fibrils, and FibJ8."
  )
}

below_expression_floor_gene_table <- analysis_gene_table %>%
  dplyr::filter(!passes_panel_expression_floor) %>%
  dplyr::arrange(quadrant, gene)

expression_floor_audit <- tibble(
  item = c(
    "Effect-qualified genes before panel expression floor",
    paste0(
      "Genes retained for panels 5a, 5b, and 5d after mean CPM >= ",
      sprintf("%.2f", minimum_mean_cpm_panels_a_d)
    ),
    "Q2 genes before panel expression floor",
    "Q2 genes retained for panel 5d",
    "Q4 genes before panel expression floor",
    "Q4 genes retained for panel 5d"
  ),
  n = c(
    nrow(analysis_gene_table),
    nrow(high_expression_analysis_gene_table),
    sum(as.character(analysis_gene_table$quadrant) == "Q2", na.rm = TRUE),
    sum(
      as.character(high_expression_analysis_gene_table$quadrant) == "Q2",
      na.rm = TRUE
    ),
    sum(as.character(analysis_gene_table$quadrant) == "Q4", na.rm = TRUE),
    sum(
      as.character(high_expression_analysis_gene_table$quadrant) == "Q4",
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  high_expression_analysis_gene_table,
  file.path(output_dir, "Figure5_panels_a_b_d_meanCPM_filtered_genes.csv")
)
readr::write_csv(
  below_expression_floor_gene_table,
  file.path(output_dir, "Figure5_genes_removed_by_panels_a_b_d_meanCPM_floor.csv")
)
readr::write_csv(
  expression_floor_audit,
  file.path(output_dir, "Figure5_panels_a_b_d_expression_floor_audit.csv")
)

quadrant_text <- c(
  Q1 = "Q1: same direction (+,+)",
  Q2 = "Q2: reversal (-,+)",
  Q3 = "Q3: same direction (-,-)",
  Q4 = "Q4: reversal (+,-)"
)

# Panel 5b uses the same expression-qualified set as panels 5a, 5b, and 5d.
panel_b_gene_table <- high_expression_analysis_gene_table

quadrant_counts <- panel_b_gene_table %>%
  dplyr::count(quadrant, name = "gene_count", .drop = FALSE) %>%
  dplyr::mutate(
    quadrant_label = factor(
      quadrant_text[as.character(quadrant)],
      levels = unname(quadrant_text)
    )
  )

quadrant_count_cpm_audit <- bind_rows(
  analysis_gene_table %>%
    dplyr::count(quadrant, name = "gene_count", .drop = FALSE) %>%
    dplyr::mutate(filter_set = "Before mean-CPM floor"),
  panel_b_gene_table %>%
    dplyr::count(quadrant, name = "gene_count", .drop = FALSE) %>%
    dplyr::mutate(
      filter_set = paste0(
        "Mean CPM >= ",
        sprintf("%.2f", minimum_mean_cpm_panels_a_d)
      )
    )
) %>%
  dplyr::mutate(
    quadrant = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))
  ) %>%
  dplyr::arrange(filter_set, quadrant)

base_valid_mask <-
  prepared_gene_table$valid_gene_identifier &
  prepared_gene_table$finite_differences
scope_mask <-
  base_valid_mask &
  !is.na(prepared_gene_table$is_in_gene_scope) &
  prepared_gene_table$is_in_gene_scope
post_exclusion_mask <-
  scope_mask & prepared_gene_table$passes_requested_exclusions
post_delta_F_mask <-
  post_exclusion_mask &
  abs(prepared_gene_table$delta_F) >= effect_threshold

filter_audit <- tibble(
  step = c(
    "Rows in input table",
    "Rows with nonempty identifiers and finite delta_F/delta_J",
    "Rows mapped to an approved HGNC symbol",
    "Unannotated ENSG-only IDs identified",
    "LINC RNAs identified",
    "Pseudogenes identified by HGNC",
    "Chromosome-Y genes identified by HGNC",
    paste0("Rows in selected scope: ", gene_scope_label),
    "Rows after requested gene-class exclusions",
    paste0("Rows after |delta_F| >= ", effect_threshold),
    paste0(
      "Final rows after |delta_F| >= ", effect_threshold,
      " and |delta_J| >= ", effect_threshold
    ),
    paste0(
      "Rows retained for panels 5a, 5b, and 5d after mean CPM >= ",
      sprintf("%.2f", minimum_mean_cpm_panels_a_d)
    )
  ),
  n = c(
    nrow(prepared_gene_table),
    sum(base_valid_mask, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$hgnc_mapped, na.rm = TRUE),
    sum(
      base_valid_mask & prepared_gene_table$is_unannotated_ensembl_id,
      na.rm = TRUE
    ),
    sum(base_valid_mask & prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(base_valid_mask & prepared_gene_table$is_chrY_gene, na.rm = TRUE),
    sum(scope_mask, na.rm = TRUE),
    sum(post_exclusion_mask, na.rm = TRUE),
    sum(post_delta_F_mask, na.rm = TRUE),
    nrow(analysis_gene_table),
    nrow(panel_b_gene_table)
  )
)

exclusion_summary <- tibble(
  exclusion = c(
    "Unannotated ENSG-only ID",
    "LINC RNA",
    "Pseudogene",
    "Chromosome Y"
  ),
  enabled = c(
    exclude_unannotated_ensembl_ids,
    exclude_lincRNAs,
    exclude_pseudogenes,
    exclude_chrY_genes
  ),
  n_identified = c(
    sum(prepared_gene_table$is_unannotated_ensembl_id, na.rm = TRUE),
    sum(prepared_gene_table$is_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$is_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$is_chrY_gene, na.rm = TRUE)
  ),
  n_flagged_for_exclusion = c(
    sum(prepared_gene_table$excluded_unannotated_ensembl, na.rm = TRUE),
    sum(prepared_gene_table$excluded_lincRNA, na.rm = TRUE),
    sum(prepared_gene_table$excluded_pseudogene, na.rm = TRUE),
    sum(prepared_gene_table$excluded_chrY, na.rm = TRUE)
  )
)

excluded_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    !is_in_gene_scope |
      !passes_requested_exclusions
  ) %>%
  dplyr::select(
    gene,
    hgnc_symbol,
    hgnc_name,
    hgnc_locus_group,
    hgnc_locus_type,
    hgnc_location,
    is_hgnc_protein_coding,
    is_unannotated_ensembl_id,
    is_lincRNA,
    is_pseudogene,
    is_chrY_gene,
    exclusion_reason,
    analysis_exclusion_reason
  ) %>%
  dplyr::arrange(gene)

readr::write_csv(
  analysis_gene_table,
  file.path(output_dir, "Figure5_sheet_only_analysis_genes.csv")
)
readr::write_csv(
  quadrant_counts,
  file.path(output_dir, "Figure5_CPM_filtered_quadrant_counts.csv")
)
readr::write_csv(
  quadrant_count_cpm_audit,
  file.path(output_dir, "Figure5_quadrant_counts_before_after_CPM_floor.csv")
)
readr::write_csv(
  filter_audit,
  file.path(output_dir, "Figure5_sheet_only_filter_audit.csv")
)
readr::write_csv(
  exclusion_summary,
  file.path(output_dir, "Figure5_gene_exclusion_summary.csv")
)
readr::write_csv(
  excluded_gene_table,
  file.path(output_dir, "Figure5_excluded_genes_and_reasons.csv")
)

message("Gene scope: ", gene_scope_label)
message("Gene-scope source: ", gene_scope_source)
message("HGNC annotation cache: ", hgnc_annotation_file)
message("Final effect-qualified analysis genes: ", scales::comma(nrow(analysis_gene_table)))
message(
  "Genes retained for panels 5a/5b/5d after mean CPM >= ",
  sprintf("%.2f", minimum_mean_cpm_panels_a_d),
  ": ",
  scales::comma(nrow(high_expression_analysis_gene_table))
)
message(
  "Quadrant counts: ",
  paste0(
    as.character(quadrant_counts$quadrant),
    "=",
    quadrant_counts$gene_count,
    collapse = ", "
  )
)

############################################################
# 6) COMMON STYLE
############################################################

delta_symbol <- "\u0394"
greater_equal_symbol <- "\u2265"
minus_symbol <- "\u2212"

reversal_blue <- "#2C7FB8"
same_direction_orange <- "#E67E39"

category_colors <- c(
  "Directional reversal (Q2/Q4)" = reversal_blue,
  "Same direction (Q1/Q3)" = same_direction_orange
)

quadrant_colors <- c(
  Q1 = "#F8766D",
  Q2 = "#7CAE00",
  Q3 = "#00BFC4",
  Q4 = "#C77CFF"
)

box_colors <- c(
  "Down in Fibrils" = "#70B96B",
  "Up in Fibrils" = "#B784C6"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5),
    plot.caption = element_text(hjust = 0, size = 7.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 8),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 7, 7, 7)
  )

############################################################
# 7) PANEL a: DESCRIPTIVE EFFECT-DIRECTION SCATTER
############################################################

# All points and labels in panel 5a pass the mean-CPM expression floor.
# This avoids ranking large log2 differences driven by very low abundance.
panel_a_plot_data <- high_expression_analysis_gene_table

if (panel_a_labels_per_quadrant > 0L) {
  panel_a_labels <- panel_a_plot_data %>%
    group_by(quadrant) %>%
    dplyr::arrange(
      desc(paired_effect),
      desc(total_effect),
      desc(mean_cpm_across_conditions),
      gene,
      .by_group = TRUE
    ) %>%
    slice_head(n = panel_a_labels_per_quadrant) %>%
    ungroup()
} else {
  panel_a_labels <- panel_a_plot_data[0, , drop = FALSE]
}

readr::write_csv(
  panel_a_plot_data,
  file.path(output_dir, "Figure5_panel_a_meanCPM_filtered_plot_data.csv")
)
readr::write_csv(
  panel_a_labels,
  file.path(output_dir, "Figure5_panel_a_data_selected_labels.csv")
)

panel_a_subtitle <- paste0(
  "Descriptive log2CPM differences per condition (single pooled library each).\n",
  "Expression floor: mean CPM across Control/Fibrils/FibJ8 ",
  greater_equal_symbol, " ",
  sprintf("%.2f", minimum_mean_cpm_panels_a_d), "."
)

panel_a_has_trend <- (
  nrow(panel_a_plot_data) >= 2L &&
    dplyr::n_distinct(panel_a_plot_data$delta_F) >= 2L
)

panel_a <- ggplot(
  panel_a_plot_data,
  aes(
    x = delta_F,
    y = delta_J,
    color = direction_class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.45
  ) +
  geom_point(size = 1.15, alpha = 0.58)

# Add the descriptive trend before labels so the line cannot obscure text.
if (panel_a_has_trend) {
  panel_a <- panel_a +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.70
    )
}

panel_a <- panel_a +
  ggrepel::geom_text_repel(
    data = panel_a_labels,
    aes(label = gene),
    color = "black",
    size = 2.45,
    box.padding = 0.25,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.linewidth = 0.3,
    max.overlaps = Inf,
    seed = 1234,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Q2",
    hjust = -0.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "Q1",
    hjust = 1.35,
    vjust = 1.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    label = "Q3",
    hjust = -0.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = "Q4",
    hjust = 1.35,
    vjust = -0.35,
    fontface = "bold.italic",
    size = 4
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Fibrils-associated expression and FibJ8 response",
    subtitle = panel_a_subtitle,
    caption = paste0(
      "Retained ", scales::comma(nrow(panel_a_plot_data)), " of ",
      scales::comma(nrow(analysis_gene_table)),
      " effect-qualified genes after the expression floor. Labels: top ",
      panel_a_labels_per_quadrant,
      " genes per quadrant by min(|", delta_symbol, "F|, |", delta_symbol,
      "J|); mean CPM breaks later ties. ",
      if (panel_a_has_trend) {
        "Black line: descriptive least-squares trend."
      } else {
        "A trend line could not be estimated from the retained values."
      }
    ),
    x = paste0(
      delta_symbol, "F = log2CPM(Fibrils) ", minus_symbol,
      " log2CPM(Control)"
    ),
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    ),
    color = "Effect direction"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.55, "cm"),
    plot.margin = margin(7, 18, 7, 7)
  )

############################################################
# 8) PANEL b: QUADRANT COUNTS FOR THE EXACT ANALYSIS SET
############################################################

panel_b_subtitle <- paste0(
  gene_scope_label, ", |", delta_symbol, "F| ",
  greater_equal_symbol, " ", sprintf("%.2f", effect_threshold),
  " & |", delta_symbol, "J| ", greater_equal_symbol, " ",
  sprintf("%.2f", effect_threshold),
  "; mean CPM ", greater_equal_symbol, " ",
  sprintf("%.2f", minimum_mean_cpm_panels_a_d)
)

panel_b <- ggplot(
  quadrant_counts,
  aes(x = quadrant_label, y = gene_count, fill = quadrant)
) +
  geom_col(width = 0.67) +
  geom_text(
    aes(label = scales::comma(gene_count)),
    vjust = -0.45,
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = quadrant_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = scales::comma
  ) +
  labs(
    title = "Effect-direction quadrant distribution",
    subtitle = panel_b_subtitle,
    caption = paste0(
      "Total genes after effect and expression filters: ",
      scales::comma(nrow(panel_b_gene_table)),
      ". Mean CPM across Control, Fibrils, and FibJ8 ",
      greater_equal_symbol, " ",
      sprintf("%.2f", minimum_mean_cpm_panels_a_d),
      ". ENSG-only IDs, LINC RNAs, pseudogenes, and chromosome-Y genes were excluded."
    ),
    x = NULL,
    y = "Gene count"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.1),
    plot.margin = margin(7, 7, 18, 7)
  )
############################################################
# 9) PANEL c: GENE-LEVEL DESCRIPTIVE DISTRIBUTION
############################################################

panel_c <- ggplot(
  analysis_gene_table,
  aes(
    x = fibrils_direction,
    y = delta_J,
    fill = fibrils_direction
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.45,
    alpha = 0.78
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 1234),
    size = 0.50,
    alpha = 0.22,
    color = "black"
  ) +
  scale_fill_manual(values = box_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "FibJ8 response stratified by Fibrils-effect direction",
    subtitle = "Gene-level description; no inferential test performed.",
    caption = paste0(
      "Boxes and points summarize genes in the same threshold-qualified set: ",
      gene_scope_label, "."
    ),
    x = NULL,
    y = paste0(
      delta_symbol, "J = log2CPM(FibJ8) ", minus_symbol,
      " log2CPM(Fibrils)"
    )
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = -18, hjust = 0, vjust = 0.5, size = 7.5),
    plot.margin = margin(7, 12, 18, 7)
  )

############################################################
# 10) PANEL d: Q2 AND Q4 RANKED BY abs(delta_J)
############################################################

# Both directional-reversal quadrants are retained:
#   Q2: delta_F < 0 and delta_J > 0
#       down in Fibrils, then up with J8
#   Q4: delta_F > 0 and delta_J < 0
#       up in Fibrils, then down with J8
#
# Ranking is performed separately within Q2 and Q4 so that both directions
# are represented in panel d whenever both contain qualifying genes.
directional_reversal_gene_table <- high_expression_analysis_gene_table %>%
  dplyr::filter(quadrant %in% c("Q2", "Q4")) %>%
  dplyr::arrange(
    quadrant,
    desc(reversal_strength),
    desc(abs(delta_F)),
    desc(mean_cpm_across_conditions),
    gene
  )

readr::write_csv(
  directional_reversal_gene_table,
  file.path(
    output_dir,
    "Figure5_all_Q2_Q4_genes_ranked_by_abs_delta_J.csv"
  )
)

panel_d_data <- directional_reversal_gene_table %>%
  group_by(quadrant) %>%
  dplyr::arrange(
    desc(reversal_strength),
    desc(abs(delta_F)),
    desc(mean_cpm_across_conditions),
    gene,
    .by_group = TRUE
  ) %>%
  slice_head(n = top_reversal_n_per_quadrant) %>%
  ungroup() %>%
  dplyr::arrange(reversal_strength, quadrant, gene) %>%
  dplyr::mutate(
    quadrant = factor(quadrant, levels = c("Q2", "Q4")),
    gene_plot_label = paste0(gene, " [", quadrant, "]"),
    gene_plot = factor(gene_plot_label, levels = gene_plot_label)
  )

readr::write_csv(
  panel_d_data,
  file.path(
    output_dir,
    "Figure5_top_Q2_Q4_by_abs_delta_J.csv"
  )
)

q2_total <- sum(
  as.character(directional_reversal_gene_table$quadrant) == "Q2",
  na.rm = TRUE
)
q4_total <- sum(
  as.character(directional_reversal_gene_table$quadrant) == "Q4",
  na.rm = TRUE
)
q2_shown <- sum(as.character(panel_d_data$quadrant) == "Q2", na.rm = TRUE)
q4_shown <- sum(as.character(panel_d_data$quadrant) == "Q4", na.rm = TRUE)

panel_d_subtitle <- paste0(
  "Reversal strength = |", delta_symbol, "J| (|FibJ8 ",
  minus_symbol, " Fibrils|)\n",
  "Q2 = down in Fibrils, up with J8; ",
  "Q4 = up in Fibrils, down with J8\n",
  "Mean CPM across conditions ", greater_equal_symbol, " ",
  sprintf("%.2f", minimum_mean_cpm_panels_a_d)
)

if (nrow(panel_d_data) == 0L) {
  panel_d <- blank_panel(
    "Top directional-reversal genes (Q2 and Q4)",
    paste0(
      "No Q2 or Q4 genes met the selected gene-scope and ",
      "effect-threshold criteria: ",
      gene_scope_label,
      "."
    )
  )
} else {
  panel_d <- ggplot(
    panel_d_data,
    aes(
      x = reversal_strength,
      y = gene_plot,
      fill = quadrant
    )
  ) +
    geom_col(width = 0.86) +
    geom_text(
      aes(label = scales::number(reversal_strength, accuracy = 0.01)),
      hjust = -0.12,
      size = 2.35
    ) +
    scale_fill_manual(
      values = quadrant_colors[c("Q2", "Q4")],
      breaks = c("Q2", "Q4"),
      labels = c(
        "Q2: down in Fibrils, up with J8",
        "Q4: up in Fibrils, down with J8"
      ),
      name = "Reversal quadrant",
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.18)),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    labs(
      title = "Top directional-reversal genes (Q2 and Q4)",
      subtitle = panel_d_subtitle,
      caption = paste0(
        "Shown: Q2 ", q2_shown, " of ", q2_total,
        "; Q4 ", q4_shown, " of ", q4_total,
        ". All displayed genes pass mean CPM ", greater_equal_symbol, " ",
        sprintf("%.2f", minimum_mean_cpm_panels_a_d),
        " across Control/Fibrils/FibJ8. Genes are ranked separately ",
        "within each quadrant by |", delta_symbol,
        "J|; ties are ordered by |", delta_symbol, "F| and mean CPM."
      ),
      x = paste0("Reversal strength = |", delta_symbol, "J|"),
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.margin = margin(7, 12, 7, 7)
    )
}

############################################################
# 11) COMBINE AND SAVE
############################################################

figure_5_abcd <- (
  panel_a + panel_b + plot_layout(widths = c(1.45, 1.00))
) / (
  panel_c + panel_d + plot_layout(widths = c(1.20, 1.00))
) +
  plot_layout(heights = c(1.00, 1.00)) +
  plot_annotation(
    title = "Figure 5",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 1, size = 13),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

combined_png <- file.path(output_dir, "Figure5_abcd_sheet_only.png")
combined_pdf <- file.path(output_dir, "Figure5_abcd_sheet_only.pdf")

ggsave(
  combined_png,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  combined_pdf,
  figure_5_abcd,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

panel_list <- list(a = panel_a, b = panel_b, c = panel_c, d = panel_d)
panel_sizes <- list(
  a = c(7.0, 5.2),
  b = c(5.4, 4.7),
  c = c(6.2, 4.9),
  d = c(6.0, 6.2)
)

for (panel_name in names(panel_list)) {
  ggsave(
    file.path(
      output_dir,
      paste0("Figure5_panel_", panel_name, "_sheet_only.pdf")
    ),
    panel_list[[panel_name]],
    width = panel_sizes[[panel_name]][1],
    height = panel_sizes[[panel_name]][2],
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )
}

message("Figure and audit files written to:")
message("  ", output_dir)
message("Combined PNG: ", combined_png)
message("Combined PDF: ", combined_pdf)

############################################################
# 12) Q2/Q4 GENE LISTS FOR SHINYGO AND ENRICHMENT SETTINGS
############################################################

# This section pools both directional-reversal quadrants for enrichment:
#   Q2: delta_F < 0 and delta_J > 0
#   Q4: delta_F > 0 and delta_J < 0
#
# The pooled test asks which pathways/processes are over-represented among
# all directional-reversal genes, regardless of whether the reversal is Q2
# or Q4. Separate Q2 and Q4 symbol lists are also exported for optional
# direction-specific analysis in ShinyGO.

enrichment_output_dir <- file.path(output_dir, "Q2_Q4_enrichment")
dir.create(enrichment_output_dir, recursive = TRUE, showWarnings = FALSE)

enrichment_fdr_cutoff <- 0.1
min_gene_set_size <- 10L
max_gene_set_size <- 500L

top_terms_kegg <- 12L
top_terms_reactome <- 12L
top_terms_go_bp <- 12L

# KEGG uses the current online KEGG annotation when FALSE. This step requires
# internet access while the script is running.
kegg_use_internal_data <- FALSE

enrichment_foreground_source <- if (
  isTRUE(apply_expression_floor_to_q2_q4_enrichment)
) {
  high_expression_analysis_gene_table
} else {
  analysis_gene_table
}

q2_q4_gene_table <- enrichment_foreground_source %>%
  dplyr::filter(quadrant %in% c("Q2", "Q4")) %>%
  dplyr::mutate(quadrant = as.character(quadrant)) %>%
  dplyr::arrange(
    quadrant,
    desc(reversal_strength),
    desc(abs(delta_F)),
    desc(mean_cpm_across_conditions),
    gene
  )

if (nrow(q2_q4_gene_table) == 0L) {
  stop(
    "No Q2 or Q4 genes are available for the enrichment analysis after ",
    "the selected gene-scope, exclusion, effect-size, and optional ",
    "expression-floor filters."
  )
}

q2_symbols <- q2_q4_gene_table %>%
  dplyr::filter(quadrant == "Q2") %>%
  pull(gene) %>%
  unique() %>%
  sort()

q4_symbols <- q2_q4_gene_table %>%
  dplyr::filter(quadrant == "Q4") %>%
  pull(gene) %>%
  unique() %>%
  sort()

q2_q4_symbols <- sort(unique(c(q2_symbols, q4_symbols)))

# Header-free, one-symbol-per-line files can be pasted directly into ShinyGO.
writeLines(
  q2_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q2_symbols.txt")
)
writeLines(
  q4_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q4_symbols.txt")
)
writeLines(
  q2_q4_symbols,
  file.path(enrichment_output_dir, "ShinyGO_Q2_Q4_combined_symbols.txt")
)

readr::write_csv(
  tibble(gene = q2_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q2_symbols.csv")
)
readr::write_csv(
  tibble(gene = q4_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q4_symbols.csv")
)
readr::write_csv(
  tibble(gene = q2_q4_symbols),
  file.path(enrichment_output_dir, "ShinyGO_Q2_Q4_combined_symbols.csv")
)

readr::write_csv(
  q2_q4_gene_table %>%
    dplyr::select(
      gene,
      quadrant,
      delta_F,
      delta_J,
      reversal_strength,
      cpm_Control,
      cpm_Fibrils,
      cpm_FibJ8,
      mean_cpm_across_conditions,
      log2cpm_Control,
      log2cpm_Fibrils,
      log2cpm_FibJ8
    ),
  file.path(enrichment_output_dir, "Q2_Q4_gene_table_with_effects.csv")
)

readr::write_csv(
  tibble(
    gene_set = c("Q2", "Q4", "Q2 + Q4 pooled"),
    n_symbols = c(
      length(q2_symbols),
      length(q4_symbols),
      length(q2_q4_symbols)
    )
  ),
  file.path(enrichment_output_dir, "Q2_Q4_gene_list_counts.csv")
)

readr::write_csv(
  tibble(
    setting = c(
      "effect_threshold_abs_delta_F",
      "effect_threshold_abs_delta_J",
      "apply_mean_CPM_floor_to_Q2_Q4_enrichment",
      "minimum_mean_CPM"
    ),
    value = c(
      as.character(effect_threshold),
      as.character(effect_threshold),
      as.character(apply_expression_floor_to_q2_q4_enrichment),
      as.character(minimum_mean_cpm_panels_a_d)
    )
  ),
  file.path(enrichment_output_dir, "Q2_Q4_enrichment_filter_settings.csv")
)

# The ShinyGO-ready symbol files above are written before this package check.
# Therefore, even when an enrichment package is missing, the exact Q2/Q4
# lists remain available in the output directory.

required_enrichment_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "org.Hs.eg.db"
)

missing_enrichment_packages <- required_enrichment_packages[
  !vapply(
    required_enrichment_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_enrichment_packages) > 0L) {
  stop(
    "Install these Bioconductor packages before running enrichment: ",
    paste(missing_enrichment_packages, collapse = ", "),
    "\n\nRun:\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"AnnotationDbi\", ",
    "\"clusterProfiler\", \"org.Hs.eg.db\"))"
  )
}

# ReactomePA is optional. This avoids stopping KEGG and GO:BP when the package
# is unavailable for the active R/Bioconductor combination.
has_reactomepa <- requireNamespace("ReactomePA", quietly = TRUE)
if (!has_reactomepa) {
  warning(
    "ReactomePA is not installed. Reactome analysis will be skipped; ",
    "KEGG and GO:BP will continue."
  )
}

############################################################
# 13) MAP SYMBOLS TO ENTREZ IDS AND DEFINE THE TESTED BACKGROUND
############################################################

# The background is not the whole human genome. It is every input-sheet gene
# that had finite delta_F/delta_J values, belonged to the selected HGNC gene
# scope, and passed the requested ENSG/LINC/pseudogene/chrY exclusions. The
# effect-size threshold is intentionally NOT applied to the background.
background_gene_table <- prepared_gene_table %>%
  dplyr::filter(
    valid_gene_identifier,
    finite_differences,
    !is.na(is_in_gene_scope),
    is_in_gene_scope,
    passes_requested_exclusions
  )

# When the foreground uses the mean-CPM floor, the enrichment background must
# use the same eligibility rule. Otherwise low-expression genes that could
# never enter the selected set would remain in the universe.
if (isTRUE(apply_expression_floor_to_q2_q4_enrichment)) {
  background_gene_table <- background_gene_table %>%
    dplyr::filter(passes_panel_expression_floor)
}

background_gene_table <- background_gene_table %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(gene)

background_symbols <- sort(unique(background_gene_table$gene))

orgdb <- org.Hs.eg.db::org.Hs.eg.db

map_symbols_to_entrez <- function(symbols, org_db) {
  symbols <- sort(unique(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  
  if (length(symbols) == 0L) {
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }
  
  suppressMessages(
    AnnotationDbi::select(
      org_db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID")
    )
  ) %>%
    as_tibble() %>%
    dplyr::transmute(
      SYMBOL = trimws(as.character(SYMBOL)),
      ENTREZID = trimws(as.character(ENTREZID))
    ) %>%
    dplyr::filter(
      !is.na(SYMBOL),
      nzchar(SYMBOL),
      !is.na(ENTREZID),
      nzchar(ENTREZID)
    ) %>%
    distinct(SYMBOL, ENTREZID)
}

all_symbol_entrez_map <- map_symbols_to_entrez(
  union(background_symbols, q2_q4_symbols),
  orgdb
)

selected_symbol_entrez_map <- all_symbol_entrez_map %>%
  dplyr::filter(SYMBOL %in% q2_q4_symbols)

background_symbol_entrez_map <- all_symbol_entrez_map %>%
  dplyr::filter(SYMBOL %in% background_symbols)

selected_entrez <- sort(unique(selected_symbol_entrez_map$ENTREZID))
background_entrez <- sort(unique(background_symbol_entrez_map$ENTREZID))

if (length(selected_entrez) == 0L) {
  stop(
    "None of the pooled Q2/Q4 symbols mapped to Entrez IDs in org.Hs.eg.db."
  )
}

if (length(background_entrez) == 0L) {
  stop(
    "None of the eligible tested background symbols mapped to Entrez IDs ",
    "in org.Hs.eg.db."
  )
}

# Ensure the selected set is contained in the supplied universe.
background_entrez <- sort(unique(c(background_entrez, selected_entrez)))

symbol_mapping_audit <- tibble(SYMBOL = union(background_symbols, q2_q4_symbols)) %>%
  left_join(
    q2_q4_gene_table %>%
      distinct(gene, quadrant) %>%
      dplyr::rename(SYMBOL = gene),
    by = "SYMBOL"
  ) %>%
  left_join(all_symbol_entrez_map, by = "SYMBOL") %>%
  dplyr::mutate(
    is_Q2_Q4_selected = SYMBOL %in% q2_q4_symbols,
    is_background_gene = SYMBOL %in% background_symbols,
    mapped_to_entrez = !is.na(ENTREZID) & nzchar(ENTREZID)
  ) %>%
  arrange(desc(is_Q2_Q4_selected), quadrant, SYMBOL, ENTREZID)

readr::write_csv(
  symbol_mapping_audit,
  file.path(enrichment_output_dir, "Q2_Q4_symbol_to_Entrez_mapping_audit.csv")
)

unmapped_selected_symbols <- setdiff(
  q2_q4_symbols,
  unique(selected_symbol_entrez_map$SYMBOL)
)
writeLines(
  unmapped_selected_symbols,
  file.path(enrichment_output_dir, "Q2_Q4_unmapped_symbols.txt")
)

mapping_summary <- tibble(
  item = c(
    "Q2 symbols",
    "Q4 symbols",
    "Pooled Q2/Q4 symbols",
    "Pooled Q2/Q4 Entrez IDs",
    "Eligible tested background symbols",
    "Eligible tested background Entrez IDs",
    "Unmapped pooled Q2/Q4 symbols"
  ),
  n = c(
    length(q2_symbols),
    length(q4_symbols),
    length(q2_q4_symbols),
    length(selected_entrez),
    length(background_symbols),
    length(background_entrez),
    length(unmapped_selected_symbols)
  )
)

readr::write_csv(
  mapping_summary,
  file.path(enrichment_output_dir, "Q2_Q4_enrichment_input_summary.csv")
)

############################################################
# 14) KEGG, REACTOME, AND GO:BP OVER-REPRESENTATION ANALYSIS
############################################################

# pvalueCutoff and qvalueCutoff are set to 1 so the full result tables are
# retained. Plotting and significant-result exports then apply the explicit
# BH-adjusted p-value threshold enrichment_fdr_cutoff.

kegg_result <- tryCatch(
  clusterProfiler::enrichKEGG(
    gene = selected_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_entrez,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = 1,
    use_internal_data = kegg_use_internal_data
  ),
  error = function(e) {
    warning("KEGG enrichment failed: ", conditionMessage(e))
    NULL
  }
)

reactome_result <- if (has_reactomepa) {
  tryCatch(
    ReactomePA::enrichPathway(
      gene = selected_entrez,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe = background_entrez,
      minGSSize = min_gene_set_size,
      maxGSSize = max_gene_set_size,
      qvalueCutoff = 1,
      organism = "human",
      readable = TRUE
    ),
    error = function(e) {
      warning("Reactome enrichment failed: ", conditionMessage(e))
      NULL
    }
  )
} else {
  NULL
}

go_bp_result <- tryCatch(
  clusterProfiler::enrichGO(
    gene = selected_entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = background_entrez,
    qvalueCutoff = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    readable = TRUE,
    pool = FALSE
  ),
  error = function(e) {
    warning("GO:BP enrichment failed: ", conditionMessage(e))
    NULL
  }
)

ratio_to_numeric <- function(x) {
  x <- as.character(x)
  
  vapply(
    strsplit(x, "/", fixed = TRUE),
    FUN.VALUE = numeric(1),
    FUN = function(parts) {
      if (length(parts) != 2L) {
        return(NA_real_)
      }
      
      numerator <- suppressWarnings(as.numeric(parts[1]))
      denominator <- suppressWarnings(as.numeric(parts[2]))
      
      if (
        !is.finite(numerator) ||
        !is.finite(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }
      
      numerator / denominator
    }
  )
}

entrez_to_symbol <- setNames(
  all_symbol_entrez_map$SYMBOL,
  all_symbol_entrez_map$ENTREZID
)

translate_gene_id_field <- function(gene_id_string) {
  if (is.na(gene_id_string) || !nzchar(gene_id_string)) {
    return(NA_character_)
  }
  
  ids <- strsplit(as.character(gene_id_string), "/", fixed = TRUE)[[1]]
  ids <- ids[nzchar(ids)]
  
  if (length(ids) == 0L) {
    return(NA_character_)
  }
  
  # readable=TRUE already gives symbols for GO and Reactome. KEGG usually
  # returns numeric Entrez IDs, which are converted here.
  if (all(grepl("^[0-9]+$", ids))) {
    translated <- unname(entrez_to_symbol[ids])
    translated[is.na(translated) | !nzchar(translated)] <-
      ids[is.na(translated) | !nzchar(translated)]
    ids <- translated
  }
  
  paste(unique(ids), collapse = "/")
}

empty_enrichment_table <- function() {
  tibble(
    database = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    GeneRatio_numeric = numeric(),
    BgRatio_numeric = numeric(),
    FoldEnrichment = numeric(),
    gene_symbols = character()
  )
}

as_enrichment_table <- function(result_object, database_name) {
  if (is.null(result_object)) {
    return(empty_enrichment_table())
  }
  
  result_table <- as.data.frame(result_object) %>%
    as_tibble()
  
  if (nrow(result_table) == 0L) {
    return(empty_enrichment_table())
  }
  
  if (!"Count" %in% names(result_table)) {
    result_table$Count <- NA_integer_
  }
  if (!"GeneRatio" %in% names(result_table)) {
    result_table$GeneRatio <- NA_character_
  }
  if (!"BgRatio" %in% names(result_table)) {
    result_table$BgRatio <- NA_character_
  }
  if (!"geneID" %in% names(result_table)) {
    result_table$geneID <- NA_character_
  }
  
  result_table %>%
    mutate(
      database = database_name,
      Count = suppressWarnings(as.integer(Count)),
      GeneRatio_numeric = ratio_to_numeric(GeneRatio),
      BgRatio_numeric = ratio_to_numeric(BgRatio),
      FoldEnrichment = if_else(
        is.finite(GeneRatio_numeric) &
          is.finite(BgRatio_numeric) &
          BgRatio_numeric > 0,
        GeneRatio_numeric / BgRatio_numeric,
        NA_real_
      ),
      gene_symbols = vapply(
        geneID,
        translate_gene_id_field,
        character(1)
      )
    ) %>%
    relocate(database, .before = 1)
}

kegg_table <- as_enrichment_table(kegg_result, "KEGG")
reactome_table <- as_enrichment_table(reactome_result, "Reactome")
go_bp_table <- as_enrichment_table(go_bp_result, "GO:BP")

write_enrichment_outputs <- function(result_table, file_prefix, fdr_cutoff) {
  all_file <- file.path(
    enrichment_output_dir,
    paste0(file_prefix, "_all_results.csv")
  )
  significant_file <- file.path(
    enrichment_output_dir,
    paste0(file_prefix, "_BH_FDR05_results.csv")
  )
  
  readr::write_csv(result_table, all_file)
  
  significant_table <- result_table %>%
    filter(
      is.finite(p.adjust),
      p.adjust <= fdr_cutoff
    ) %>%
    arrange(p.adjust, desc(Count), Description)
  
  readr::write_csv(significant_table, significant_file)
  invisible(significant_table)
}

kegg_significant <- write_enrichment_outputs(
  kegg_table,
  "KEGG_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)
reactome_significant <- write_enrichment_outputs(
  reactome_table,
  "Reactome_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)
go_bp_significant <- write_enrichment_outputs(
  go_bp_table,
  "GO_BP_Q2_Q4_pooled",
  enrichment_fdr_cutoff
)

############################################################
# 15) FIGURE 5e-f STYLE BAR PLOTS FOR THE POOLED Q2 + Q4 SET
############################################################

wrap_plot_labels <- function(x, width = 34L) {
  vapply(
    as.character(x),
    FUN.VALUE = character(1),
    FUN = function(label) {
      paste(strwrap(label, width = width), collapse = "\n")
    }
  )
}

empty_enrichment_plot <- function(title, subtitle, database_name) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0(
        "No ", database_name,
        " terms met BH FDR ", less_than_or_equal_symbol,
        " ", sprintf("%.2f", enrichment_fdr_cutoff), "."
      ),
      size = 3.7,
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = element_text(
        hjust = 0,
        size = 8,
        face = "italic",
        color = "grey35"
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Define the symbol here because it was not needed in panels a-d.
less_than_or_equal_symbol <- "\u2264"

make_enrichment_barplot <- function(
    result_table,
    title,
    database_name,
    x_variable = c("Count", "FoldEnrichment"),
    x_axis_title,
    top_n = 12L,
    bar_fill = "#2B8CBE"
) {
  x_variable <- match.arg(x_variable)
  
  expression_floor_note <- if (
    isTRUE(apply_expression_floor_to_q2_q4_enrichment)
  ) {
    paste0(
      "; mean CPM ", greater_equal_symbol, " ",
      sprintf("%.2f", minimum_mean_cpm_panels_a_d)
    )
  } else {
    ""
  }
  
  plot_subtitle <- paste0(
    "ORA on ", scales::comma(length(selected_entrez)),
    " mapped Q2 + Q4 genes; background = ",
    scales::comma(length(background_entrez)),
    " eligible tested genes", expression_floor_note,
    "; BH FDR ", less_than_or_equal_symbol, " ",
    sprintf("%.2f", enrichment_fdr_cutoff)
  )
  
  plot_table <- result_table %>%
    filter(
      is.finite(p.adjust),
      p.adjust <= enrichment_fdr_cutoff,
      is.finite(.data[[x_variable]]),
      .data[[x_variable]] > 0
    ) %>%
    arrange(p.adjust, desc(Count), Description) %>%
    slice_head(n = top_n)
  
  if (nrow(plot_table) == 0L) {
    return(
      empty_enrichment_plot(
        title = title,
        subtitle = plot_subtitle,
        database_name = database_name
      )
    )
  }
  
  plot_table <- plot_table %>%
    arrange(.data[[x_variable]], desc(Count), Description) %>%
    mutate(
      plot_key = make.unique(paste0(Description, " [", ID, "]")),
      term_label = wrap_plot_labels(Description, width = 34L),
      plot_key = factor(plot_key, levels = plot_key)
    )
  
  y_label_lookup <- setNames(
    plot_table$term_label,
    as.character(plot_table$plot_key)
  )
  
  ggplot(
    plot_table,
    aes(x = .data[[x_variable]], y = plot_key)
  ) +
    geom_col(width = 0.72, fill = bar_fill) +
    scale_y_discrete(labels = y_label_lookup) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = scales::label_number(accuracy = if (x_variable == "Count") 1 else 0.1)
    ) +
    labs(
      title = title,
      subtitle = plot_subtitle,
      caption = paste0(
        "Top ", min(top_n, nrow(plot_table)),
        " adjusted-significant terms ordered by BH-adjusted p-value; ",
        "ties use overlapping-gene count."
      ),
      x = x_axis_title,
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 11),
      plot.subtitle = element_text(
        hjust = 0,
        size = 7.8,
        face = "italic",
        color = "grey35",
        margin = margin(b = 7)
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 6.8,
        color = "grey35",
        margin = margin(t = 7)
      ),
      axis.title.x = element_text(size = 8.5, margin = margin(t = 8)),
      axis.text.x = element_text(size = 7.6, color = "grey30"),
      axis.text.y = element_text(size = 7.2, color = "grey30"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.35,
        linetype = "dashed"
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

panel_e_kegg <- make_enrichment_barplot(
  result_table = kegg_table,
  title = "KEGG pathways (Q2 + Q4 pooled set)",
  database_name = "KEGG",
  x_variable = "Count",
  x_axis_title = "Overlapping genes",
  top_n = top_terms_kegg,
  bar_fill = "#2B8CBE"
) +
  labs(tag = "e") +
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

panel_f_go_bp <- make_enrichment_barplot(
  result_table = go_bp_table,
  title = "GO biological processes (Q2 + Q4 pooled set)",
  database_name = "GO:BP",
  x_variable = "FoldEnrichment",
  x_axis_title = "Fold enrichment",
  top_n = top_terms_go_bp,
  bar_fill = "#6A3D9A"
) +
  labs(tag = "f") +
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

panel_reactome <- if (has_reactomepa) {
  make_enrichment_barplot(
    result_table = reactome_table,
    title = "Reactome pathways (Q2 + Q4 pooled set)",
    database_name = "Reactome",
    x_variable = "Count",
    x_axis_title = "Overlapping genes",
    top_n = top_terms_reactome,
    bar_fill = "#3A9D5D"
  )
} else {
  blank_panel(
    title = "Reactome pathways (Q2 + Q4 pooled set)",
    message = paste0(
      "Reactome analysis was skipped because ReactomePA is not installed.\n",
      "KEGG and GO:BP outputs were still generated."
    )
  )
}

figure_5_ef_q2_q4 <- panel_e_kegg + panel_f_go_bp +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1))

figure_5_three_database_q2_q4 <-
  make_enrichment_barplot(
    result_table = kegg_table,
    title = "KEGG pathways",
    database_name = "KEGG",
    x_variable = "Count",
    x_axis_title = "Overlapping genes",
    top_n = top_terms_kegg,
    bar_fill = "#2B8CBE"
  ) +
  panel_reactome +
  make_enrichment_barplot(
    result_table = go_bp_table,
    title = "GO biological processes",
    database_name = "GO:BP",
    x_variable = "FoldEnrichment",
    x_axis_title = "Fold enrichment",
    top_n = top_terms_go_bp,
    bar_fill = "#6A3D9A"
  ) +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(
    title = "Q2 + Q4 pooled directional-reversal enrichment",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )
  )

# Two-panel Figure 5e-f analogue: KEGG plus GO:BP.
ggsave(
  file.path(
    enrichment_output_dir,
    "Figure5_ef_Q2_Q4_pooled_KEGG_GO_BP.png"
  ),
  figure_5_ef_q2_q4,
  width = 13.5,
  height = 5.6,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    enrichment_output_dir,
    "Figure5_ef_Q2_Q4_pooled_KEGG_GO_BP.pdf"
  ),
  figure_5_ef_q2_q4,
  width = 13.5,
  height = 5.6,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# Reactome is saved separately because the requested two-panel figure mirrors
# the attached KEGG + GO:BP structure.
ggsave(
  file.path(enrichment_output_dir, "Reactome_Q2_Q4_pooled.png"),
  panel_reactome,
  width = 6.8,
  height = 5.6,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(enrichment_output_dir, "Reactome_Q2_Q4_pooled.pdf"),
  panel_reactome,
  width = 6.8,
  height = 5.6,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# A three-panel version is also exported so KEGG, Reactome, and GO:BP can be
# inspected together without changing the requested two-panel Figure 5e-f.
ggsave(
  file.path(
    enrichment_output_dir,
    "Q2_Q4_pooled_KEGG_Reactome_GO_BP_three_panel.png"
  ),
  figure_5_three_database_q2_q4,
  width = 19.5,
  height = 5.8,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    enrichment_output_dir,
    "Q2_Q4_pooled_KEGG_Reactome_GO_BP_three_panel.pdf"
  ),
  figure_5_three_database_q2_q4,
  width = 19.5,
  height = 5.8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

message("Q2/Q4 ShinyGO lists and enrichment outputs written to:")
message("  ", enrichment_output_dir)
message(
  "Pooled Q2/Q4 symbols: ", length(q2_q4_symbols),
  "; mapped Entrez IDs: ", length(selected_entrez)
)
message(
  "Q2/Q4 enrichment expression floor applied: ",
  apply_expression_floor_to_q2_q4_enrichment,
  if (isTRUE(apply_expression_floor_to_q2_q4_enrichment)) {
    paste0(" (mean CPM >= ", sprintf("%.2f", minimum_mean_cpm_panels_a_d), ")")
  } else {
    ""
  }
)
message(
  "Eligible background symbols: ", length(background_symbols),
  "; mapped Entrez IDs: ", length(background_entrez)
)
message(
  "BH FDR <= ", enrichment_fdr_cutoff,
  " terms: KEGG=", nrow(kegg_significant),
  ", Reactome=",
  if (has_reactomepa) nrow(reactome_significant) else "skipped",
  ", GO:BP=", nrow(go_bp_significant)
)
#########################################################
############################################################
# FIGURE 5e-f FROM SHINYGO v0.85 EXPORTS
# Pooled Q2 + Q4 directional-reversal gene set
#
# Input files:
#   KEGG_pathways.csv (or .xlsx)
#   GO_BP.csv         (or .xlsx)
#
# Uses the Enrichment FDR values reported by ShinyGO and
# retains terms with FDR <= 0.05.
#
# Two versions are exported:
#   1) automatic: no biological keyword selection
#   2) neuron-focused: significant neuronal terms only,
#      with redundant GO terms reduced using gene-set Jaccard
#      overlap calculated from the ShinyGO "Genes" column.
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 1) PACKAGES
############################################################

required_packages <- c(
  "dplyr",
  "forcats",
  "ggplot2",
  "patchwork",
  "readr",
  "readxl",
  "scales",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\n\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(readxl)
  library(scales)
  library(stringr)
  library(tibble)
})

############################################################
# 2) USER SETTINGS
############################################################

project_dir <- "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq_2"

# Change these extensions to .xlsx when using Excel workbooks.
kegg_file <- file.path(project_dir, "KEGG_pathways_1.csv")
go_bp_file <- file.path(project_dir, "GO_BP_1.csv")

# Used only for .xlsx or .xls inputs.
kegg_excel_sheet <- 1
go_bp_excel_sheet <- 1

output_dir <- file.path(project_dir, "Figure5_ShinyGO_Q2_Q4_FDR005_1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
top_n_kegg <- 15L
top_n_go_bp <- 15L

# Automatic plot settings:
# KEGG is ranked by overlapping genes.
# GO:BP is ranked by fold enrichment. Requiring at least 10
# overlapping genes prevents very small 4- to 6-gene categories
# from dominating only because their fold enrichment is large.
go_minimum_overlap_automatic <- 15L

# Neuron-focused plot settings. The keyword filtering is applied
# only to pathway names already significant at FDR <= 0.05.
go_jaccard_cutoff <- 0.50

kegg_color <- "#2E8FAF"
go_color <- "#5B3C93"

############################################################
# 3) READ AND STANDARDIZE SHINYGO TABLES
############################################################

read_shinygo_file <- function(path, sheet = 1) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
  
  extension <- tolower(tools::file_ext(path))
  
  if (extension == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else if (extension %in% c("xlsx", "xls")) {
    readxl::read_excel(path, sheet = sheet)
  } else {
    stop(
      "Unsupported input extension for ", basename(path),
      ". Use .csv, .xlsx, or .xls.",
      call. = FALSE
    )
  }
}

standardize_shinygo <- function(x, database) {
  required_columns <- c(
    "Enrichment FDR",
    "nGenes",
    "Pathway Genes",
    "Fold Enrichment",
    "Pathway",
    "URL",
    "Genes"
  )
  
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      database,
      " input is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- x %>%
    dplyr::transmute(
      database = database,
      enrichment_fdr = suppressWarnings(as.numeric(.data[["Enrichment FDR"]])),
      overlapping_genes = suppressWarnings(as.integer(.data[["nGenes"]])),
      pathway_genes = suppressWarnings(as.integer(.data[["Pathway Genes"]])),
      fold_enrichment = suppressWarnings(as.numeric(.data[["Fold Enrichment"]])),
      pathway_raw = stringr::str_squish(as.character(.data[["Pathway"]])),
      url = as.character(.data[["URL"]]),
      genes = stringr::str_squish(as.character(.data[["Genes"]]))
    ) %>%
    dplyr::filter(
      is.finite(enrichment_fdr),
      is.finite(overlapping_genes),
      is.finite(pathway_genes),
      is.finite(fold_enrichment),
      nzchar(pathway_raw)
    )
  
  if (identical(database, "KEGG")) {
    out <- out %>%
      dplyr::mutate(
        term_id = stringr::str_extract(pathway_raw, "(?<=^Path:)hsa[0-9]+"),
        term = stringr::str_remove(pathway_raw, "^Path:hsa[0-9]+\\s+")
      )
  } else {
    out <- out %>%
      dplyr::mutate(
        term_id = stringr::str_extract(pathway_raw, "^GO:[0-9]+"),
        term = stringr::str_remove(pathway_raw, "^GO:[0-9]+\\s+")
      )
  }
  
  out %>%
    dplyr::mutate(term = stringr::str_squish(term)) %>%
    distinct(database, term_id, term, .keep_all = TRUE)
}

kegg_all <- read_shinygo_file(
  kegg_file,
  sheet = kegg_excel_sheet
) %>%
  standardize_shinygo(database = "KEGG")

go_bp_all <- read_shinygo_file(
  go_bp_file,
  sheet = go_bp_excel_sheet
) %>%
  standardize_shinygo(database = "GO:BP")

kegg_fdr05 <- kegg_all %>%
  dplyr::filter(enrichment_fdr <= fdr_cutoff)

go_bp_fdr05 <- go_bp_all %>%
  dplyr::filter(enrichment_fdr <= fdr_cutoff)

############################################################
# 4) OBJECTIVE AUTOMATIC TERM SELECTION
############################################################

# No biological keyword filter is used here.
kegg_automatic <- kegg_fdr05 %>%
  dplyr::arrange(
    desc(overlapping_genes),
    enrichment_fdr,
    desc(fold_enrichment),
    term
  ) %>%
  slice_head(n = top_n_kegg)

go_bp_automatic <- go_bp_fdr05 %>%
  dplyr::filter(overlapping_genes >= go_minimum_overlap_automatic) %>%
  dplyr::arrange(
    desc(fold_enrichment),
    enrichment_fdr,
    desc(overlapping_genes),
    term
  ) %>%
  slice_head(n = top_n_go_bp)

############################################################
# 5) OPTIONAL NEURON-FOCUSED, REDUNDANCY-REDUCED SELECTION
############################################################

kegg_neuron_pattern <- stringr::regex(
  paste(
    c(
      "neuro",
      "alzheimer",
      "parkinson",
      "huntington",
      "amyotrophic",
      "prion",
      "motor",
      "synap",
      "axon",
      "guidance",
      "signaling",
      "alcoholism"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

go_neuron_pattern <- stringr::regex(
  paste(
    c(
      "neuron",
      "neural",
      "nervous system",
      "neurogenesis",
      "axon",
      "dendrit",
      "synap",
      "brain",
      "forebrain",
      "cell projection",
      "neurotrophin"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

split_gene_set <- function(x) {
  x <- stringr::str_squish(x)
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  unique(strsplit(x, "\\s+")[[1]])
}

gene_set_jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  union_set <- union(a, b)
  if (length(union_set) == 0L) {
    return(0)
  }
  length(intersect(a, b)) / length(union_set)
}

reduce_gene_set_redundancy <- function(x, n, cutoff = 0.50) {
  if (nrow(x) == 0L) {
    return(x)
  }
  
  selected_rows <- integer(0)
  selected_gene_sets <- list()
  
  for (i in seq_len(nrow(x))) {
    current_genes <- split_gene_set(x$genes[[i]])
    
    retain_current <- if (length(selected_gene_sets) == 0L) {
      TRUE
    } else {
      all(
        vapply(
          selected_gene_sets,
          function(previous_genes) {
            gene_set_jaccard(current_genes, previous_genes) < cutoff
          },
          logical(1)
        )
      )
    }
    
    if (retain_current) {
      selected_rows <- c(selected_rows, i)
      selected_gene_sets[[length(selected_gene_sets) + 1L]] <- current_genes
    }
    
    if (length(selected_rows) >= n) {
      break
    }
  }
  
  x[selected_rows, , drop = FALSE]
}

kegg_neuron_focused <- kegg_fdr05 %>%
  dplyr::filter(
    stringr::str_detect(term, kegg_neuron_pattern),
    !stringr::str_detect(
      term,
      stringr::regex("non-alcoholic", ignore_case = TRUE)
    )
  ) %>%
  dplyr::arrange(
    desc(overlapping_genes),
    enrichment_fdr,
    desc(fold_enrichment),
    term
  ) %>%
  slice_head(n = top_n_kegg)

go_bp_neuron_ranked <- go_bp_fdr05 %>%
  dplyr::filter(stringr::str_detect(term, go_neuron_pattern)) %>%
  dplyr::arrange(
    enrichment_fdr,
    desc(fold_enrichment),
    desc(overlapping_genes),
    term
  )

go_bp_neuron_focused <- reduce_gene_set_redundancy(
  go_bp_neuron_ranked,
  n = top_n_go_bp,
  cutoff = go_jaccard_cutoff
)

############################################################
# 6) PLOTTING FUNCTIONS
############################################################

figure_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      color = "grey88",
      linewidth = 0.45,
      linetype = "dashed"
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_text(
      color = "grey35",
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(
      color = "grey30",
      size = 8.1,
      lineheight = 0.92
    ),
    axis.text.x = element_text(color = "grey35", size = 8),
    plot.title = element_text(
      size = 12,
      face = "plain",
      color = "grey12",
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 8.5,
      face = "italic",
      color = "grey42",
      margin = margin(b = 10)
    ),
    plot.tag = element_text(
      face = "bold",
      size = 12,
      color = "grey12"
    ),
    plot.tag.position = c(-0.07, 1.08),
    plot.background = element_rect(
      fill = "white",
      color = "grey84",
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 12, 10, 14)
  )

make_enrichment_panel <- function(
    x,
    value_column,
    title,
    subtitle,
    x_label,
    fill_color,
    panel_tag,
    wrap_width = 34) {
  
  if (nrow(x) == 0L) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No terms met the plotting criteria",
          size = 4
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
          title = title,
          subtitle = subtitle,
          tag = panel_tag
        ) +
        theme_void() +
        theme(
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8.5, face = "italic"),
          plot.tag = element_text(face = "bold", size = 12),
          plot.background = element_rect(
            fill = "white",
            color = "grey84"
          )
        )
    )
  }
  
  plot_table <- x %>%
    dplyr::mutate(
      plotting_value = .data[[value_column]],
      term_wrapped = stringr::str_wrap(term, width = wrap_width)
    ) %>%
    dplyr::arrange(plotting_value, term) %>%
    dplyr::mutate(
      term_wrapped = factor(
        term_wrapped,
        levels = unique(term_wrapped)
      )
    )
  
  ggplot(
    plot_table,
    aes(x = plotting_value, y = term_wrapped)
  ) +
    geom_col(
      width = 0.72,
      fill = fill_color
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      tag = panel_tag
    ) +
    figure_theme
}

make_two_panel_figure <- function(
    kegg_table,
    go_table,
    kegg_subtitle,
    go_subtitle) {
  
  panel_e <- make_enrichment_panel(
    kegg_table,
    value_column = "overlapping_genes",
    title = "KEGG pathways (Q2 + Q4 pooled set)",
    subtitle = kegg_subtitle,
    x_label = "Overlapping genes",
    fill_color = kegg_color,
    panel_tag = "e"
  )
  
  panel_f <- make_enrichment_panel(
    go_table,
    value_column = "fold_enrichment",
    title = "GO biological processes (Q2 + Q4 pooled set)",
    subtitle = go_subtitle,
    x_label = "Fold enrichment",
    fill_color = go_color,
    panel_tag = "f"
  )
  
  panel_e + panel_f +
    patchwork::plot_layout(
      ncol = 2,
      widths = c(1, 1)
    )
}

############################################################
# 7) CREATE AND SAVE BOTH VERSIONS
############################################################

automatic_figure <- make_two_panel_figure(
  kegg_automatic,
  go_bp_automatic,
  kegg_subtitle = paste0(
    "ShinyGO v0.85; FDR <= ",
    fdr_cutoff,
    "; top terms ranked by overlapping genes"
  ),
  go_subtitle = paste0(
    "ShinyGO v0.85; FDR <= ",
    fdr_cutoff,
    "; overlap >= ",
    go_minimum_overlap_automatic,
    "; ranked by fold enrichment"
  )
)

neuron_focused_figure <- make_two_panel_figure(
  kegg_neuron_focused,
  go_bp_neuron_focused,
  kegg_subtitle = paste0(
    "ShinyGO v0.85; FDR <= ",
    fdr_cutoff,
    "; significant neuron-related terms"
  ),
  go_subtitle = paste0(
    "ShinyGO v0.85; FDR <= ",
    fdr_cutoff,
    "; neuron-focused; Jaccard < ",
    sprintf("%.2f", go_jaccard_cutoff)
  )
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure5_ef_ShinyGO_Q2_Q4_FDR005_automatic.png"
  ),
  plot = automatic_figure,
  width = 13.8,
  height = 5.5,
  dpi = 400,
  bg = "#F7F8FA"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure5_ef_ShinyGO_Q2_Q4_FDR005_automatic.pdf"
  ),
  plot = automatic_figure,
  width = 13.8,
  height = 5.5,
  bg = "#F7F8FA"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure5_ef_ShinyGO_Q2_Q4_FDR005_neuron_focused.png"
  ),
  plot = neuron_focused_figure,
  width = 13.8,
  height = 5.5,
  dpi = 400,
  bg = "#F7F8FA"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure5_ef_ShinyGO_Q2_Q4_FDR005_neuron_focused.pdf"
  ),
  plot = neuron_focused_figure,
  width = 13.8,
  height = 5.5,
  bg = "#F7F8FA"
)

############################################################
# 8) EXPORT AUDIT TABLES
############################################################

all_significant_terms <- bind_rows(
  kegg_fdr05,
  go_bp_fdr05
) %>%
  dplyr::select(
    database,
    term_id,
    term,
    enrichment_fdr,
    overlapping_genes,
    pathway_genes,
    fold_enrichment,
    url,
    genes
  ) %>%
  dplyr::arrange(database, enrichment_fdr, desc(fold_enrichment))

terms_used_in_plots <- bind_rows(
  kegg_automatic %>% dplyr::mutate(plot_mode = "automatic"),
  go_bp_automatic %>% dplyr::mutate(plot_mode = "automatic"),
  kegg_neuron_focused %>% dplyr::mutate(plot_mode = "neuron_focused"),
  go_bp_neuron_focused %>% dplyr::mutate(plot_mode = "neuron_focused")
) %>%
  dplyr::select(
    plot_mode,
    database,
    term_id,
    term,
    enrichment_fdr,
    overlapping_genes,
    pathway_genes,
    fold_enrichment,
    url,
    genes
  )

readr::write_csv(
  all_significant_terms,
  file.path(output_dir, "ShinyGO_Q2_Q4_all_terms_FDR_le_0.05.csv")
)

readr::write_csv(
  terms_used_in_plots,
  file.path(output_dir, "ShinyGO_Q2_Q4_terms_used_in_plots.csv")
)

message("ShinyGO plotting complete.")
message("Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message(
  "Significant terms at FDR <= ",
  fdr_cutoff,
  ": KEGG = ",
  nrow(kegg_fdr05),
  "; GO:BP = ",
  nrow(go_bp_fdr05)
)
message(
  "Automatic plot terms: KEGG = ",
  nrow(kegg_automatic),
  "; GO:BP = ",
  nrow(go_bp_automatic)
)
message(
  "Neuron-focused plot terms: KEGG = ",
  nrow(kegg_neuron_focused),
  "; GO:BP = ",
  nrow(go_bp_neuron_focused)
)
#############################################################################
