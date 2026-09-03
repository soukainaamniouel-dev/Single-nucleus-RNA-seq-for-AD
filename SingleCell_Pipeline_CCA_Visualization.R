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
# PUBLICATION-QUALITY CELL-TYPE COMPOSITION PLOTS
#
# Inputs:
#   1) CCA object with clusters 4 and 5 merged
#   2) CCA object with clusters 4 and 5 separate
#
# Required metadata in each Seurat object:
#   sample
#   cell_type
#
# Outputs for each object:
#   - Long and wide count/proportion tables
#   - 100% stacked composition plot
#   - Grouped proportion plot
#   - Grouped count plot
#   - Dot plot: size = count, color = proportion
#   - Main two-panel figure
#   - Vector PDF and 1000-dpi TIFF versions
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

############################################################
# 1) INPUT FILES
############################################################

input_dir <- paste0(
  "H:/Documents/Qi_Projects/AD_scRNAseq/AD_Qi_scRNAseq/",
  "output_new_process_scRNAseq/"
)

merged_file <- file.path(
  input_dir,
  paste0(
    "integrated_CCA_removed_8_9_10_12_14_merged_4_5_",
    "celltypes_annotated.qs"
  )
)

separate_file <- file.path(
  input_dir,
  paste0(
    "integrated_CCA_removed_clusters_8_9_10_12_14_",
    "celltypes_annotated.qs"
  )
)

if (!file.exists(merged_file)) {
  message(
    "Select the annotated QS object with clusters 4 and 5 merged."
  )
  merged_file <- file.choose()
}

if (!file.exists(separate_file)) {
  message(
    "Select the annotated QS object with clusters 4 and 5 separate."
  )
  separate_file <- file.choose()
}

############################################################
# 2) DISPLAY SETTINGS
############################################################

# Sample order used in all figures.
sample_order <- c(
  "Ctrl",
  "Fibrils",
  "FibJ8"
)

# Publication labels displayed on figure axes and legends.
sample_display_labels <- c(
  "Ctrl" = "Control",
  "Fibrils" = "Fibrils",
  "FibJ8" = "Fibrils + J8"
)

# Biological display order for cell types.
# Any unexpected labels found in an object are appended automatically.
cell_type_order <- c(
  "L6 IT ExN",
  "Migrating cortical ExN",
  "Mid-layer IT ExN",
  "Immature upper-layer ExN",
  "L6 CT-like ExN",
  "L5 deep-layer ExN",
  "GABAergic neurons",
  "Developing GABAergic neurons",
  "Immature GABAergic neurons",
  "Early OPC",
  "Developmental astrocytes",
  "Perivascular mesenchymal"
)

# Consistent cell-type colors across both annotation versions.
cell_type_palette <- c(
  "L6 IT ExN" = "#3C5488",
  "Migrating cortical ExN" = "#4DBBD5",
  "Mid-layer IT ExN" = "#8491B4",
  "Immature upper-layer ExN" = "#91D1C2",
  "L6 CT-like ExN" = "#6A3D9A",
  "L5 deep-layer ExN" = "#1F78B4",
  "GABAergic neurons" = "#E64B35",
  "Developing GABAergic neurons" = "#F39B7F",
  "Immature GABAergic neurons" = "#DC0000",
  "Early OPC" = "#00A087",
  "Developmental astrocytes" = "#FDBF6F",
  "Perivascular mesenchymal" = "#7E6148"
)

# Three high-contrast sample colors.
sample_palette <- c(
  "Ctrl" = "#4D4D4D",
  "Fibrils" = "#E64B35",
  "FibJ8" = "#3C5488"
)

# Arial is commonly accepted for journal figures on Windows.
# Change to "sans" if Arial is unavailable on your system.
figure_font <- "Arial"

# Raster line-art resolution. The vector PDF is the preferred
# submission file; the high-resolution TIFF is an alternate.
raster_dpi <- 1000L

############################################################
# 3) JOURNAL-STYLE THEME AND SAVE HELPER
############################################################

theme_cell_reports <- function(
    base_size = 8,
    base_family = figure_font
) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 10,
        face = "bold",
        hjust = 0,
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8,
        color = "#333333",
        margin = ggplot2::margin(b = 7)
      ),
      axis.title = ggplot2::element_text(
        size = 8.5,
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        size = 7.5,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.35,
        color = "black"
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.4,
        color = "black"
      ),
      legend.title = ggplot2::element_text(
        size = 8,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 7
      ),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width = grid::unit(0.34, "cm"),
      legend.spacing.y = grid::unit(0.04, "cm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = 8,
        face = "bold"
      ),
      plot.margin = ggplot2::margin(7, 9, 7, 7)
    )
}

save_publication_plot <- function(
    plot,
    file_stub,
    width,
    height
) {
  pdf_file <- paste0(file_stub, ".pdf")
  tiff_file <- paste0(file_stub, ".tiff")
  
  # Vector PDF. Use Cairo when available for cleaner font rendering.
  tryCatch(
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      limitsize = FALSE
    ),
    error = function(e) {
      warning(
        "Cairo PDF failed; saving with the standard PDF device: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      ggplot2::ggsave(
        filename = pdf_file,
        plot = plot,
        device = "pdf",
        width = width,
        height = height,
        units = "in",
        bg = "white",
        limitsize = FALSE
      )
    }
  )
  
  # High-resolution raster copy for journal submission systems.
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    units = "in",
    dpi = raster_dpi,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

############################################################
# 4) PALETTE HELPER
############################################################

get_cell_type_palette <- function(cell_types) {
  palette_out <- cell_type_palette[cell_types]
  
  missing_color <- is.na(palette_out)
  
  if (any(missing_color)) {
    fallback_colors <- grDevices::hcl.colors(
      n = sum(missing_color),
      palette = "Dark 3"
    )
    
    palette_out[missing_color] <- fallback_colors
  }
  
  names(palette_out) <- cell_types
  palette_out
}

############################################################
# 5) MAKE TABLES AND PLOTS FOR ONE ANNOTATED OBJECT
############################################################

make_composition_plots <- function(
    input_file,
    output_stub,
    dataset_title
) {
  message("\nLoading: ", input_file)
  
  obj <- qs::qread(input_file)
  
  if (!inherits(obj, "Seurat")) {
    stop(
      "The selected QS file does not contain a Seurat object: ",
      input_file
    )
  }
  
  metadata <- obj[[]]
  
  required_columns <- c(
    "sample",
    "cell_type"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required metadata column(s) in ",
      basename(input_file),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  sample_values <- as.character(
    metadata$sample
  )
  
  cell_type_values <- as.character(
    metadata$cell_type
  )
  
  if (
    anyNA(sample_values) ||
    any(sample_values == "")
  ) {
    stop(
      "The sample column contains missing or empty values in ",
      basename(input_file)
    )
  }
  
  if (
    anyNA(cell_type_values) ||
    any(cell_type_values == "")
  ) {
    stop(
      "The cell_type column contains missing or empty values in ",
      basename(input_file)
    )
  }
  
  observed_samples <- unique(
    sample_values
  )
  
  sample_levels <- c(
    sample_order[
      sample_order %in% observed_samples
    ],
    setdiff(
      observed_samples,
      sample_order
    )
  )
  
  observed_cell_types <- unique(
    cell_type_values
  )
  
  cell_type_levels <- c(
    cell_type_order[
      cell_type_order %in% observed_cell_types
    ],
    setdiff(
      observed_cell_types,
      cell_type_order
    )
  )
  
  plot_cell_type_palette <- get_cell_type_palette(
    cell_type_levels
  )
  
  plot_sample_palette <- sample_palette[
    sample_levels
  ]
  
  missing_sample_colors <- is.na(
    plot_sample_palette
  )
  
  if (any(missing_sample_colors)) {
    plot_sample_palette[missing_sample_colors] <- grDevices::hcl.colors(
      n = sum(missing_sample_colors),
      palette = "Set 2"
    )
  }
  
  names(plot_sample_palette) <- sample_levels
  
  plot_sample_labels <- sample_display_labels[
    sample_levels
  ]
  
  missing_sample_labels <- is.na(
    plot_sample_labels
  )
  
  plot_sample_labels[missing_sample_labels] <- sample_levels[
    missing_sample_labels
  ]
  
  names(plot_sample_labels) <- sample_levels
  
  raw_counts <- data.frame(
    sample = sample_values,
    cell_type = cell_type_values,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::count(
      sample,
      cell_type,
      name = "n_cells"
    )
  
  complete_grid <- tidyr::expand_grid(
    sample = sample_levels,
    cell_type = cell_type_levels
  )
  
  composition <- complete_grid %>%
    dplyr::left_join(
      raw_counts,
      by = c(
        "sample",
        "cell_type"
      )
    ) %>%
    dplyr::mutate(
      n_cells = tidyr::replace_na(
        n_cells,
        0L
      )
    ) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      total_cells = sum(n_cells),
      proportion = n_cells / total_cells,
      percent = 100 * proportion
    ) %>%
    dplyr::ungroup()
  
  composition$sample <- factor(
    composition$sample,
    levels = sample_levels
  )
  
  composition$cell_type <- factor(
    composition$cell_type,
    levels = cell_type_levels
  )
  
  sample_totals <- composition %>%
    dplyr::distinct(
      sample,
      total_cells
    )
  
  output_dir <- file.path(
    dirname(input_file),
    "cell_type_composition_publication",
    output_stub
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ##########################################################
  # Export count and proportion tables
  ##########################################################
  
  write.csv(
    composition,
    file.path(
      output_dir,
      paste0(output_stub, "_cell_type_counts_and_proportions_long.csv")
    ),
    row.names = FALSE
  )
  
  counts_wide <- composition %>%
    dplyr::select(
      cell_type,
      sample,
      n_cells
    ) %>%
    tidyr::pivot_wider(
      names_from = sample,
      values_from = n_cells
    )
  
  proportions_wide <- composition %>%
    dplyr::select(
      cell_type,
      sample,
      percent
    ) %>%
    tidyr::pivot_wider(
      names_from = sample,
      values_from = percent
    )
  
  write.csv(
    counts_wide,
    file.path(
      output_dir,
      paste0(output_stub, "_cell_type_counts_wide.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    proportions_wide,
    file.path(
      output_dir,
      paste0(output_stub, "_cell_type_percentages_wide.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    sample_totals,
    file.path(
      output_dir,
      paste0(output_stub, "_sample_total_cells.csv")
    ),
    row.names = FALSE
  )
  
  ##########################################################
  # Plot 1: 100% stacked composition bars
  ##########################################################
  
  p_stacked <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = sample,
      y = proportion,
      fill = cell_type
    )
  ) +
    ggplot2::geom_col(
      width = 0.72,
      color = "white",
      linewidth = 0.22
    ) +
    ggplot2::geom_text(
      data = sample_totals,
      mapping = ggplot2::aes(
        x = sample,
        y = 1.04,
        label = paste0(
          "n = ",
          scales::comma(total_cells)
        )
      ),
      inherit.aes = FALSE,
      family = figure_font,
      size = 2.5
    ) +
    ggplot2::scale_x_discrete(
      labels = plot_sample_labels
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.09),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = ggplot2::expansion(
        mult = c(0, 0)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = plot_cell_type_palette,
      breaks = cell_type_levels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = dataset_title,
      subtitle = "Cell-type composition within each sample",
      x = NULL,
      y = "Cells within sample (%)",
      fill = "Cell type"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(
          color = NA
        )
      )
    ) +
    theme_cell_reports() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(
        face = "bold"
      )
    )
  
  ##########################################################
  # Plot 2: grouped proportions by cell type and sample
  ##########################################################
  
  p_proportion <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = cell_type,
      y = proportion,
      fill = sample
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(
        width = 0.82,
        preserve = "single",
        padding = 0.08
      ),
      width = 0.72,
      color = "black",
      linewidth = 0.14
    ) +
    ggplot2::coord_flip(
      clip = "off"
    ) +
    ggplot2::scale_x_discrete(
      limits = rev(cell_type_levels),
      labels = scales::wrap_format(
        width = 30
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = plot_sample_palette,
      breaks = sample_levels,
      labels = plot_sample_labels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = dataset_title,
      subtitle = "Within-sample cell-type proportions",
      x = NULL,
      y = "Cells within sample (%)",
      fill = "Sample"
    ) +
    theme_cell_reports() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )
  
  ##########################################################
  # Plot 3: grouped absolute cell counts
  ##########################################################
  
  p_count <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = cell_type,
      y = n_cells,
      fill = sample
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(
        width = 0.82,
        preserve = "single",
        padding = 0.08
      ),
      width = 0.72,
      color = "black",
      linewidth = 0.14
    ) +
    ggplot2::coord_flip(
      clip = "off"
    ) +
    ggplot2::scale_x_discrete(
      limits = rev(cell_type_levels),
      labels = scales::wrap_format(
        width = 30
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma_format(),
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = plot_sample_palette,
      breaks = sample_levels,
      labels = plot_sample_labels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = dataset_title,
      subtitle = "Absolute retained-cell counts; not normalized for library size",
      x = NULL,
      y = "Number of cells",
      fill = "Sample"
    ) +
    theme_cell_reports() +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )
  
  ##########################################################
  # Plot 4: count/proportion dot plot
  ##########################################################
  
  max_count <- max(
    composition$n_cells,
    na.rm = TRUE
  )
  
  size_breaks <- pretty(
    c(0, max_count),
    n = 4
  )
  
  size_breaks <- size_breaks[
    size_breaks > 0 &
      size_breaks <= max_count
  ]
  
  if (length(size_breaks) == 0L) {
    size_breaks <- max_count
  }
  
  p_dot <- ggplot2::ggplot(
    composition,
    ggplot2::aes(
      x = sample,
      y = cell_type
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = n_cells,
        color = proportion
      ),
      alpha = 0.95
    ) +
    ggplot2::scale_x_discrete(
      labels = plot_sample_labels
    ) +
    ggplot2::scale_y_discrete(
      limits = rev(cell_type_levels),
      labels = scales::wrap_format(
        width = 30
      )
    ) +
    ggplot2::scale_size_area(
      max_size = 10,
      breaks = size_breaks,
      labels = scales::comma_format(),
      name = "Cell count"
    ) +
    ggplot2::scale_color_viridis_c(
      option = "C",
      begin = 0.10,
      end = 0.92,
      labels = scales::percent_format(
        accuracy = 1
      ),
      name = "Proportion\nwithin sample"
    ) +
    ggplot2::labs(
      title = dataset_title,
      subtitle = "Dot size represents cell count; color represents within-sample proportion",
      x = NULL,
      y = NULL
    ) +
    theme_cell_reports() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "right"
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        order = 1
      ),
      color = ggplot2::guide_colorbar(
        order = 2,
        barheight = grid::unit(2.4, "cm"),
        barwidth = grid::unit(0.35, "cm")
      )
    )
  
  ##########################################################
  # Main two-panel composition figure
  ##########################################################
  
  # Remove internal titles from the assembled journal figure.
  # The figure title and explanation belong in the manuscript legend.
  p_stacked_journal <- p_stacked +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL
    )
  
  p_dot_journal <- p_dot +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL
    )
  
  p_main <- (
    p_stacked_journal /
      p_dot_journal
  ) +
    patchwork::plot_layout(
      heights = c(0.90, 1.35)
    ) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(
          family = figure_font,
          face = "bold",
          size = 11
        )
      )
    )
  
  ##########################################################
  # Save publication files
  ##########################################################
  
  save_publication_plot(
    plot = p_stacked,
    file_stub = file.path(
      output_dir,
      paste0(output_stub, "_stacked_cell_type_proportions")
    ),
    width = 7.2,
    height = 5.3
  )
  
  save_publication_plot(
    plot = p_proportion,
    file_stub = file.path(
      output_dir,
      paste0(output_stub, "_grouped_cell_type_proportions")
    ),
    width = 7.2,
    height = 6.4
  )
  
  save_publication_plot(
    plot = p_count,
    file_stub = file.path(
      output_dir,
      paste0(output_stub, "_grouped_cell_type_counts")
    ),
    width = 7.2,
    height = 6.4
  )
  
  save_publication_plot(
    plot = p_dot,
    file_stub = file.path(
      output_dir,
      paste0(output_stub, "_cell_type_count_proportion_dotplot")
    ),
    width = 7.2,
    height = 6.5
  )
  
  save_publication_plot(
    plot = p_main,
    file_stub = file.path(
      output_dir,
      paste0(output_stub, "_main_cell_type_composition_figure")
    ),
    width = 8.4,
    height = 10.0
  )
  
  message(
    "Saved composition plots and tables to: ",
    normalizePath(
      output_dir,
      winslash = "/",
      mustWork = FALSE
    )
  )
  
  # Keep plot data and plots, but release the large Seurat object.
  rm(obj)
  invisible(gc())
  
  list(
    composition = composition,
    sample_totals = sample_totals,
    cell_type_levels = cell_type_levels,
    cell_type_palette = plot_cell_type_palette,
    p_stacked = p_stacked,
    p_dot = p_dot,
    p_main = p_main,
    output_dir = output_dir
  )
}

############################################################
# 6) MERGED 4/5 OBJECT
############################################################

merged_results <- make_composition_plots(
  input_file = merged_file,
  output_stub = "CCA_removed_8_9_10_12_14_merged_4_5",
  dataset_title = "Cell-type composition: clusters 4 and 5 merged"
)

############################################################
# 7) SEPARATE 4/5 OBJECT
############################################################

separate_results <- make_composition_plots(
  input_file = separate_file,
  output_stub = "CCA_removed_8_9_10_12_14_clusters_4_5_separate",
  dataset_title = "Cell-type composition: clusters 4 and 5 separate"
)

############################################################
# 8) OPTIONAL SENSITIVITY FIGURE: BOTH ANNOTATION VERSIONS
############################################################

comparison_data <- dplyr::bind_rows(
  merged_results$composition %>%
    dplyr::mutate(
      annotation_version = "Clusters 4 and 5 merged"
    ),
  separate_results$composition %>%
    dplyr::mutate(
      annotation_version = "Clusters 4 and 5 separate"
    )
)

comparison_cell_types <- unique(
  c(
    merged_results$cell_type_levels,
    separate_results$cell_type_levels
  )
)

comparison_palette <- get_cell_type_palette(
  comparison_cell_types
)

comparison_data$annotation_version <- factor(
  comparison_data$annotation_version,
  levels = c(
    "Clusters 4 and 5 separate",
    "Clusters 4 and 5 merged"
  )
)

comparison_data$cell_type <- factor(
  as.character(comparison_data$cell_type),
  levels = comparison_cell_types
)

p_annotation_sensitivity <- ggplot2::ggplot(
  comparison_data,
  ggplot2::aes(
    x = sample,
    y = proportion,
    fill = cell_type
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.20
  ) +
  ggplot2::facet_wrap(
    ~annotation_version,
    ncol = 1
  ) +
  ggplot2::scale_x_discrete(
    labels = sample_display_labels
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.01)
    )
  ) +
  ggplot2::scale_fill_manual(
    values = comparison_palette,
    breaks = comparison_cell_types,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Cell-type composition under alternative cluster 4/5 annotations",
    subtitle = "The same curated cells are shown; only the treatment of clusters 4 and 5 differs",
    x = NULL,
    y = "Cells within sample (%)",
    fill = "Cell type"
  ) +
  theme_cell_reports() +
  ggplot2::theme(
    legend.position = "right",
    axis.text.x = ggplot2::element_text(
      face = "bold"
    )
  )

comparison_output_dir <- file.path(
  dirname(merged_file),
  "cell_type_composition_publication",
  "annotation_sensitivity_comparison"
)

dir.create(
  comparison_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_publication_plot(
  plot = p_annotation_sensitivity,
  file_stub = file.path(
    comparison_output_dir,
    "CCA_cell_type_composition_4_5_annotation_sensitivity"
  ),
  width = 8.4,
  height = 8.2
)

write.csv(
  comparison_data,
  file.path(
    comparison_output_dir,
    "CCA_cell_type_composition_4_5_annotation_sensitivity_data.csv"
  ),
  row.names = FALSE
)

message("\nAll publication-quality composition plots are complete.")
message("Merged-object output: ", merged_results$output_dir)
message("Separate-object output: ", separate_results$output_dir)
message("Comparison output: ", comparison_output_dir)
################################################################################
############################################################
# CELL REPORTS-STYLE t-SNE, STANDALONE LEGEND, AND HEATMAP
#
# Input:
#   A cell-type-annotated CCA Seurat object saved as .qs.
#   The object must contain:
#     - cell_type metadata
#     - a stored PCA reduction
#     - an RNA assay
#
# The script works with either of these annotation versions:
#   1) clusters 4 and 5 merged as "GABAergic neurons"
#   2) clusters 4 and 5 kept separate as developing and immature
#      GABAergic neurons
#
# Primary outputs are saved separately as PDF and TIFF:
#   1) t-SNE
#   2) cell-type legend/key
#   3) compact marker heatmap
#
# An optional all-marker heatmap is also produced.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(ggplot2)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

############################################################
# 1) USER SETTINGS
############################################################

set.seed(1234)

# Select the annotated QS object when the file browser opens.
# Examples:
#   integrated_CCA_removed_8_9_10_12_14_merged_4_5_celltypes_annotated.qs
#   integrated_CCA_removed_clusters_8_9_10_12_14_celltypes_annotated.qs
qs_file <- file.choose()

# t-SNE settings. The original CCA workflow used PCs 1:20.
pcs_requested <- 1:20
tsne_seed <- 1234

# Number of marker genes displayed per cell type in the compact heatmap.
# These are selected from the user-supplied marker list using the data.
top_markers_per_cell_type <- 6L

# Also create a large heatmap containing all available marker genes.
make_all_marker_heatmap <- TRUE

# Save a new large QS object containing the new t-SNE reduction.
# FALSE avoids writing another approximately gigabyte-sized object.
save_object_with_tsne <- FALSE

# Use a generic sans-serif family, which maps to Arial/Helvetica on most systems.
font_family <- "sans"

# Main raster output resolution. Elsevier treats color heatmaps and
# dimensional-reduction plots as combination artwork; 600 dpi exceeds
# the usual 500 dpi minimum.
main_tiff_dpi <- 600
legend_tiff_dpi <- 1000

############################################################
# 2) USER-SUPPLIED MARKER LIST
############################################################

neuron_markers <- list(
  L6_IT_ExN = c(
    "STK32B", "CDH13", "MCTP1", "DOK5", "NEO1"
  ),
  
  Migrating_Cortical_ExN = c(
    "DPY19L1", "PRSS12", "PLXNA2", "SEMA3C", "PTPN4",
    "FRMD4B", "ST3GAL6"
  ),
  
  Perivascular_Mesenchymal = c(
    "PDGFRB", "PDE3A", "TNC", "FBN1", "COL11A1", "PTN",
    "PTPRM", "IL33", "GLI3", "LIPG", "SLCO1C1", "INTU",
    "DYNC2H1"
  ),
  
  Mid_layer_IT_ExN = c(
    "LRRC4C", "RORB", "GRM7", "KCNH5", "NTNG1", "CDH18",
    "TENM1", "LUZP2"
  ),
  
  Developing_GABAergic_neurons = c(
    "GAD1", "PDZRN3", "NRXN3", "GAD2", "DLX6-AS1", "ST18",
    "DCLK2"
  ),
  
  Immature_GABAergic_neurons = c(
    "COPG2IT1", "PAX6", "MEST", "ST18", "DLX6-AS1", "GAD2",
    "ARX", "NRXN3", "GAD1", "SOX2-OT"
  ),
  
  Early_OPC = c(
    "OLIG1", "FERMT1", "SMOC1", "PDGFRA", "PLP1", "SOX10",
    "CLDN11", "MBP", "DLL3", "CSPG4", "SLC24A3", "CEROX1",
    "MAP3K1", "POLR2F"
  ),
  
  Immature_upper_layer_ExN = c(
    "UNC5D", "SOX11", "DCC", "NRP1", "LRP8", "EPHA3", "DOK6",
    "KCNQ3", "CNR1", "CDH4"
  ),
  
  L6_CT_like_ExN = c(
    "HS3ST4", "TRPM3", "TLE4", "FOXP2", "ZFPM2", "SOX5",
    "ADAMTSL1", "SEMA3E"
  ),
  
  Developmental_astrocytes = c(
    "AQP4", "SLC4A4", "GFAP", "MEGF10", "SPARCL1", "CD44",
    "GLIS3", "RFX4", "ID4", "EDNRB", "LIFR", "LAMA1",
    "ITGA6", "NID1"
  ),
  
  L5_deep_layer_ExN = c(
    "SGCZ", "PEX5L", "RALYL", "CD36", "NWD2", "ST6GALNAC3",
    "RELN", "HCN1", "RYR2", "KCNIP4", "CACNA2D3", "CADPS2",
    "GRIK2", "GRM8", "GRM5"
  )
)

############################################################
# 3) HELPER FUNCTIONS
############################################################

wrap_text <- function(x, width = 43L) {
  vapply(
    x,
    function(txt) {
      paste(strwrap(txt, width = width), collapse = "\n")
    },
    character(1)
  )
}

save_ggplot_publication <- function(
    plot,
    file_stem,
    width,
    height,
    tiff_dpi = main_tiff_dpi
) {
  pdf_file <- paste0(file_stem, ".pdf")
  tiff_file <- paste0(file_stem, ".tiff")
  
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      family = font_family,
      bg = "white",
      limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = "pdf",
      useDingbats = FALSE,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "tiff",
    dpi = tiff_dpi,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

save_complex_heatmap <- function(
    heatmap_object,
    file_stem,
    width,
    height,
    tiff_dpi = main_tiff_dpi
) {
  pdf_file <- paste0(file_stem, ".pdf")
  tiff_file <- paste0(file_stem, ".tiff")
  
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = font_family,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      file = pdf_file,
      width = width,
      height = height,
      family = font_family,
      useDingbats = FALSE,
      bg = "white"
    )
  }
  
  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(3, 3, 3, 3), "mm")
  )
  grDevices::dev.off()
  
  grDevices::tiff(
    filename = tiff_file,
    width = width,
    height = height,
    units = "in",
    res = tiff_dpi,
    compression = "lzw",
    bg = "white"
  )
  
  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(3, 3, 3, 3), "mm")
  )
  grDevices::dev.off()
}

broad_class_from_cell_type <- function(cell_type) {
  if (grepl("GABAergic", cell_type, ignore.case = TRUE)) {
    return("InN")
  }
  
  if (grepl("OPC", cell_type, ignore.case = TRUE)) {
    return("OPC")
  }
  
  if (grepl("astro", cell_type, ignore.case = TRUE)) {
    return("Astro")
  }
  
  if (grepl(
    "Perivascular|mesenchymal",
    cell_type,
    ignore.case = TRUE
  )) {
    return("Perivascular")
  }
  
  "ExN"
}

marker_key_from_cell_type <- function(cell_type) {
  aliases <- c(
    "L6 IT ExN" = "L6_IT_ExN",
    "L6 IT excitatory neurons" = "L6_IT_ExN",
    
    "Migrating cortical ExN" = "Migrating_Cortical_ExN",
    "Migrating cortical excitatory neurons" = "Migrating_Cortical_ExN",
    
    "Perivascular mesenchymal" = "Perivascular_Mesenchymal",
    "Perivascular mesenchymal cells" = "Perivascular_Mesenchymal",
    
    "Mid-layer IT ExN" = "Mid_layer_IT_ExN",
    "Mid-layer IT-like excitatory neurons" = "Mid_layer_IT_ExN",
    
    "Developing GABAergic neurons" = "Developing_GABAergic_neurons",
    
    "Immature GABAergic neurons" = "Immature_GABAergic_neurons",
    "Immature local GABAergic neurons" = "Immature_GABAergic_neurons",
    
    "Early OPC" = "Early_OPC",
    "Early OPCs" = "Early_OPC",
    
    "Immature upper-layer ExN" = "Immature_upper_layer_ExN",
    "Immature upper-layer excitatory neurons" = "Immature_upper_layer_ExN",
    
    "L6 CT-like ExN" = "L6_CT_like_ExN",
    "L6 corticothalamic-like excitatory neurons" = "L6_CT_like_ExN",
    
    "Developmental astrocytes" = "Developmental_astrocytes",
    
    "L5 deep-layer ExN" = "L5_deep_layer_ExN",
    "L5 deep-layer excitatory neurons" = "L5_deep_layer_ExN"
  )
  
  if (cell_type == "GABAergic neurons") {
    return("Merged_GABAergic_neurons")
  }
  
  unname(aliases[cell_type])
}

fallback_marker_label <- function(cell_type) {
  labels <- c(
    "L6 IT ExN" = "STK32B+/CDH13+ L6 IT ExN",
    "L6 IT excitatory neurons" = "STK32B+/CDH13+ L6 IT ExN",
    
    "Migrating cortical ExN" =
      "DPY19L1+/SEMA3C+ migrating cortical ExN",
    "Migrating cortical excitatory neurons" =
      "DPY19L1+/SEMA3C+ migrating cortical ExN",
    
    "Perivascular mesenchymal" =
      "PDGFRB+/PDE3A+ perivascular mesenchymal cells",
    "Perivascular mesenchymal cells" =
      "PDGFRB+/PDE3A+ perivascular mesenchymal cells",
    
    "Mid-layer IT ExN" = "RORB+ mid-layer IT-like ExN",
    "Mid-layer IT-like excitatory neurons" =
      "RORB+ mid-layer IT-like ExN",
    
    "GABAergic neurons" = "GAD1+/DLX6-AS1+ GABAergic neurons",
    "Developing GABAergic neurons" =
      "DLX6-AS1+/ST18+ developing GABAergic neurons",
    "Immature GABAergic neurons" =
      "PAX6+/ARX+ immature GABAergic neurons",
    "Immature local GABAergic neurons" =
      "PAX6+/ARX+ immature GABAergic neurons",
    
    "Early OPC" = "OLIG1+ early OPC",
    "Early OPCs" = "OLIG1+ early OPC",
    
    "Immature upper-layer ExN" =
      "UNC5D+ immature upper-layer ExN",
    "Immature upper-layer excitatory neurons" =
      "UNC5D+ immature upper-layer ExN",
    
    "L6 CT-like ExN" =
      "FOXP2+/TLE4+ L6 corticothalamic-like ExN",
    "L6 corticothalamic-like excitatory neurons" =
      "FOXP2+/TLE4+ L6 corticothalamic-like ExN",
    
    "Developmental astrocytes" =
      "AQP4+/SLC4A4+ developmental astrocytes",
    
    "L5 deep-layer ExN" = "RELN+ L5 deep-layer ExN",
    "L5 deep-layer excitatory neurons" =
      "RELN+ L5 deep-layer ExN"
  )
  
  value <- unname(labels[cell_type])
  
  if (is.na(value)) {
    value <- cell_type
  }
  
  value
}

############################################################
# 4) LOAD AND VALIDATE THE OBJECT
############################################################

message("Loading annotated object: ", qs_file)
obj <- qs::qread(qs_file)

if (!inherits(obj, "Seurat")) {
  stop("The selected QS file does not contain a Seurat object.")
}

if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop(
    "The selected object does not contain the metadata column ",
    "'cell_type'. Load one of the cell-type-annotated QS files."
  )
}

if (!"pca" %in% names(obj@reductions)) {
  stop(
    "The selected object does not contain a stored PCA reduction. ",
    "The t-SNE in this script is calculated from the CCA-derived PCA."
  )
}

if (!"RNA" %in% names(obj@assays)) {
  stop("The selected object does not contain an RNA assay.")
}

cell_type_values <- as.character(obj$cell_type)

if (anyNA(cell_type_values) || any(cell_type_values == "")) {
  stop("The cell_type metadata contains missing or empty labels.")
}

if (is.factor(obj$cell_type)) {
  cell_type_levels <- levels(droplevels(obj$cell_type))
} else {
  cell_type_levels <- unique(cell_type_values)
}

obj$cell_type <- factor(
  cell_type_values,
  levels = cell_type_levels
)

SeuratObject::Idents(obj) <- "cell_type"

message(
  "Cell types in the selected object: ",
  paste(cell_type_levels, collapse = ", ")
)

############################################################
# 4b) CLEAN BIOLOGICAL ORDER FOR THE HEATMAP ONLY
#
# This order keeps all excitatory-neuron populations together,
# followed by GABAergic populations, then glial/mesenchymal types.
# It does not change the Seurat identities, t-SNE coordinates, or
# cell-type assignments; it only controls heatmap rows and columns.
############################################################

preferred_heatmap_order <- c(
  # Excitatory neurons (ExN)
  "L6 IT ExN",
  "L6 IT excitatory neurons",
  "Migrating cortical ExN",
  "Migrating cortical excitatory neurons",
  "Mid-layer IT ExN",
  "Mid-layer IT-like excitatory neurons",
  "Immature upper-layer ExN",
  "Immature upper-layer excitatory neurons",
  "L6 CT-like ExN",
  "L6 corticothalamic-like excitatory neurons",
  "L5 deep-layer ExN",
  "L5 deep-layer excitatory neurons",
  
  # GABAergic neurons (InN)
  "GABAergic neurons",
  "Developing GABAergic neurons",
  "Immature GABAergic neurons",
  "Immature local GABAergic neurons",
  
  # Other populations
  "Early OPC",
  "Early OPCs",
  "Developmental astrocytes",
  "Perivascular mesenchymal",
  "Perivascular mesenchymal cells"
)

heatmap_cell_type_order <- preferred_heatmap_order[
  preferred_heatmap_order %in% cell_type_levels
]

# Preserve any unexpected labels rather than dropping them.
unlisted_heatmap_types <- setdiff(
  cell_type_levels,
  heatmap_cell_type_order
)

if (length(unlisted_heatmap_types) > 0L) {
  warning(
    "The following cell types were not in preferred_heatmap_order and ",
    "were appended at the end: ",
    paste(unlisted_heatmap_types, collapse = ", ")
  )
  
  heatmap_cell_type_order <- c(
    heatmap_cell_type_order,
    unlisted_heatmap_types
  )
}

if (!setequal(heatmap_cell_type_order, cell_type_levels)) {
  stop("The heatmap cell-type order does not match the observed cell types.")
}

message(
  "Heatmap order: ",
  paste(heatmap_cell_type_order, collapse = " -> ")
)

############################################################
# 5) OUTPUT DIRECTORY
############################################################

object_stub <- tools::file_path_sans_ext(
  basename(qs_file)
)

output_dir <- file.path(
  dirname(qs_file),
  paste0(object_stub, "_tsne_legend_heatmap")
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Output directory: ", output_dir)

############################################################
# 6) CELL-TYPE COLORS AND FIGURE KEY
############################################################

color_dictionary <- c(
  "L6 IT ExN" = "#1B9E77",
  "L6 IT excitatory neurons" = "#1B9E77",
  
  "Migrating cortical ExN" = "#66C2A5",
  "Migrating cortical excitatory neurons" = "#66C2A5",
  
  "Perivascular mesenchymal" = "#A676B7",
  "Perivascular mesenchymal cells" = "#A676B7",
  
  "Mid-layer IT ExN" = "#66A61E",
  "Mid-layer IT-like excitatory neurons" = "#66A61E",
  
  "GABAergic neurons" = "#D95F02",
  "Developing GABAergic neurons" = "#E69F00",
  "Immature GABAergic neurons" = "#C44E52",
  "Immature local GABAergic neurons" = "#C44E52",
  
  "Early OPC" = "#00A6D6",
  "Early OPCs" = "#00A6D6",
  
  "Immature upper-layer ExN" = "#4C78A8",
  "Immature upper-layer excitatory neurons" = "#4C78A8",
  
  "L6 CT-like ExN" = "#2F5597",
  "L6 corticothalamic-like excitatory neurons" = "#2F5597",
  
  "Developmental astrocytes" = "#E377C2",
  
  "L5 deep-layer ExN" = "#7A5195",
  "L5 deep-layer excitatory neurons" = "#7A5195"
)

cell_type_colors <- unname(
  color_dictionary[cell_type_levels]
)

names(cell_type_colors) <- cell_type_levels

missing_color <- is.na(cell_type_colors)

if (any(missing_color)) {
  fallback_colors <- grDevices::hcl.colors(
    sum(missing_color),
    palette = "Dark 3"
  )
  
  cell_type_colors[missing_color] <- fallback_colors
  
  warning(
    "Fallback colors were assigned to: ",
    paste(names(cell_type_colors)[missing_color], collapse = ", ")
  )
}

# These are display keys for the figure, not Seurat cluster numbers.
figure_key <- data.frame(
  plot_id = seq_along(cell_type_levels),
  cell_type = cell_type_levels,
  broad_class = vapply(
    cell_type_levels,
    broad_class_from_cell_type,
    character(1)
  ),
  color = unname(cell_type_colors[cell_type_levels]),
  stringsAsFactors = FALSE
)

# Prefer the marker-rich labels stored in the annotated object.
if ("cell_type_markers" %in% colnames(obj@meta.data)) {
  marker_label_lookup <- tapply(
    as.character(obj$cell_type_markers),
    as.character(obj$cell_type),
    function(x) unique(x)[1]
  )
  
  figure_key$marker_label <- unname(
    marker_label_lookup[figure_key$cell_type]
  )
} else {
  figure_key$marker_label <- NA_character_
}

missing_marker_label <- is.na(figure_key$marker_label) |
  figure_key$marker_label == ""

if (any(missing_marker_label)) {
  figure_key$marker_label[missing_marker_label] <- vapply(
    figure_key$cell_type[missing_marker_label],
    fallback_marker_label,
    character(1)
  )
}

write.csv(
  figure_key,
  file.path(output_dir, "cell_type_figure_key.csv"),
  row.names = FALSE
)

############################################################
# 7) RUN t-SNE FROM THE SAVED CCA PCA REDUCTION
############################################################

pca_dimensions_available <- ncol(
  SeuratObject::Embeddings(obj, reduction = "pca")
)

pcs_use <- pcs_requested[
  pcs_requested <= pca_dimensions_available
]

if (length(pcs_use) < 2L) {
  stop("Fewer than two PCA dimensions are available for t-SNE.")
}

if (ncol(obj) < 16L) {
  stop("Too few cells are available for a stable t-SNE plot.")
}

perplexity_use <- min(
  30,
  floor((ncol(obj) - 1) / 3)
)

perplexity_use <- max(
  5,
  perplexity_use
)

reduction_name <- "tsne_celltype"

if (reduction_name %in% names(obj@reductions)) {
  obj[[reduction_name]] <- NULL
}

message(
  "Running t-SNE with PCA dimensions ",
  min(pcs_use),
  ":",
  max(pcs_use),
  " and perplexity ",
  perplexity_use,
  "."
)

obj <- Seurat::RunTSNE(
  object = obj,
  reduction = "pca",
  dims = pcs_use,
  seed.use = tsne_seed,
  tsne.method = "Rtsne",
  reduction.name = reduction_name,
  reduction.key = "tSNEct_",
  perplexity = perplexity_use,
  check_duplicates = FALSE,
  verbose = TRUE
)

SeuratObject::Idents(obj) <- "cell_type"

embedding <- as.data.frame(
  SeuratObject::Embeddings(
    obj,
    reduction = reduction_name
  )
)

colnames(embedding)[1:2] <- c(
  "tSNE_1",
  "tSNE_2"
)

embedding$cell <- rownames(embedding)
embedding$cell_type <- as.character(obj$cell_type)
embedding$plot_id <- figure_key$plot_id[
  match(
    embedding$cell_type,
    figure_key$cell_type
  )
]

centers <- stats::aggregate(
  cbind(tSNE_1, tSNE_2) ~ cell_type,
  data = embedding,
  FUN = stats::median
)

centers$plot_id <- figure_key$plot_id[
  match(
    centers$cell_type,
    figure_key$cell_type
  )
]

centers$cell_type <- factor(
  centers$cell_type,
  levels = cell_type_levels
)

write.csv(
  centers,
  file.path(output_dir, "tsne_cell_type_centers.csv"),
  row.names = FALSE
)

# Rasterize the point layer inside the PDF when scattermore is present.
use_raster <- requireNamespace(
  "scattermore",
  quietly = TRUE
)

p_tsne <- Seurat::DimPlot(
  object = obj,
  reduction = reduction_name,
  group.by = "cell_type",
  cols = unname(cell_type_colors[cell_type_levels]),
  pt.size = 0.22,
  shuffle = TRUE,
  seed = tsne_seed,
  label = FALSE,
  alpha = 0.88,
  raster = use_raster,
  raster.dpi = c(1200, 1200)
) +
  ggplot2::geom_point(
    data = centers,
    mapping = ggplot2::aes(
      x = tSNE_1,
      y = tSNE_2,
      fill = cell_type
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 5.6,
    stroke = 0.45,
    color = "black"
  ) +
  ggplot2::geom_text(
    data = centers,
    mapping = ggplot2::aes(
      x = tSNE_1,
      y = tSNE_2,
      label = plot_id
    ),
    inherit.aes = FALSE,
    family = font_family,
    fontface = "bold",
    size = 3.1,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = cell_type_colors,
    guide = "none"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    x = "t-SNE 1",
    y = "t-SNE 2"
  ) +
  ggplot2::theme_classic(
    base_size = 8,
    base_family = font_family
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(
      linewidth = 0.45,
      color = "black"
    ),
    axis.title = ggplot2::element_text(
      size = 8,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4,
      unit = "mm"
    )
  )

save_ggplot_publication(
  plot = p_tsne,
  file_stem = file.path(
    output_dir,
    "Figure_TSNE_cell_types"
  ),
  width = 7.2,
  height = 5.8,
  tiff_dpi = main_tiff_dpi
)

############################################################
# 8) CREATE A STANDALONE CELL-TYPE LEGEND
############################################################

broad_order <- c(
  "ExN",
  "InN",
  "OPC",
  "Astro",
  "Perivascular"
)

broad_order <- broad_order[
  broad_order %in% figure_key$broad_class
]

legend_df <- figure_key
legend_df$cell_type <- factor(
  legend_df$cell_type,
  levels = cell_type_levels
)
legend_df$broad_class <- factor(
  legend_df$broad_class,
  levels = broad_order
)

legend_df <- legend_df[
  order(
    legend_df$broad_class,
    legend_df$plot_id
  ),
  ,
  drop = FALSE
]

legend_df$item_y <- ave(
  seq_len(nrow(legend_df)),
  legend_df$broad_class,
  FUN = function(x) rev(seq_along(x))
)

legend_df$wrapped_label <- wrap_text(
  legend_df$marker_label,
  width = 43L
)

p_legend <- ggplot2::ggplot(
  legend_df,
  ggplot2::aes(
    x = 0,
    y = item_y
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(fill = cell_type),
    shape = 21,
    size = 5.5,
    stroke = 0.45,
    color = "black"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = plot_id),
    family = font_family,
    fontface = "bold",
    size = 3.0,
    color = "black"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0.55,
      label = wrapped_label
    ),
    hjust = 0,
    vjust = 0.5,
    family = font_family,
    size = 3.2,
    lineheight = 0.94,
    color = "black"
  ) +
  ggplot2::facet_wrap(
    ~ broad_class,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::scale_fill_manual(
    values = cell_type_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-0.35, 10.5),
    expand = c(0, 0)
  ) +
  ggplot2::theme_void(
    base_family = font_family,
    base_size = 9
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 10,
      hjust = 0
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F1F1F1",
      color = "black",
      linewidth = 0.45
    ),
    panel.spacing = grid::unit(5, "mm"),
    plot.margin = ggplot2::margin(
      t = 4,
      r = 5,
      b = 4,
      l = 5,
      unit = "mm"
    )
  )

legend_height <- max(
  4.2,
  0.55 * max(table(legend_df$broad_class)) + 1.8
)

save_ggplot_publication(
  plot = p_legend,
  file_stem = file.path(
    output_dir,
    "Figure_cell_type_legend"
  ),
  width = 7.6,
  height = legend_height,
  tiff_dpi = legend_tiff_dpi
)

############################################################
# 9) PREPARE RNA EXPRESSION FOR THE MARKER HEATMAP
############################################################

SeuratObject::DefaultAssay(obj) <- "RNA"

if (inherits(obj[["RNA"]], "Assay5")) {
  rna_layers <- SeuratObject::Layers(obj[["RNA"]])
  
  multiple_count_layers <- sum(grepl("^counts", rna_layers)) > 1L
  multiple_data_layers <- sum(grepl("^data", rna_layers)) > 1L
  
  if (multiple_count_layers || multiple_data_layers) {
    message("Joining RNA layers before calculating the heatmap.")
    obj <- SeuratObject::JoinLayers(
      object = obj,
      assay = "RNA"
    )
  }
  
  rna_layers <- SeuratObject::Layers(obj[["RNA"]])
  
  if (!"data" %in% rna_layers) {
    message("RNA data layer is absent; running NormalizeData().")
    obj <- Seurat::NormalizeData(
      object = obj,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
} else {
  message("Normalizing the legacy RNA assay for heatmap expression.")
  obj <- Seurat::NormalizeData(
    object = obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

SeuratObject::Idents(obj) <- "cell_type"

############################################################
# 10) MAP MARKER PANELS TO THE CELL TYPES IN THIS OBJECT
############################################################

marker_panels <- list()

for (cell_type in heatmap_cell_type_order) {
  marker_key <- marker_key_from_cell_type(cell_type)
  
  if (is.na(marker_key)) {
    stop(
      "No marker-list mapping was defined for cell type: ",
      cell_type
    )
  }
  
  if (marker_key == "Merged_GABAergic_neurons") {
    marker_panels[[cell_type]] <- unique(
      c(
        neuron_markers$Developing_GABAergic_neurons,
        neuron_markers$Immature_GABAergic_neurons
      )
    )
  } else {
    marker_panels[[cell_type]] <- neuron_markers[[marker_key]]
  }
}

all_candidate_genes <- unique(
  unlist(marker_panels, use.names = FALSE)
)

rna_features <- rownames(obj[["RNA"]])

marker_availability <- do.call(
  rbind,
  lapply(
    names(marker_panels),
    function(cell_type) {
      data.frame(
        cell_type = cell_type,
        gene = marker_panels[[cell_type]],
        present_in_RNA = marker_panels[[cell_type]] %in% rna_features,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  marker_availability,
  file.path(output_dir, "marker_gene_availability.csv"),
  row.names = FALSE
)

available_genes <- all_candidate_genes[
  all_candidate_genes %in% rna_features
]

missing_genes <- setdiff(
  all_candidate_genes,
  available_genes
)

if (length(missing_genes) > 0L) {
  warning(
    length(missing_genes),
    " requested marker genes were not present in the RNA assay. ",
    "See marker_gene_availability.csv."
  )
}

if (length(available_genes) < 2L) {
  stop("Fewer than two requested marker genes were found in the RNA assay.")
}

rna_data <- SeuratObject::GetAssayData(
  object = obj,
  assay = "RNA",
  layer = "data"
)

rna_data <- rna_data[
  available_genes,
  ,
  drop = FALSE
]

############################################################
# 11) CALCULATE AVERAGE LOG-NORMALIZED EXPRESSION
############################################################

average_expression <- vapply(
  heatmap_cell_type_order,
  function(cell_type) {
    cells <- colnames(obj)[
      as.character(obj$cell_type) == cell_type
    ]
    
    if (length(cells) == 0L) {
      stop("No cells were found for cell type: ", cell_type)
    }
    
    Matrix::rowMeans(
      rna_data[
        ,
        cells,
        drop = FALSE
      ]
    )
  },
  numeric(length(available_genes))
)

rownames(average_expression) <- available_genes
colnames(average_expression) <- heatmap_cell_type_order

write.csv(
  average_expression,
  file.path(
    output_dir,
    "average_log_normalized_expression_all_requested_markers.csv"
  ),
  row.names = TRUE
)

############################################################
# 12) RANK USER-SUPPLIED MARKERS BY CELL-TYPE SPECIFICITY
############################################################

candidate_table <- do.call(
  rbind,
  lapply(
    names(marker_panels),
    function(cell_type) {
      genes <- marker_panels[[cell_type]]
      genes <- genes[genes %in% rownames(average_expression)]
      
      data.frame(
        cell_type = cell_type,
        gene = genes,
        supplied_order = seq_along(genes),
        stringsAsFactors = FALSE
      )
    }
  )
)

candidate_table <- unique(candidate_table)

candidate_table$target_average <- mapply(
  function(gene, cell_type) {
    average_expression[gene, cell_type]
  },
  candidate_table$gene,
  candidate_table$cell_type
)

candidate_table$maximum_other_average <- mapply(
  function(gene, cell_type) {
    other_types <- setdiff(heatmap_cell_type_order, cell_type)
    
    if (length(other_types) == 0L) {
      return(0)
    }
    
    max(
      average_expression[gene, other_types],
      na.rm = TRUE
    )
  },
  candidate_table$gene,
  candidate_table$cell_type
)

candidate_table$specificity_score <- (
  candidate_table$target_average -
    candidate_table$maximum_other_average
)

# A duplicated gene is shown only once in the heatmap. Assign it to
# the marker panel in which it is most specific in this dataset.
candidate_table <- candidate_table[
  order(
    candidate_table$gene,
    -candidate_table$specificity_score,
    candidate_table$supplied_order
  ),
  ,
  drop = FALSE
]

candidate_unique <- candidate_table[
  !duplicated(candidate_table$gene),
  ,
  drop = FALSE
]

# Restore cell-type order before selecting top genes.
candidate_unique$cell_type <- factor(
  candidate_unique$cell_type,
  levels = heatmap_cell_type_order
)

selected_compact <- do.call(
  rbind,
  lapply(
    heatmap_cell_type_order,
    function(cell_type) {
      x <- candidate_unique[
        candidate_unique$cell_type == cell_type,
        ,
        drop = FALSE
      ]
      
      x <- x[
        order(
          -x$specificity_score,
          x$supplied_order
        ),
        ,
        drop = FALSE
      ]
      
      head(
        x,
        n = min(
          top_markers_per_cell_type,
          nrow(x)
        )
      )
    }
  )
)

selected_compact$cell_type <- factor(
  selected_compact$cell_type,
  levels = heatmap_cell_type_order
)

selected_compact <- selected_compact[
  order(
    selected_compact$cell_type,
    -selected_compact$specificity_score
  ),
  ,
  drop = FALSE
]

write.csv(
  selected_compact,
  file.path(output_dir, "compact_heatmap_selected_markers.csv"),
  row.names = FALSE
)

write.csv(
  candidate_unique,
  file.path(output_dir, "all_unique_marker_candidates_ranked.csv"),
  row.names = FALSE
)

############################################################
# 13) HEATMAP-BUILDING FUNCTION
############################################################

build_marker_heatmap <- function(
    marker_table,
    heatmap_name = "Expression"
) {
  marker_table <- marker_table[
    !duplicated(marker_table$gene),
    ,
    drop = FALSE
  ]
  
  selected_genes <- marker_table$gene
  marker_groups <- as.character(marker_table$cell_type)
  
  average_selected <- average_expression[
    selected_genes,
    heatmap_cell_type_order,
    drop = FALSE
  ]
  
  gene_means <- rowMeans(
    average_selected,
    na.rm = TRUE
  )
  
  gene_sd <- apply(
    average_selected,
    MARGIN = 1,
    FUN = stats::sd,
    na.rm = TRUE
  )
  
  gene_sd[
    is.na(gene_sd) |
      !is.finite(gene_sd) |
      gene_sd == 0
  ] <- 1
  
  z_matrix <- sweep(
    average_selected,
    MARGIN = 1,
    STATS = gene_means,
    FUN = "-"
  )
  
  z_matrix <- sweep(
    z_matrix,
    MARGIN = 1,
    STATS = gene_sd,
    FUN = "/"
  )
  
  z_matrix[!is.finite(z_matrix)] <- 0
  z_matrix[z_matrix > 2.5] <- 2.5
  z_matrix[z_matrix < -2.5] <- -2.5
  
  # Rows are cell types and columns are genes, matching the reference style.
  heatmap_matrix <- t(z_matrix)
  
  rownames(heatmap_matrix) <- heatmap_cell_type_order
  colnames(heatmap_matrix) <- selected_genes
  
  marker_groups <- factor(
    marker_groups,
    levels = heatmap_cell_type_order
  )
  
  row_broad_class <- factor(
    vapply(
      heatmap_cell_type_order,
      broad_class_from_cell_type,
      character(1)
    ),
    levels = broad_order
  )
  
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Marker_panel = marker_groups,
    col = list(
      Marker_panel = cell_type_colors
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    annotation_height = grid::unit(3.0, "mm"),
    simple_anno_size = grid::unit(3.0, "mm")
  )
  
  left_annotation <- ComplexHeatmap::rowAnnotation(
    Cell_type = factor(
      heatmap_cell_type_order,
      levels = heatmap_cell_type_order
    ),
    col = list(
      Cell_type = cell_type_colors
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = grid::unit(3.2, "mm"),
    simple_anno_size = grid::unit(3.2, "mm")
  )
  
  heatmap_colors <- circlize::colorRamp2(
    c(-2.5, 0, 2.5),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )
  
  ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = heatmap_name,
    col = heatmap_colors,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    row_split = row_broad_class,
    row_gap = grid::unit(1.2, "mm"),
    row_title = NULL,
    
    column_split = marker_groups,
    column_gap = grid::unit(1.0, "mm"),
    column_title = NULL,
    
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(
      fontsize = 8.5,
      fontface = "bold",
      fontfamily = font_family
    ),
    row_names_max_width = grid::unit(48, "mm"),
    
    show_column_names = TRUE,
    column_names_rot = 90,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(
      fontsize = 6.5,
      fontface = "italic",
      fontfamily = font_family,
      col = unname(
        cell_type_colors[as.character(marker_groups)]
      )
    ),
    
    rect_gp = grid::gpar(
      col = "#D9D9D9",
      lwd = 0.35
    ),
    border = TRUE,
    
    heatmap_legend_param = list(
      title = "Z-scored\nexpression",
      title_gp = grid::gpar(
        fontsize = 8,
        fontface = "bold",
        fontfamily = font_family
      ),
      labels_gp = grid::gpar(
        fontsize = 7,
        fontfamily = font_family
      ),
      at = c(-2, 0, 2),
      legend_height = grid::unit(28, "mm"),
      grid_width = grid::unit(3.5, "mm")
    ),
    
    use_raster = FALSE
  )
}

############################################################
# 14) COMPACT PUBLICATION HEATMAP
############################################################

if (nrow(selected_compact) == 0L) {
  stop("No markers were available for the compact heatmap.")
}

compact_heatmap <- build_marker_heatmap(
  marker_table = selected_compact
)

compact_gene_n <- nrow(selected_compact)
compact_width <- max(
  8.0,
  3.5 + 0.115 * compact_gene_n
)
compact_height <- max(
  5.1,
  2.2 + 0.34 * length(heatmap_cell_type_order)
)

save_complex_heatmap(
  heatmap_object = compact_heatmap,
  file_stem = file.path(
    output_dir,
    "Figure_marker_heatmap_compact"
  ),
  width = compact_width,
  height = compact_height,
  tiff_dpi = main_tiff_dpi
)

# Save the exact matrix underlying the compact figure.
compact_genes <- selected_compact$gene
compact_average <- average_expression[
  compact_genes,
  heatmap_cell_type_order,
  drop = FALSE
]

compact_z <- t(
  scale(
    t(compact_average)
  )
)
compact_z[!is.finite(compact_z)] <- 0

write.csv(
  t(compact_z),
  file.path(
    output_dir,
    "compact_heatmap_z_scored_expression_matrix.csv"
  ),
  row.names = TRUE
)

############################################################
# 15) OPTIONAL ALL-MARKER HEATMAP
############################################################

if (isTRUE(make_all_marker_heatmap)) {
  all_marker_table <- candidate_unique
  
  all_marker_table <- all_marker_table[
    order(
      all_marker_table$cell_type,
      -all_marker_table$specificity_score,
      all_marker_table$supplied_order
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(all_marker_table) > 0L) {
    all_marker_heatmap <- build_marker_heatmap(
      marker_table = all_marker_table
    )
    
    all_gene_n <- nrow(all_marker_table)
    all_width <- max(
      10.0,
      3.8 + 0.105 * all_gene_n
    )
    all_height <- compact_height
    
    save_complex_heatmap(
      heatmap_object = all_marker_heatmap,
      file_stem = file.path(
        output_dir,
        "Supplement_marker_heatmap_all_available_markers"
      ),
      width = all_width,
      height = all_height,
      tiff_dpi = main_tiff_dpi
    )
  }
}

############################################################
# 16) OPTIONAL: SAVE OBJECT WITH THE NEW t-SNE REDUCTION
############################################################

if (isTRUE(save_object_with_tsne)) {
  output_qs <- file.path(
    output_dir,
    paste0(object_stub, "_with_celltype_tsne.qs")
  )
  
  qs::qsave(
    obj,
    output_qs
  )
  
  message("Saved object containing t-SNE: ", output_qs)
}

############################################################
# 17) COMPLETION MESSAGE
############################################################

message("\nFigure generation complete.")
message("Primary t-SNE: Figure_TSNE_cell_types.pdf/.tiff")
message("Standalone legend: Figure_cell_type_legend.pdf/.tiff")
message("Primary heatmap: Figure_marker_heatmap_compact.pdf/.tiff")
message("Outputs saved to: ", output_dir)
##################################################################

############################################################
# PUBLICATION-QUALITY NUMBERED t-SNE
# Same cell-type colors and biological order as the heatmap
#
# Input:
#   A cell-type-annotated CCA Seurat object saved as .qs.
#   The object should contain:
#     - cell_type metadata (preferred), or celltype_final
#     - a PCA reduction
#
# Works with either annotation version:
#   1) clusters 4 and 5 merged as "GABAergic neurons"
#   2) clusters 4 and 5 kept separate
#
# Main outputs:
#   Figure_tSNE_cell_types_numbered.pdf
#   Figure_tSNE_cell_types_numbered.tiff
#   Figure_tSNE_cell_types_numbered.png
#   Figure_tSNE_number_key.pdf
#   Figure_tSNE_number_key.tiff
#   tSNE_cell_type_figure_key.csv
#   tSNE_cell_coordinates.csv
#   tSNE_label_positions.csv
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

############################################################
# 1) USER SETTINGS
############################################################

set.seed(1234)

# Select one annotated object when the file browser opens, for example:
#   integrated_CCA_removed_8_9_10_12_14_merged_4_5_celltype_annotated.qs
# or:
#   integrated_CCA_removed_clusters_8_9_10_12_14_celltype_annotated.qs
qs_file <- file.choose()

# Match the PCA input in the example code supplied by the user.
# The script automatically uses only dimensions available in the object.
pcs_requested <- 1:30

tsne_seed <- 1234
reduction_name <- "tsne_celltype"

# FALSE reuses reduction_name when it already exists in the selected object.
# TRUE deletes it and calculates a new t-SNE.
force_recompute_tsne <- FALSE

# NULL lets the script use min(30, floor((n_cells - 1) / 3)).
# Set a number, such as 30, to specify it directly.
perplexity_user <- NULL

# t-SNE orientation is mathematically arbitrary. These settings rotate or
# mirror the plot without changing distances or biological structure.
# They are useful only when you want a preferred visual orientation.
rotation_degrees <- 0
swap_axes <- FALSE
flip_x <- FALSE
flip_y <- FALSE

# Point and label appearance. The white point layer creates the clean
# separated-dot appearance in the supplied reference figure.
cell_point_size <- 0.48
cell_halo_size <- 0.92
cell_alpha <- 0.96
number_circle_size <- 7.0
number_text_size <- 4.15
number_circle_stroke <- 0.85

# Figure dimensions chosen to resemble the wide reference image.
figure_width <- 7.25
figure_height <- 4.35
figure_tiff_dpi <- 600
preview_png_dpi <- 300

font_family <- "sans"

# Save a new large QS object containing the t-SNE reduction.
# The coordinate CSV is always saved, so FALSE is usually sufficient.
save_object_with_tsne <- FALSE

# Also save a separate number-to-cell-type key.
make_number_key <- TRUE

############################################################
# 2) LOAD AND VALIDATE THE ANNOTATED OBJECT
############################################################

if (!file.exists(qs_file)) {
  stop("Selected QS file does not exist: ", qs_file)
}

obj <- qs::qread(qs_file)

if (!inherits(obj, "Seurat")) {
  stop("The selected QS file does not contain a Seurat object.")
}

metadata_columns <- colnames(obj@meta.data)

cell_type_column <- if ("cell_type" %in% metadata_columns) {
  "cell_type"
} else if ("celltype_final" %in% metadata_columns) {
  "celltype_final"
} else {
  stop(
    "The object contains neither 'cell_type' nor 'celltype_final' metadata. ",
    "Available metadata columns: ",
    paste(metadata_columns, collapse = ", ")
  )
}

cell_type_values <- as.character(obj@meta.data[[cell_type_column]])

if (anyNA(cell_type_values) || any(cell_type_values == "")) {
  stop("The selected cell-type metadata contains missing or empty labels.")
}

if (!("pca" %in% names(obj@reductions))) {
  stop("The selected object does not contain a PCA reduction.")
}

message("Loaded: ", normalizePath(qs_file, winslash = "/", mustWork = FALSE))
message("Cells: ", ncol(obj))
message("Cell-type column: ", cell_type_column)


############################################################
# 3) BIOLOGICAL ORDER USED BY BOTH HEATMAP AND t-SNE KEY
#
# All excitatory-neuron populations are together first,
# followed by GABAergic populations, then other populations.
############################################################

preferred_cell_type_order <- c(
  # Excitatory neurons
  "L6 IT ExN",
  "L6 IT excitatory neurons",
  "Migrating cortical ExN",
  "Migrating cortical excitatory neurons",
  "Mid-layer IT ExN",
  "Mid-layer IT-like excitatory neurons",
  "Immature upper-layer ExN",
  "Immature upper-layer excitatory neurons",
  "L6 CT-like ExN",
  "L6 corticothalamic-like excitatory neurons",
  "L5 deep-layer ExN",
  "L5 deep-layer excitatory neurons",
  
  # GABAergic neurons
  "GABAergic neurons",
  "Developing GABAergic neurons",
  "Immature GABAergic neurons",
  "Immature local GABAergic neurons",
  
  # Other populations
  "Early OPC",
  "Early OPCs",
  "Developmental astrocytes",
  "Perivascular mesenchymal",
  "Perivascular mesenchymal cells"
)

observed_cell_types <- unique(cell_type_values)

cell_type_order <- preferred_cell_type_order[
  preferred_cell_type_order %in% observed_cell_types
]

unexpected_cell_types <- setdiff(
  observed_cell_types,
  cell_type_order
)

if (length(unexpected_cell_types) > 0L) {
  warning(
    "These cell types were not in the preferred heatmap order and were ",
    "appended at the end: ",
    paste(unexpected_cell_types, collapse = ", "),
    call. = FALSE
  )
  
  cell_type_order <- c(
    cell_type_order,
    unexpected_cell_types
  )
}

if (!setequal(cell_type_order, observed_cell_types)) {
  stop("The ordered cell-type list does not match the observed cell types.")
}

obj$cell_type_plot <- factor(
  cell_type_values,
  levels = cell_type_order
)

SeuratObject::Idents(obj) <- "cell_type_plot"

message(
  "t-SNE number order: ",
  paste(
    paste0(seq_along(cell_type_order), " = ", cell_type_order),
    collapse = "; "
  )
)

############################################################
# 4) EXACT COLOR DICTIONARY USED FOR THE ORDERED HEATMAP
############################################################

color_dictionary <- c(
  "L6 IT ExN" = "#1B9E77",
  "L6 IT excitatory neurons" = "#1B9E77",
  
  "Migrating cortical ExN" = "#66C2A5",
  "Migrating cortical excitatory neurons" = "#66C2A5",
  
  "Mid-layer IT ExN" = "#66A61E",
  "Mid-layer IT-like excitatory neurons" = "#66A61E",
  
  "Immature upper-layer ExN" = "#4C78A8",
  "Immature upper-layer excitatory neurons" = "#4C78A8",
  
  "L6 CT-like ExN" = "#2F5597",
  "L6 corticothalamic-like excitatory neurons" = "#2F5597",
  
  "L5 deep-layer ExN" = "#7A5195",
  "L5 deep-layer excitatory neurons" = "#7A5195",
  
  "GABAergic neurons" = "#D95F02",
  "Developing GABAergic neurons" = "#E69F00",
  "Immature GABAergic neurons" = "#C44E52",
  "Immature local GABAergic neurons" = "#C44E52",
  
  "Early OPC" = "#00A6D6",
  "Early OPCs" = "#00A6D6",
  
  "Developmental astrocytes" = "#E377C2",
  
  "Perivascular mesenchymal" = "#A676B7",
  "Perivascular mesenchymal cells" = "#A676B7"
)

cell_type_colors <- color_dictionary[cell_type_order]

missing_colors <- is.na(cell_type_colors)

if (any(missing_colors)) {
  fallback_colors <- grDevices::hcl.colors(
    sum(missing_colors),
    palette = "Dark 3"
  )
  
  cell_type_colors[missing_colors] <- fallback_colors
  
  warning(
    "Fallback colors were assigned to: ",
    paste(names(cell_type_colors)[missing_colors], collapse = ", "),
    call. = FALSE
  )
}

############################################################
# 5) OUTPUT DIRECTORY AND FIGURE KEY
############################################################

object_stub <- tools::file_path_sans_ext(
  basename(qs_file)
)

output_dir <- file.path(
  dirname(qs_file),
  paste0(object_stub, "_numbered_tsne_publication")
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

marker_labels <- rep(NA_character_, length(cell_type_order))
names(marker_labels) <- cell_type_order

if ("cell_type_markers" %in% metadata_columns) {
  marker_lookup <- tapply(
    as.character(obj$cell_type_markers),
    as.character(obj$cell_type_plot),
    function(x) unique(x)[1]
  )
  
  marker_labels[names(marker_lookup)] <- unname(marker_lookup)
}

marker_labels[
  is.na(marker_labels) | marker_labels == ""
] <- names(marker_labels)[
  is.na(marker_labels) | marker_labels == ""
]

figure_key <- data.frame(
  plot_number = seq_along(cell_type_order),
  cell_type = cell_type_order,
  marker_label = unname(marker_labels[cell_type_order]),
  color = unname(cell_type_colors[cell_type_order]),
  stringsAsFactors = FALSE
)

write.csv(
  figure_key,
  file.path(output_dir, "tSNE_cell_type_figure_key.csv"),
  row.names = FALSE
)

############################################################
# 6) RUN OR REUSE t-SNE
############################################################

pca_embeddings <- SeuratObject::Embeddings(
  obj,
  reduction = "pca"
)

n_pcs_available <- ncol(pca_embeddings)

pcs_use <- pcs_requested[
  pcs_requested <= n_pcs_available
]

if (length(pcs_use) < 2L) {
  stop(
    "Fewer than two requested PCA dimensions are available. ",
    "PCA dimensions found: ",
    n_pcs_available
  )
}

if (!is.null(perplexity_user)) {
  perplexity_use <- as.numeric(perplexity_user)
} else {
  perplexity_use <- min(
    30,
    floor((ncol(obj) - 1) / 3)
  )
}

perplexity_use <- max(
  5,
  perplexity_use
)

if (perplexity_use >= (ncol(obj) - 1) / 3) {
  perplexity_use <- max(
    2,
    floor((ncol(obj) - 1) / 3) - 1
  )
}

if (force_recompute_tsne && reduction_name %in% names(obj@reductions)) {
  obj[[reduction_name]] <- NULL
}

if (!(reduction_name %in% names(obj@reductions))) {
  message(
    "Running t-SNE from PCA dimensions ",
    min(pcs_use),
    ":",
    max(pcs_use),
    " with perplexity ",
    perplexity_use,
    " and seed ",
    tsne_seed,
    "."
  )
  
  obj <- Seurat::RunTSNE(
    object = obj,
    reduction = "pca",
    dims = pcs_use,
    seed.use = tsne_seed,
    tsne.method = "Rtsne",
    reduction.name = reduction_name,
    reduction.key = "tSNEct_",
    perplexity = perplexity_use,
    check_duplicates = FALSE,
    verbose = TRUE
  )
} else {
  message("Reusing existing reduction: ", reduction_name)
}

############################################################
# 7) BUILD EMBEDDING TABLE AND OPTIONAL ORIENTATION
############################################################

embedding <- as.data.frame(
  SeuratObject::Embeddings(
    obj,
    reduction = reduction_name
  )
)

if (ncol(embedding) < 2L) {
  stop("The t-SNE reduction has fewer than two dimensions.")
}

colnames(embedding)[1:2] <- c(
  "tSNE_1",
  "tSNE_2"
)

embedding <- embedding[, c("tSNE_1", "tSNE_2"), drop = FALSE]
embedding$cell <- rownames(embedding)
embedding$cell_type <- factor(
  as.character(obj@meta.data[embedding$cell, "cell_type_plot"]),
  levels = cell_type_order
)

# Optional axis swap, rotation, and reflection. These preserve all pairwise
# relationships and are purely presentational.
if (isTRUE(swap_axes)) {
  old_x <- embedding$tSNE_1
  embedding$tSNE_1 <- embedding$tSNE_2
  embedding$tSNE_2 <- old_x
}

if (!isTRUE(all.equal(rotation_degrees, 0))) {
  theta <- rotation_degrees * pi / 180
  
  x_rot <- embedding$tSNE_1 * cos(theta) -
    embedding$tSNE_2 * sin(theta)
  
  y_rot <- embedding$tSNE_1 * sin(theta) +
    embedding$tSNE_2 * cos(theta)
  
  embedding$tSNE_1 <- x_rot
  embedding$tSNE_2 <- y_rot
}

if (isTRUE(flip_x)) {
  embedding$tSNE_1 <- -embedding$tSNE_1
}

if (isTRUE(flip_y)) {
  embedding$tSNE_2 <- -embedding$tSNE_2
}

# Draw large groups first and small groups last so small populations remain
# visible. ggplot draws later rows on top of earlier rows.
cell_type_n <- table(embedding$cell_type)

embedding$n_cell_type <- as.integer(
  cell_type_n[as.character(embedding$cell_type)]
)

embedding <- embedding %>%
  dplyr::arrange(
    dplyr::desc(.data$n_cell_type)
  )

write.csv(
  embedding[, c("cell", "tSNE_1", "tSNE_2", "cell_type")],
  file.path(output_dir, "tSNE_cell_coordinates.csv"),
  row.names = FALSE
)

############################################################
# 8) NUMBERED CIRCLE POSITIONS
#
# For each cell type, calculate the 2-D median and then choose the actual
# cell closest to that median. This keeps the numbered circle inside the
# population even when the population is curved or contains an empty center.
############################################################

median_positions <- embedding %>%
  dplyr::group_by(.data$cell_type) %>%
  dplyr::summarise(
    median_x = stats::median(.data$tSNE_1),
    median_y = stats::median(.data$tSNE_2),
    .groups = "drop"
  )

label_positions <- embedding %>%
  dplyr::left_join(
    median_positions,
    by = "cell_type"
  ) %>%
  dplyr::mutate(
    distance_to_median =
      (.data$tSNE_1 - .data$median_x)^2 +
      (.data$tSNE_2 - .data$median_y)^2
  ) %>%
  dplyr::group_by(.data$cell_type) %>%
  dplyr::slice_min(
    order_by = .data$distance_to_median,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    cell_type,
    tSNE_1,
    tSNE_2
  )

# Optional manual offsets. Add rows here only when a numbered circle needs
# a small adjustment after inspecting the first plot. Units are t-SNE units.
manual_label_offsets <- data.frame(
  cell_type = character(),
  dx = numeric(),
  dy = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(manual_label_offsets) > 0L) {
  label_positions <- label_positions %>%
    dplyr::left_join(
      manual_label_offsets,
      by = "cell_type"
    ) %>%
    dplyr::mutate(
      dx = dplyr::coalesce(.data$dx, 0),
      dy = dplyr::coalesce(.data$dy, 0),
      tSNE_1 = .data$tSNE_1 + .data$dx,
      tSNE_2 = .data$tSNE_2 + .data$dy
    ) %>%
    dplyr::select(
      -dplyr::all_of(c("dx", "dy"))
    )
}

label_positions$plot_number <- figure_key$plot_number[
  match(
    as.character(label_positions$cell_type),
    figure_key$cell_type
  )
]

label_positions$cell_type <- factor(
  as.character(label_positions$cell_type),
  levels = cell_type_order
)

write.csv(
  label_positions,
  file.path(output_dir, "tSNE_label_positions.csv"),
  row.names = FALSE
)

############################################################
# 9) CUSTOM BOTTOM-LEFT ARROW AXES
############################################################

x_range <- range(
  embedding$tSNE_1,
  na.rm = TRUE
)

y_range <- range(
  embedding$tSNE_2,
  na.rm = TRUE
)

x_span <- diff(x_range)
y_span <- diff(y_range)

if (x_span <= 0 || y_span <= 0) {
  stop("Invalid t-SNE coordinate range.")
}

# Add controlled white space for the custom arrow axes.
plot_xlim <- c(
  x_range[1] - 0.20 * x_span,
  x_range[2] + 0.035 * x_span
)

plot_ylim <- c(
  y_range[1] - 0.17 * y_span,
  y_range[2] + 0.035 * y_span
)

x0 <- x_range[1] - 0.145 * x_span
y0 <- y_range[1] - 0.105 * y_span

x1 <- x0 + 0.145 * x_span
y1 <- y0 + 0.165 * y_span

############################################################
# 10) BUILD THE t-SNE PLOT
############################################################

# ggrastr keeps the point cloud sharp while preventing a very large PDF.
use_ggrastr <- requireNamespace(
  "ggrastr",
  quietly = TRUE
)

p_tsne <- ggplot2::ggplot(
  embedding,
  ggplot2::aes(
    x = .data$tSNE_1,
    y = .data$tSNE_2
  )
)

if (use_ggrastr) {
  p_tsne <- p_tsne +
    ggrastr::geom_point_rast(
      color = "white",
      size = cell_halo_size,
      alpha = 0.98,
      raster.dpi = 1200
    ) +
    ggrastr::geom_point_rast(
      ggplot2::aes(
        color = .data$cell_type
      ),
      size = cell_point_size,
      alpha = cell_alpha,
      raster.dpi = 1200
    )
} else {
  p_tsne <- p_tsne +
    ggplot2::geom_point(
      color = "white",
      size = cell_halo_size,
      alpha = 0.98,
      stroke = 0
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = .data$cell_type
      ),
      size = cell_point_size,
      alpha = cell_alpha,
      stroke = 0
    )
}

p_tsne <- p_tsne +
  ggplot2::scale_color_manual(
    values = cell_type_colors,
    limits = cell_type_order,
    drop = FALSE,
    guide = "none"
  ) +
  
  # Filled numbered circles, matching the supplied reference style.
  ggplot2::geom_point(
    data = label_positions,
    mapping = ggplot2::aes(
      x = .data$tSNE_1,
      y = .data$tSNE_2,
      fill = .data$cell_type
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = number_circle_size,
    stroke = number_circle_stroke,
    color = "black",
    alpha = 0.98
  ) +
  ggplot2::geom_text(
    data = label_positions,
    mapping = ggplot2::aes(
      x = .data$tSNE_1,
      y = .data$tSNE_2,
      label = .data$plot_number
    ),
    inherit.aes = FALSE,
    family = font_family,
    fontface = "bold",
    size = number_text_size,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = cell_type_colors,
    limits = cell_type_order,
    drop = FALSE,
    guide = "none"
  ) +
  
  # Custom bottom-left arrow axes.
  ggplot2::annotate(
    geom = "segment",
    x = x0,
    xend = x1,
    y = y0,
    yend = y0,
    linewidth = 0.78,
    color = "black",
    lineend = "butt",
    arrow = grid::arrow(
      length = grid::unit(3.2, "mm"),
      type = "closed"
    )
  ) +
  ggplot2::annotate(
    geom = "segment",
    x = x0,
    xend = x0,
    y = y0,
    yend = y1,
    linewidth = 0.78,
    color = "black",
    lineend = "butt",
    arrow = grid::arrow(
      length = grid::unit(3.2, "mm"),
      type = "closed"
    )
  ) +
  ggplot2::annotate(
    geom = "text",
    x = (x0 + x1) / 2,
    y = y0 - 0.055 * y_span,
    label = "tSNE_1",
    family = font_family,
    size = 4.0,
    color = "black"
  ) +
  ggplot2::annotate(
    geom = "text",
    x = x0 - 0.050 * x_span,
    y = (y0 + y1) / 2,
    label = "tSNE_2",
    angle = 90,
    family = font_family,
    size = 4.0,
    color = "black"
  ) +
  ggplot2::coord_fixed(
    xlim = plot_xlim,
    ylim = plot_ylim,
    expand = FALSE,
    clip = "off"
  ) +
  ggplot2::theme_void(
    base_size = 14,
    base_family = font_family
  ) +
  ggplot2::theme(
    legend.position = "none",
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = ggplot2::margin(
      t = 2.5,
      r = 2.5,
      b = 2.5,
      l = 2.5,
      unit = "mm"
    )
  )

print(p_tsne)

############################################################
# 11) SAVE PUBLICATION-QUALITY t-SNE OUTPUTS
############################################################

pdf_file <- file.path(
  output_dir,
  "Figure_tSNE_cell_types_numbered.pdf"
)

tiff_file <- file.path(
  output_dir,
  "Figure_tSNE_cell_types_numbered.tiff"
)

png_file <- file.path(
  output_dir,
  "Figure_tSNE_cell_types_numbered.png"
)

if (capabilities("cairo")) {
  ggplot2::ggsave(
    filename = pdf_file,
    plot = p_tsne,
    device = grDevices::cairo_pdf,
    width = figure_width,
    height = figure_height,
    units = "in",
    bg = "white"
  )
} else {
  grDevices::pdf(
    file = pdf_file,
    width = figure_width,
    height = figure_height,
    family = font_family,
    useDingbats = FALSE
  )
  
  print(p_tsne)
  grDevices::dev.off()
}

ggplot2::ggsave(
  filename = tiff_file,
  plot = p_tsne,
  device = "tiff",
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = figure_tiff_dpi,
  compression = "lzw",
  bg = "white"
)

ggplot2::ggsave(
  filename = png_file,
  plot = p_tsne,
  device = "png",
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = preview_png_dpi,
  bg = "white"
)

############################################################
# 12) OPTIONAL MATCHING NUMBER KEY
############################################################

if (isTRUE(make_number_key)) {
  key_plot_data <- figure_key
  
  # Top item first in the vertical key.
  key_plot_data$cell_type <- factor(
    key_plot_data$cell_type,
    levels = rev(cell_type_order)
  )
  
  key_plot_data$display_label <- paste0(
    key_plot_data$plot_number,
    "  ",
    key_plot_data$marker_label
  )
  
  p_key <- ggplot2::ggplot(
    key_plot_data,
    ggplot2::aes(
      x = 1,
      y = .data$cell_type
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = .data$cell_type
      ),
      shape = 21,
      size = 5.6,
      stroke = 0.75,
      color = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = .data$plot_number
      ),
      family = font_family,
      fontface = "bold",
      size = 3.25,
      color = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 1.12,
        label = .data$marker_label
      ),
      hjust = 0,
      family = font_family,
      size = 3.45,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = cell_type_colors,
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.92, 2.55),
      clip = "off"
    ) +
    ggplot2::theme_void(
      base_family = font_family
    ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      plot.margin = ggplot2::margin(
        4,
        10,
        4,
        4,
        unit = "mm"
      )
    )
  
  key_height <- max(
    3.0,
    0.42 * nrow(key_plot_data) + 0.5
  )
  
  key_pdf <- file.path(
    output_dir,
    "Figure_tSNE_number_key.pdf"
  )
  
  key_tiff <- file.path(
    output_dir,
    "Figure_tSNE_number_key.tiff"
  )
  
  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = key_pdf,
      plot = p_key,
      device = grDevices::cairo_pdf,
      width = 5.8,
      height = key_height,
      units = "in",
      bg = "white"
    )
  } else {
    grDevices::pdf(
      file = key_pdf,
      width = 5.8,
      height = key_height,
      family = font_family,
      useDingbats = FALSE
    )
    
    print(p_key)
    grDevices::dev.off()
  }
  
  ggplot2::ggsave(
    filename = key_tiff,
    plot = p_key,
    device = "tiff",
    width = 5.8,
    height = key_height,
    units = "in",
    dpi = 1000,
    compression = "lzw",
    bg = "white"
  )
}

############################################################
# 13) OPTIONAL SAVE OF THE OBJECT WITH t-SNE
############################################################

if (isTRUE(save_object_with_tsne)) {
  object_out <- file.path(
    output_dir,
    paste0(object_stub, "_with_numbered_tsne.qs")
  )
  
  qs::qsave(
    obj,
    object_out
  )
  
  message("Object with t-SNE saved to: ", object_out)
}

message("Finished.")
message("Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
message("Main PDF: ", normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
message("Main TIFF: ", normalizePath(tiff_file, winslash = "/", mustWork = FALSE))
###################################################################

############################################################
# PUBLICATION-QUALITY CELL-TYPE LEGEND
#
# This legend matches the ordered heatmap and numbered t-SNE:
#   1-6  = excitatory neurons
#   7-8  = GABAergic neurons
#   9    = early OPC
#   10   = developmental astrocytes
#   11   = perivascular mesenchymal cells
#
# IMPORTANT: This 11-entry legend matches the object in which
# clusters 4 and 5 remain separate.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(grid)
  library(patchwork)
})

############################################################
# 1) OUTPUT SETTINGS
############################################################

output_dir <- file.path(
  getwd(),
  "cell_type_legend_publication"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

font_family <- if (.Platform$OS.type == "windows") {
  "Arial"
} else {
  "sans"
}

############################################################
# 2) EXACT ORDER, NAMES, MARKERS, AND COLORS
#
# The cell-type order and colors are identical to the ordered
# heatmap and numbered t-SNE scripts.
############################################################

legend_map <- tibble::tribble(
  ~broad_group, ~cell_type,                           ~plot_number, ~gene1,   ~gene2,      ~color,
  "ExN",        "L6 IT ExN",                                  1L, "STK32B", "MCTP1",     "#1B9E77",
  "ExN",        "Migrating cortical ExN",                     2L, "PLXNA2", "PTPN4",     "#66C2A5",
  "ExN",        "Mid-layer IT ExN",                           3L, "LRRC4C", "RORB",      "#66A61E",
  "ExN",        "Immature upper-layer ExN",                   4L, "UNC5D",  "KCNQ3",     "#4C78A8",
  "ExN",        "L6 CT-like ExN",                             5L, "TLE4",   "FOXP2",     "#2F5597",
  "ExN",        "L5 deep-layer ExN",                          6L, "RELN",   "GRIK2",     "#7A5195",
  "InN",        "Developing GABAergic neurons",               7L, "GAD1",   "ST18",      "#E69F00",
  "InN",        "Immature GABAergic neurons",                 8L, "GAD2",   "DLX6-AS1",  "#C44E52",
  "OPC",        "Early OPC",                                  9L, "OLIG1",  "PDGFRA",    "#00A6D6",
  "Astro",      "Developmental astrocytes",                  10L, "GFAP",   "AQP4",      "#E377C2",
  "Perivascular", "Perivascular mesenchymal",                11L, "PDGFRB", "PTN",       "#A676B7"
) %>%
  mutate(
    marker_pair = paste(gene1, gene2, sep = "–"),
    cell_type = factor(
      cell_type,
      levels = cell_type
    )
  )

# Header gradients use the same broad-class anchor colors as the
# corresponding cell-type circles.
header_colors <- c(
  "ExN" = "#1B9E77",
  "InN" = "#E69F00",
  "OPC" = "#00A6D6",
  "Astro" = "#E377C2",
  "Perivascular" = "#A676B7"
)

write.csv(
  legend_map,
  file.path(output_dir, "cell_type_legend_key.csv"),
  row.names = FALSE
)

############################################################
# 3) HELPER: ONE DASHED GROUP PANEL
############################################################

make_legend_panel <- function(
    group_name,
    x_max = 10.4,
    gradient_steps = 260L
) {
  
  panel_df <- legend_map %>%
    filter(broad_group == group_name) %>%
    mutate(
      row_id = row_number(),
      y = rev(seq_len(n()))
    )
  
  n_rows <- nrow(panel_df)
  
  if (n_rows == 0L) {
    stop("No rows found for legend group: ", group_name)
  }
  
  header_ymin <- n_rows + 0.40
  header_ymax <- n_rows + 1.10
  border_ymin <- 0.30
  border_ymax <- n_rows + 1.35
  
  gradient_x <- seq(
    0.35,
    x_max - 0.35,
    length.out = gradient_steps + 1L
  )
  
  gradient_cols <- grDevices::colorRampPalette(
    c(header_colors[[group_name]], "white")
  )(gradient_steps)
  
  gradient_df <- tibble(
    xmin = gradient_x[-length(gradient_x)],
    xmax = gradient_x[-1L],
    ymin = header_ymin,
    ymax = header_ymax,
    fill_hex = gradient_cols
  )
  
  ggplot() +
    # Broad-class gradient strip.
    geom_rect(
      data = gradient_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill_hex
      ),
      color = NA
    ) +
    
    # Exact heatmap-color circles.
    geom_point(
      data = panel_df,
      aes(x = 0.98, y = y, fill = color),
      shape = 21,
      size = 10.2,
      stroke = 1.05,
      color = "black"
    ) +
    
    # t-SNE figure number.
    geom_text(
      data = panel_df,
      aes(x = 0.98, y = y, label = plot_number),
      family = font_family,
      fontface = "bold",
      size = 5.0,
      color = "black"
    ) +
    
    # Exact heatmap cell-type name.
    geom_text(
      data = panel_df,
      aes(x = 1.78, y = y + 0.16, label = cell_type),
      hjust = 0,
      vjust = 0.5,
      family = font_family,
      fontface = "bold",
      size = 4.15,
      color = "black"
    ) +
    
    # User-selected marker pair.
    geom_text(
      data = panel_df,
      aes(x = 1.78, y = y - 0.23, label = marker_pair),
      hjust = 0,
      vjust = 0.5,
      family = font_family,
      size = 3.75,
      color = "#222222"
    ) +
    
    # Broad-class header.
    annotate(
      "text",
      x = 0.55,
      y = (header_ymin + header_ymax) / 2,
      label = group_name,
      hjust = 0,
      vjust = 0.5,
      family = font_family,
      fontface = "bold",
      size = 6.2,
      color = "black"
    ) +
    
    # Dashed panel border, matching the reference design.
    annotate(
      "rect",
      xmin = 0.06,
      xmax = x_max - 0.06,
      ymin = border_ymin,
      ymax = border_ymax,
      fill = NA,
      color = "black",
      linewidth = 0.82,
      linetype = "22"
    ) +
    
    scale_fill_identity() +
    coord_cartesian(
      xlim = c(0, x_max),
      ylim = c(0.20, border_ymax + 0.06),
      clip = "off",
      expand = FALSE
    ) +
    theme_void(base_family = font_family) +
    theme(
      plot.margin = margin(2, 2, 2, 2, unit = "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

############################################################
# 4) BUILD PANELS
############################################################

p_exn <- make_legend_panel("ExN")
p_inn <- make_legend_panel("InN")
p_opc <- make_legend_panel("OPC")
p_astro <- make_legend_panel("Astro")
p_perivascular <- make_legend_panel("Perivascular")

# Layout modeled after the supplied reference:
# left column = Astro, InN, Perivascular
# right column = ExN, OPC
left_column <- (
  p_astro /
    p_inn /
    p_perivascular /
    patchwork::plot_spacer()
) +
  patchwork::plot_layout(
    heights = c(2.0, 3.0, 2.0, 2.15)
  )

right_column <- (
  p_exn /
    p_opc
) +
  patchwork::plot_layout(
    heights = c(7.15, 2.0)
  )

p_legend <- ((
  left_column |
    right_column
) +
  patchwork::plot_layout(
    widths = c(1.03, 1.42)
  )) &
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p_legend)

############################################################
# 5) SAVE PUBLICATION FILES
############################################################

pdf_file <- file.path(
  output_dir,
  "Figure_cell_type_number_legend_exact_heatmap_colors.pdf"
)

tiff_file <- file.path(
  output_dir,
  "Figure_cell_type_number_legend_exact_heatmap_colors.tiff"
)

png_file <- file.path(
  output_dir,
  "Figure_cell_type_number_legend_exact_heatmap_colors.png"
)

# Vector PDF; use Cairo when available to reduce font problems.
pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  "pdf"
}

ggsave(
  filename = pdf_file,
  plot = p_legend,
  width = 11.2,
  height = 7.4,
  units = "in",
  device = pdf_device,
  bg = "white"
)

# High-resolution TIFF for journal submission.
ggsave(
  filename = tiff_file,
  plot = p_legend,
  width = 11.2,
  height = 7.4,
  units = "in",
  dpi = 1000,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

# Convenient preview.
ggsave(
  filename = png_file,
  plot = p_legend,
  width = 11.2,
  height = 7.4,
  units = "in",
  dpi = 400,
  device = "png",
  bg = "white"
)

message("Legend complete.")
message("PDF: ", normalizePath(pdf_file, winslash = "/", mustWork = FALSE))
message("TIFF: ", normalizePath(tiff_file, winslash = "/", mustWork = FALSE))
message("PNG: ", normalizePath(png_file, winslash = "/", mustWork = FALSE))

################################################################################
