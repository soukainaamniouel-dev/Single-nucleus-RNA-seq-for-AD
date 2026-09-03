
############################################################
# Human brain organoid scRNA-seq pipeline
# Seurat v5 | SCT integration | per-sample DoubletFinder
# Design: 2 iPSC lines (206, 114) x genotype x 2 timepoints
#         *** n = 1 per condition -> DESCRIPTIVE analysis ***
############################################################

##########################################################
# ) PACKAGES
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

###########################################################
# 1) PATHS + HELPERS
############################################################
out_folder  <- "output_new_process_scRNAseq"
out_name    <- "v11_output"
plot_folder <- "plots"

out_path  <- file.path(out_folder, out_name)
plot_path <- file.path(out_path, plot_folder)

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(out_path,   showWarnings = FALSE, recursive = TRUE)
dir.create(plot_path,  showWarnings = FALSE, recursive = TRUE)

save_plot_pdf <- function(p, filename, width = 7, height = 5) {
  pdf(filename, width = width, height = height); print(p); dev.off()
}

###########################################################
# 2) READ 10X + CREATE SEURAT OBJECTS (per sample)
############################################################

############################################################
# READ 10X DATA + CREATE SEURAT OBJECTS
# Samples: Ctrl, Fibrils, FibJ8
############################################################

base_dir <- paste0(
  "H:/Documents/Qi_Projects/AD_scRNAseq/",
  "scRNAseq_data/Cell_Ranger/"
)

sample_dirs <- list(
  Ctrl = file.path(
    base_dir,
    "Ctrl_count/outs/filtered_feature_bc_matrix"
  ),
  Fibrils = file.path(
    base_dir,
    "Fibrils_count/outs/filtered_feature_bc_matrix"
  ),
  FibJ8 = file.path(
    base_dir,
    "FibJ8_count/outs/filtered_feature_bc_matrix"
  )
)

############################################################
# Experimental design metadata
############################################################

sample_metadata <- data.frame(
  sample = c(
    "Ctrl",
    "Fibrils",
    "FibJ8"
  ),
  condition = c(
    "Ctrl",
    "Fibrils",
    "FibJ8"
  ),
  fibrils = c(
    "No",
    "Yes",
    "Yes"
  ),
  J8_treatment = c(
    "No",
    "No",
    "Yes"
  ),
  stringsAsFactors = FALSE
)

rownames(sample_metadata) <- sample_metadata$sample

# Confirm that metadata and directory names match
if (!setequal(names(sample_dirs), rownames(sample_metadata))) {
  stop(
    "The names in sample_dirs do not match the samples ",
    "in sample_metadata."
  )
}

############################################################
# Read each Cell Ranger output and create one Seurat object
############################################################

seurat_list <- lapply(names(sample_dirs), function(nm) {
  
  message("Reading sample: ", nm)
  
  counts <- Read10X(
    data.dir = sample_dirs[[nm]]
  )
  
  # Read10X returns a list when multiple feature types are present,
  # for example Gene Expression and Antibody Capture.
  if (is.list(counts)) {
    
    if ("Gene Expression" %in% names(counts)) {
      
      counts <- counts[["Gene Expression"]]
      
    } else {
      
      stop(
        "No 'Gene Expression' matrix found for sample: ",
        nm,
        ". Available feature types: ",
        paste(names(counts), collapse = ", ")
      )
    }
  }
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = nm,
    assay = "RNA"
  )
  
  # Add sample-level metadata
  obj$sample <- nm
  
  obj$condition <- sample_metadata[
    nm,
    "condition"
  ]
  
  obj$fibrils <- sample_metadata[
    nm,
    "fibrils"
  ]
  
  obj$J8_treatment <- sample_metadata[
    nm,
    "J8_treatment"
  ]
  
  # Give each cell a sample-specific cell name.
  # Example: Ctrl_AAAC..., Fibrils_AAAC..., FibJ8_AAAC...
  obj <- RenameCells(
    object = obj,
    add.cell.id = nm
  )
  
  obj
})

names(seurat_list) <- names(sample_dirs)

############################################################
# Set consistent factor levels
############################################################

seurat_list <- lapply(seurat_list, function(obj) {
  
  obj$sample <- factor(
    obj$sample,
    levels = c(
      "Ctrl",
      "Fibrils",
      "FibJ8"
    )
  )
  
  obj$condition <- factor(
    obj$condition,
    levels = c(
      "Ctrl",
      "Fibrils",
      "FibJ8"
    )
  )
  
  obj$fibrils <- factor(
    obj$fibrils,
    levels = c(
      "No",
      "Yes"
    )
  )
  
  obj$J8_treatment <- factor(
    obj$J8_treatment,
    levels = c(
      "No",
      "Yes"
    )
  )
  
  obj
})

############################################################
# Verify sample metadata and cell counts
############################################################

sample_check <- do.call(
  rbind,
  lapply(names(seurat_list), function(nm) {
    
    obj <- seurat_list[[nm]]
    
    data.frame(
      sample = nm,
      condition = as.character(unique(obj$condition)),
      fibrils = as.character(unique(obj$fibrils)),
      J8_treatment = as.character(unique(obj$J8_treatment)),
      n_cells = ncol(obj),
      stringsAsFactors = FALSE
    )
  })
)

print(sample_check)
###########################################################
# 3) MERGE + QC METRICS
############################################################
dataset.combined <- merge(
  seurat_list[[1]], y = seurat_list[-1],
  add.cell.ids = names(seurat_list), project = "scRNAseq"
)

# Mito % (human -> ^MT-)
dataset.combined <- PercentageFeatureSet(dataset.combined, pattern = "^MT-",
                                         col.name = "percent_mt")

lvls <- names(sample_dirs)
dataset.combined$sample    <- factor(dataset.combined$sample,    levels = lvls)
dataset.combined$condition <- factor(dataset.combined$condition, levels = lvls)
Idents(dataset.combined) <- "sample"

# Pre-filter QC plots
p_vln_qc <- VlnPlot(dataset.combined,
                    features = c("nFeature_RNA","nCount_RNA","percent_mt"),
                    group.by = "sample", ncol = 3, pt.size = 0)
ggsave(file.path(plot_path, "qc_violin_prefilter_by_sample.pdf"),
       p_vln_qc, width = 12, height = 6)

p_s1 <- FeatureScatter(dataset.combined, "nCount_RNA", "percent_mt",   group.by = "sample")
p_s2 <- FeatureScatter(dataset.combined, "nCount_RNA", "nFeature_RNA", group.by = "sample")
save_plot_pdf(p_s1 + p_s2, file.path(plot_path, "qc_scatter_prefilter.pdf"), 10, 5)

###########################################################
# 4) FILTERING
#    NOTE for organoids: inspect the violins before trusting
#    percent_mt < 20 (hypoxic cores can carry real high-mito cells).
############################################################
dataset.filtered <- subset(
  dataset.combined,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mt < 20
)
write.csv(table(dataset.filtered$sample),
          file.path(out_path, "metadata_cells_per_sample_filtered.csv"))


###########################################################
# 5) DOUBLET DETECTION with scDblFinder (per sample, pre-integration)
#    scDblFinder is matrix-native (works on a SingleCellExperiment),
#    so it sidesteps ALL the Seurat v5 assay-layer issues that make
#    DoubletFinder fragile on v5. It runs per-sample via samples=,
#    does its own normalization/PCA internally, and is the current
#    Bioconductor-standard doublet caller.
############################################################
plan(sequential)

s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
pcs_use   <- 1:20

# --- ensure a single clean 'counts' layer (v5 merge() splits them) ---
dataset.filtered[["RNA"]] <- JoinLayers(dataset.filtered[["RNA"]])
message("RNA layers after join: ",
        paste(Layers(dataset.filtered[["RNA"]]), collapse = ", "))

# --- run scDblFinder on raw counts, per sample ---
cts <- GetAssayData(dataset.filtered, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays = list(counts = cts))
sce$sample <- dataset.filtered$sample          # per-sample detection

set.seed(1234)

sce <- scDblFinder(sce, samples = "sample")    # detects doublets within each sample

# --- harvest calls back onto the Seurat object (matched by barcode) ---
m <- match(colnames(dataset.filtered), colnames(sce))
dataset.filtered$scDblFinder.class <- sce$scDblFinder.class[m]
dataset.filtered$scDblFinder.score <- sce$scDblFinder.score[m]
# unify column name used downstream (labels are lowercase: singlet/doublet)
dataset.filtered$DF_class <- ifelse(
  dataset.filtered$scDblFinder.class == "doublet", "Doublet", "Singlet"
)

# --- summary table (QC check: expect ~5-15% doublets per sample) ---
doublet_summary <- as.data.frame.matrix(
  table(dataset.filtered$sample, dataset.filtered$DF_class)
)
doublet_summary$sample <- rownames(doublet_summary)
doublet_summary$total  <- doublet_summary$Singlet + doublet_summary$Doublet
doublet_summary$pct_doublet <- round(100 * doublet_summary$Doublet /
                                       doublet_summary$total, 2)
print(doublet_summary)
write.csv(doublet_summary,
          file.path(out_path, "doublet_summary_per_sample.csv"), row.names = FALSE)

# --- doublet QC UMAPs per sample (BEFORE removal) ---
# quick per-sample embedding just for the QC figure
qc_umap <- dataset.filtered
DefaultAssay(qc_umap) <- "RNA"
qc_umap <- NormalizeData(qc_umap, verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(dims = pcs_use, verbose = FALSE)
p_db <- DimPlot(qc_umap, group.by = "DF_class", split.by = "sample", ncol = 4,
                order = "Doublet",
                cols = c("Singlet" = "grey80", "Doublet" = "red"), pt.size = 0.3) +
  ggtitle("scDblFinder calls (pre-removal)")
ggsave(file.path(plot_path, "doublets_umap_by_sample.pdf"), p_db, width = 16, height = 8)
rm(qc_umap); gc()

# --- REMOVE doublets ---
dataset.filtered <- subset(dataset.filtered, subset = DF_class == "Singlet")
message("After doublet removal: ", ncol(dataset.filtered), " cells remain.")

dataset.filtered <- subset(
  dataset.filtered,
  subset = DF_class == "Singlet"
)
############################################################
# 5b) SPLIT + PER-SAMPLE SCT (on singlets) for integration
############################################################
options(future.globals.maxSize = 8 * 1024^3)
obj.list <- SplitObject(dataset.filtered, split.by = "sample")

obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x, assay = "RNA", verbose = FALSE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes,
                        set.ident = FALSE)
  x <- SCTransform(x, assay = "RNA", new.assay.name = "SCT",
                   variable.features.n = 3000,
                   vars.to.regress = c("percent_mt","S.Score","G2M.Score"),
                   method = "glmGamPoi", verbose = FALSE)
  DefaultAssay(x) <- "SCT"
  x
})


###########################################################
# 6) (optional) PRE-INTEGRATION LOOK (doublet-free)
############################################################
pre <- merge(obj.list[[1]], y = obj.list[-1])
DefaultAssay(pre) <- "RNA"
pre <- NormalizeData(pre, verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)
p_pre <- DimPlot(pre, group.by = "sample") + ggtitle("Before integration (singlets)")
ggsave(file.path(plot_path, "umap_before_integration.pdf"), p_pre, width = 8, height = 6)
rm(pre); gc()

###########################################################
# CCA (SCT) vs Harmony (log-norm) integration benchmark
# Run AFTER section 5b (obj.list = 8 singlet-only SCT objects,
# dataset.filtered = merged singlets with joined RNA counts).
#
# Goal: compare integration QUALITY + SPEED on YOUR data, then
# pick one method for the expensive downstream (annotation/DE).
############################################################
library(harmony)   # install.packages("harmony") if needed
n_pcs_use <- 20
res_use   <- 0.5

timings <- list()   # collect wall-clock times

###########################################################
# PATH A -- SCT + CCA  (your reference method)
############################################################
t0 <- Sys.time()

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = obj.list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
cca <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(cca) <- "integrated"
cca <- ScaleData(cca, verbose = FALSE)
cca <- RunPCA(cca, npcs = 50, verbose = FALSE)
cca <- FindNeighbors(cca, dims = 1:n_pcs_use)
cca <- FindClusters(cca, resolution = res_use)
cca <- RunUMAP(cca, dims = 1:n_pcs_use, reduction = "pca")

timings$CCA <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("CCA done in %.1f min", as.numeric(timings$CCA)))

ggsave(file.path(plot_path, "elbowplot_cca.pdf"),
       ElbowPlot(cca, ndims = 50), width = 7, height = 5)

###########################################################
# PATH B -- log-normalize + Harmony  (fast candidate)
#   Harmony is standardly run on log-normalized PCA, not SCT.
#   This compares the two COMPLETE workflows you'd choose between.
############################################################
t0 <- Sys.time()

# start from the merged singlet object (RNA counts already joined in 5b)
harm <- dataset.filtered
DefaultAssay(harm) <- "RNA"
harm <- NormalizeData(harm, verbose = FALSE)
harm <- FindVariableFeatures(harm, nfeatures = 3000, verbose = FALSE)
harm <- ScaleData(harm, vars.to.regress = c("percent_mt","S.Score","G2M.Score"),
                  verbose = FALSE)  # match CCA's regression for fair comparison
harm <- RunPCA(harm, npcs = 50, verbose = FALSE)

# the one-line integration: correct the PCA embedding by sample
harm <- harmony::RunHarmony(harm, group.by.vars = "sample",
                            reduction.use = "pca",
                            reduction.save = "harmony")

harm <- FindNeighbors(harm, reduction = "harmony", dims = 1:n_pcs_use)
harm <- FindClusters(harm, resolution = res_use)
harm <- RunUMAP(harm, reduction = "harmony", dims = 1:n_pcs_use)

timings$Harmony <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("Harmony done in %.1f min", as.numeric(timings$Harmony)))

###########################################################
# COMPARISON 1 -- SPEED
############################################################
speed_tab <- data.frame(
  method   = c("CCA (SCT)", "Harmony (lognorm)"),
  minutes  = c(as.numeric(timings$CCA), as.numeric(timings$Harmony))
)
print(speed_tab)
write.csv(speed_tab, file.path(out_path, "cmp_speed.csv"), row.names = FALSE)

###########################################################
# COMPARISON 2 -- MIXING (are samples well integrated?)
#   iLISI: higher = better sample mixing. Compare distributions.
############################################################
lisi_cca  <- compute_lisi(Embeddings(cca,  "pca")[, 1:n_pcs_use],
                          cca@meta.data,  "sample")$sample
lisi_harm <- compute_lisi(Embeddings(harm, "harmony")[, 1:n_pcs_use],
                          harm@meta.data, "sample")$sample

mix_tab <- data.frame(
  method    = c("CCA", "Harmony"),
  mean_iLISI = c(mean(lisi_cca), mean(lisi_harm)),
  med_iLISI  = c(median(lisi_cca), median(lisi_harm))
)
print(mix_tab)
write.csv(mix_tab, file.path(out_path, "cmp_mixing_iLISI.csv"), row.names = FALSE)

pdf(file.path(plot_path, "cmp_iLISI_hist.pdf"), width = 10, height = 4)
par(mfrow = c(1, 2))
hist(lisi_cca,  breaks = 30, main = "CCA iLISI",     xlab = "iLISI")
hist(lisi_harm, breaks = 30, main = "Harmony iLISI", xlab = "iLISI")
dev.off()

###########################################################
# COMPARISON 3 -- STRUCTURE PRESERVED (not over-mixed?)
#   n clusters + side-by-side UMAPs. Over-integration shows up as
#   too few clusters / biology collapsed together.
############################################################
struct_tab <- data.frame(
  method     = c("CCA", "Harmony"),
  n_clusters = c(length(unique(cca$seurat_clusters)),
                 length(unique(harm$seurat_clusters))),
  n_cells    = c(ncol(cca), ncol(harm))
)
print(struct_tab)
write.csv(struct_tab, file.path(out_path, "cmp_structure.csv"), row.names = FALSE)

p_cca  <- DimPlot(cca,  group.by = "seurat_clusters", label = TRUE) +
  NoLegend() + ggtitle("CCA: clusters")
p_harm <- DimPlot(harm, group.by = "seurat_clusters", label = TRUE) +
  NoLegend() + ggtitle("Harmony: clusters")
p_cca_s  <- DimPlot(cca,  group.by = "sample") + ggtitle("CCA: sample")
p_harm_s <- DimPlot(harm, group.by = "sample") + ggtitle("Harmony: sample")

ggsave(file.path(plot_path, "cmp_umap_clusters.pdf"),
       p_cca | p_harm, width = 14, height = 6)
ggsave(file.path(plot_path, "cmp_umap_sample.pdf"),
       p_cca_s | p_harm_s, width = 16, height = 6)

###########################################################
# COMPARISON 4 -- KEY MARKERS look sane in BOTH?
#   The real test: do known organoid populations resolve cleanly
#   in each? Score expression on RNA (join layers first).
############################################################
check_markers <- c("SOX2","PAX6","EOMES","NEUROD6","RBFOX3","DCX",
                   "GFAP","AQP4","MKI67","DLX2")

for (obj_name in c("cca","harm")) {
  o <- get(obj_name)
  DefaultAssay(o) <- if ("integrated" %in% names(o@assays)) "RNA" else "RNA"
  o[["RNA"]] <- JoinLayers(o[["RNA"]])
  o <- NormalizeData(o, verbose = FALSE)
  p <- FeaturePlot(o, features = check_markers, order = TRUE, ncol = 5)
  ggsave(file.path(plot_path, paste0("cmp_markers_", obj_name, ".pdf")),
         p, width = 20, height = 8)
  assign(obj_name, o)
}

###########################################################
# SAVE both so you don't have to rerun; annotate ONLY the winner
############################################################
qsave(cca,  file.path(out_path, "integrated_CCA.qs"))
qsave(harm, file.path(out_path, "integrated_Harmony.qs"))


############################################################
# CCA CLUSTER MARKERS
#
# Cluster identities: CCA-integrated clustering
# Marker expression: joined, log-normalized RNA assay
############################################################
# ----------------------------------------------------------
# 1) Preserve and verify the CCA cluster assignments
# ----------------------------------------------------------

if (!"seurat_clusters" %in% colnames(cca[[]])) {
  stop(
    "'seurat_clusters' is missing from the CCA object. ",
    "Run FindNeighbors() and FindClusters() first."
  )
}

# Save an explicitly named CCA cluster column
cca$cca_cluster <- factor(
  as.character(cca$seurat_clusters),
  levels = levels(cca$seurat_clusters)
)

Idents(cca) <- "cca_cluster"

message(
  "Finding markers for ",
  length(levels(Idents(cca))),
  " CCA-derived clusters: ",
  paste(levels(Idents(cca)), collapse = ", ")
)

# Number of cells in each cluster
cca_cluster_sizes <- data.frame(
  cluster = names(table(Idents(cca))),
  n_cells = as.integer(table(Idents(cca))),
  row.names = NULL
)

print(cca_cluster_sizes)

write.csv(
  cca_cluster_sizes,
  file.path(out_path, "cca_cluster_cell_counts.csv"),
  row.names = FALSE
)

# Warn about very small clusters
small_clusters <- cca_cluster_sizes$cluster[
  cca_cluster_sizes$n_cells < 20
]

if (length(small_clusters) > 0) {
  warning(
    "These CCA clusters contain fewer than 20 cells and their ",
    "marker results should be interpreted cautiously: ",
    paste(small_clusters, collapse = ", ")
  )
}

# ----------------------------------------------------------
# 2) Check cluster representation across samples
# ----------------------------------------------------------

cca_cluster_by_sample <- table(
  cluster = cca$cca_cluster,
  sample = cca$sample
)

print(cca_cluster_by_sample)

write.csv(
  as.data.frame.matrix(cca_cluster_by_sample),
  file.path(out_path, "cca_cluster_cell_counts_by_sample.csv"),
  row.names = TRUE
)

cca_cluster_pct_by_sample <- round(
  100 * prop.table(cca_cluster_by_sample, margin = 1),
  digits = 2
)

write.csv(
  as.data.frame.matrix(cca_cluster_pct_by_sample),
  file.path(out_path, "cca_cluster_percent_by_sample.csv"),
  row.names = TRUE
)

# ----------------------------------------------------------
# 3) Prepare unintegrated RNA expression for marker testing
# ----------------------------------------------------------

DefaultAssay(cca) <- "RNA"

# Required when RNA data are stored as multiple Seurat v5 layers
cca[["RNA"]] <- JoinLayers(cca[["RNA"]])

message(
  "CCA RNA layers used for marker testing: ",
  paste(Layers(cca[["RNA"]]), collapse = ", ")
)

if (!"counts" %in% Layers(cca[["RNA"]])) {
  stop("The joined RNA assay does not contain a counts layer.")
}

# Create/recreate the log-normalized RNA data layer
cca <- NormalizeData(
  object = cca,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

if (!"data" %in% Layers(cca[["RNA"]])) {
  stop("NormalizeData() did not create an RNA data layer.")
}

# Restore CCA clusters as identities after assay preparation
Idents(cca) <- "cca_cluster"

# ----------------------------------------------------------
# 4) Find positive markers for every CCA cluster
#
# Each cluster is compared against all other cells.
# ----------------------------------------------------------

set.seed(1234)

cca_markers <- FindAllMarkers(
  object = cca,
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05,
  random.seed = 1234,
  verbose = TRUE
)

if (nrow(cca_markers) == 0) {
  stop(
    "FindAllMarkers() returned no markers. Consider lowering ",
    "min.pct or logfc.threshold."
  )
}

# Make sure gene names are stored in an explicit column
if (!"gene" %in% colnames(cca_markers)) {
  cca_markers$gene <- rownames(cca_markers)
}

# Handle Seurat versions that use different FC column names
fc_col <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(cca_markers)
)[1]

if (is.na(fc_col)) {
  stop(
    "Could not find avg_log2FC or avg_logFC in the ",
    "FindAllMarkers output."
  )
}

# Detection difference is useful for ranking marker specificity
cca_markers$delta_pct <- (
  cca_markers$pct.1 -
    cca_markers$pct.2
)

cca_markers <- cca_markers %>%
  arrange(
    cluster,
    desc(.data[[fc_col]]),
    desc(delta_pct)
  )

# Save complete marker table
write.csv(
  cca_markers,
  file.path(out_path, "cca_all_positive_cluster_markers.csv"),
  row.names = FALSE
)

# Also save an R object preserving column types
saveRDS(
  cca_markers,
  file.path(out_path, "cca_all_positive_cluster_markers.rds")
)

# ----------------------------------------------------------
# 5) Select top 10 markers per cluster
# ----------------------------------------------------------

cca_top10_markers <- cca_markers %>%
  filter(
    p_val_adj < 0.05,
    .data[[fc_col]] >= 0.25
  ) %>%
  group_by(cluster) %>%
  arrange(
    desc(.data[[fc_col]]),
    desc(delta_pct),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  cca_top10_markers,
  file.path(out_path, "cca_top10_markers_per_cluster.csv"),
  row.names = FALSE
)

print(
  cca_top10_markers %>%
    select(
      cluster,
      gene,
      all_of(fc_col),
      pct.1,
      pct.2,
      delta_pct,
      p_val_adj
    ),
  n = Inf
)

############################################################
# Number of markers recovered for each CCA cluster
############################################################

# Confirm that cca_markers is still the FindAllMarkers table
if (!is.data.frame(cca_markers)) {
  stop(
    "'cca_markers' is not a data frame. Current class: ",
    paste(class(cca_markers), collapse = ", "),
    ". It may have been overwritten by a list."
  )
}

required_marker_columns <- c("cluster", "gene")

missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(cca_markers)
)

if (length(missing_marker_columns) > 0) {
  stop(
    "Missing required column(s) from cca_markers: ",
    paste(missing_marker_columns, collapse = ", ")
  )
}

# Make cluster a standard character vector
cca_markers$cluster <- as.character(
  cca_markers[["cluster"]]
)

# Count distinct marker genes per cluster
cca_marker_counts <- cca_markers %>%
  dplyr::group_by(.data$cluster) %>%
  dplyr::summarise(
    n_markers = dplyr::n_distinct(.data$gene),
    .groups = "drop"
  )

print(cca_marker_counts)

write.csv(
  cca_marker_counts,
  file.path(
    out_path,
    "cca_number_of_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

# ----------------------------------------------------------
# 6) Create an annotation-friendly marker table
#
# Keep the complete table above. This second table only removes
# commonly uninformative mitochondrial/ribosomal/housekeeping
# genes from the top-marker ranking.
# ----------------------------------------------------------

cca_annotation_markers <- cca_markers %>%
  filter(
    p_val_adj < 0.05,
    .data[[fc_col]] >= 0.25,
    delta_pct >= 0.10,
    !grepl("^(MT-|RPL|RPS|HIST)", gene),
    gene != "MALAT1"
  )

cca_top10_annotation <- cca_annotation_markers %>%
  group_by(cluster) %>%
  arrange(
    desc(.data[[fc_col]]),
    desc(delta_pct),
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  cca_annotation_markers,
  file.path(out_path, "cca_annotation_candidate_markers.csv"),
  row.names = FALSE
)

write.csv(
  cca_top10_annotation,
  file.path(out_path, "cca_top10_annotation_markers_per_cluster.csv"),
  row.names = FALSE
)
# Marker table split into one data frame per cluster
cca_markers_by_cluster <- split(
  cca_markers,
  cca_markers$cluster
)

saveRDS(
  cca_markers_by_cluster,
  file.path(out_path, "cca_markers_split_by_cluster.rds")
)

############################################################
# 7) Dot plot of top annotation markers
############################################################

cca_top5_genes <- cca_top10_annotation %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

if (length(cca_top5_genes) > 0) {
  
  p_cca_markers_dot <- DotPlot(
    object = cca,
    features = cca_top5_genes,
    assay = "RNA",
    group.by = "cca_cluster",
    dot.scale = 6
  ) +
    RotatedAxis() +
    ggtitle("Top markers of CCA-derived clusters") +
    xlab(NULL) +
    ylab("CCA cluster")
  
  ggsave(
    filename = file.path(
      plot_path,
      "cca_top_cluster_markers_dotplot.pdf"
    ),
    plot = p_cca_markers_dot,
    width = max(12, 0.28 * length(cca_top5_genes)),
    height = max(6, 0.40 * length(levels(Idents(cca))))
  )
}

############################################################
# 8) Heatmap of top annotation markers
############################################################

if (length(cca_top5_genes) > 0) {
  
  # Scale only the marker genes required for this heatmap
  cca <- ScaleData(
    object = cca,
    assay = "RNA",
    features = cca_top5_genes,
    verbose = FALSE
  )
  
  # Downsample cells for a readable heatmap
  set.seed(1234)
  
  heatmap_cells <- unlist(
    lapply(levels(Idents(cca)), function(cl) {
      
      cluster_cells <- WhichCells(
        object = cca,
        idents = cl
      )
      
      sample(
        cluster_cells,
        size = min(100, length(cluster_cells))
      )
    }),
    use.names = FALSE
  )
  
  p_cca_markers_heatmap <- DoHeatmap(
    object = cca,
    features = cca_top5_genes,
    cells = heatmap_cells,
    assay = "RNA",
    group.by = "cca_cluster",
    slot = "scale.data",
    raster = TRUE
  ) +
    NoLegend() +
    ggtitle("Top markers of CCA-derived clusters")
  
  ggsave(
    filename = file.path(
      plot_path,
      "cca_top_cluster_markers_heatmap.pdf"
    ),
    plot = p_cca_markers_heatmap,
    width = 14,
    height = max(8, 0.22 * length(cca_top5_genes))
  )
}

############################################################
# 9) Save updated CCA object
############################################################

qsave(
  cca,
  file.path(out_path, "integrated_CCA_with_cluster_markers.qs")
)

message(
  "CCA cluster-marker analysis complete. Marker tables were ",
  "saved to: ",
  normalizePath(out_path, winslash = "/", mustWork = FALSE)
)

############################################################
# 10) For sparse genes or smaller developmental populations, rerun with: 
############################################################
cca_markers_sensitive <- FindAllMarkers(
  object = cca,
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  return.thresh = 0.05,
  verbose = TRUE
)

write.csv(
  cca_markers_sensitive,
  file.path(out_path, "cca_cluster_markers_sensitive.csv"),
  row.names = FALSE
)

############################################################
# OPTIONAL: CONSERVED MARKERS ACROSS SAMPLES
############################################################

if (!requireNamespace("metap", quietly = TRUE)) {
  message(
    "Skipping conserved markers because the 'metap' package ",
    "is not installed."
  )
} else {
  
  Idents(cca) <- "cca_cluster"
  DefaultAssay(cca) <- "RNA"
  
  cluster_sample_n <- table(
    cca$cca_cluster,
    cca$sample
  )
  
  # Require at least 10 cells from every sample in a cluster
  eligible_clusters <- rownames(cluster_sample_n)[
    apply(
      cluster_sample_n,
      MARGIN = 1,
      FUN = function(x) all(x >= 10)
    )
  ]
  
  message(
    "Clusters eligible for conserved-marker testing: ",
    paste(eligible_clusters, collapse = ", ")
  )
  
  conserved_marker_list <- lapply(
    eligible_clusters,
    function(cl) {
      
      message(
        "Finding conserved markers for CCA cluster ",
        cl
      )
      
      result <- FindConservedMarkers(
        object = cca,
        ident.1 = cl,
        grouping.var = "sample",
        assay = "RNA",
        slot = "data",
        only.pos = TRUE,
        min.pct = 0.10,
        logfc.threshold = 0.25,
        verbose = FALSE
      )
      
      result$gene <- rownames(result)
      result$cluster <- cl
      rownames(result) <- NULL
      
      result
    }
  )
  
  names(conserved_marker_list) <- eligible_clusters
  
  cca_conserved_markers <- bind_rows(
    conserved_marker_list
  )
  
  write.csv(
    cca_conserved_markers,
    file.path(
      out_path,
      "cca_conserved_markers_across_samples.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    conserved_marker_list,
    file.path(
      out_path,
      "cca_conserved_markers_by_cluster.rds"
    )
  )
}

############################################################
# STRICT CONSERVED CCA MARKERS ACROSS ALL THREE SAMPLES
#
# A gene must satisfy, in Ctrl, Fibrils, and FibJ8:
#   1. avg_log2FC > 1
#   2. Seurat expression p_val_adj < 0.05
#   3. pct.1 - pct.2 >= 0.10
#   4. adjusted Fisher detection-rate p-value < 0.05
############################################################

# ----------------------------------------------------------
# User-defined thresholds
# ----------------------------------------------------------

required_samples <- c(
  "Ctrl",
  "Fibrils",
  "FibJ8"
)

min_cells_per_group <- 10L

min_log2fc <- 1

# Fractions, not percentages:
# 0.10 means a difference of at least 10 percentage points
min_delta_pct <- 0.10

alpha <- 0.05

# Adjustment for the separate detection-rate tests
pct_p_adjust_method <- "BH"


# ----------------------------------------------------------
# Package and metadata checks
# ----------------------------------------------------------

if (!requireNamespace("metap", quietly = TRUE)) {
  stop(
    "The 'metap' package is required for FindConservedMarkers()."
  )
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop(
    "The 'Matrix' package is required for sparse-matrix calculations."
  )
}

if (!"cca_cluster" %in% colnames(cca[[]])) {
  stop(
    "Metadata column 'cca_cluster' is missing."
  )
}

if (!"sample" %in% colnames(cca[[]])) {
  stop(
    "Metadata column 'sample' is missing."
  )
}

samples_present <- unique(
  as.character(cca$sample)
)

if (!setequal(samples_present, required_samples)) {
  stop(
    "Expected samples: ",
    paste(required_samples, collapse = ", "),
    ". Samples found: ",
    paste(samples_present, collapse = ", ")
  )
}


# ----------------------------------------------------------
# Prepare the unintegrated RNA assay
#
# CCA clusters define the cell groups, but RNA expression
# should be used for marker testing.
# ----------------------------------------------------------

Idents(cca) <- "cca_cluster"
DefaultAssay(cca) <- "RNA"

cca[["RNA"]] <- JoinLayers(
  cca[["RNA"]]
)

cca <- NormalizeData(
  object = cca,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

if (!"data" %in% Layers(cca[["RNA"]])) {
  stop(
    "The normalized RNA data layer is missing."
  )
}

Idents(cca) <- "cca_cluster"

# Log-normalized expression matrix.
# Expression > 0 is used to define detection.
rna_data <- GetAssayData(
  object = cca,
  assay = "RNA",
  layer = "data"
)


# ----------------------------------------------------------
# Determine which clusters have enough cells in every sample
#
# Require at least min_cells_per_group cells in:
#   - the cluster
#   - all other clusters combined
# within every sample.
# ----------------------------------------------------------

cluster_ids <- levels(
  droplevels(
    factor(
      as.character(cca$cca_cluster)
    )
  )
)

cluster_sample_n <- table(
  cluster = factor(
    as.character(cca$cca_cluster),
    levels = cluster_ids
  ),
  sample = factor(
    as.character(cca$sample),
    levels = required_samples
  )
)

sample_totals <- colSums(
  cluster_sample_n
)

other_sample_n <- matrix(
  sample_totals,
  nrow = nrow(cluster_sample_n),
  ncol = ncol(cluster_sample_n),
  byrow = TRUE,
  dimnames = dimnames(cluster_sample_n)
) - cluster_sample_n

eligible_clusters <- rownames(cluster_sample_n)[
  apply(
    cluster_sample_n >= min_cells_per_group &
      other_sample_n >= min_cells_per_group,
    MARGIN = 1,
    FUN = all
  )
]

message(
  "Clusters eligible for conserved-marker testing: ",
  paste(eligible_clusters, collapse = ", ")
)

write.csv(
  as.data.frame.matrix(cluster_sample_n),
  file.path(
    out_path,
    "cca_conserved_marker_cluster_cells_by_sample.csv"
  ),
  row.names = TRUE
)


# ----------------------------------------------------------
# Helper: add exact detection-rate statistics
#
# For each sample, compare:
#   cells in the selected cluster
# versus
#   all other cells from the same sample
#
# A one-sided Fisher exact test asks whether the gene is
# detected in a greater fraction of cluster cells.
# ----------------------------------------------------------

############################################################
# Add exact detection-rate tests for each sample
############################################################

add_detection_rate_tests <- function(
    marker_table,
    cluster_id
) {
  
  if (
    is.null(marker_table) ||
    nrow(marker_table) == 0
  ) {
    return(NULL)
  }
  
  marker_table$gene <- rownames(marker_table)
  marker_table$cluster <- as.character(cluster_id)
  rownames(marker_table) <- NULL
  
  missing_genes <- setdiff(
    marker_table$gene,
    rownames(rna_data)
  )
  
  if (length(missing_genes) > 0) {
    stop(
      "Some marker genes were not found in the RNA data layer: ",
      paste(
        head(missing_genes),
        collapse = ", "
      )
    )
  }
  
  # Preserve the same gene order as marker_table
  expression_subset <- rna_data[
    marker_table$gene,
    ,
    drop = FALSE
  ]
  
  cluster_vector <- as.character(
    cca$cca_cluster
  )
  
  sample_vector <- as.character(
    cca$sample
  )
  
  for (sample_id in required_samples) {
    
    cluster_cells <- colnames(cca)[
      cluster_vector == as.character(cluster_id) &
        sample_vector == sample_id
    ]
    
    other_cells <- colnames(cca)[
      cluster_vector != as.character(cluster_id) &
        sample_vector == sample_id
    ]
    
    n_cluster <- length(cluster_cells)
    n_other <- length(other_cells)
    
    if (
      n_cluster < min_cells_per_group ||
      n_other < min_cells_per_group
    ) {
      stop(
        "Insufficient cells for cluster ",
        cluster_id,
        " in sample ",
        sample_id,
        ". Cluster cells: ",
        n_cluster,
        "; other cells: ",
        n_other
      )
    }
    
    # Number of cluster cells expressing each gene
    detected_cluster <- Matrix::rowSums(
      expression_subset[
        ,
        cluster_cells,
        drop = FALSE
      ] > 0
    )
    
    # Number of other cells expressing each gene
    detected_other <- Matrix::rowSums(
      expression_subset[
        ,
        other_cells,
        drop = FALSE
      ] > 0
    )
    
    # Exact detection proportions
    pct_cluster <- detected_cluster / n_cluster
    pct_other <- detected_other / n_other
    
    # Positive values mean greater detection in the cluster
    delta_pct <- pct_cluster - pct_other
    
    # One-sided Fisher exact test:
    # detection in cluster > detection outside cluster
    #
    # This vectorized hypergeometric calculation is equivalent
    # to performing a one-sided Fisher test for every gene.
    fisher_p <- stats::phyper(
      q = detected_cluster - 1,
      m = detected_cluster + detected_other,
      n = (
        n_cluster +
          n_other -
          detected_cluster -
          detected_other
      ),
      k = n_cluster,
      lower.tail = FALSE
    )
    
    # Replace any numerical NA values conservatively
    fisher_p[is.na(fisher_p)] <- 1
    
    # Adjust across genes within this cluster/sample comparison
    fisher_p_adj <- stats::p.adjust(
      fisher_p,
      method = pct_p_adjust_method
    )
    
    # Add the statistics as dynamically named columns
    marker_table[[paste0(
      sample_id,
      "_n_cluster"
    )]] <- rep(
      n_cluster,
      nrow(marker_table)
    )
    
    marker_table[[paste0(
      sample_id,
      "_n_other"
    )]] <- rep(
      n_other,
      nrow(marker_table)
    )
    
    marker_table[[paste0(
      sample_id,
      "_detected_cluster"
    )]] <- as.integer(
      detected_cluster
    )
    
    marker_table[[paste0(
      sample_id,
      "_detected_other"
    )]] <- as.integer(
      detected_other
    )
    
    marker_table[[paste0(
      sample_id,
      "_pct_cluster_exact"
    )]] <- as.numeric(
      pct_cluster
    )
    
    marker_table[[paste0(
      sample_id,
      "_pct_other_exact"
    )]] <- as.numeric(
      pct_other
    )
    
    marker_table[[paste0(
      sample_id,
      "_delta_pct_exact"
    )]] <- as.numeric(
      delta_pct
    )
    
    marker_table[[paste0(
      sample_id,
      "_pct_fisher_p_val"
    )]] <- as.numeric(
      fisher_p
    )
    
    marker_table[[paste0(
      sample_id,
      "_pct_fisher_p_val_adj"
    )]] <- as.numeric(
      fisher_p_adj
    )
  }
  
  return(marker_table)
}
# ----------------------------------------------------------
# Run FindConservedMarkers for every eligible cluster
# ----------------------------------------------------------

set.seed(1234)

conserved_marker_list_all <- setNames(
  lapply(
    eligible_clusters,
    function(cluster_id) {
      
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
          only.pos = TRUE,
          
          # Prefilter candidates within every sample
          logfc.threshold = min_log2fc,
          min.pct = 0.10,
          min.diff.pct = min_delta_pct,
          
          min.cells.group = min_cells_per_group,
          base = 2,
          random.seed = 1234,
          meta.method = metap::minimump,
          verbose = FALSE
        ),
        error = function(e) {
          
          warning(
            "Conserved-marker analysis failed for cluster ",
            cluster_id,
            ": ",
            conditionMessage(e)
          )
          
          NULL
        }
      )
      
      if (is.null(result) || nrow(result) == 0) {
        
        warning(
          "No conserved-marker candidates were found for cluster ",
          cluster_id
        )
        
        return(NULL)
      }
      
      result <- add_detection_rate_tests(
        marker_table = result,
        cluster_id = cluster_id
      )
      
      result <- apply_strict_conserved_filter(
        marker_table = result
      )
      
      result
    }
  ),
  eligible_clusters
)

# Remove clusters for which no table was returned
conserved_marker_list_all <- Filter(
  f = Negate(is.null),
  x = conserved_marker_list_all
)


# ----------------------------------------------------------
# Retain only strict markers
# ----------------------------------------------------------

conserved_marker_list_strict <- lapply(
  conserved_marker_list_all,
  function(marker_table) {
    
    marker_table[
      marker_table$pass_strict_all_samples,
      ,
      drop = FALSE
    ]
  }
)

conserved_marker_list_strict <- Filter(
  f = function(x) {
    !is.null(x) && nrow(x) > 0
  },
  x = conserved_marker_list_strict
)


# ----------------------------------------------------------
# Combine and rank results
# ----------------------------------------------------------

cca_conserved_markers_all <- dplyr::bind_rows(
  conserved_marker_list_all
)

cca_conserved_markers_strict <- dplyr::bind_rows(
  conserved_marker_list_strict
)

if (nrow(cca_conserved_markers_strict) > 0) {
  
  cca_conserved_markers_strict <-
    cca_conserved_markers_strict %>%
    dplyr::arrange(
      cluster,
      dplyr::desc(minimum_log2FC),
      dplyr::desc(minimum_delta_pct),
      maximum_pct_fisher_p_adj
    )
}


# ----------------------------------------------------------
# Save both unfiltered and strict tables
# ----------------------------------------------------------

write.csv(
  cca_conserved_markers_all,
  file.path(
    out_path,
    "cca_conserved_markers_all_candidates.csv"
  ),
  row.names = FALSE
)

write.csv(
  cca_conserved_markers_strict,
  file.path(
    out_path,
    paste0(
      "cca_conserved_markers_",
      "log2FC_gt1_",
      "deltaPct_ge10_",
      "significant_all_samples.csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  conserved_marker_list_all,
  file.path(
    out_path,
    "cca_conserved_marker_candidates_by_cluster.rds"
  )
)

saveRDS(
  conserved_marker_list_strict,
  file.path(
    out_path,
    "cca_strict_conserved_markers_by_cluster.rds"
  )
)


# ----------------------------------------------------------
# Summary: number of strict markers per cluster
# ----------------------------------------------------------

if (nrow(cca_conserved_markers_strict) > 0) {
  
  strict_marker_counts <-
    cca_conserved_markers_strict %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      n_strict_markers = dplyr::n_distinct(gene),
      .groups = "drop"
    )
  
} else {
  
  strict_marker_counts <- data.frame(
    cluster = character(),
    n_strict_markers = integer()
  )
  
  warning(
    "No genes passed all strict conserved-marker criteria."
  )
}

print(strict_marker_counts)

write.csv(
  strict_marker_counts,
  file.path(
    out_path,
    "cca_strict_conserved_marker_counts.csv"
  ),
  row.names = FALSE
)


# ----------------------------------------------------------
# Print the most useful columns
# ----------------------------------------------------------

if (nrow(cca_conserved_markers_strict) > 0) {
  
  columns_to_print <- c(
    "cluster",
    "gene",
    
    "Ctrl_avg_log2FC",
    "Fibrils_avg_log2FC",
    "FibJ8_avg_log2FC",
    
    "Ctrl_delta_pct_exact",
    "Fibrils_delta_pct_exact",
    "FibJ8_delta_pct_exact",
    
    "Ctrl_pct_fisher_p_val_adj",
    "Fibrils_pct_fisher_p_val_adj",
    "FibJ8_pct_fisher_p_val_adj",
    
    "minimum_log2FC",
    "minimum_delta_pct"
  )
  
  columns_to_print <- intersect(
    columns_to_print,
    colnames(cca_conserved_markers_strict)
  )
  
  print(
    cca_conserved_markers_strict[
      ,
      columns_to_print,
      drop = FALSE
    ]
  )
}

