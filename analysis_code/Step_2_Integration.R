# Script 2: Multi-Batch Integration, Dimensionality Reduction, and Cell Line Annotation
# Reference: https://satijalab.org/seurat/articles/merge_vignette.html

library(Signac)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(pheatmap) # Added for correlation visualization
library(dplyr)    # For data manipulation
library(cowplot)  # For plotting

# --- USER INPUT & SETUP 🛠️ ---
# Define paths relative to the script location.
# USERS: Pre-processed h5Seurat objects must be placed in the 'data/processed/' directory.
# USERS: Bulk RNA-seq annotation file (GSE89413 file) must be placed in the 'data/annotation/' directory.
data_dir_base <- "data/"
processed_dir <- paste0(data_dir_base, "processed/")
annotation_dir <- paste0(data_dir_base, "annotation/")
output_dir <- "results/integration/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(1234)

# --- LOAD PRE-PROCESSED BATCH OBJECTS ---
NBcell_0601 <- LoadH5Seurat(paste0(processed_dir, "NBcell_0601.h5Seurat"))
NBcell_0602 <- LoadH5Seurat(paste0(processed_dir, "NBcell_0602.h5Seurat"))

# --- PREPARE BATCH METADATA ---
# Prefix cluster IDs with batch IDs to track batch contribution after merging
NBcell_0601@meta.data[["clusters"]] <- paste0('0601_', NBcell_0601@meta.data[['seurat_clusters']])
NBcell_0602@meta.data[["clusters"]] <- paste0('0602_', NBcell_0602@meta.data[['seurat_clusters']])

# --- MERGE BATCHES ---
NB.combined <- merge(
    NBcell_0601,
    y = c(NBcell_0602),
    add.cell.ids = c("0601", "0602"),
    project = "NBsinglecell"
)

# Add combined cluster metadata
NB.combined <- AddMetaData(NB.combined, as.factor(NB.combined@meta.data[["clusters"]]),
                           col.name = 'combined_clusters')


# --- RE-PROCESS GENE EXPRESSION DATA (RNA) ---
DefaultAssay(NB.combined) <- "RNA"
NB.combined <- SCTransform(NB.combined, verbose = FALSE)
NB.combined <- RunPCA(NB.combined, verbose = FALSE)

# --- RE-PROCESS DNA ACCESSIBILITY DATA (PEAKS) ---
DefaultAssay(NB.combined) <- "peaks"
NB.combined <- FindTopFeatures(NB.combined, min.cutoff = 5)
NB.combined <- RunTFIDF(NB.combined)
NB.combined <- RunSVD(NB.combined)


# --- JOINT MULTI-MODAL INTEGRATION (WNN) ---
# Build a joint neighbor graph using both assays
NB.combined <- FindMultiModalNeighbors(
    object = NB.combined,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:40), # Parameters matching preprocessing
    modality.weight.name = "RNA.weight",
    verbose = TRUE
)

# Find clusters on the integrated WNN graph
NB.combined <- FindClusters(NB.combined, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = FALSE)

# Build a joint UMAP visualization
NB.combined <- RunUMAP(
    object = NB.combined,
    nn.name = "weighted.nn",
    assay = "RNA", # Assay used for default color scaling, but embedding comes from WNN
    verbose = TRUE
)

# --- SAVE COMBINED OBJECT ---
save_path <- paste0(output_dir, "NBcell_combined.rds")
saveRDS(NB.combined, save_path)
