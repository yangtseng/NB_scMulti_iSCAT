# Script 1: Single-Batch Multiome Preprocessing (Batch 0601)
# Reference: https://satijalab.org/signac/articles/pbmc_multiomic.html

library(Signac)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(EnsDb.Hsapiens.v86)
library(ggplot2)

# --- USER INPUT & SETUP 🛠️ ---
# Define paths relative to the script location.
# USERS: Place 10x output in 'data/batch_0601/outs/' and Souporcell output in 'data/souporcell_0601/'.
# USERS: Set your desired output directory.
data_dir_base <- "data/"
output_dir <- "results/preprocessing/" # Public output directory
batch_id <- "0601"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(1234) # Set random seed for reproducibility

# --- LOAD 10X MULTIOME DATA ---
# Construct paths to 10x files
data_path_10x <- paste0(data_dir_base, "batch_", batch_id, "/outs/")
counts_path <- paste0(data_path_10x, "filtered_feature_bc_matrix/")
fragpath <- paste0(data_path_10x, 'atac_fragments.tsv.gz')

# Load 10x multiome data (RNA and ATAC)
counts <- Read10X(data.dir = counts_path)

# Prepare genomic annotation for ATAC assay
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# Create Seurat object with RNA data
NBcell_0601 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  project = "NB_singlecell"
)

# Create ATAC assay and add it to the object
NBcell_0601[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# --- DOUBLET FILTERING (SOUPORCELL) ---
# Load Souporcell result (assumes file is in data_dir_base/souporcell_0601/clusters.tsv)
souporcell_path <- paste0(data_dir_base, 'souporcell_', batch_id, '/clusters.tsv')
cluster0601 <- read.table(souporcell_path, header = TRUE)

# Assign cell status metadata (singlet/doublet/unassigned)
NBcell_0601 <- AddMetaData(NBcell_0601, cluster0601$status, col.name = 'status')

# Filter cells: keep only singlets
NBcell_0601 <- subset(x = NBcell_0601, subset = status == 'singlet')

# Assign Souporcell cluster result as metadata
NBcell_0601 <- AddMetaData(NBcell_0601,
                            cluster0601$assignment[cluster0601$status == 'singlet'],
                            col.name = 'cluster')

# --- QUALITY CONTROL (QC) ---
DefaultAssay(NBcell_0601) <- "ATAC"

# Compute ATAC quality metrics
NBcell_0601 <- NucleosomeSignal(NBcell_0601)
NBcell_0601 <- TSSEnrichment(NBcell_0601)

# Plot QC metrics before filtering
png(filename = paste0(output_dir, batch_id, "_nucleosome_tss_prefilter.png"), width = 16, height = 9, units = "in", res = 300)
VlnPlot(object = NBcell_0601, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0)
dev.off()

# Filter out low quality cells based on established thresholds
NBcell_0601 <- subset(
  x = NBcell_0601,
  subset = nCount_RNA > 1000 &
           nCount_RNA < 25000 &
           nCount_ATAC > 1000 &
           nCount_ATAC < 100000 &
           nucleosome_signal < 1.1 &
           TSS.enrichment > 5
)

# Plot QC metrics after filtering
png(filename = paste0(output_dir, batch_id, "_nucleosome_tss_filtered.png"), width = 16, height = 9, units = "in", res = 300)
VlnPlot(object = NBcell_0601, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0)
dev.off()

# --- PEAK CALLING AND REFINEMENT ---
# NOTE: macs2.path should be NULL if MACS2 is in the system PATH.
peaks <- CallPeaks(NBcell_0601, macs2.path = NULL) 

# Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
# NOTE: blacklist_hg38 must be loaded/defined previously (e.g., from Signac's helper functions)
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38, invert = TRUE)

# Quantify counts in each new MACS2 peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(NBcell_0601),
  features = peaks,
  cells = colnames(NBcell_0601)
)

# Create a new assay using the MACS2 peak set
NBcell_0601[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# --- SAVE PROCESSED OBJECT ---
save_path <- paste0(output_dir, "NBcell_", batch_id, ".h5Seurat")
SaveH5Seurat(NBcell_0601, filename = save_path, overwrite = TRUE)
