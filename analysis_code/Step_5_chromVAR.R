# Script 4: Motif Activity Analysis (ChromVAR) and Initial Visualization

library(Signac)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38) # Required for genome-based analysis
library(TFBSTools) # Required for motif matrix handling
library(JASPAR2020) # Required for motif database
library(ggplot2)

# --- USER INPUT & SETUP 🛠️ ---
# Define paths relative to the script location.
# USERS: Place the integrated Seurat object in the 'data/processed_integrated/' directory.
# USERS: Set your desired output directory.
input_dir <- "data/processed_integrated/"
output_dir <- "results/motif_analysis/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(1234)

# Load integrated Seurat object (NB_combined_processed.rds)
# NOTE: Replace the filename if your integrated object has a different name
NB_combined <- readRDS(paste0(input_dir, "NB_combined_signatures.rds")) 

# --- STEP 1: PEAK-TO-GENE LINKING ---
DefaultAssay(NB_combined) <- 'peaks'

# Compute GC content and sequence depth information for peaks
NB_combined <- RegionStats(NB_combined, genome = BSgenome.Hsapiens.UCSC.hg38)

# Link accessible peaks to genes based on correlation with gene expression (SCT assay)
NB_combined <- LinkPeaks(
  object = NB_combined, 
  peak.assay = "peaks", 
  expression.assay = "SCT"
)

# --- STEP 2: MOTIF IDENTIFICATION ---
# Get TF motif position frequency matrices (PFM) from JASPAR 2020
pfm <- getMatrixSet(
  x = JASPAR2020, 
  opts = list(collection = "CORE", tax_group = "vertebrates", all_version = FALSE)
)

# NOTE: Ensure the genome metadata is correctly set for motif scanning
NB_combined@assays[["peaks"]]@ranges@seqinfo@genome[] <- "hg38"

# Add motif information to the Seurat object
NB_combined <- AddMotifs(
  object = NB_combined, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  pfm = pfm
)


# --- STEP 3: MOTIF ACTIVITY CALCULATION (ChromVAR) ---
# Calculate per-cell motif accessibility deviation scores
NB_combined <- RunChromVAR(
  object = NB_combined, 
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# --- GLOBAL ACCESSIBILITY INDEX CALCULATION ---
# ChromVAR results are stored in a new 'chromvar' assay (deviations and counts)
DefaultAssay(NB_combined) <- 'chromvar'

# Calculate Total Motif-Associated Fragment Count (nCount_chromvar)
NB_combined$nCount_chromvar <- Matrix::colSums(NB_combined@assays[["chromvar"]]$counts)

# Calculate Total Accessible Motifs Count (nFeature_chromvar)
# The 'counts' slot is where the raw counts are stored. Counting features with >0 fragments.
NB_combined$nFeature_chromvar <- Matrix::colSums(NB_combined@assays[["chromvar"]]$counts > 0)


# --- INITIAL VISUALIZATION OF CHROMVAR INDICES ---

# 1. Visualize Total Fragment Count (nCount_chromvar) on UMAP
png(paste0(output_dir, "chromvar_nCount_umap.png"), width = 7, height = 6, units = 'in', res = 300)
FeaturePlot(NB_combined, features = c('nCount_chromvar'), order = TRUE)
dev.off()

# 2. Visualize Total Fragment Count (nCount_chromvar) by cell line (VlnPlot)
png(paste0(output_dir, "chromvar_nCount_vlnplot.png"), width = 10, height = 6, units = 'in', res = 300)
VlnPlot(NB_combined, features = c('nCount_chromvar'), pt.size = 0, group.by = 'cell_line4') + NoLegend()
dev.off()

# 3. Scatter Plot: ChromVAR Count vs. MES/ADRN Signature Scores (Using expanded PMID:286 list)
# Note: Assuming 'signature_1PMID286_all' and 'signature_1PMID378' scores were calculated in a previous step
DefaultAssay(NB_combined) <- 'RNA' # Switch back for FeatureScatter axis labels if needed

# Scatter: ChromVAR Count vs. Expanded ADRN Score (PMID:286)
png(paste0(output_dir, "scatter_chromvar_ADRN286all.png"), width = 7, height = 6, units = 'in', res = 300)
FeatureScatter(NB_combined, feature1 = 'nCount_chromvar', feature2 = "signature_1PMID286_ADRN_all1", group.by = 'cell_line5') 
dev.off()

# Scatter: ChromVAR Count vs. Expanded MES Score (PMID:286)
png(paste0(output_dir, "scatter_chromvar_MES286all.png"), width = 7, height = 6, units = 'in', res = 300)
FeatureScatter(NB_combined, feature1 = 'nCount_chromvar', feature2 = "signature_1PMID286_MES_all1", group.by = 'cell_line5')
dev.off()


# --- SAVE OBJECT ---
save_path_rds <- paste0(output_dir, 'NB_combined_motif_ready.rds')
saveRDS(NB_combined, save_path_rds)

print(paste0("ChromVAR analysis complete. Object saved to: ", save_path_rds))
