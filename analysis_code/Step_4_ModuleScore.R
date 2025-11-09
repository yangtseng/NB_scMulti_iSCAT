# Script 4: MES/ADRN Signature Score Calculation (PMID:286 Focus)

library(Seurat)
library(Signac)
library(SeuratData)
library(ggplot2)
library(UCell) # Used for UCell score calculation

# --- USER INPUT & SETUP 🛠️ ---
# Define paths relative to the script location.
# USERS: Processed Seurat object must be in the current directory or specified path.
# USERS: Gene signature files must be placed in the 'data/signatures/' directory.
signatures_dir <- "data/signatures/"
output_dir <- "results/mes_adrn/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
set.seed(1234)

# Load integrated Seurat object (saved from Script 2)
# NOTE: Replace 'NB_combined_processed.rds' with your actual file name
NB_combined <- readRDS("results/integration/NB_combined_new.rds") 


# --- PREPARE CELL LINE METADATA ---
# Convert current active identity to a factor with a fixed order (cell_line4)
NB_combined$cell_line4 <- factor(NB_combined@active.ident, 
                                 levels = c("SHSY5Y_ctrl","SHSY5Y_treat","SKNBE2C_ctrl","SKNBE2C_treat",
                                            "SKNDZ_ctrl",'SKNDZ_treat','SKNSH_ctrl','SKNSH_treat'))

# Create complex cell line metadata (cell_line5) separating SHSY5Y subtypes
cell_line <- as.data.frame(NB_combined$cell_line4)
cell_line$current_ident <- as.character(cell_line$`NB_combined_250530$cell_line4`) # Renamed variable for clarity
cell_line$clusters <- NB_combined$clusters

cell_line_list5 <- ifelse(cell_line$current_ident == 'SHSY5Y_treat',
                          # NOTE: cluster IDs ('0602_9', '0601_3') are hardcoded based on previous clustering
                          ifelse(cell_line$clusters == '0602_9', 'SHSY5Y_S_treat', 'SHSY5Y_N_treat'),
                          ifelse(cell_line$current_ident == 'SHSY5Y_ctrl',
                                 ifelse(cell_line$clusters == '0601_3', 'SHSY5Y_S_ctrl', 'SHSY5Y_N_ctrl'),
                                 cell_line$current_ident))

NB_combined$cell_line5 <- cell_line_list5
DimPlot(NB_combined, group.by = 'cell_line5') # Visualization step remains

# --- LOAD PMID:286 GENE SIGNATURES ---
# ADRN/MES signatures from PMID:28650485 (Nature Genetics, 2017)
PMID286 <- read.table(paste0(signatures_dir, 'ADRN_MES_sig_PMID28650485.txt'), header = TRUE)
PMID286_all <- read.table(paste0(signatures_dir, 'ADRN_MES_sig_PMID28650485_all.txt'), header = TRUE)

# Define feature lists
ADRN_286 <- PMID286$gene[PMID286$term %in% "ADRN"]      # Original (Small) List
MES_286 <- PMID286$gene[PMID286$term %in% "MES"]        # Original (Small) List
ADRN_286_all <- PMID286_all$gene[PMID286_all$term %in% "ADRN"] # Expanded List
MES_286_all <- PMID286_all$gene[PMID286_all$term %in% "MES"]   # Expanded List

# --- CALCULATE MODULE SCORES (Seurat::AddModuleScore) ---
DefaultAssay(NB_combined) <- 'RNA'
# 1. Expanded Gene Set ("all")
NB_combined <- AddModuleScore(NB_combined, features = list(ADRN_286_all), assay = 'RNA', name = 'ADRN_PMID286_all', seed = 1234)
NB_combined <- AddModuleScore(NB_combined, features = list(MES_286_all), assay = 'RNA', name = 'MES_PMID286_all', seed = 1234)
# 2. Original Gene Set
NB_combined <- AddModuleScore(NB_combined, features = list(ADRN_286), assay = 'RNA', name = 'ADRN_PMID286', seed = 1234)
NB_combined <- AddModuleScore(NB_combined, features = list(MES_286), assay = 'RNA', name = 'MES_PMID286', seed = 1234)

# --- CALCULATE MODULE SCORES (UCell) ---
# UCell is often preferred for signature scoring
DefaultAssay(NB_combined) <- 'RNA'
# 1. Expanded Gene Set ("all")
gene.sets.286.ADRN.all <- list(c(ADRN_286_all)) 
gene.sets.286.MES.all <- list(c(MES_286_all))   
NB_combined <- AddModuleScore_UCell(NB_combined, features = gene.sets.286.ADRN.all, name="PMID286_ADRN_all")
NB_combined <- AddModuleScore_UCell(NB_combined, features = gene.sets.286.MES.all, name="PMID286_MES_all")
# 2. Original Gene Set
gene.sets.286.ADRN <- list(c(ADRN_286)) 
gene.sets.286.MES <- list(c(MES_286)) 
NB_combined <- AddModuleScore_UCell(NB_combined, features = gene.sets.286.ADRN, name="PMID286_ADRN")
NB_combined <- AddModuleScore_UCell(NB_combined, features = gene.sets.286.MES, name="PMID286_MES")


# --- CORRELATION ANALYSIS (PMID:286 vs PMID:286) ---
# Correlate the two lists (original vs. expanded) for verification
print("Correlation between original and expanded PMID:286 ADRN scores:")
cor.test(NB_combined@meta.data[["ADRN_PMID286_all1"]], NB_combined@meta.data[["ADRN_PMID2861"]])

print("Correlation between original and expanded PMID:286 MES scores:")
cor.test(NB_combined@meta.data[["MES_PMID286_all1"]], NB_combined@meta.data[["MES_PMID2861"]])

print("Correlation between ADRN and MES (Expanded) scores:")
cor.test(NB_combined@meta.data[["ADRN_PMID286_all1"]], NB_combined@meta.data[["MES_PMID286_all1"]])


# --- SAVE OBJECT ---
save_path_rds <- paste0(output_dir, 'NB_combined_signatures.rds')
saveRDS(NB_combined, save_path_rds)
print(paste0("Signature score calculation complete. Object saved to: ", save_path_rds))
