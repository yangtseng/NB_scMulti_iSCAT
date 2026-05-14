library(Signac)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(pheatmap) # Added for correlation visualization
library(dplyr)    # For data manipulation
library(cowplot)  # For plotting

input_dir <- "results/integration/"
# --- CELL LINE ANNOTATION VIA BULK CORRELATION ---
# NOTE: Bulk data file is assumed to be accessible via the GEO accession used in the original code.
# The user needs to download this file and place it in the annotation directory.
bulk_file_path <- paste0(annotation_dir, "GSE89413_2016-10-30-NBL-cell-line-STAR-fpkm.txt.gz")
NB_cellline <- read.table(gzfile(bulk_file_path), sep = "\t", header = TRUE)

# Prepare bulk data
rownames(NB_cellline) <- NB_cellline$GeneID
NB_cellline <- NB_cellline %>% 
    dplyr::select(all_of(c("SKNDZ", "SKNBE2C", "SHSY5Y", "SKNSH"))) # Select relevant cell lines
NB_cellline <- log2(NB_cellline + 1) # Log2 transformation

# Calculate average expression for each scRNA cluster (using SCT assay)
NB.combined <- readRDS(paste0(input_dir, "NBcell_combined.rds"))
DefaultAssay(NB.combined) <- "SCT"
cellLine.avg <- AverageExpression(NB.combined, assays = 'SCT', slot = 'data')[['SCT']]

# Find overlapping genes and order them
corroverlapgenes <- intersect(rownames(NB_cellline), rownames(cellLine.avg))
cellLine.avg <- cellLine.avg[corroverlapgenes, ]
NB_cellline <- NB_cellline[corroverlapgenes, ]
cellLine.avg <- cellLine.avg[order(rownames(cellLine.avg)),]
NB_cellline <- NB_cellline[order(rownames(NB_cellline)),]

# Calculate Spearman correlation
corr <- cor(x = NB_cellline, y = cellLine.avg, method = c('spearman'))

# Visualize correlation results
p <- pheatmap(scale(corr), colorRampPalette(c('#8E9AB1','white','#D99388'))(100),
              cellwidth = 35, cellheight = 20,
              angle_col = 45, fontsize = 20, border_color = "white") 
p$gtable$grobs[[3]]$gp <- gpar(lwd = 4)
p

png(
  filename  = paste("fig3b.png", sep=""),
  width     = 13,
  height    = 7,
  units     = "in",
  res       = 300,
)

grid.newpage()
grid.draw(p$gtable)

dev.off()



# NOTE: The manual assignment of cell line predictions based on the heatmap
# requires manual review of the heatmap. The resulting vector is provided here
# as an example for downstream code execution, but users MUST verify the order
# against their own heatmap output (the columns of the 'corr' matrix).

# ********** MANUAL STEP START **********
# USERS: Replace this with the order derived from your Pheatmap output.
# This assumes the combined cluster identity is the active identity/levels.
# Example: NB.combined@active.ident (or the levels of the 'seurat_clusters' used for correlation)
# The length of this vector MUST match the number of clusters found (levels(NB.combined@active.ident)).
NB_Prediction_new <- c(
    "SKNDZ_treat","SKNBE2C_treat","SKNSH_ctrl","SHSY5Y_ctrl","SHSY5Y_ctrl","SKNDZ_treat",
    "SKNDZ_treat","SKNSH_ctrl","SKNSH_ctrl","SKNBE2C_treat","SKNDZ_ctrl","SKNSH_treat",
    "SHSY5Y_treat","SKNDZ_ctrl","SKNSH_treat","SKNSH_treat","SKNDZ_ctrl","SKNBE2C_ctrl",
    "SHSY5Y_treat","SKNSH_treat","SKNBE2C_ctrl" # This example vector assumes 21 clusters
)
names(NB_Prediction_new) <- levels(NB.combined) # Set names to match cluster IDs
# ********** MANUAL STEP END **********

# Assign predicted cell line to the Seurat object
NB.combined <- RenameIdents(NB.combined, NB_Prediction_new)

# --- SAVE COMBINED OBJECT ---
save_path_h5 <- paste0(output_dir, "NBcell_Combined.h5Seurat")
SaveH5Seurat(NB.combined, filename = save_path_h5, overwrite = TRUE)

save_path_rds <- paste0(output_dir, "NB_combined_new.rds")
saveRDS(NB.combined, save_path_rds)

print(paste0("Integration and Cell Line Annotation complete. Objects saved to: ", output_dir))
