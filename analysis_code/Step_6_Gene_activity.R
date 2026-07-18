library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)
### Load data
NB_combined_250617 <- readRDS("~/R/OtherLabMember/kjh/TYT/NB_combined_250617.rds")
genename <- as.data.frame(NB_combined_250617@assays[["RNA"]]@data@Dimnames[[1]])
batch <- as.character(str_extract(NB_combined_250617@meta.data[["clusters"]], "[^_]+"))
treatment <- as.character(sub(".*_", "",NB_combined_250617@meta.data[["cell_line4"]]))

NB_combined_250617[['batch']] <- as.factor(batch)
NB_combined_250617[['treatment']] <- as.factor(treatment)
library(Seurat)
library(dplyr)

# 1. Set the active identity to the treatment column
# Replace 'treatment' with the actual name of your metadata column
Idents(NB_combined_250617) <- "treatment"

# 2. Find markers upregulated in 'ATRA' compared to 'Control'
# ident.1 = the group of interest (Treated)
# ident.2 = the reference group (Control)
atra_vs_control_markers <- FindMarkers(
  object = NB_combined_250617,
  ident.1 = "treat",
  ident.2 = "ctrl",
  only.pos = TRUE,          # Only return genes upregulated in ATRA
  min.pct = 0.25,           # Gene must be detected in at least 25% of cells in either group
  logfc.threshold = 0.25,   # Minimum log-fold change threshold
  test.use = "wilcox"       # Default Wilcoxon Rank Sum test
)
rownames(atra_vs_control_markers)
# 3. Filter for highly expressed and significant genes
# We filter for a strict p-value and a higher log-fold change
highly_expressed_genes <- atra_vs_control_markers %>%
  filter(p_val_adj < 0.01) %>%
  arrange(desc(avg_log2FC))

# View the top genes
head(highly_expressed_genes, n = 20)

# Optional: Save the results to a CSV file
# write.csv(highly_expressed_genes, "ATRA_upregulated_genes.csv")


DefaultAssay(NB_combined_250617) <- 'RNA'


DefaultAssay(NB_combined_260514) <- 'SCT'
VlnPlot(NB_combined_260514, features = c('NTRK2'), slot = 'scale.data', pt.size = 0, group.by = 'cell_line4')
DefaultAssay(NB_combined_260514) <- 'RNA'
VlnPlot(NB_combined_260514, features = c('NTRK2'), pt.size = 0, group.by = 'cell_line4')

DefaultAssay(NB_combined_250617) <- "SCT"
png('fig3d.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c('DHRS3'), slot = 'scale.data', pt.size = .5, order = T, reduction = 'umap',cols = c('grey80',"#D4524E")) + 
  theme_classic() + xlab('UMAP-1') + ylab("UMAP-2") + ggtitle("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

### Gene activity
library(Signac)
library(Seurat)

# 2. Calculate gene activity
# This counts fragments in the gene body + 2kb upstream (default)
DefaultAssay(NB_combined_250617) <- "peaks"
gene_activities <- GeneActivity(NB_combined_250617)

# 3. Add to the Seurat object as a new assay
NB_combined_250617[['gene_activity']] <- CreateAssayObject(counts = gene_activities)

# 4. Normalize the activity scores
# This makes the scores comparable across cells
DefaultAssay(NB_combined_250617) <- "gene_activity"
NB_combined_250617 <- NormalizeData(
  object = NB_combined_250617,
  assay = "gene_activity",
  normalization.method = "LogNormalize")

png('fig3e.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c('DHRS3'), pt.size = .5, order = F, 
            reduction = 'umap',cols = c('grey80',"#D4524E")) + 
  theme_classic() + xlab('UMAP-1') + ylab("UMAP-2") + ggtitle("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()


library(ggplot2)
library(Seurat)
library(ggpubr)
library(ggsignif)
setwd('~/new_fig5')
### new_fig. 5A and related supplementary new_figures
png('new_fig5a.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("nCount_gene_activity"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("Gene activity nCount")
dev.off()

averages <- NB_combined_250617@meta.data %>%
  group_by(cell_line4) %>%
  summarise(mean_gene_activity = mean(nCount_gene_activity, na.rm = TRUE))
### new_fig. 5B
# A tibble: 8 × 2
# cell_line4    mean_gene_activity
# <fct>                      <dbl>
#   1 SHSY5Y_ctrl               14373.
# 2 SHSY5Y_treat              12028.
# 3 SKNBE2C_ctrl              11290.
# 4 SKNBE2C_treat             13897.
# 5 SKNDZ_ctrl                15960.
# 6 SKNDZ_treat               15953.
# 7 SKNSH_ctrl                13634.
# 8 SKNSH_treat                9820.

### new_fig. 5B
png('new_fig5b.png',
    width = 10,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nCount_gene_activity"), pt.size = 0, group.by = c("cell_line4"),
        cols = c('#D38681','#F6DB88','#8CADCF','#85B69E','#C44943','#F3C064','#4874A7','#4B7C59')) +
  geom_signif(comparisons = list(c("SHSY5Y_ctrl","SHSY5Y_treat"), c("SKNBE2C_ctrl","SKNBE2C_treat"),
                                 c("SKNDZ_ctrl",'SKNDZ_treat'),c('SKNSH_ctrl','SKNSH_treat')), 
              y_position = c(68000,68000,68000,68000), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(0,81000) +
  theme_classic() + ylab("Gene activity\nnCount") + xlab("") + ggtitle("") + NoLegend() +
  annotate("text", x = 1.5, y = 81000, label = "Δ = -2345", size = 5) +
  annotate("text", x = 3.5, y = 81000, label = "Δ = 2607", size = 5) +
  annotate("text", x = 5.5, y = 81000, label = "Δ = -7", size = 5) +
  annotate("text", x = 7.5, y = 81000, label = "Δ = -3794", size = 5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial')) 
dev.off()

png('new_fig5c.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("nFeature_gene_activity"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("Gene activity nFeature")
dev.off()

averages <- NB_combined_250617@meta.data %>%
  group_by(cell_line4) %>%
  summarise(mean_gene_activity = mean(nFeature_gene_activity, na.rm = TRUE))
# A tibble: 8 × 2
# cell_line4    mean_gene_activity
# <fct>                      <dbl>
#   1 SHSY5Y_ctrl                7145.
# 2 SHSY5Y_treat               6509.
# 3 SKNBE2C_ctrl               6313.
# 4 SKNBE2C_treat              7157.
# 5 SKNDZ_ctrl                 7838.
# 6 SKNDZ_treat                7856.
# 7 SKNSH_ctrl                 7281.
# 8 SKNSH_treat                5964.
png('new_fig5d.png',
    width = 10,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nFeature_gene_activity"), pt.size = 0, group.by = c("cell_line4"),
        cols = c('#D38681','#F6DB88','#8CADCF','#85B69E','#C44943','#F3C064','#4874A7','#4B7C59')) +
  geom_signif(comparisons = list(c("SHSY5Y_ctrl","SHSY5Y_treat"), c("SKNBE2C_ctrl","SKNBE2C_treat"),
                                 c("SKNDZ_ctrl",'SKNDZ_treat'),c('SKNSH_ctrl','SKNSH_treat')), 
              y_position = c(15000,15000,15000,15000), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(0,18000) +
  theme_classic() + ylab("Gene activity\nnFeature") + xlab("") + ggtitle("") + NoLegend() +
  annotate("text", x = 1.5, y = 18000, label = "Δ = -636", size = 5) +
  annotate("text", x = 3.5, y = 18000, label = "Δ = 681", size = 5) +
  annotate("text", x = 5.5, y = 18000, label = "Δ = 18", size = 5) +
  annotate("text", x = 7.5, y = 18000, label = "Δ = -1317", size = 5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial')) 
dev.off()

### new_fig. 5E
png('new_fig5e_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SHSY5Y_ctrl") |NB_combined_250617$cell_line4 %in% c("SHSY5Y_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nCount") + ylab("MES signature gene score") + ggtitle("Cell line: SH-SY5Y (Corr. = 0.08)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5e_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNBE2C_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNBE2C_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nCount") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-BE(2)C (Corr. = 0.13)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5e_3.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNDZ_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNDZ_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nCount") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-DZ (Corr. = -0.02)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5e_4.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNSH_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNSH_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nCount") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-SH (Corr. = 0.01)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

### new_fig. 5F
png('new_fig5f_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SHSY5Y_ctrl") |NB_combined_250617$cell_line4 %in% c("SHSY5Y_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nFeature") + ylab("MES signature gene score") + ggtitle("Cell line: SH-SY5Y (Corr. = 0.05)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5f_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNBE2C_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNBE2C_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nFeature") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-BE(2)C (Corr. = 0.13)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5f_3.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNDZ_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNDZ_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nFeature") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-DZ (Corr. = -0.04)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

png('new_fig5f_4.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_gene_activity",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNSH_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNSH_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("Gene activity nFeature") + ylab("MES signature gene score") + ggtitle("Cell line: SK-N-SH (Corr. = -0.01)") + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()





#############################################

# MES score: vector of length 14284
mes_score <- NB_combined_250617$signature_1PMID286_MES_all  # UPDATE column name

# Gene Activity matrix: 19607 genes x 14284 cells
ga_mat <- GetAssayData(NB_combined_250617, assay = "gene_activity", slot = "data")  # UPDATE assay name
chrom_mat <- GetAssayData(NB_combined_250617, assay = 'chromvar', slot = 'data')
# Spearman correlation: each gene's activity vs MES score
cor_results <- apply(ga_mat, 1, function(x) {
  cor(as.numeric(x), mes_score, method = "spearman", use = "complete.obs")
})

# Rank high to low
cor_df <- data.frame(
  gene = names(cor_results),
  rho = cor_results
) %>% arrange(desc(rho))

top20_cor_df <- as.data.frame(head(cor_df, 20))

png('fig6d.png',
    width = 7,
    height = 3.5,
    res = 300,
    unit = 'in')
ggplot(top20_cor_df, aes(x =  reorder(gene, -rho), y = rho)) + 
  geom_col(fill = "#A8D0C7", width = 0.6) +
  theme_classic() + ggtitle("Correlation between MES score & gene activity") +
  ylab('SCC') + 
  xlab('Gene') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 15, color = 'black', hjust = 1, angle = 45), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = 0), 
        axis.title.y = element_text(size = 20, vjust = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(0, 0.35), expand = c(0,1))
dev.off()
  

cor_results_2 <- apply(chrom_mat, 1, function(x) {
  cor(as.numeric(x), mes_score, method = "spearman", use = "complete.obs")
})

# Rank high to low
cor_df_2 <- data.frame(
  gene = names(cor_results_2),
  rho = cor_results_2
) %>% arrange(desc(rho))

head(cor_df_2, 20)

# Get motif-to-TF mapping from your Seurat object
library(JASPAR2020)
library(TFBSTools)

pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", species = "Homo sapiens"))
motif_names <- sapply(pfm, function(x) x@name)
motif_ids <- sapply(pfm, function(x) x@ID)
motif_map <- data.frame(id = motif_ids, tf_name = motif_names)

# Map your top hits
top_motifs <- head(cor_df_2, 20)$gene
motif_map[motif_map$id %in% top_motifs, ]

top20_cor_df_2 <- as.data.frame(head(cor_df_2, 20))

png('fig6e.png',
    width = 7,
    height = 3.5,
    res = 300,
    unit = 'in')
ggplot(top20_cor_df_2, aes(x =  reorder(tf, -rho), y = rho)) + 
  geom_col(fill = "#A8D0C7", width = 0.6) +
  theme_classic() + ggtitle("Correlation between MES score & chromVAR") +
  ylab('SCC') + 
  xlab('TF motif') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 15, color = 'black', hjust = 1, angle = 45), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = 0), 
        axis.title.y = element_text(size = 20, vjust = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(0, 0.55), expand = c(0,1))
dev.off()
# MES score variance per cell line per treatment
cat("=== MES Score Variance ===\n")
mes_var <- tapply(NB_combined_250617$signature_1PMID286_MES_all,
                  paste(NB_combined_250617$cell_line4, NB_combined_250617$treatment, sep = "_"),
                  var)
print(sort(mes_var, decreasing = TRUE))

# Also compute coefficient of variation (CV) for fair comparison across different means
mes_cv <- tapply(NB_combined_250617$signature_1PMID286_MES_all,
                 paste(NB_combined_250617$cell_line4, NB_combined_250617$treatment, sep = "_"),
                 function(x) sd(x) / abs(mean(x)))
cat("\n=== MES Score CV ===\n")
print(sort(mes_cv, decreasing = TRUE))

# AP-1 motif variance — pick one representative: FOS::JUN (MA0099.3)
cat("\n=== AP-1 (FOS::JUN) Motif Variance ===\n")
ap1_scores <- GetAssayData(NB_combined_250617, assay = "chromvar")["MA0099.3", ]
ap1_var <- tapply(ap1_scores,
                  paste(NB_combined_250617$cell_line4, NB_combined_250617$treatment, sep = "_"),
                  var)
print(sort(ap1_var, decreasing = TRUE))

ap1_cv <- tapply(ap1_scores,
                 paste(NB_combined_250617$cell_line4, NB_combined_250617$treatment, sep = "_"),
                 function(x) sd(x) / abs(mean(x)))
cat("\n=== AP-1 (FOS::JUN) Motif CV ===\n")
print(sort(ap1_cv, decreasing = TRUE))

