library(ggplot2)
library(Seurat)
library(ggpubr)
library(ggsignif)

### Fig. 5A and related supplementary figures
png('fig5a.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("nCount_chromvar"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("ChromVAR nCount")
dev.off()

averages <- NB_combined_250617@meta.data %>%
  group_by(cell_line4) %>%
  summarise(mean_chromvar = mean(nCount_chromvar, na.rm = TRUE))
### Fig. 5B
# A tibble: 8 × 2
# cell_line4    mean_chromvar
# <fct>                 <dbl>
#   1 SHSY5Y_ctrl           622.
# 2 SHSY5Y_treat          231.
# 3 SKNBE2C_ctrl          -23.7
# 4 SKNBE2C_treat         323.
# 5 SKNDZ_ctrl           -199.
# 6 SKNDZ_treat           -78.9
# 7 SKNSH_ctrl            -45.7
# 8 SKNSH_treat          -166.
### Fig. 5B
png('fig5b.png',
    width = 10,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nCount_chromvar"), pt.size = 0, group.by = c("cell_line_4"),
        cols = c('#D38681','#F6DB88','#8CADCF','#85B69E','#C44943','#F3C064','#4874A7','#4B7C59')) +
  geom_signif(comparisons = list(c("SHSY5Y_ctrl","SHSY5Y_treat"), c("SKNBE2C_ctrl","SKNBE2C_treat"),
                                 c("SKNDZ_ctrl",'SKNDZ_treat'),c('SKNSH_ctrl','SKNSH_treat')), 
              y_position = c(2500, 1800,1500,1500), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(-800,3200) +
  theme_classic() + ylab("ChromVAR\nnCount") + xlab("") + ggtitle("") + NoLegend() +
  annotate("text", x = 1.5, y = 3150, label = "Δ = -390.4", size = 5) +
  annotate("text", x = 3.5, y = 2450, label = "Δ = 347.1", size = 5) +
  annotate("text", x = 5.5, y = 2150, label = "Δ = 120.5", size = 5) +
  annotate("text", x = 7.5, y = 2150, label = "Δ = -120.7", size = 5) +
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

png('fig5c.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("nFeature_chromvar"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("ChromVAR nFeature")
dev.off()

averages <- NB_combined_250617@meta.data %>%
  group_by(cell_line4) %>%
  summarise(mean_chromvar = mean(nFeature_chromvar, na.rm = TRUE))
# # A tibble: 8 × 2
# cell_line4    mean_chromvar
# <fct>                 <dbl>
#   1 SHSY5Y_ctrl            469.
# 2 SHSY5Y_treat           402.
# 3 SKNBE2C_ctrl           342.
# 4 SKNBE2C_treat          437.
# 5 SKNDZ_ctrl             324.
# 6 SKNDZ_treat            356.
# 7 SKNSH_ctrl             364.
# 8 SKNSH_treat            316.
png('fig5d.png',
    width = 10,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nFeature_chromvar"), pt.size = 0, group.by = c("cell_line_4"),
        cols = c('#D38681','#F6DB88','#8CADCF','#85B69E','#C44943','#F3C064','#4874A7','#4B7C59')) +
  geom_signif(comparisons = list(c("SHSY5Y_ctrl","SHSY5Y_treat"), c("SKNBE2C_ctrl","SKNBE2C_treat"),
                                 c("SKNDZ_ctrl",'SKNDZ_treat'),c('SKNSH_ctrl','SKNSH_treat')), 
              y_position = c(700, 700,700,700), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(100,900) +
  theme_classic() + ylab("ChromVAR\nnFeature") + xlab("") + ggtitle("") + NoLegend() +
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

### Fig. 5E
png('fig5e_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SHSY5Y_ctrl") |NB_combined_250617$cell_line4 %in% c("SHSY5Y_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("Cell line: SH-SY5Y (Corr. = 0.67)") + NoLegend() +
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

png('fig5e_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNBE2C_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNBE2C_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("Cell line: SK-N-BE(2)C (Corr. = 0.3)") + NoLegend() +
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

png('fig5e_3.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNDZ_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNDZ_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("Cell line: SK-N-DZ (Corr. = 0.1)") + NoLegend() +
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

png('fig5e_4.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNSH_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNSH_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("Cell line: SK-N-SH (Corr. = 0.05)") + NoLegend() +
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

### Fig. 5F
png('fig5f_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SHSY5Y_ctrl") |NB_combined_250617$cell_line4 %in% c("SHSY5Y_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR Feature") + ggtitle("Cell line: SH-SY5Y (Corr. = 0.51)") + NoLegend() +
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

png('fig5f_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNBE2C_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNBE2C_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nFeature") + ggtitle("Cell line: SK-N-BE(2)C (Corr. = 0.28)") + NoLegend() +
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

png('fig5f_3.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNDZ_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNDZ_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nFeature") + ggtitle("Cell line: SK-N-DZ (Corr. = 0.07)") + NoLegend() +
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

png('fig5f_4.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line4 %in% c("SKNSH_ctrl") |NB_combined_250617$cell_line4 %in% c("SKNSH_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nFeature") + ggtitle("Cell line: SK-N-SH (Corr. = 0.04)") + NoLegend() +
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
