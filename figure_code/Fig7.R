### Fig. 7
library(Seurat)
library(ggplot2)
library(Seurat)
library(ggpubr)
library(ggsignif)

setwd("~/R/iSCAT/")
NB_combined_250617 <- readRDS("~/R/OtherLabMember/kjh/TYT/NB_combined_250617.rds")

### Fig 7A
DefaultAssay(NB_combined_250617) <- "SCT"
png('fig7a1.png',
    width = 7,
    height = 11,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, pt.size = .2, features = c("CD44", "FN1"), 
            ncol = 1, cols = c('grey90',"#D4524E"), reduction = 'umap') &
  theme_classic() & xlab('UMAP-1') & ylab("UMAP-2") &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, family = 'Arial'),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()

### Fig 7B
png('fig7b.png',
    width = 8,
    height = 5,
    res = 300,
    unit = 'in')
DimPlot(NB_combined_250617, pt.size = .5, group.by = 'cell_line5', order = T,
        cols = c('#D38681','#C44943','#EAB389','#B86029','grey90','grey90','grey90','grey90','grey90','grey90')) + 
  theme_classic() + xlab('UMAP-1') + ylab("UMAP-2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial'))
dev.off()


SHSY5Y_cellname <- NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617@active.ident %in% c("SHSY5Y_ctrl","SHSY5Y_treat")]

# c('#D38681','#EAB389','#F6DB88','#8CADCF','#85B69E','#C44943','#B86029','#F3C064','#4874A7','#4B7C59')
png('fig7c.png',
    width = 4,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("signature_1PMID286_MES_all"), pt.size = 0, group.by = c("cell_line5"),idents = c("SHSY5Y_ctrl","SHSY5Y_treat"),
        cols = c('#D38681','#EAB389','#C44943','#B86029','#F6DB88','#8CADCF','#85B69E','#F3C064','#4874A7','#4B7C59')) +
  ylim(0.17,0.38) +
  geom_signif(comparisons = list(c("SHSY5Y_N_ctrl","SHSY5Y_N_treat"), c("SHSY5Y_S_ctrl","SHSY5Y_S_treat")),
              y_position = c(0.3, 0.35), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  theme_classic() + ylab("MES signature \n gene score") + xlab("") + ggtitle("") + NoLegend() +
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

CDF_df_generate <- function(cell1, cell2,type){
  df <- as.data.frame(NB_combined_250617@meta.data[[type]][NB_combined_250617$cell_line5 %in% cell1])
  df$cell <-NB_combined_250617$cell_line4[NB_combined_250617$cell_line5 %in% cell1]
  colnames(df) <- c("CDF","Cell")
  df2 <- as.data.frame(NB_combined_250617@meta.data[[type]][NB_combined_250617$cell_line5 %in% cell2])
  df2$cell <-NB_combined_250617$cell_line4[NB_combined_250617$cell_line5 %in% cell2]
  colnames(df2) <- c("CDF","Cell")
  df <- rbind(df, df2)
  
  return(df)
}

df <- CDF_df_generate('SHSY5Y_S_treat','SHSY5Y_S_ctrl', "signature_1PMID286_MES_all")

png('fig7e.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
ggplot(df, aes(CDF, colour = Cell)) + ylab("CDF") + xlab("MES signature gene score") +
  stat_ecdf(size = 1.5) + ggtitle("SH-SY5Y-S") + theme_classic() + NoLegend() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title = element_text(size = 20, family = "Arial"),
        axis.title.x = element_text(size = 20, family = 'Arial', vjust = 0),
        axis.title.y = element_text(size = 20, family = 'Arial'),
        axis.text.x = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 15, family = 'Arial', colour = 'black'),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15, family = 'Arial')) +
  scale_color_manual(values = c('#A8D0C7','#D78D82'))
dev.off()


png('fig7g.png',
    width = 4,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nCount_chromvar"), pt.size = 0, group.by = c("cell_line5"),idents = c("SHSY5Y_ctrl","SHSY5Y_treat"),
        cols = c('#D38681','#EAB389','#C44943','#B86029','#F6DB88','#8CADCF','#85B69E','#F3C064','#4874A7','#4B7C59')) +
  ylim(0.17,0.38) +
  geom_signif(comparisons = list(c("SHSY5Y_N_ctrl","SHSY5Y_N_treat"), c("SHSY5Y_S_ctrl","SHSY5Y_S_treat")),
              y_position = c(1500,2500), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(-800,2800) +
  theme_classic() + ylab("ChromVAR\nnCount") + xlab("") + ggtitle("") + NoLegend() +
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


png('fig7h.png',
    width = 4,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("nFeature_chromvar"), pt.size = 0, group.by = c("cell_line5"),idents = c("SHSY5Y_ctrl","SHSY5Y_treat"),
        cols = c('#D38681','#EAB389','#C44943','#B86029','#F6DB88','#8CADCF','#85B69E','#F3C064','#4874A7','#4B7C59')) +
  ylim(0.17,0.38) +
  geom_signif(comparisons = list(c("SHSY5Y_N_ctrl","SHSY5Y_N_treat"), c("SHSY5Y_S_ctrl","SHSY5Y_S_treat")),
              y_position = c(640,630), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  ylim(150,700) +
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

### Fig. 7i
png('fig7i_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line5', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line5 %in% c("SHSY5Y_N_ctrl") |NB_combined_250617$cell_line5 %in% c("SHSY5Y_N_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("SH-SY5Y-N (Corr. = 0.44)") + NoLegend() +
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

png('fig7i_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nCount_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line5 %in% c("SHSY5Y_S_ctrl") |NB_combined_250617$cell_line5 %in% c("SHSY5Y_S_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nCount") + ggtitle("SH-SY5Y-S (Corr. = 0.22)") + NoLegend() +
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

### Fig. 7j
png('fig7j_1.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line5', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line5 %in% c("SHSY5Y_N_ctrl") |NB_combined_250617$cell_line5 %in% c("SHSY5Y_N_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D83')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nFeature") + ggtitle("SH-SY5Y-N (Corr. = 0.28)") + NoLegend() +
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

png('fig7j_2.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
FeatureScatter(object = NB_combined_250617, 
               feature1 = "nFeature_chromvar",
               feature2 = 'signature_1PMID286_MES_all', 
               group.by = 'cell_line4', 
               cells = NB_combined_250617@assays[["RNA"]]@data@Dimnames[[2]][NB_combined_250617$cell_line5 %in% c("SHSY5Y_S_ctrl") |NB_combined_250617$cell_line5 %in% c("SHSY5Y_S_treat")],
               shuffle = T,
               cols = c('#A8D0C7','#D78D82')) + 
  xlab("MES signature gene score") + ylab("ChromVAR nFeature") + ggtitle("SH-SY5Y-S (Corr. = 0.19)") + NoLegend() +
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
