library(ggplot2)
library(Seurat)
library(ggpubr)
library(ggsignif)

### Load most recent data
NB_combined_250617 <- readRDS("~/R/OtherLabMember/kjh/TYT/NB_combined_250617.rds")

### Fig. 4A and related supplementary figures
png('fig4a.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("signature_1PMID286_MES_all"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("MES signature gene score")
dev.off()

png('supp_fig4a.png',
    width = 6,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, features = c("signature_1PMID286_ADRN_all"), order = T, cols = c('grey80',"#D4524E")) +
  ggtitle("ADRN signature gene score")
dev.off()

### Fig. 4B
png('fig4b.png',
    width = 10,
    height = 5,
    res = 300,
    unit = 'in')
VlnPlot(NB_combined_250617, features = c("signature_1PMID286_MES_all"), pt.size = 0, group.by = c("cell_line_4"),
        cols = c('#D38681','#F6DB88','#8CADCF','#85B69E','#C44943','#F3C064','#4874A7','#4B7C59')) +
  ylim(0.17,0.38) +
  geom_signif(comparisons = list(c("SHSY5Y_ctrl","SHSY5Y_treat"), c("SKNBE2C_ctrl","SKNBE2C_treat"),
                                 c("SKNDZ_ctrl",'SKNDZ_treat'),c('SKNSH_ctrl','SKNSH_treat')), 
              y_position = c(0.35, 0.28,0.28,0.25), test = 'wilcox.test',
              textsize = 5, size = 0.8, vjust = -.1) + 
  theme_classic() + ylab("MES signature \n gene score") + xlab("") + ggtitle("") + NoLegend() +
  annotate("text", x = 1.5, y = .39, label = "Δ = -0.0176", size = 5) +
  annotate("text", x = 3.5, y = .32, label = "Δ = 0.0074", size = 5) +
  annotate("text", x = 5.5, y = .32, label = "Δ = -0.0003", size = 5) +
  annotate("text", x = 7.5, y = .29, label = "Δ = -0.0002", size = 5) +
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

### Fig. 4C
###CDF
CDF_df_generate <- function(cell1, cell2,type){
  df <- as.data.frame(NB_combined_250617@meta.data[[type]][NB_combined_250617$cell_line4 %in% cell1])
  df$cell <-NB_combined_250617$cell_line4[NB_combined_250617$cell_line4 %in% cell1]
  colnames(df) <- c("CDF","Cell")
  df2 <- as.data.frame(NB_combined_250617@meta.data[[type]][NB_combined_250617$cell_line4 %in% cell2])
  df2$cell <-NB_combined_250617$cell_line4[NB_combined_250617$cell_line4 %in% cell2]
  colnames(df2) <- c("CDF","Cell")
  df <- rbind(df, df2)
  
  return(df)
}

df <- CDF_df_generate('SKNSH_treat','SKNSH_ctrl', "signature_1PMID286_MES_all")

png('fig4f.png',
    width = 5,
    height = 4.5,
    res = 300,
    unit = 'in')
ggplot(df, aes(CDF, colour = Cell)) + ylab("CDF") + xlab("MES signature gene score") +
  stat_ecdf(size = 1.5) + ggtitle("Cell line: SK-N-SH") + theme_classic() + NoLegend() +
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
