### iSCAT paper
### visualization

### Fig. 3
# UMAP-leiden clustering (3A)
# Corr. heatmap (3B)
# Cell line annotation (3C)

# Libraries
library(ggplot2)
library(Seurat)
library(Signac)

setwd("~/R/iSCAT/")

png('fig3a.png',
    width = 8,
    height = 5,
    res = 300,
    unit = 'in')
DimPlot(NB_combined_250617, group.by = 'clusters', pt.size = .5, cols = c('stepped')) + 
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



png('fig3c.png',
    width = 8,
    height = 5,
    res = 300,
    unit = 'in')
DimPlot(NB_combined_250617, pt.size = .5, 
        cols = c('#4B7C59','#85B69E','#4874A7','#8CADCF','#F3C064','#F6DB88','#C44943','#D38681')) + 
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


png('fig3d.png',
    width = 8,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, pt.size = .5, features = c("GAP43"), order = T) +
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

png('fig3e.png',
    width = 8,
    height = 5,
    res = 300,
    unit = 'in')
FeaturePlot(NB_combined_250617, pt.size = .5, features = c("CALR"), order = T, max.cutoff = 10) +
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
