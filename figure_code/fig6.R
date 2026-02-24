### iSCAT ### 
### Fig. 6
library(ggplot2)
library(ggpubr)
library(ggsignif)

set.seed(1234)
### Boxplot
SY5Y_iSCAT <- read.table("~/Desktop/SY5Y_iSCAT.txt", header = T)
df <- data.frame(data1 = c(SY5Y_iSCAT$control, SY5Y_iSCAT$ATRA_20uM),
                 data2 = c(rep("Control\n(n=100)", 156), rep("ATRA\n(n=156)", 156)))
df <- na.omit(df)
colnames(df) <- c('data','group')

p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=100)','ATRA\n(n=156)')), y=data)) + 
  geom_boxplot(colour = c('#D38681', '#C44943'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.6,
             position=position_dodge2(width=0.7), size = 1.5) + 
  scale_fill_manual(values = c( '#C44943', '#D38681')) +
  geom_signif(comparisons = list(c('Control\n(n=100)','ATRA\n(n=156)')), 
              annotations = "**", textsize = 10, size = 0.8, vjust = .5, 
              y_position = 2.54) +
  annotate("text", x = 1, y = 2.5, label = "2.01 ± 0.09", size = 5) +
  annotate("text", x = 2, y = 2.5, label = "2.13 ± 0.12", size = 5) +
  theme_classic() + ylab('Condensation level') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(1.7, 2.6)) + 
  scale_y_continuous(breaks = c(1.8, 2, 2.2,2.4))

png('fig6a.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
BE_iSCAT <- read.table("~/Desktop/BE_iSCAT.txt", header = T)
df <- data.frame(data1 = c(BE_iSCAT$control, BE_iSCAT$ATRA_20uM),
                 data2 = c(rep("Control\n(n=32)", 34), rep("ATRA\n(n=34)", 34)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=32)','ATRA\n(n=34)')), y=data)) + 
  geom_boxplot(colour = c('#F6DB88', '#F3C064'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#F3C064','#F6DB88')) +
  geom_signif(comparisons = list(c("Control\n(n=32)", "ATRA\n(n=34)")), 
              annotations = "**", textsize = 10, size = 0.8, vjust = .5, 
              y_position = 2.255) +
  annotate("text", x = 1, y = 2.23, label = "2.02 ± 0.08", size = 5) +
  annotate("text", x = 2, y = 2.23, label = "1.89 ± 0.08", size = 5) +
  theme_classic() + ylab('Condensation level') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(1.7, 2.3)) + 
  scale_y_continuous(breaks = c(1.8, 2, 2.2))

png('fig6b.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
DZ_iSCAT <- read.table("~/Desktop/DZ_iSCAT.txt", header = T)
df <- data.frame(data1 = c(DZ_iSCAT$control, DZ_iSCAT$ATRA_20uM),
                 data2 = c(rep("Control\n(n=40)", 40), rep("ATRA\n(n=40)", 40)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=40)','ATRA\n(n=40)')), y=data)) + 
  geom_boxplot(colour = c('#8CADCF','#4874A7'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#4874A7','#8CADCF')) +
  geom_signif(comparisons = list(c("Control\n(n=40)", "ATRA\n(n=40)")), 
              annotations = "n.s.", textsize = 5, size = 0.8, vjust = -.2, 
              y_position = 2.265) +
  annotate("text", x = 1, y = 2.25, label = "2.04 ± 0.05", size = 5) +
  annotate("text", x = 2, y = 2.25, label = "2.08 ± 0.06", size = 5) +
  theme_classic() + ylab('Condensation level') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(1.95, 2.29)) + 
  scale_y_continuous(breaks = c(2,2.1, 2.2))

png('fig6c.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
SH_iSCAT <- read.table("~/Desktop/SH_iSCAT.txt", header = T)
df <- data.frame(data1 = c(SH_iSCAT$control, SH_iSCAT$ATRA_20uM),
                 data2 = c(rep("Control\n(n=45)", 45), rep("ATRA\n(n=45)", 45)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=45)','ATRA\n(n=45)')), y=data)) + 
  geom_boxplot(colour = c('#85B69E','#4B7C59'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#4B7C59','#85B69E')) +
  geom_signif(comparisons = list(c("Control\n(n=45)", "ATRA\n(n=45)")), 
              annotations = "n.s.", textsize = 5, size = 0.8, vjust = -.2, 
              y_position = 2.2) +
  annotate("text", x = 1, y = 2.18, label = "2.01 ± 0.07", size = 5) +
  annotate("text", x = 2, y = 2.18, label = "1.98 ± 0.07", size = 5) +
  theme_classic() + ylab('Condensation level') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(1.83, 2.23)) + 
  scale_y_continuous(breaks = c(1.9,2, 2.1))

png('fig6d.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()
