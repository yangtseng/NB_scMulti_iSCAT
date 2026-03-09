### CV ### 
### Fig. Supp. 5B
library(ggplot2)
library(ggpubr)
library(ggsignif)

set.seed(1234)
### Boxplot
SY5Y_CV <- read.table("~/Downloads/SY5Y_cv.txt", header = T)
wilcox.test(SY5Y_CV$ATRA, SY5Y_CV$ctrl)
df <- data.frame(data1 = c(SY5Y_CV$ctrl, SY5Y_CV$ATRA),
                 data2 = c(rep("Control\n(n=56)", 56), rep("ATRA\n(n=52)", 56)))
df <- na.omit(df)
colnames(df) <- c('data','group')

p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=56)','ATRA\n(n=52)')), y=data)) + 
  geom_boxplot(colour = c('#D38681', '#C44943'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.6,
             position=position_dodge2(width=0.7), size = 1.5) + 
  scale_fill_manual(values = c( '#C44943', '#D38681')) +
  geom_signif(comparisons = list(c('Control\n(n=56)','ATRA\n(n=52)')), 
              annotations = "***", textsize = 10, size = 0.8, vjust = .5, 
              y_position = 1.75) +

  theme_classic() + ylab('Coefficient of variation') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(0.65, 1.85)) + 
  scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4 ,1.6, 1.8))
p

png('supp_fig5b_a.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
BE_CV <- read.table("~/Downloads/BE2C_CV.txt", header = T)
wilcox.test(BE_CV$ATRA, BE_CV$ctrl)
df <- data.frame(data1 = c(BE_CV$ctrl, BE_CV$ATRA),
                 data2 = c(rep("Control\n(n=50)", 50), rep("ATRA\n(n=50)", 50)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=50)','ATRA\n(n=50)')), y=data)) + 
  geom_boxplot(colour = c('#F6DB88', '#F3C064'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#F3C064','#F6DB88')) +
  geom_signif(comparisons = list(c("Control\n(n=50)", "ATRA\n(n=50)")), 
              annotations = "***", textsize = 10, size = 0.8, vjust = .5, 
              y_position = 1.35) +
  theme_classic() + ylab('Coefficient of variation') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .4),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(0.7, 1.4)) + 
  scale_y_continuous(breaks = c(0.8, 1, 1.2))
p
png('supp_fig5b_b.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
DZ_CV <- read.table("~/Downloads/DZ_CV.txt", header = T)
wilcox.test(DZ_CV$ATRA, DZ_CV$ctrl)
df <- data.frame(data1 = c(DZ_CV$ctrl, DZ_CV$ATRA),
                 data2 = c(rep("Control\n(n=50)", 50), rep("ATRA\n(n=50)", 50)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=50)','ATRA\n(n=50)')), y=data)) + 
  geom_boxplot(colour = c('#8CADCF','#4874A7'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#4874A7','#8CADCF')) +
  geom_signif(comparisons = list(c("Control\n(n=50)", "ATRA\n(n=50)")), 
              annotations = "***", textsize = 10, size = 0.8, vjust = .5, 
              y_position = 1.05) +
  theme_classic() + ylab('Coefficient of variation') + xlab('') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
        axis.text.y = element_text(size = 20, color = 'black'),
        axis.title.x = element_text(size = 20, vjust = -1), 
        axis.title.y = element_text(size = 20, vjust = 2, hjust = .3),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.9),
        axis.ticks = element_line(size = 0.9),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  coord_cartesian(ylim=c(0.5, 1.1)) + 
  scale_y_continuous(breaks = c(0.6, 0.8, 1.0))
p
png('supp_fig5b_c.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()

### Boxplot
SH_CV <- read.table("~/Downloads/SH_CV.txt", header = T)
wilcox.test(SH_CV$ATRA, SH_CV$ctrl)
df <- data.frame(data1 = c(SH_CV$ctrl, SH_CV$ATRA),
                 data2 = c(rep("Control\n(n=50)", 50), rep("ATRA\n(n=50)", 50)))

colnames(df) <- c('data','group')
df <- na.omit(df)
p <- ggplot(df, aes(x=factor(group, levels = c('Control\n(n=50)','ATRA\n(n=50)')), y=data)) + 
  geom_boxplot(colour = c('#85B69E','#4B7C59'), outlier.shape = NA, size = .8) + 
  geom_point(aes(fill = group), colour = 'black', pch = 21, stroke = 0.8,
             position=position_dodge2(width=0.7), size = 2) + 
  scale_fill_manual(values = c('#4B7C59','#85B69E')) +
  geom_signif(comparisons = list(c("Control\n(n=50)", "ATRA\n(n=50)")), 
              annotations = "n.s.", textsize = 6, size = 0.8, vjust = 0, 
              y_position = 0.9) +
  theme_classic() + ylab('Coefficient of variation') + xlab('') +
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
  coord_cartesian(ylim=c(0.4, 0.95)) + 
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9))
p
png('supp_fig5b_d.png',
    width = 4,
    height = 4,
    res = 300,
    unit = 'in')
p
dev.off()
