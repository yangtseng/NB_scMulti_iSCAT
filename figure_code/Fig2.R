### Number of neurites
library(ggplot2)
BE_neurites <- data.frame(data1 = c(2.5, 41.5, 43.5), 
                          data2 = c(1.63299, 5.24934, 2.86744),
                          data3 = c(0, 10, 20))

colnames(BE_neurites) <- c('no','sd', 'conc')

library(ggplot2)
library(ggsignif)

# Define significance annotations (stars) - modify these as needed
star_for_second_bar <- "**"  # Significance for comparison: group 0 vs group 1
star_for_third_bar <- "***"    # Significance for comparison: group 0 vs group 5

# Define y positions for stars - adjust based on your data heights
y_pos_second_bar <- 48  # Y position for star above second bar
y_pos_third_bar <- 48   # Y position for star above third bar

p <- ggplot(BE_neurites, aes(x=factor(conc), y=no)) + 
  geom_bar(stat = 'identity', fill = c('#F6DB88'), width = .5, linewidth = 1, color = 'black') +
  geom_errorbar(aes(ymin=no-sd, ymax=no+sd), 
                width=.1,
                size=1.5,  # Thicker error bars
                position=position_dodge(.9)) +
  annotate("text", x = 2, y = y_pos_second_bar, 
           label = star_for_second_bar, size = 10) +  # Star on second bar
  annotate("text", x = 3, y = y_pos_third_bar, 
           label = star_for_third_bar, size = 10) +   # Star on third bar
  theme_classic() + 
  ylab('Number of neurites') + 
  xlab('ATRA Conc. (μM)') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
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
  coord_cartesian(ylim=c(0, 52), expand = c(0,1))

p

png('fig2a.png',
    width = 3.5,
    height = 3,
    res = 300,
    unit = 'in')
p
dev.off()

### Number of neurites
library(ggplot2)
DZ_neurites <- data.frame(data1 = c(14.25, 19.25, 20), 
                          data2 = c(4.112988, 2.061553, 12.49),
                          data3 = c(0, 10, 20))

colnames(DZ_neurites) <- c('no','sd', 'conc')

library(ggplot2)
library(ggsignif)

# Define significance annotations (stars) - modify these as needed
star_for_second_bar <- "*"  # Significance for comparison: group 0 vs group 1
star_for_third_bar <- ""    # Significance for comparison: group 0 vs group 5

# Define y positions for stars - adjust based on your data heights
y_pos_second_bar <- 25  # Y position for star above second bar
y_pos_third_bar <- 20   # Y position for star above third bar

p <- ggplot(DZ_neurites, aes(x=factor(conc), y=no)) + 
  geom_bar(stat = 'identity', fill = c('#8CADCF'), width = .5, linewidth = 1, color = 'black') +
  geom_errorbar(aes(ymin=no-sd, ymax=no+sd), 
                width=.1,
                size=1.5,  # Thicker error bars
                position=position_dodge(.9)) +
  annotate("text", x = 2, y = y_pos_second_bar, 
           label = star_for_second_bar, size = 10) +  # Star on second bar
  annotate("text", x = 3, y = y_pos_third_bar, 
           label = star_for_third_bar, size = 10) +   # Star on third bar
  theme_classic() + 
  ylab('Number of neurites') + 
  xlab('ATRA Conc. (μM)') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
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
  coord_cartesian(ylim=c(0, 35), expand = c(0,1))

p

png('fig2b.png',
    width = 3.5,
    height = 3,
    res = 300,
    unit = 'in')
p
dev.off()

### Number of neurites
library(ggplot2)
SH_neurites <- data.frame(data1 = c(18.67, 58.67, 57.67), 
                            data2 = c(3.05505, 8.0829, 15.5671),
                            data3 = c(0, 10, 20))

colnames(SH_neurites) <- c('no','sd', 'conc')

library(ggplot2)
library(ggsignif)

# Define significance annotations (stars) - modify these as needed
star_for_second_bar <- "**"  # Significance for comparison: group 0 vs group 1
star_for_third_bar <- "**"    # Significance for comparison: group 0 vs group 5

# Define y positions for stars - adjust based on your data heights
y_pos_second_bar <- 69  # Y position for star above second bar
y_pos_third_bar <- 77   # Y position for star above third bar

p <- ggplot(SH_neurites, aes(x=factor(conc), y=no)) + 
  geom_bar(stat = 'identity', fill = c('#85B69E'), width = .5, linewidth = 1, color = 'black') +
  geom_errorbar(aes(ymin=no-sd, ymax=no+sd), 
                width=.1,
                size=1.5,  # Thicker error bars
                position=position_dodge(.9)) +
  annotate("text", x = 2, y = y_pos_second_bar, 
           label = star_for_second_bar, size = 10) +  # Star on second bar
  annotate("text", x = 3, y = y_pos_third_bar, 
           label = star_for_third_bar, size = 10) +   # Star on third bar
  theme_classic() + 
  ylab('Number of neurites') + 
  xlab('ATRA Conc. (μM)') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
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
  coord_cartesian(ylim=c(0, 84), expand = c(0,1))

p

png('fig2c.png',
    width = 3.5,
    height = 3,
    res = 300,
    unit = 'in')
p
dev.off()


library(ggplot2)
SY5Y_neurites <- data.frame(data1 = c(14.67, 29.50, 46.00), 
                            data2 = c(7.094599, 7.778175, 8.888194),
                            data3 = c(0, 10, 20))

colnames(SY5Y_neurites) <- c('no','sd', 'conc')

library(ggplot2)
library(ggsignif)

# Define significance annotations (stars) - modify these as needed
star_for_second_bar <- ""  # Significance for comparison: group 0 vs group 1
star_for_third_bar <- "**"    # Significance for comparison: group 0 vs group 5

# Define y positions for stars - adjust based on your data heights
y_pos_second_bar <- 40  # Y position for star above second bar
y_pos_third_bar <- 58   # Y position for star above third bar

p <- ggplot(SY5Y_neurites, aes(x=factor(conc), y=no)) + 
  geom_bar(stat = 'identity', fill = c('#D38681'), width = .5, linewidth = 1, color = 'black') +
  geom_errorbar(aes(ymin=no-sd, ymax=no+sd), 
                width=.1,
                size=1.5,  # Thicker error bars
                position=position_dodge(.9)) +
  annotate("text", x = 2, y = y_pos_second_bar, 
           label = star_for_second_bar, size = 10) +  # Star on second bar
  annotate("text", x = 3, y = y_pos_third_bar, 
           label = star_for_third_bar, size = 10) +   # Star on third bar
  theme_classic() + 
  ylab('Number of neurites') + 
  xlab('ATRA Conc. (μM)') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(size = 18, color = 'black', vjust = 1), 
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
  coord_cartesian(ylim=c(0, 62), expand = c(0,1))

p

png('fig2d.png',
    width = 3.5,
    height = 3,
    res = 300,
    unit = 'in')
p
dev.off()
