library(ggplot2)
library(ggtree)
library(ggthemes)
library(gginnards)
library(treeio)
library(scales)
library(deeptime)
library(Cairo)

sa.tree = read.nhx("sb2-total-evidence.mcc.tree")

sa.tree@data$sa = as.logical(sa.tree@data$sa)

mcc_plot_revts = revts(ggtree(sa.tree))

p = mcc_plot_revts +
  geom_nodepoint(fill = "black", size = 2, aes(shape = ifelse(sa, 16, 32))) +
  annotate("rect", xmin = -2.58, xmax = -5.333, ymin = 0, ymax = Inf, fill = "grey80") +
  annotate("rect", xmin = -23.03, xmax = -33.9, ymin = 0, ymax = Inf, fill = "grey80") +
  coord_geo(xlim = c(-32,0), ylim = c(0,Ntip(sa.tree) + 1), pos = "bottom",
            skip = "Holocene", dat = list("epochs"), abbrv = FALSE,
            size = 4, neg = TRUE) +
  scale_shape_identity() +
  theme_tree2() +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank())

CairoSVG(file = "cover", width = 15, height = 15)
move_layers(p, "GeomRect", position = "bottom")
dev.off()
