library(ggplot2)
library(ggtree)
library(ggthemes)
library(treeio)
library(scales)
library(deeptime)
library(Cairo)

mcc.tree = read.nhx("sb2-extant-only.mcc.tree")

mcc.tree@phylo$tip.label = gsub("_", " ", mcc.tree@phylo$tip.label)

mcc_plot_revts = revts(ggtree(mcc.tree))

CairoPDF(file = "sb2_ultrametric.pdf", width = 15, height = 8)
mcc_plot_revts +
  geom_tiplab(fontface = 3) +
  coord_geo(xlim = c(-32,5), ylim = c(0,Ntip(mcc.tree) + 1), pos = "bottom",
            skip = "Holocene", dat = list("epochs"), abbrv = FALSE,
            size = 4, neg = TRUE) +
  scale_x_continuous(breaks = seq(-50,0,5), labels = abs(seq(-50,0,5))) +
  scale_shape_identity() +
  labs(tag = "Million years ago") +
  theme_tree2() +
  theme(plot.tag.position = c(0.45, 0), plot.tag = element_text(size = 12, vjust = 1), plot.margin = unit(c(0,0,1,0), "lines"))
dev.off()
