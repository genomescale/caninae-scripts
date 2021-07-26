library(ggplot2)
library(ggtree)
library(ggthemes)
library(treeio)
library(scales)
library(Cairo)

tas2r.tree = read.tree("tas2r.tree")

tas2r.tree$tip.label = gsub("_", " ", tas2r.tree$tip.label)

CairoPDF(file = "tas2r.pdf", width = 15, height = 17)
ggtree(tas2r.tree) +
  theme_tree() +
  geom_tiplab(size = 2.5) +
  geom_treescale(x = 0) +
  scale_x_continuous(limits = c(0, 2.2))
dev.off()
