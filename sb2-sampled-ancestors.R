library(ggplot2)
library(ggtree)
library(ggthemes)
library(treeio)
library(scales)
library(deeptime)
library(Cairo)

sa.tree = read.nhx("sb2-total-evidence.mcc.tree")

sa.tree@phylo$tip.label = gsub("_", " ", sa.tree@phylo$tip.label)
sa.tree@phylo$node.label = gsub("_", " ", sa.tree@phylo$node.label)
sa.tree@data$sa = as.logical(sa.tree@data$sa)

draw.label.above = c(65, 67, 99, 100, 109, 115)

mcc_plot_revts = revts(ggtree(sa.tree))

CairoPDF(file = "sb2_sampled_ancestors.pdf", width = 15, height = 15)
mcc_plot_revts %>% flip(124, 127) %>% flip(109, 102) %>% flip(97, 101) %>% flip(19, 99) +
  geom_nodepoint(fill = "black", size = 2, aes(shape = ifelse(sa, ifelse(node %in% draw.label.above, 24, 25), 32))) +
  geom_nodelab(geom = "label", label.size = 0, label.padding = unit(0.35, "lines"), fill = NA, hjust = 1, nudge_x = 0.6, aes(vjust = ifelse(node %in% draw.label.above, 0, 0.9)), fontface = 3) +# geom_text(aes(label = node)) +
  geom_tiplab(fontface = 3) +
  coord_geo(xlim = c(-32,5), ylim = c(0,Ntip(sa.tree) + 1), pos = "bottom",
            skip = "Holocene", dat = list("epochs"), abbrv = FALSE,
            size = 4, neg = TRUE) +
  scale_x_continuous(breaks = seq(-50,0,5), labels = abs(seq(-50,0,5))) +
  scale_shape_identity() +
  labs(tag = "Million years ago") +
  theme_tree2() +
  theme(plot.tag.position = c(0.45, 0), plot.tag = element_text(size = 12, vjust = 1), plot.margin = unit(c(0,0,1,0), "lines"))
dev.off()
