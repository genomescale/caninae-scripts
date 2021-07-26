library(ggplot2)
library(ggtree)
library(ggthemes)
library(treeio)
library(scales)
library(cowplot)
library(Cairo)

mcc_tree = read.nhx("sb2-total-evidence-pruned-concat-total-evidence-pruned.tree")
mcc_tree@phylo$tip.label = gsub("_", " ", mcc_tree@phylo$tip.label)

mcc_plot_revts = revts(ggtree(mcc_tree, aes(linetype = probability < 0.005)) + scale_linetype_discrete(guide = "none"))

diff_table = read.csv("sb2-total-evidence-pruned-concat-total-evidence-pruned.csv")
diff_table = subset(merge(mcc_plot_revts$data, diff_table, by = "id"), probability >= 0.005)
diff_table$direction = ifelse(diff_table$ldiff < 0, "Shorter", "Longer")

CairoPDF(file = "fbd_length_differences.pdf", width = 7, height = 4)
mcc_plot_revts +
  geom_segment(data = diff_table, size = 0.7, lineend = "round", aes(x = start, y = y, xend = end, yend = y, color = direction)) +
  geom_tiplab(color = "black", size = 2.5, fontface = 3) +
  geom_label(color = "black", size = 2.0, label.size = 0.0, nudge_x = 0.12, fill = NA, aes(label = tag)) +
  scale_x_continuous(breaks = seq(-14, 0), labels = abs(seq(-14, 0)), limits = c(-13.5, 2.5)) +
  labs(color = "FBD-concatenation branch length relative to FBD-MSC", tag = "Million years ago") +
  scale_colour_tableau() +
  theme_tree2() +
  background_grid(major = "x") +
  theme(legend.position = "top", legend.justification = "left", plot.tag.position = c(0.5, 0), plot.tag = element_text(size = 11, vjust = 1), plot.margin = unit(c(0,0,1,0), "lines"))
dev.off()
