library(ggplot2)
library(ggthemes)
library(ggstance)
library(cowplot)
library(scales)
library(Cairo)

locus_sources = read.csv("whitelist.csv", header = FALSE, col.names = c("locus"))
sequence_sources = read.csv("sequence-sources.csv")
label_mapping = read.csv("label-mapping.csv", header = FALSE, col.names = c("organism", "binomial"))
filtered_sources = subset(merge(sequence_sources, label_mapping, by = "organism"), locus %in% locus_sources$locus)
filtered_sources$binomial = sub("_", " ", filtered_sources$binomial)
filtered_sources$source = factor(filtered_sources$source, levels = c("Bardeleben2005a", "Bardeleben2005b", "Koepfli2015", "LindbladToh2005", "Shang2017"))
levels(filtered_sources$source) = c("Bardeleben et al. (2005a)", "Bardeleben et al. (2005b)", "Koepfli et al. (2015)", "Lindblad-Toh et al. (2005)", "Shang et al. (2017)")

# msa_lengths = read.csv("msa-lengths.csv")

sources_plot = ggplot(filtered_sources, aes(x = binomial, y = locus, fill = source)) +
  geom_tile() +
  scale_fill_tableau(guide = guide_legend(nrow = 3)) +
  labs(x = "Species", y = "Locus", fill = "Source") +
  theme_cowplot(11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.text = element_text(margin = margin(r = 10)))

# msa_plot = ggplot(msa_lengths, aes(y = locus, x = length)) +
#   geom_barh(stat = "identity") +
#   labs(x = "Length") +
#   theme_cowplot(11) +
#   background_grid(major = "x") +
#   theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5))

CairoPDF(file = "filtered_sources.pdf", width = 6, height = 11)
sources_plot
# plot_grid(sources_plot, msa_plot, ncol = 2, rel_widths = c(0.66, 0.33), align = "h", axis = "b")
dev.off()

write.csv(filtered_sources, file = "accession-source-table.csv", row.names = FALSE)

length(filtered_sources$binomial)
