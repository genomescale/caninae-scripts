library(ggplot2)
library(ggstance)
library(scales)
library(ggthemes)
library(cowplot)
library(Cairo)

sb2_log = read.table("sb2-total-evidence.combined.log", header = TRUE)
concat_log = read.table("concat-total-evidence.combined.log", header = TRUE)

sb2_clock = read.csv("sb2-total-evidence-pruned.clock.csv")
concat_clock = read.csv("concat-total-evidence-pruned.clock.csv")

sb2_log$Method = "FBD-MSC"
concat_log$Method = "FBD-concatenation"

sb2_clock$Method = "FBD-MSC"
concat_clock$Method = "FBD-concatenation"

sb2_clock$ChildLeafCount = as.factor(sb2_clock$ChildLeafCount)
concat_clock$ChildLeafCount = as.factor(concat_clock$ChildLeafCount)

sb2_molecular = subset(sb2_clock, Rate == "molecular" & Age > 0)
sb2_morphology = subset(sb2_clock, Rate == "morphology" & Age > 0)

concat_molecular = subset(concat_clock, Rate == "molecular" & Age > 0)
concat_morphology = subset(concat_clock, Rate == "morphology" & Age > 0)

molecular_times = rbind(sb2_molecular, concat_molecular)
morphology_times = rbind(sb2_morphology, concat_morphology)

combined_log = data.frame(MolecularClockRate = c(sb2_log$strictClockRate.c.molecular, concat_log$strictClockRate.c.molecular),
                          MorphologyClockRate = c(sb2_log$strictClockRate.c.morphology, concat_log$strictClockRate.c.morphology),
                          Method = c(sb2_log$Method, concat_log$Method))

molecular_times = ggplot(molecular_times, aes(x = Age, y = SubsPerSite, color = Method)) +
  geom_point(size = 0.1, alpha = 0.1) +
  scale_y_log10(limits = c(min(molecular_times$SubsPerSite), max(molecular_times$SubsPerSite))) +
  scale_x_log10(limits = c(min(molecular_times$Age), max(molecular_times$Age)), breaks = 2^(seq(-10, 10, 2))) +
  scale_color_tableau() +
  labs(x = "Million years ago", y = "Expected substitutions/site") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme_half_open() +
  background_grid() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

morphology_times = ggplot(morphology_times, aes(x = Age, y = SubsPerSite, color = Method)) +
  geom_point(size = 0.1, alpha = 0.1) +
  scale_y_log10(limits = c(min(morphology_times$SubsPerSite), max(morphology_times$SubsPerSite))) +
  scale_x_log10(limits = c(min(morphology_times$Age), max(morphology_times$Age)), breaks = 2^(seq(-10, 10, 2))) +
  scale_color_tableau() +
  labs(x = "Million years ago", y = "Expected substitutions/trait") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme_half_open() +
  background_grid() +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

legend_b = get_legend(molecular_times + theme(legend.justification="left"))
prow = plot_grid(molecular_times + theme(legend.position = "none"),
                 morphology_times + theme(legend.position = "none"),
                 ncol = 2, align = "v", labels = c("A", "B"))

CairoPNG(filename = "substitutions.png", width = 2600, height = 1300, dpi = 300, bg = "white")
plot_grid(legend_b, prow, ncol = 1, rel_heights = c(0.05, 0.95))
dev.off()
