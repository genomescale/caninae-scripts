library(ggplot2)
library(ggstance)
library(scales)
library(ggthemes)
library(cowplot)
library(Cairo)

sb2_extant_log = read.table("sb2-extant-only.combined.log", header = TRUE)
sb2_full_log = read.table("sb2-total-evidence.combined.log", header = TRUE)
concat_full_log = read.table("concat-total-evidence.combined.log", header = TRUE)

sb2_extant_clock = read.csv("sb2-extant-only-pruned.clock.csv")
sb2_full_clock = read.csv("sb2-total-evidence-pruned.clock.csv")
concat_full_clock = read.csv("concat-total-evidence-pruned.clock.csv")

sb2_extant_log$Method = "BD-MSC"
sb2_full_log$Method = "FBD-MSC"
concat_full_log$Method = "FBD-concatenation"

sb2_extant_clock$Method = "BD-MSC"
sb2_full_clock$Method = "FBD-MSC"
concat_full_clock$Method = "FBD-concatenation"

sb2_extant_clock$ChildLeafCount = as.factor(sb2_extant_clock$ChildLeafCount)
sb2_full_clock$ChildLeafCount = as.factor(sb2_full_clock$ChildLeafCount)
concat_full_clock$ChildLeafCount = as.factor(concat_full_clock$ChildLeafCount)

sb2_extant_molecular = subset(sb2_extant_clock, Rate == "molecular" & Age > 0)
sb2_extant_morphology = subset(sb2_extant_clock, Rate == "morphology" & Age > 0)

sb2_full_molecular = subset(sb2_full_clock, Rate == "molecular" & Age > 0)
sb2_full_morphology = subset(sb2_full_clock, Rate == "morphology" & Age > 0)

concat_full_molecular = subset(concat_full_clock, Rate == "molecular" & Age > 0)
concat_full_morphology = subset(concat_full_clock, Rate == "morphology" & Age > 0)

combined_log = data.frame(MolecularClockRate = c(sb2_extant_log$strictClockRate.c.molecular, sb2_full_log$strictClockRate.c.molecular, concat_full_log$strictClockRate.c.molecular),
                          MorphologyClockRate = c(sb2_extant_log$strictClockRate.c.morphology, sb2_full_log$strictClockRate.c.morphology, concat_full_log$strictClockRate.c.morphology),
                          Method = c(sb2_extant_log$Method, sb2_full_log$Method, concat_full_log$Method))

fbd_log = subset(combined_log, startsWith(Method, "FBD"))
msc_log = subset(combined_log, endsWith(Method, "MSC"))

fbd_molecular_rates = ggplot(fbd_log, aes(x = MolecularClockRate, y = (..count..)/1984, fill = Method)) +
  geom_histogram(binwidth = 0.00002, position = "identity", alpha = 0.7, boundary = 0) +
  scale_x_continuous(breaks = seq(0.0004, 0.0012, 0.0002), labels = scientific) +
  scale_y_continuous(labels = percent, limits = c(0, 0.2)) +
  scale_fill_tableau(palette = "Tableau 10") +
  theme_half_open() +
  background_grid() +
  labs(y = "Posterior probability", x = "Substitutions/site/million years") +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

fbd_morphology_rates = ggplot(fbd_log, aes(x = MorphologyClockRate, y = (..count..)/1984, fill = Method)) +
  geom_histogram(binwidth = 0.0005, position = "identity", alpha = 0.7, boundary = 0) +
  scale_x_continuous(labels = scientific) +
  scale_y_continuous(labels = percent, limits = c(0, 0.2)) +
  scale_fill_tableau(palette = "Tableau 10") + 
  theme_half_open() +
  background_grid() +
  labs(y = "", x = "Substitutions/trait/million years") +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

fbd_legend = get_legend(fbd_molecular_rates + theme(legend.justification="left"))
fbd_prow = plot_grid(fbd_molecular_rates + theme(legend.position = "none"),
                 fbd_morphology_rates + theme(legend.position = "none"),
                 ncol = 2, align = "v", labels = c("A", "B"))

cairo_pdf(filename = "fbd_clocks.pdf", width = 7, height = 3.75)
plot_grid(fbd_legend, fbd_prow, ncol = 1, rel_heights = c(0.05, 0.95))
dev.off()

msc_molecular_rates = ggplot(msc_log, aes(x = MolecularClockRate, y = (..count..)/1984, fill = Method)) +
  geom_histogram(binwidth = 0.00002, position = "identity", alpha = 0.7, boundary = 0) +
  scale_x_continuous(breaks = seq(0.0004, 0.0012, 0.0002), labels = scientific) +
  scale_y_continuous(labels = percent, limits = c(0, 0.2)) +
  scale_fill_tableau(palette = "Tableau 10") +
  theme_half_open() +
  background_grid() +
  labs(y = "Posterior probability", x = "Substitutions/site/million years") +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

msc_morphology_rates = ggplot(msc_log, aes(x = MorphologyClockRate, y = (..count..)/1984, fill = Method)) +
  geom_histogram(binwidth = 0.0005, position = "identity", alpha = 0.7, boundary = 0) +
  scale_x_continuous(labels = scientific) +
  scale_y_continuous(labels = percent, limits = c(0, 0.2)) +
  scale_fill_tableau(palette = "Tableau 10") + 
  theme_half_open() +
  background_grid() +
  labs(y = "", x = "Substitutions/trait/million years") +
  theme(legend.title = element_blank(), legend.direction = "horizontal", legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))

msc_legend = get_legend(msc_molecular_rates + theme(legend.justification="left"))
msc_prow = plot_grid(msc_molecular_rates + theme(legend.position = "none"),
                     msc_morphology_rates + theme(legend.position = "none"),
                     ncol = 2, align = "v", labels = c("A", "B"))

cairo_pdf(filename = "msc_clocks.pdf", width = 7, height = 3.75)
plot_grid(msc_legend, msc_prow, ncol = 1, rel_heights = c(0.05, 0.95))
dev.off()
