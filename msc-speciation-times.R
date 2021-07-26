library(ggplot2)
library(ggstance)
library(ggthemes)
library(scales)
library(cowplot)

sb2_extant_nodes = read.csv("sb2-extant-only-pruned.nodes.csv")
sb2_nodes = read.csv("sb2-total-evidence-pruned.nodes.csv")

sb2_extant_nodes$label = merge(sb2_extant_nodes, sb2_nodes, by = c("bitmask"), all.x = TRUE)$label.y

sb2_extant_nodes$method = "BD-MSC"
sb2_nodes$method = "FBD-MSC"

nodes = rbind(sb2_nodes, sb2_extant_nodes)

mcc_nodes = subset(nodes, label != "" & !is.na(label) & age_mean > 0.001 & probability >= 0.005)
mcc_nodes$label = as.factor(mcc_nodes$label)
mcc_nodes$label = factor(mcc_nodes$label, levels = rev(levels(mcc_nodes$label)))

times_plot = ggplot(mcc_nodes, aes(y = label, x = age_mean, color = method)) +
	geom_hline(yintercept = seq(1.5, length(unique(mcc_nodes$label))-0.5, 1), lwd = 0.5, color = "grey85") +
	geom_pointrangeh(aes(xmin = age_low, xmax = age_high), position = position_dodgev(height = 0.65), fatten = 0.25) +
	labs(x = "Million years ago", y = "Caninae tree node") +
	scale_x_reverse(breaks = seq(0, 30, 5), limits = c(32, 0)) +
	scale_y_discrete(position = "right") +
	scale_color_tableau(palette = "Tableau 10") +
	theme_half_open() +
	background_grid(major = "x") +
	theme(axis.title.x = element_blank(), legend.position = "top", legend.title = element_blank(), legend.direction = "horizontal", axis.text.y = element_text(size = 9))

cairo_pdf(filename = "msc_speciation_times.pdf", width = 7, height = 7)
times_plot
dev.off()
