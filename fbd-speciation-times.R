library(ggplot2)
library(ggstance)
library(ggthemes)
library(scales)
library(cowplot)

sb2_nodes = read.csv("sb2-total-evidence-pruned.nodes.csv")
concat_nodes = read.csv("concat-total-evidence-pruned.nodes.csv")

concat_nodes$label = merge(concat_nodes, sb2_nodes, by = c("bitmask"), all.x = TRUE)$label.y

sb2_nodes$method = "FBD-MSC"
concat_nodes$method = "FBD-concatenation"

nodes = rbind(sb2_nodes, concat_nodes)

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

sb2_ltt = read.csv("sb2-total-evidence.combined.ltt.csv")
concat_ltt = read.csv("concat-total-evidence.combined.ltt.csv")

sb2_ltt$method = "FBD-MSC"
concat_ltt$method = "FBD-concatenation"

ltt = rbind(sb2_ltt, concat_ltt)

ltt_plot = ggplot(ltt, aes(x = mya, y = mean, ymin = hpd_low, ymax = hpd_high)) +
	geom_step(aes(color = method)) +
	geom_ribbon(alpha = 0.3, aes(fill = method)) +
	# geom_abline(slope = 1, intercept = 0) +
	labs(x = "Million years ago", y = "Number of lineages") +
	scale_x_continuous(breaks = seq(-50, 0, 5), labels = seq(50, 0, -5)) +
	scale_y_continuous(position = "right") +
	# scale_y_continuous(limits = c(0, 1)) +
	scale_color_tableau(palette = "Tableau 10") +
	scale_fill_tableau(palette = "Tableau 10") +
	theme_half_open() +
	background_grid() +
	theme(legend.position = "top", legend.title = element_blank(), legend.direction = "horizontal", axis.text.y = element_text(size = 9))

legend_b = get_legend(times_plot + theme(legend.justification="left"))
prow = plot_grid(times_plot + theme(legend.position = "none"),
                 ltt_plot + theme(legend.position = "none"),
                 rel_heights = c(3, 2), ncol = 1, align = "v", labels = c("A", "B"))

cairo_pdf(filename = "fbd_speciation_times.pdf", width = 7, height = 7)
plot_grid(legend_b, prow, ncol = 1, rel_heights = c(0.05, 0.95))
dev.off()
