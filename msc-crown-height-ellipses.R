library(ggplot2)
library(ggstance)
library(ggforce)
library(scales)
library(cowplot)

sb2_extant_nodes = read.csv("sb2-extant-only-pruned.nodes.csv")
sb2_full_nodes = read.csv("sb2-total-evidence-pruned.nodes.csv")

ticks = 2^seq(-3, 4)

fbd_v_bd_nodes = merge(sb2_full_nodes, sb2_extant_nodes, by = c("bitmask"), suffixes = c(".fbd", ".bd"))

fbd_v_bd_mcc_nodes = subset(fbd_v_bd_nodes, label.fbd != "" & !is.na(label.fbd) & age_mean.fbd > 0.001 & probability.bd >= 0.005)
fbd_v_bd_mcc_nodes$label = as.factor(fbd_v_bd_mcc_nodes$label.fbd)
fbd_v_bd_mcc_nodes$label = factor(fbd_v_bd_mcc_nodes$label.fbd, levels = rev(levels(fbd_v_bd_mcc_nodes$label.fbd)))

fbd_v_bd_ellipse_plot = ggplot(fbd_v_bd_mcc_nodes, aes(x = log_age_mean.fbd, y = log_age_mean.bd)) +
  geom_smooth(formula = y ~ poly(x, 2), method = lm, fullrange = TRUE, se = FALSE, linetype = 3, color = "black") +
  geom_ellipse(aes(x0 = log_age_mean.fbd, y0 = log_age_mean.bd, a = log_age_sd.fbd, b = log_age_sd.bd, angle = 0)) +
  geom_text(aes(label = label.fbd)) +
  scale_x_continuous(limits = c(log(min(ticks)), log(max(ticks))), breaks = log(ticks), labels = ticks) +
  scale_y_continuous(limits = c(log(0.125), log(max(ticks))), breaks = log(ticks), labels = ticks) +
  geom_abline(linetype = 5) +
  labs(x = "FBD-MSC clade age (million years ago)", y = "BD-MSC clade age") +
  theme_cowplot()

cairo_pdf(filename = "msc_crown_height_ellipses.pdf", width = 5, height = 5)
fbd_v_bd_ellipse_plot
dev.off()
