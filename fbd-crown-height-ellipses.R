library(ggplot2)
library(ggstance)
library(ggforce)
library(scales)
library(cowplot)

sb2_full_nodes = read.csv("sb2-total-evidence-pruned.nodes.csv")
concat_full_nodes = read.csv("concat-total-evidence-pruned.nodes.csv")

ticks = 2^seq(-3, 4)

sb2_v_concat_nodes = merge(sb2_full_nodes, concat_full_nodes, by = c("bitmask"), suffixes = c(".sb2", ".concat"))

sb2_v_concat_mcc_nodes = subset(sb2_v_concat_nodes, label.sb2 != "" & !is.na(label.sb2) & age_mean.sb2 > 0.001 & probability.concat >= 0.005)
sb2_v_concat_mcc_nodes$label = as.factor(sb2_v_concat_mcc_nodes$label.sb2)
sb2_v_concat_mcc_nodes$label = factor(sb2_v_concat_mcc_nodes$label.sb2, levels = rev(levels(sb2_v_concat_mcc_nodes$label.sb2)))

sb2_v_concat_ellipse_plot = ggplot(sb2_v_concat_mcc_nodes, aes(x = log_age_mean.sb2, y = log_age_mean.concat)) +
  geom_smooth(formula = y ~ poly(x, 2), method = lm, fullrange = TRUE, se = FALSE, linetype = 3, color = "black") +
  geom_ellipse(aes(x0 = log_age_mean.sb2, y0 = log_age_mean.concat, a = log_age_sd.sb2, b = log_age_sd.concat, angle = 0)) +
  geom_text(aes(label = label.sb2)) +
  scale_x_continuous(limits = c(log(min(ticks)), log(max(ticks))), breaks = log(ticks), labels = ticks) +
  scale_y_continuous(limits = c(log(0.125), log(max(ticks))), breaks = log(ticks), labels = ticks) +
  geom_abline(linetype = 5) +
  labs(x = "FBD-MSC clade age (million years ago)", y = "FBD-concatenation clade age") +
  theme_cowplot()

cairo_pdf(filename = "fbd_crown_height_ellipses.pdf", width = 5, height = 5)
sb2_v_concat_ellipse_plot
dev.off()
