rm(list = ls())
library(tidyverse)
load("./output/Figure4/prop_diabt_glucose_raw.RData")

p_prop_diabt_glucose = ggplot(data = df_prop_diabt_glucose_trans) +
  theme_light() +
  geom_line(aes(x = Time, y = W1, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W2, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W3, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W4, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W5, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W6, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W7, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W8, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W9, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W10, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W11, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W12, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W13, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W14, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W15, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W16, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W17, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W18, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W19, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W20, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
  labs(x = "Days of follow-up", y = "Cumulative statistic", title = "(a)") +
  scale_x_continuous(limits = c(900,9700), breaks = seq(900,9700,800), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-1,1.5), breaks = c(seq(-1,1.5,0.5)), expand = expansion(0,0)) +
  scale_color_manual("", values = c(alpha(c("orange", "blue"), 1)), breaks = c("Simulated", "Observed")) +
  scale_linetype_manual("", values = c("dashed", "solid"), breaks = c("Simulated", "Observed")) +
  scale_linewidth_manual("", values = c(0.5,1), breaks = c("Simulated", "Observed")) +
  annotate("text", x = 6500, y = -0.75, label = "P-value < 0.001", hjust = 0, size = 6, color = "red") +
  theme(legend.position = "none",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_blank(), 
        legend.key.height = unit(17, "pt"),
        legend.key.width = unit(30, "pt"),
        axis.title.x = element_text(margin = margin(t=8)),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        axis.text = element_text(size = 12),
        plot.margin = margin(t=0.5,0.5,0.5,0.5, "lines"),
        aspect.ratio=8/8,
        title = element_text(size = 18))
p_prop_diabt_glucose
save(p_prop_diabt_glucose, file = "./output/Figure4/prop_diabt_glucose_plot.RData")

sessionInfo()