rm(list = ls())
library(tidyverse)
load("./output/Figure5/ff_hyper_DiaBP_raw.RData")

p_ff_hyper_DiaBP = ggplot(data = df_ff_hyper_DiaBP) +
  theme_light() +
  geom_line(aes(x = DiaBP, y = W1, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W2, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W3, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W4, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W5, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W6, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W7, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W8, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W9, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W10, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W11, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W12, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W13, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W14, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W15, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W16, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W17, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W18, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W19, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W20, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = DiaBP, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
  labs(x = "Diastolic blood pressure", y = "Cumulative statistic", title = "(b)") +
  scale_x_continuous(limits = c(40,90), breaks = seq(40,90,10), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-0.8,0.8), breaks = c(seq(-0.8,0,0.2), seq(0.2,0.8,0.2)), expand = expansion(0,0)) +
  scale_color_manual("", values = c(alpha(c("orange", "blue"), 1)), breaks = c("Simulated", "Observed")) +
  scale_linewidth_manual("", values = c(0.5,1), breaks = c("Simulated", "Observed")) +
  scale_linetype_manual("", values = c("dashed", "solid"), breaks = c("Simulated", "Observed")) +
  annotate("text", x = 45, y = -0.7, label = "P-value = 0.005", hjust = 0, size = 7, color = "red") +
  theme(legend.position = "none",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_blank(), 
        legend.key.height = unit(17, "pt"),
        legend.key.width = unit(30, "pt"),
        axis.title.x = element_text(margin = margin(t=8)),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        axis.text = element_text(size = 15),
        plot.margin = margin(t=0.5,0.5,0.5,0.5, "lines"),
        aspect.ratio=8/8,
        title = element_text(size = 18))
p_ff_hyper_DiaBP
save(p_ff_hyper_DiaBP, file = "./output/Figure5/ff_hyper_DiaBP_plot.RData")

sessionInfo()