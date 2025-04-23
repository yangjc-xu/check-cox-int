rm(list = ls())
library(tidyverse)
load("./output/Figure4/ff_diabt_bmi_logtGlucose_raw.RData")

p_ff_diabt_bmi_logtGlucose = ggplot(data = df_ff_diabt_bmi) +
  theme_light() +
  geom_line(aes(x = BMI, y = W1, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W2, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W3, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W4, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W5, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W6, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W7, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W8, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W9, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W10, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W11, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W12, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W13, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W14, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W15, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W16, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W17, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W18, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W19, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W20, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = BMI, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
  labs(x = "Body mass index", y = "Cumulative statistic", title = "(b)") +
  scale_x_continuous(limits = c(15,55), breaks = seq(15,55,5), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-0.6,0.6), breaks = c(seq(-0.6,0,0.2), seq(0.2,0.6,0.2)), expand = expansion(0,0)) +
  scale_color_manual("", values = c(alpha(c("orange", "blue"), 1)), breaks = c("Simulated", "Observed")) +
  scale_linewidth_manual("", values = c(0.5,1), breaks = c("Simulated", "Observed")) +
  scale_linetype_manual("", values = c("dashed", "solid"), breaks = c("Simulated", "Observed")) +
  annotate("text", x = 40, y = -0.5, label = "P-value = 0.009", hjust = 0, size = 6, color = "red") +
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
p_ff_diabt_bmi_logtGlucose
save(p_ff_diabt_bmi_logtGlucose, file = "./output/Figure4/ff_diabt_bmi_logtGlucose_plot.RData")
