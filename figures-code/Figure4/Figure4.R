rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./output/Figure4/prop_diabt_glucose_plot.RData")
load("./output/Figure4/ff_diabt_bmi_logtGlucose_plot.RData")

plot = ggarrange(p_prop_diabt_glucose, p_ff_diabt_bmi_logtGlucose,
                 nrow = 1, ncol = 2,
                 label.y = 0.98, label.x = 0.04,
                 font.label = list(size = 18, color = "black", face = "bold", family = NULL))

ggsave(plot, file = "./figures/Figure4.pdf", width = 14, height = 6)
ggsave(plot, file = "./figures/Figure4.eps", width = 14, height = 6)

sessionInfo()