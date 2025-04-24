rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./output/Figure5/ff_hyper_SysBP_plot.RData")
load("./output/Figure5/ff_hyper_DiaBP_plot.RData")


plot = ggarrange(p_ff_hyper_SysBP, p_ff_hyper_DiaBP,
                 nrow = 1, ncol = 2,
                 label.y = 0.98, label.x = 0.04,
                 font.label = list(size = 18, color = "black", face = "bold", family = NULL))

ggsave(plot, file = "./figures/Figure5.pdf", width = 14, height = 6)
ggsave(plot, file = "./figures/Figure5.eps", width = 14, height = 6)

sessionInfo()
