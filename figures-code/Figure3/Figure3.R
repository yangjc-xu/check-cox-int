rm(list = ls())
library(tidyverse)
library(ggpubr)
load("./output/Figure3/prop_mono_n200_plot.RData")
load("./output/Figure3/prop_mono_n400_plot.RData")
load("./output/Figure3/prop_mono_n800_plot.RData")
load("./output/Figure3/prop_quad_n200_plot.RData")
load("./output/Figure3/prop_quad_n400_plot.RData")
load("./output/Figure3/prop_quad_n800_plot.RData")

p_ff = ggarrange(p_prop_optimal_n200, p_prop_optimal_n400, p_prop_optimal_n800,
                 p_prop_quad_n200, p_prop_quad_n400, p_prop_quad_n800,
                 ncol = 3, nrow = 2)
ggsave(p_ff, file = "./figures/Figure3.pdf", width = 16, height = 9, units = "in")
ggsave(p_ff, file = "./figures/Figure3.eps", width = 16, height = 9, units = "in")

sessionInfo()