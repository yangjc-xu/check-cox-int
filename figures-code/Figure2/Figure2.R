rm(list = ls())
library(tidyverse)
library(ggpubr)
load("./output/Figure2/ff_quad_n200_plot.RData")
load("./output/Figure2/ff_quad_n400_plot.RData")
load("./output/Figure2/ff_quad_n800_plot.RData")
load("./output/Figure2/ff_twisted_n200_plot.RData")
load("./output/Figure2/ff_twisted_n400_plot.RData")
load("./output/Figure2/ff_twisted_n800_plot.RData")


p_ff = ggarrange(p_ff_optimal_n200,p_ff_optimal_n400,p_ff_optimal_n800,
                 p_ff_cates_n200, p_ff_cates_n400, p_ff_cates_n800,
                 ncol = 3, nrow = 2)
ggsave(p_ff, file = "./figures/Figure2.pdf", width = 16, height = 9, units = "in")
ggsave(p_ff, file = "./figures/Figure2.eps", width = 16, height = 9, units = "in")

sessionInfo()