rm(list = ls())
library(tidyverse)
library(ggpubr)
load("./output/Figure1/quad_plot.RData")
load("./output/Figure1/log_plot.RData")
load("./output/Figure1/ind_plot.RData")
load("./output/Figure1/cubic_plot.RData")
load("./output/Figure1/hr_plot.RData")
load("./output/Figure1/hr_div_plot.RData")
load("./output/Figure1/exp_plot.RData")
load("./output/Figure1/sqrt_plot.RData")

p = ggarrange(p_quad, p_log, 
              p_exp, p_sqrt,
              p_ind, p_cubic, 
              p_hr, p_hr_div,
              ncol = 4, nrow = 2)
ggsave(p, width = 15, height = 7.5, file = "./figures/Figure1.pdf", , units = "in")
ggsave(p, width = 15, height = 7.5, file = "./figures/Figure1.eps", units = "in")


