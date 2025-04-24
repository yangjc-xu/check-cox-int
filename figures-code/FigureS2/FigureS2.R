rm(list = ls())
library(tidyverse)
library(ggpubr)
load("./output/FigureS2/prop_quad_n200_graph_plot.RData")
load("./output/FigureS2/prop_quad_n400_graph_plot.RData")
load("./output/FigureS2/prop_quad_n800_graph_plot.RData")
load("./output/FigureS2/prop_mono_n200_graph_plot.RData")
load("./output/FigureS2/prop_mono_n400_graph_plot.RData")
load("./output/FigureS2/prop_mono_n800_graph_plot.RData")


p_prop = ggarrange(p_prop_mono_n200, p_prop_mono_n400, p_prop_mono_n800,
                   p_prop_quad_n200,p_prop_quad_n400,p_prop_quad_n800,
                   ncol = 3, nrow = 2)
ggsave(p_prop, file = "./figures/FigureS2.pdf", width = 16, height = 9, units = "in")
ggsave(p_prop, file = "./figures/FigureS2.eps", width = 16, height = 9, units = "in")

sessionInfo()