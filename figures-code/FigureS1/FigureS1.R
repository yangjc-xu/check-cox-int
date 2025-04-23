rm(list = ls())
library(tidyverse)
library(ggpubr)
load("./output/FigureS1/ff_quad_n200_graph_plot.RData")
load("./output/FigureS1/ff_quad_n400_graph_plot.RData")
load("./output/FigureS1/ff_quad_n800_graph_plot.RData")
load("./output/FigureS1/ff_twisted_n200_graph_plot.RData")
load("./output/FigureS1/ff_twisted_n400_graph_plot.RData")
load("./output/FigureS1/ff_twisted_n800_graph_plot.RData")


p_ff = ggarrange(p_ff_quad_n200,p_ff_quad_n400,p_ff_quad_n800,
                 p_ff_twisted_n200, p_ff_twisted_n400, p_ff_twisted_n800,
                 ncol = 3, nrow = 2)
ggsave(p_ff, file = "./figures/FigureS1.pdf", width = 16, height = 9, units = "in")
ggsave(p_ff, file = "./figures/FigureS1.eps", width = 16, height = 9, units = "in")
