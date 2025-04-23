rm(list = ls())
library(tidyverse)
library(latex2exp)
load("./output/Figure2/ff_twisted_n800_raw.RData")

p_ff_cates_n800 = ggplot(data = df, aes(x = r, y = power, color = Test, shape = Test)) +
  geom_point(size = 2) +
  geom_line() + 
  labs(title = "Twisted functional form, n = 800") +
  scale_x_continuous(limits = c(0,1), breaks =  seq(0,1,0.2)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,length.out=6), name = "Empirical rejection rate")+
  scale_color_discrete(name = "Method", labels = c("Proposed_ff" = "Functional form test", 
                                                   "Score" = TeX(r'(Score test for $X_2^2$)'), 
                                                   "Proposed_omn" = "Omnibus test"))+
  scale_shape_discrete(name = "Method", labels = c("Proposed_ff" = "Functional form test", 
                                                   "Score" = TeX(r'(Score test for $X_2^2$)'), 
                                                   "Proposed_omn" = "Omnibus test")) +
  geom_hline(yintercept = 0.05, color = "red",linetype="dashed")+
  theme_light() +
  theme(legend.position = c(0.2, 0.8),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.title = element_blank(),
        legend.background = element_blank(), 
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(20, "pt"),
        axis.title.x = element_text(size = 15, margin = margin(t=8)),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        plot.margin = margin(t=1,1,1,1, "lines"),
        aspect.ratio=6/8,
        title = element_text(size = 15))
p_ff_cates_n800

save(p_ff_cates_n800, file = "./output/Figure2/ff_twisted_n800_plot.RData")



