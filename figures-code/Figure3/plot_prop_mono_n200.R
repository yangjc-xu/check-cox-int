rm(list = ls())
library(tidyverse)
load("./output/Figure3/prop_mono_n200_raw.RData")

p_prop_optimal_n200 = ggplot(data = df, aes(x = r, y = power, color = Test, shape = Test)) +
  geom_point(size = 2) +
  geom_line() + 
  geom_line(data = df, aes(x=r, y=power)) +
  labs(title = "Monotone hazard ratio, n = 200") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,length.out=6), name = "Empirical rejection rate")+
  scale_color_discrete(name = "Method", labels = c("Proposed" = "Proportionality test", 
                                                   "Score_logt" = TeX(r'(Score test for $X_1\log t$)'), 
                                                   "Score_t" = TeX(r'(Score test for $X_1t$)')))+
  scale_shape_discrete(name = "Method", labels = c("Proposed" = "Proportionality test", 
                                                   "Score_logt" = TeX(r'(Score test for $X_1\log t$)'), 
                                                   "Score_t" = TeX(r'(Score test for $X_1t$)')))+
  geom_hline(yintercept = 0.05, color = "red",linetype="dashed")+
  theme_light() +
  theme(legend.position = c(0.75, 0.2),
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
p_prop_optimal_n200

save(p_prop_optimal_n200, file = "./output/Figure3/prop_mono_n200_plot.RData")

sessionInfo()
  