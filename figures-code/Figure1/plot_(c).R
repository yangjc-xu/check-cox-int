rm(list = ls())
library(tidyverse)
library(latex2exp)
load("./output/Figure1/exp_raw.RData")
df = data.frame(y = colMeans(W_score_all), x = x)
intervals <- list(c(0, 3))
smooth_piecewise <- function(data, intervals, span = 0.4) {
  smoothed_data <- data.frame(x = numeric(), y = numeric())
  
  for (interval in intervals) {
    subset_data <- data %>% filter(x >= interval[1] & x < interval[2])
    loess_model <- loess(y ~ x, data = subset_data, span = span)
    smoothed_values <- predict(loess_model, subset_data$x)
    smoothed_data <- rbind(smoothed_data, data.frame(x = subset_data$x, y = smoothed_values))
  }
  
  return(smoothed_data)
}
smoothed_data <- smooth_piecewise(df, intervals)

p_exp = ggplot() +
  theme_light() +
  geom_line(data = smoothed_data, aes(x = x, y = y))+
  labs(x = "Covariate", y = "Cumulative statistic", title = "(c)")+
  scale_x_continuous(limits = c(0,3), breaks = seq(0,3,0.5), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-0.5,0.4), breaks = c(seq(-0.4,0,0.2), seq(0.2,0.4,0.2)), expand = expansion(0,0))+
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
        title = element_text(size = 15))
p_exp
save(p_exp, file = "./output/Figure1/exp_plot.RData")

sessionInfo()