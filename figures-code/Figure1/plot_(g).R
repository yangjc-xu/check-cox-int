rm(list = ls())
library(tidyverse)
library(latex2exp)
load("./output/Figure1/hr_raw.RData")
df = data.frame(y = colMeans(W_score_all), x = x)
intervals <- list(c(0, 5))
smooth_piecewise <- function(data, intervals, span = 0.5) {
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
p_hr = ggplot() +
  theme_light() +
  geom_line(data = smoothed_data,aes(x = x, y = y))+
  labs(x = "Time", y = "Cumulative statistic", title = "(g)")+
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,1), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(0,0.7), breaks = seq(0,0.8,0.2), expand = expansion(0,0))+
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
p_hr
save(p_hr, file = "./output/Figure1/hr_plot.RData")
