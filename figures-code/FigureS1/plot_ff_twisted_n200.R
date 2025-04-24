rm(list = ls())
set.seed(66324)
library(dplyr)
library(Rcpp)
sourceCpp("./simulation-code/ff_n200_twisted/unireg_indep.cpp")
sourceCpp("./simulation-code/ff_n200_twisted/profile_new.cpp")
##################simulation##############################################

find_time = function(simpdat){
  Time_unique = unique(c(simpdat[which(simpdat[,2] != 0),2], simpdat[which(simpdat[,3] != Inf),3])) %>%
    sort(decreasing = TRUE)
}

rtruncnorm <- function(n, a, b, mu, s) {
  F.a <- pnorm(a, mean = mu, sd = s)
  F.b <- pnorm(b, mean = mu, sd = s)
  u <- runif(n, min = F.a, max = F.b)
  return(qnorm(u, mean = mu, sd = s))
}

powerx = function(x, r){
  trans = sign(x)*abs(x)^r
  return(trans)
}

SimIntData = function(n_sample, r){
  
  N = n_sample
  beta = c(0.5, -0.5)
  Subject_ID = 1:N
  X1 = rbinom(n = N, size = 1, prob = c(0.5))
  X2 = runif(n = N, min = -1, max = 1)
  linear_term = beta[1] * X1 + beta[2] * powerx(X2, r)
  
  EventTime = rep(NA, N)
  for (i in 1:N) {
    Time = rexp(n = 1, rate = exp(linear_term[i]))
    EventTime[i] = 3*sqrt(Time)
  }
  
  ExTime = matrix(NA, nrow = N, ncol = 6)
  Left_Time = rep(NA, N)
  Right_Time = rep(NA, N)
  for (i in 1:N) {
    tmp_ExTime = sort(runif(n = 6, min = 0, max = 5), decreasing = F)
    ExTime[i,1] = tmp_ExTime[1] + runif(1, min = 0, max = 0.1)
    ExTime[i,2] = tmp_ExTime[2] + runif(1, min = 0.1, max = 0.2)
    ExTime[i,3] = tmp_ExTime[3] + runif(1, min = 0.2, max = 0.3)
    ExTime[i,4] = tmp_ExTime[4] + runif(1, min = 0.3, max = 0.4)
    ExTime[i,5] = tmp_ExTime[5] + runif(1, min = 0.4, max = 0.5)
    ExTime[i,6] = tmp_ExTime[6] + runif(1, min = 0.5, max = 0.6)
    
    if(EventTime[i] > ExTime[i,6]){
      Left_Time[i] = ExTime[i,6]
      Right_Time[i] = Inf
    }else if(EventTime[i] < ExTime[i,1]){
      Left_Time[i] = 0
      Right_Time[i] = ExTime[i,1]
    }else{
      index = which(ExTime[i,] >= EventTime[i])[1]
      Left_Time[i] = ExTime[i,index-1]
      Right_Time[i] = ExTime[i,index]
    }
  }
  SimData = data.frame(Subject_ID = Subject_ID,
                       Left_Time = Left_Time,
                       Right_Time = Right_Time,
                       X1 = X1, X2 = X2, Event_Time = EventTime)
  return(SimData)
}

get_est = function(SimData){
  ## Data processing
  n = dim(SimData)[1]
  simpdat_simulated = data.frame(Subject_ID = 1:n, 
                                 SimData[,c(2,3)],
                                 R_star = ifelse(SimData$Right_Time==Inf, SimData$Left_Time, SimData$Right_Time))
  simpdat_simulated$new_ID = 1:nrow(simpdat_simulated)
  simpdat_simulated = simpdat_simulated[order(simpdat_simulated$R_star, decreasing = TRUE),]
  Time_simulated = find_time(simpdat_simulated)
  simpdat_simulated = as.matrix(simpdat_simulated)
  X_simulated = model.matrix(~ X1 + X2, SimData)[,-1]
  ## Model fitting
  startT = Sys.time()
  result_simulated = unireg_indep_EM(X_simulated, simpdat_simulated, 0.001, Time_simulated)
  cat("Time difference = ",Sys.time()-startT,"\n\n")
  beta = result_simulated$beta
  time = Time_simulated
  lambda = result_simulated$lambda
  eff_score = result_simulated$eff_score
  Cov_first = result_simulated$Cov1
  I_first = result_simulated$I1
  
  return(list(beta = beta, lambda = lambda, time = time, X = X_simulated, 
              simpdat = simpdat_simulated, eff_score = eff_score,
              Cov_first = Cov_first, I_first = I_first))
}

n = 200
rep = 10000
n_simulation = 1
r_seq = 0
power1 = numeric(length(r_seq))
for (l in 1:length(r_seq)) {
  r = r_seq[l]
  cat("r =", r, "\n")
  p1 = numeric(n_simulation)
  for (k in 1:n_simulation) {
    cat("k =", k, "\n")
    SimData = SimIntData(n_sample = n, r = r)
    EstRes = get_est(SimData = SimData)
    p = length(EstRes$beta)
    constant = 1
    beta_new = c(EstRes$beta, 0)
    x = seq(-1, 1, length.out = 51)
    #x = sort(unique(SimData$X2))
    W = numeric(length(x))
    W_simu1 = matrix(data = 0, nrow = rep, ncol = length(x))
    G = matrix(data = rnorm(n = rep*n), nrow = rep, ncol = n)
    
    for (i in 1:length(x)) {
      cat("i =", i, "\n")
      cat("x =", x[i], "\n")
      X_extra = ifelse(EstRes$X[,2] <= x[i], 1, 0)
      X_new = cbind(EstRes$X, X_extra)
      
      res = get_W(beta_new, X_new, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
      W[i] = res$W/sqrt(n)
      
      sum_unit1 = as.vector(res$Dpl[,p+1] - 
                              res$Info[p+1,1:p]%*%
                              solve(EstRes$I_first[1:p,1:p])%*%
                              t(EstRes$eff_score))
      for (j in 1:rep) { W_simu1[j,i] = sum(G[j,] * sum_unit1)/sqrt(n) }
    }
    p1[k] = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))
  }
  
  power1[l] = mean(p1 < 0.05)
}
###################################################################################

df = data.frame(Covariate = x, W = W)
for (i in 1:20) {
  cname = paste("W", i, sep = "")
  df[[cname]] = W_simu1[i,]
}

p_ff_twisted_n200 = ggplot(data = df) +
  theme_light() +
  geom_line(aes(x = Covariate, y = W1, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W2, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W3, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W4, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W5, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W6, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W7, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W8, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W9, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W10, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W11, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W12, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W13, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W14, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W15, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W16, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W17, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W18, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W19, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W20, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Covariate, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
  labs(x = "Covariate", y = "Cumulative statistic", title = "Twisted functional form, n = 200") +
  scale_x_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-1.6,1.6), breaks = c(seq(-1.6,0,0.4), seq(0.4,1.6,0.4)), expand = expansion(0,0)) +
  scale_color_manual("", values = c(alpha(c("orange", "blue"), 1)), breaks = c("Simulated", "Observed")) +
  scale_linewidth_manual("", values = c(0.5,1), breaks = c("Simulated", "Observed")) +
  scale_linetype_manual("", values = c("dashed", "solid"), breaks = c("Simulated", "Observed")) +
  theme(legend.position = c(0.85, 0.85),
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
p_ff_twisted_n200
save(p_ff_twisted_n200, file = "./output/FigureS1/ff_twisted_n200_graph_plot.RData")  

sessionInfo()