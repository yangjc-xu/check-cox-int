rm(list = ls())
set.seed(13279)
library(dplyr)
library(Rcpp)
sourceCpp("./simulation-code/prop_n800_mono/unireg_indep.cpp")
sourceCpp("./simulation-code/prop_n800_mono/profile_dep_fast.cpp")
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

time_dep_coef = function(t){
  coe = log(t)
  return(coe)
}

hazard = function(beta1, beta2, x1, x2, t){
  h = 0.2*sqrt(t)*exp(beta1 * x1 + beta2 * time_dep_coef(t) * x2)
  return(h)
}

cumhazard = function(Beta1, Beta2, X1, X2, Time){
  H = integrate(hazard, lower = 0, upper = Time, beta1 = Beta1, beta2 = Beta2, x1 = X1, x2 = X2)
  return(H$value)
}

generate_T = function(beta1, beta2, x1, x2){
  eps = 1e-3
  U = runif(n = 1, min = 0, max = 1)
  goal = -log(1-U)
  T_right = 6
  T_left = 0
  if(cumhazard(Beta1 = beta1, Beta2 = beta2, X1 = x1, X2 = x2, Time = T_right) <= goal){
    return(T_right)
  } else{
    while ((T_right - T_left) > eps) {
      T_tmp = (T_left + T_right)/2
      if(cumhazard(Beta1 = beta1, Beta2 = beta2, X1 = x1, X2 = x2, Time = T_tmp) <= goal){
        T_left = T_tmp
      } else{
        T_right = T_tmp
      }
    }
    return(T_right)
  }
}

SimIntData = function(n_sample, r){
  
  N = n_sample
  beta = c(0.5, r)
  Subject_ID = 1:N
  X1 = runif(n = N, min = 0, max = 1)
  X2 = rbinom(n = N, size = 1, prob = c(0.5))
  
  EventTime = rep(NA, N)
  for (i in 1:N) {
    Time = generate_T(beta1 = beta[1], beta2 = beta[2], x1 = X1[i], x2 = X2[i])
    EventTime[i] = Time
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
  ## Model fitting
  X_simulated = model.matrix(~ X1 + X2, SimData)[,-1]
  startT = Sys.time()
  result_simulated = unireg_indep_EM(X_simulated, simpdat_simulated, 0.001, Time_simulated)
  cat("Time difference = ",Sys.time()-startT,"\n")
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

n = 800
rep = 10000
n_simulation = 1
r_seq = 1
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
    t = seq(0, 5, length.out = 51)
    W = numeric(length(t))
    W_simu1 = matrix(data = 0, nrow = rep, ncol = length(t))
    G = matrix(data = rnorm(n = rep*n), nrow = rep, ncol = n)
    
    for (i in 1:length(t)) {
      cat("i =", i, "\n")
      cat("t =", t[i], "\n")
      X_dep_new = SimData %>%
        select(Subject_ID, X1, X2) %>%
        slice(rep(1:n(), each = length(EstRes$time))) %>%
        mutate(Time = rep(EstRes$time, n),
               X3 = ifelse(Time <= t[i], X2, 0)) %>%
        select(X1,X2,X3) %>%
        as.matrix()
      
      res = get_W(beta_new, X_dep_new, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
      W[i] = res$W/sqrt(n)
      sum_unit1 = as.vector(res$Dpl[,3] - 
                              res$Info[3,1:2]%*%
                              solve(EstRes$I_first[1:2,1:2])%*%
                              t(EstRes$eff_score))
      for (j in 1:rep) { W_simu1[j,i] = sum(G[j,] * sum_unit1)/sqrt(n) }
    }
    p1[k] = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))
  }
  
  power1[l] = mean(p1 < 0.05)
}
###################################################################################


df = data.frame(Time = t, W = W)
for (i in 1:20) {
  cname = paste("W", i, sep = "")
  df[[cname]] = W_simu1[i,]
}

p_prop_mono_n800 = ggplot(data = df) +
  theme_light() +
  geom_line(aes(x = Time, y = W1, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W2, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W3, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W4, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W5, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W6, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W7, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W8, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W9, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W10, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W11, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W12, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W13, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W14, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W15, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W16, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W17, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W18, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W19, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W20, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
  geom_line(aes(x = Time, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
  labs(x = "Time", y = "Cumulative statistic", title = "Monotone hazard ratio, n = 800") +
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,0.5), expand = expansion(0,0)) +
  scale_y_continuous(limits = c(-1.2,1.2), breaks = c(seq(-1.2,0,0.4), seq(0.4,1.2,0.4)), expand = expansion(0,0)) +
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
p_prop_mono_n800
save(p_prop_mono_n800, file = "./output/FigureS2/prop_mono_n800_graph_plot.RData")
