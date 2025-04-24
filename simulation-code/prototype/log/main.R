rm(list = ls())
seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
library(dplyr)
library(Rcpp)
sourceCpp("./simulation-code/prototype/log/unireg_indep.cpp")
sourceCpp("./simulation-code/prototype/log/profile_new.cpp")
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

SimIntData = function(n_sample, r){
  
  N = n_sample
  beta = c(0.5, -0.5)
  Subject_ID = 1:N
  X1 = rbinom(n = N, size = 1, prob = c(0.5))
  X2 = runif(n = N, min = 0.5, max = 1.5)
  linear_term = r * log(X2)       
  
  EventTime = rep(NA, N)
  for (i in 1:N) {
    Time = rexp(n = 1, rate = exp(linear_term[i]))
    EventTime[i] = 4*sqrt(Time)
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
  X_simulated = as.matrix(model.matrix(~ X2, SimData)[,-1])
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
  
  return(list(beta = beta, lambda = lambda, cumulative_lambda = result_simulated$cumulative_lambda, time = time, X = X_simulated, 
              simpdat = simpdat_simulated, eff_score = eff_score,
              Cov_first = Cov_first, I_first = I_first))
}

n = 200
r = 1
n_simulation = 100
constant = 1
x = seq(0.5, 1.5, length.out = 51)
t = seq(0, 5, length.out = 1000)
W_score_all = matrix(0, nrow = n_simulation, ncol = length(x))
W_wald_all = matrix(0, nrow = n_simulation, ncol = length(x))
beta_est = numeric(n_simulation)
cumulative_lambda_est = matrix(0, nrow = n_simulation, ncol = length(t))

for (k in 1:n_simulation) {
  SimData = SimIntData(n_sample = n, r = r)
  EstRes = get_est(SimData = SimData)
  beta_est[k] = EstRes$beta
  for (q in 1:length(t)) {
    cumulative_lambda_est[k,q] = sum(EstRes$lambda[which(EstRes$time<=t[q])])
  }
  
  ## Score
  beta_new = c(EstRes$beta, 0)
  for (i in 1:length(x)) {
    cat("i =", i, "\n")
    cat("x =", x[i], "\n")
    X_extra = ifelse(EstRes$X[,1] <= x[i], 1, 0)
    X_new = cbind(EstRes$X, X_extra)
    
    res_score = get_W(beta_new, X_new, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
    W_score_all[k,i] = res_score$W/sqrt(n)
  }
  
  ## Wald
  for (i in 1:length(x)) {
    cat("i =", i, "\n")
    cat("x =", x[i], "\n")
    X_extra = ifelse(EstRes$X[,1] <= x[i], 1, 0)
    
    res_wald = get_wald_W(X_extra, EstRes$beta, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
    W_wald_all[k,i] = res_wald$gamma
  }
}

save(W_score_all, x, t, file = "./output/Figure1/log_raw.RData")

sessionInfo()