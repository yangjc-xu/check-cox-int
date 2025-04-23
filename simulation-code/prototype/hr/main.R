rm(list = ls())
seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
library(dplyr)
library(Rcpp)
sourceCpp("./simulation-code/prototype/hr/unireg_indep.cpp")
sourceCpp("./simulation-code/prototype/hr/profile_dep_fast.cpp")
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
  beta = c(0, r)
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
r = -0.8
n_simulation = 100
constant = 1
x = seq(0, 5, length.out = 51)
t = seq(0, 5, length.out = 1000)
W_score_all = matrix(0, nrow = n_simulation, ncol = length(x))
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
    
    X_dep_new = SimData %>%
      select(Subject_ID, X2) %>%
      slice(rep(1:n(), each = length(EstRes$time))) %>%
      mutate(Time = rep(EstRes$time, n),
             X3 = ifelse(Time <= x[i], X2, 0)) %>%
      select(X2,X3) %>%
      as.matrix()
    
    res_score = get_W(beta_new, X_dep_new, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
    W_score_all[k,i] = res_score$W/sqrt(n)
  }
}

save(W_score_all, x, t, file = "./output/Figure1/hr_raw.RData")
