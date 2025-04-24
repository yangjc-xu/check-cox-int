seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
filename <- paste("simu", seed, ".RData", sep = "")
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
n_simulation = 100
r_seq = seq(0, 1, 0.2)
power1 = numeric(length(r_seq))
power_score_t = numeric(length(r_seq))
power_score_logt = numeric(length(r_seq))
for (l in 1:length(r_seq)) {
  r = r_seq[l]
  cat("r =", r, "\n")
  p1 = numeric(n_simulation)
  p_score_t = numeric(n_simulation)
  p_score_logt = numeric(n_simulation)
  
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
    
    
    e = rep(0,p+1)
    e[p+1] = 1
    H = matrix(c(1,0,0,1,0,0), nrow = p, ncol = p+1)
    ## Score t
    X_dep_score_t = SimData %>%
      select(Subject_ID, X1, X2) %>%
      slice(rep(1:n(), each = length(EstRes$time))) %>%
      mutate(Time = rep(EstRes$time, n),
             X3 = Time * X2) %>%
      select(X1,X2,X3) %>%
      as.matrix()
    
    res_score_t = get_W(beta_new, X_dep_score_t, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
    Stat_t = res_score_t$W / sqrt(n)
    A_t = e - res_score_t$Info[p+1,1:p] %*% solve(EstRes$I_first[1:p,1:p]) %*% H
    Var_t = A_t %*% res_score_t$Info %*% t(A_t)
    Test_t = Stat_t/sqrt(Var_t)
    p_score_t[k] = 2*(1-pnorm(abs(Test_t)))
    
    ## Score logt
    X_dep_score_logt = SimData %>%
      select(Subject_ID, X1, X2) %>%
      slice(rep(1:n(), each = length(EstRes$time))) %>%
      mutate(Time = rep(EstRes$time, n),
             X3 = log(Time) * X2) %>%
      select(X1,X2,X3) %>%
      as.matrix()
    
    res_score_logt = get_W(beta_new, X_dep_score_logt, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
    Stat_logt = res_score_logt$W / sqrt(n)
    A_logt = e - res_score_logt$Info[p+1,1:p] %*% solve(EstRes$I_first[1:p,1:p]) %*% H
    Var_logt = A_logt %*% res_score_logt$Info %*% t(A_logt)
    Test_logt = Stat_logt/sqrt(Var_logt)
    p_score_logt[k] = 2*(1-pnorm(abs(Test_logt)))
  }
  
  power1[l] = mean(p1 < 0.05)
  power_score_t[l] = mean(p_score_t < 0.05)
  power_score_logt[l] = mean(p_score_logt < 0.05)
}
###################################################################################

save(power1, power_score_t, power_score_logt, file = filename)

sessionInfo()

