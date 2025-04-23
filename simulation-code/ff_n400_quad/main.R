seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
filename <- paste("simu", seed, ".RData", sep = "")
library(dplyr)
library(Rcpp)
sourceCpp("./unireg_indep.cpp")
sourceCpp("./profile_new.cpp")
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
  X2 = rtruncnorm(n = N, a = 0, b = 3, mu = 1.5, s = 1) #True
  linear_term = beta[1] * X1 + beta[2] * X2 + r * X2^2       
  
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

score_test = function(EstRes, X_extra, e, H, c){
  
  p = length(EstRes$beta)
  X_new_score = cbind(EstRes$X, X_extra)
  beta_initial = c(EstRes$beta, 0)
  lambda_initial = EstRes$lambda
  res_score = get_W(beta_initial, X_new_score, lambda_initial, EstRes$simpdat, 0.001, EstRes$time, c)
  
  Stat = res_score$W/sqrt(n)
  A = e - res_score$Info[p+1,1:p] %*% solve(EstRes$I_first[1:p,1:p]) %*% H
  Var = A %*% res_score$Info %*% t(A)
  Test = Stat/sqrt(Var)
  p = 2*(1-pnorm(abs(Test)))
  p_chi = 1 - pchisq(Test^2, df = 1)
  
  return(list(Stat = Stat, Var = Var, p = p, p_chi = p_chi))
}

n = 400
rep = 10000
n_simulation = 100
r_seq = seq(0,1,0.1)
constant = 1
power1 = numeric(length(r_seq))
power_score = numeric(length(r_seq))
power_score_log = numeric(length(r_seq))
for (l in 1:length(r_seq)) {
  r = r_seq[l]
  cat("r =", r, "\n")
  p1 = numeric(n_simulation)
  p_score = numeric(n_simulation)
  p_score_log = numeric(n_simulation)
  for (k in 1:n_simulation) {
    cat("k =", k, "\n")
    SimData = SimIntData(n_sample = n, r = r)
    EstRes = get_est(SimData = SimData)
    p = length(EstRes$beta)
    
    beta_new = c(EstRes$beta, 0)
    x = seq(0, 3, length.out = 51)
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
    
    
    e = rep(0,p+1)
    e[p+1] = 1
    H = matrix(c(1,0,0,1,0,0), nrow = p, ncol = p+1)
    ## Score test quad
    test_quad = score_test(EstRes = EstRes, X_extra = SimData$X2^2, e = e, H = H, c = constant)
    p_score[k] = test_quad$p
    
    ## Score test log
    test_log = score_test(EstRes = EstRes, X_extra = log(SimData$X2), e = e, H = H, c = constant)
    p_score_log[k] = test_log$p
  }
  
  power1[l] = mean(p1 < 0.05)
  power_score[l] = mean(p_score < 0.05)
  power_score_log[l] = mean(p_score_log < 0.05)
}
###################################################################################

save(power1, power_score, power_score_log, file = filename)



