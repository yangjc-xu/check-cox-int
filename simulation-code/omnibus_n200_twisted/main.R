seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
filename <- paste("simu", seed, ".RData", sep = "")
library(dplyr)
library(Rcpp)
sourceCpp("./unireg_indep.cpp")
sourceCpp("./profile_dep_fast.cpp")
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
n_simulation = 100
r_seq = seq(1,0,-0.2)
constant = 1
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
    
    beta_new = c(EstRes$beta, 0)
    
    x1 = c(0.5, 1)
    x2 = seq(-0.6, 0.6, length.out = 4)
    t = seq(1, 5, length.out = 5)
    
    M = length(x1)*length(x2)*length(t)
    W = numeric(M)
    W_simu1 = matrix(data = 0, nrow = rep, ncol = M)
    G = matrix(data = rnorm(n = rep*n), nrow = rep, ncol = n)
    
    for (i in 1:length(x1)) {
      for (j in 1:length(x2)) {
        for (q in 1:length(t)) {
          cat("i =", i, "\n")
          cat("x1 =", x1[i], "\n")
          
          cat("j =", j, "\n")
          cat("x2 =", x2[j], "\n")
          
          cat("q =", q, "\n")
          cat("t =", t[q], "\n")
          
          X_dep_new = SimData %>%
            select(Subject_ID, X1, X2) %>%
            slice(rep(1:n(), each = length(EstRes$time))) %>%
            mutate(Time = rep(EstRes$time, n),
                   X3 = ifelse(X1 <= x1[i] & X2 <= x2[j] & Time <= t[q], 1, 0)) %>%
            select(X1,X2,X3) %>%
            as.matrix()
          
          res = get_W(beta_new, X_dep_new, EstRes$X, EstRes$lambda, EstRes$simpdat, 0.001, EstRes$time, constant)
          
          ind_W = (i-1)*length(x2)*length(t) + (j-1)*length(t) + q
          W[ind_W] = res$W/sqrt(n)
          sum_unit1 = as.vector(res$Dpl[,3] - 
                                  res$Info[3,1:2]%*%
                                  solve(EstRes$I_first[1:2,1:2])%*%
                                  t(EstRes$eff_score))
          for (q in 1:rep) { W_simu1[q,ind_W] = sum(G[q,] * sum_unit1)/sqrt(n) }
        }
      }
    }
    p1[k] = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))
  }
  
  power1[l] = mean(p1 < 0.05)
}
###################################################################################

save(power1, file = filename)
