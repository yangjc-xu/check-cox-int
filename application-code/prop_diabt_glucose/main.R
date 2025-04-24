rm(list = ls())
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("../unireg_indep.cpp")
sourceCpp("../profile_dep_fast.cpp")
load("../result_estimation_diabt_select.RData")

n = nrow(X_aric)
constant = 1
rep = 10000
p = length(result_aric$beta)
beta_new = c(result_aric$beta, 0)

## Age
args = commandArgs(trailingOnly = TRUE)
x = as.numeric(args[1])
counter = as.numeric(args[2])
W = numeric(length(x))
W_simu1 = matrix(data = 0, nrow = rep, ncol = length(x))
G = matrix(data = rnorm(n = rep*n), nrow = rep, ncol = n)

for (i in 1:length(x)) {
  cat("i =", i, "\n")
  cat("x =", x[i], "\n")
  
  X_dep_new = as.data.frame(X_aric) %>%
    slice(rep(1:n(), each = length(Time_aric))) %>%
    mutate(Time = rep(Time_aric, n),
           X_extra = ifelse(Time <= x[i], Glucose, 0)) %>%
    select(-Time) %>%
    as.matrix()
  
  res = get_W(beta_new, X_dep_new, X_aric, result_aric$lambda, simpdat_aric, 0.001, Time_aric, constant)
  W[i] = res$W/sqrt(n)
  
  sum_unit1 = as.vector(res$Dpl[,p+1] - 
                          res$Info[p+1,1:p]%*%
                          solve(result_aric$I1[1:p,1:p])%*%
                          t(result_aric$eff_score))
  for (j in 1:rep) { W_simu1[j,i] = sum(G[j,] * sum_unit1)/sqrt(n) }
}

#pval = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))

save(W, W_simu1, file = paste("prop_diabt_glucose_", counter, ".RData", sep = ""))

sessionInfo()