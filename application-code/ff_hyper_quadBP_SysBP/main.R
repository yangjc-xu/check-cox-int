rm(list = ls())
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("../unireg_indep.cpp")
sourceCpp("../profile_new.cpp")
load("../result_estimation_hyper_quadBP_select.RData")

n = nrow(X_aric)
constant = 1
rep = 10000
p = length(result_aric$beta)
beta_new = c(result_aric$beta, 0)

## SysBP
args = commandArgs(trailingOnly = TRUE)
x = as.numeric(args[1])
counter = as.numeric(args[2])
W = numeric(length(x))
W_simu1 = matrix(data = 0, nrow = rep, ncol = length(x))
G = matrix(data = rnorm(n = rep*n), nrow = rep, ncol = n)

for (i in 1:length(x)) {
  cat("i =", i, "\n")
  cat("x =", x[i], "\n")
  X_extra = ifelse(X_aric[,9] <= x[i], 1, 0)
  X_new = cbind(X_aric, X_extra)
  
  res = get_W(beta_new, X_new, result_aric$lambda, simpdat_aric, 0.001, Time_aric, constant)
  W[i] = res$W/sqrt(n)
  
  sum_unit1 = as.vector(res$Dpl[,p+1] - 
                          res$Info[p+1,1:p]%*%
                          solve(result_aric$I1[1:p,1:p])%*%
                          t(result_aric$eff_score))
  for (j in 1:rep) { W_simu1[j,i] = sum(G[j,] * sum_unit1)/sqrt(n) }
}

#pval = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))

save(W, W_simu1, file = paste("ff_hyper_SysBP_", counter, ".RData", sep = ""))
