rm(list = ls())
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../unireg_indep.cpp")

load("../dat_hyper_splitted.RData")
find_time = function(simpdat){
  Time_unique = unique(c(simpdat[which(simpdat[,2] != 0),2], simpdat[which(simpdat[,3] != Inf),3])) %>%
    sort(decreasing = TRUE)
}
## Data processing
n = dim(dat_select)[1]
simpdat_aric = data.frame(Subject_ID = dat_select$Subject_ID, 
                          dat_select[,c(3,4)],
                          R_star = ifelse(dat_select$Right_Time==Inf, dat_select$Left_Time, dat_select$Right_Time))
simpdat_aric$new_ID = dat_select$Subject_ID
simpdat_aric = simpdat_aric[order(simpdat_aric$R_star, decreasing = TRUE),]

Time_aric = find_time(simpdat_aric)
simpdat_aric = as.matrix(simpdat_aric)

## Standardization
X_aric = model.matrix(~ Center+Age+Sex+Race+BMI+Glucose+SysBP+DiaBP+I(SysBP^2)+I(DiaBP^2), dat_select)[,-1]
### Age
std_Age = scale(X_aric[,4])
X_aric[,4] = std_Age
### BMI
std_BMI = scale(X_aric[,7])
X_aric[,7] = std_BMI
### Glucose
std_Glucose = scale(X_aric[,8])
X_aric[,8] = std_Glucose
### SysBP
std_SysBP = scale(X_aric[,9])
X_aric[,9] = std_SysBP
### DiaBP
std_DiaBP = scale(X_aric[,10])
X_aric[,10] = std_DiaBP
### SysBP^2
std_SysBP_sq = scale(X_aric[,11])
X_aric[,11] = std_SysBP_sq
### DiaBP^2
std_DiaBP_sq = scale(X_aric[,12])
X_aric[,12] = std_DiaBP_sq

save(std_Age, std_BMI, std_Glucose, std_SysBP, std_DiaBP, std_SysBP_sq, std_DiaBP_sq, file = "../scale_hyper_quadBP_select.RData")

## Model fitting
cat("Model fitting!\n")
startT = Sys.time()
result_aric = unireg_indep_EM(X_aric, simpdat_aric, 0.001, Time_aric)
cat("Time difference = ",Sys.time()-startT,"\n\n")

# index = c(4,7,8,9,10,11,12)
# attr = c(attr(std_Age,"scaled:scale"),
#          attr(std_BMI,"scaled:scale"),
#          attr(std_Glucose,"scaled:scale"),
#          attr(std_SysBP,"scaled:scale"),
#          attr(std_DiaBP,"scaled:scale"),
#          attr(std_SysBP_sq,"scaled:scale"),
#          attr(std_DiaBP_sq,"scaled:scale"))
# ## beta estimate
# for (i in 1:length(index)) {
#   result_aric$beta[index[i],1] = result_aric$beta[index[i],1]/attr[i]
# }
# ## covariance
# for (i in 1:length(index)) {
#   result_aric$Cov1[index[i],] = result_aric$Cov1[index[i],]/attr[i]
#   result_aric$Cov1[,index[i]] = result_aric$Cov1[,index[i]]/attr[i]
# }

save(result_aric, X_aric, simpdat_aric, Time_aric, file = "result_estimation_hyper_quadBP_select.RData")
