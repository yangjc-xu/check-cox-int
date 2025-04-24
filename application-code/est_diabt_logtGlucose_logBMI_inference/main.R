rm(list = ls())
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../unireg_dep.cpp")

load("../dat_diabt_splitted.RData")
find_time = function(simpdat){
  Time_unique = unique(c(simpdat[which(simpdat[,2] != 0),2], simpdat[which(simpdat[,3] != Inf),3])) %>%
    sort(decreasing = TRUE)
}
## Data processing
n = dim(dat_inference)[1]
simpdat_aric = data.frame(Subject_ID = dat_inference$Subject_ID, 
                          dat_inference[,c(3,4)],
                          R_star = ifelse(dat_inference$Right_Time==Inf, dat_inference$Left_Time, dat_inference$Right_Time))
simpdat_aric$new_ID = dat_inference$Subject_ID
simpdat_aric = simpdat_aric[order(simpdat_aric$R_star, decreasing = TRUE),]

Time_aric = find_time(simpdat_aric)
simpdat_aric = as.matrix(simpdat_aric)

## Log Transformation
dat_inference$BMI = log(dat_inference$BMI)

## Standardization
X_aric = model.matrix(~ Center+Age+Sex+Race+BMI+Glucose+SysBP+DiaBP, dat_inference)[,-1]
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

save(std_Age, std_BMI, std_Glucose, std_SysBP, std_DiaBP, file = "../scale_diabt_logtGlucose_logBMI_inference.RData")

X_dep_new = as.data.frame(X_aric) %>%
  slice(rep(1:n(), each = length(Time_aric))) %>%
  mutate(Time = rep(Time_aric, n),
         X_extra1 = Glucose * log(Time)) %>%
  select(-Time) %>%
  as.matrix()

## Model fitting
cat("Model fitting!\n")
startT = Sys.time()
result_aric = unireg_dep_EM(X_dep_new, simpdat_aric, 0.001, Time_aric)
cat("Time difference = ",Sys.time()-startT,"\n\n")

load("../result_estimation_diabt_logtGlucose_logBMI_inference.RData")
index = c(4,7,8,9,10,11)
attr = c(attr(std_Age,"scaled:scale"),
         attr(std_BMI,"scaled:scale"),
         attr(std_Glucose,"scaled:scale"),
         attr(std_SysBP,"scaled:scale"),
         attr(std_DiaBP,"scaled:scale"),
         attr(std_Glucose,"scaled:scale"))
## beta estimate
for (i in 1:length(index)) {
  result_aric$beta[index[i],1] = result_aric$beta[index[i],1]/attr[i]
}
## covariance
for (i in 1:length(index)) {
  result_aric$Cov1[index[i],] = result_aric$Cov1[index[i],]/attr[i]
  result_aric$Cov1[,index[i]] = result_aric$Cov1[,index[i]]/attr[i]
}
result_aric$beta
sqrt(diag(result_aric$Cov1))
save(result_aric, X_aric, X_dep_new, simpdat_aric, Time_aric, file = "result_estimation_diabt_logtGlucose_logBMI_inference.RData")

sessionInfo()