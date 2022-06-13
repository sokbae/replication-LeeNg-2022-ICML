# Replication file: empirical example
# Lee and Ng (2022, ICML)
# Date: June 13, 2022

rm(list = ls())

library("sketching")
setwd("~/Dropbox/Github/replication-LeeNg-2022-ICML/empirical_example")

start_time <- proc.time()    

seed <- 98796  
set.seed(seed)  

###########################
### loading the dataset ###
###########################

akdata <- AK # The dataset is included in the "sketching" package.

Y <- akdata$LWKLYWGE
intercept <- akdata$CNST
X_end <- akdata$EDUC
X_exg <- akdata[,3:11]
X <- cbind(X_exg, X_end)
Z_inst <- akdata[,12:(ncol(akdata)-1)]
Z <- cbind(X_exg, Z_inst)

###########################
###      OLS            ###
###########################

fullsample <- cbind(Y,intercept,X)
n <- nrow(fullsample)
d <- ncol(X)

# choice of m (data-oblivious sketch size)
target_size <- 0.05
target_power <- 0.8
S_constant <- (qnorm(1-target_size) + qnorm(target_power))^2
tau_limit <- 10
m_ols <- n*S_constant/tau_limit^2 
m_ols <- floor(m_ols)

# sketching methods for OLS
methods <- c("fullsample","bernoulli","unif","leverage","countsketch","srht","fft")
results_ols <- array(NA, dim = c(length(methods),3))

for (met in 1:length(methods)){
  
  method <- methods[met]
  
  # generate a sketch
  
  if (method == "fullsample"){
    ys <- fullsample[,1]
    reg <- as.matrix(fullsample[,-1])
  } else if (method == "leverage"){
    subsample <- sketch_leverage(fullsample, m_ols, method = method)
    ys <- subsample$subsample[,1]
    reg <- subsample$subsample[,-1]
    prb <- subsample$prob
  } else{       
    subsample <- sketch(fullsample, m_ols, method = method)
    ys <- subsample[,1]
    reg <- subsample[,-1]
  }
  
  if (method == "leverage"){
    submodel_weight <- 1/(m_ols*prb)
    submodel <- lm(ys ~ reg - 1, weights = submodel_weight) 
  } else{
    submodel <- lm(ys ~ reg - 1) 
  }
  
  # print(summary(submodel))
  
  ztest_tilde <- lmtest::coeftest(submodel, df = Inf)
    est_tilde <- ztest_tilde[(d+1),1] 
     se_tilde <- ztest_tilde[(d+1),2]

  # use HC estimator
  ztest_robust_tilde <- lmtest::coeftest(submodel, df = Inf, vcov = sandwich::vcovHC, type = "HC0")
         se_hc_tilde <- ztest_robust_tilde[(d+1),2]
  
  results_ols[met,] <- c(est_tilde, se_tilde, se_hc_tilde)
  
}

rownames(results_ols) <- methods
colnames(results_ols) <- c("est", "non-robust se","robust se")

###########################
###       2SLS          ###
###########################

fullsample <- cbind(Y,intercept,X,intercept,Z)
n <- nrow(fullsample)
p <- ncol(X)
q <- ncol(Z)

# choice of m (data-oblivious sketch size)
target_size <- 0.05
target_power <- 0.8
S_constant <- (qnorm(1-target_size) + qnorm(target_power))^2
tau_limit <- 5
m_2sls <- n*S_constant/tau_limit^2 
m_2sls <- floor(m_2sls)

# sketching methods for 2SLS
methods <- c("fullsample","bernoulli","unif","countsketch","srht","fft")
results_2sls <- array(NA, dim = c(length(methods),3))

for (met in 1:length(methods)){
  
  method <- methods[met]
  
  # generate a sketch
  
  if (method == "fullsample"){
      ys <- fullsample[,1]
      reg <- as.matrix(fullsample[,2:(p+2)])
      inst <- as.matrix(fullsample[,(p+3):ncol(fullsample)]) 
  } else{
    subsample <- sketch(fullsample, m_2sls, method = method)
      ys <- subsample[,1]
      reg <- subsample[,2:(p+2)]
      inst <- subsample[,(p+3):ncol(subsample)] 
  }

  submodel <- ivreg::ivreg(ys ~ reg - 1 | inst - 1) 
  
  # print(summary(submodel))
  
  ztest_tilde <- lmtest::coeftest(submodel, df = Inf)
  est_tilde <- ztest_tilde[(d+1),1] 
  se_tilde <- ztest_tilde[(d+1),2]
  
  # use HC estimator
  ztest_robust_tilde <- lmtest::coeftest(submodel, df = Inf, vcov = sandwich::vcovHC, type = "HC0")
  se_hc_tilde <- ztest_robust_tilde[(d+1),2]
  
  results_2sls[met,] <- c(est_tilde, se_tilde, se_hc_tilde)
  
}

rownames(results_2sls) <- methods
colnames(results_2sls) <- c("est", "non-robust se","robust se")

########################################
###  Saving and Printing the results ###
########################################

# setup
output_text_file_name <- paste("empirical_example","summary","June2022.txt", sep = "_")
save_file_name <- paste("empirical_example","output","June2022.Rdata", sep = "_")

save(results_ols, results_2sls, file = save_file_name)

sink(file = output_text_file_name, append = FALSE)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")
ols_xtable = xtable::xtable(results_ols, digits=5)
tsls_xtable = xtable::xtable(results_2sls, digits=3)

print("m for OLS")
print(m_ols)

print("m for 2SLS")
print(m_2sls)

print("OLS")
print(ols_xtable)

print("2SLS")
print(tsls_xtable)

time_taken <- proc.time() - start_time 
print("Time Taken (seconds)")
print(time_taken)
sink()

closeAllConnections()
