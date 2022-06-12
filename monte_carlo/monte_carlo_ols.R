rm(list = ls())

library("sketching")
setwd("~/Dropbox/Github/replication-LeeNg-2022-ICML/monte_carlo")

seed <- 53426  
set.seed(seed)  

# Setup

n <- 1e+6
m <- 500
d <- 5
nrep <- 5000
alpha <- 0.05
cv <- qnorm(1 - alpha/2)

mu <- rep(0,d)
rho <- 0.5
Sigma_first_row <- rho^(0:(d-1)) 
Sigma <- stats::toeplitz(Sigma_first_row)
beta <- matrix(rep(1,d), nrow = d, ncol = 1)

true_beta <- beta[d]

heteroset <- c("homo","hetero")

for (i_het in 1:length(heteroset)){ 
  
start_time <- proc.time()    
  
hetero <- heteroset[i_het]  

if (hetero == "hetero"){
  false_beta <- true_beta + 0.4
} else{
  false_beta <- true_beta + 0.1
}   

output_text_file_name <- paste("mc",hetero,"summary","June2022.txt", sep = "_")
save_file_name <- paste("mc",hetero,"output","June2022.Rdata", sep = "_")

sink(file = output_text_file_name, append = FALSE)
options(digits=3)

print("Number of Simulations")
print(nrep)

if (hetero == "hetero"){
  print("heteroskedastic design")
} else{
  print("homoskedastic design")
}   

print("n")
print(n)

print("m")
print(m)

print("Confidence level")
print(1-alpha)

sink()

methods <- c("bernoulli","unif","leverage","countsketch","srht","fft")

test_summary <- {}

results <- array(NA, dim = c(nrep,4,length(methods)))

for (i in 1:nrep){
  
  # generate data
  X <- MASS::mvrnorm(n, mu, Sigma)
  eps <- matrix(rnorm(n), nrow = n, ncol = 1)

  if (hetero == "hetero"){
    het <- exp(X[,d])
  } else{
    het <- 1
  }
  
  Y <- X %*% beta + het*eps

  intercept <- matrix(rep(1,n), nrow = n, ncol = 1)
  fullsample <- cbind(Y,intercept,X)
  
  for (met in 1:length(methods)){
    
    method <- methods[met]
    
  # generate a sketch
    
    if (method == "leverage"){
      subsample <- sketch_leverage(fullsample, m, method = method)
      ys <- subsample$subsample[,1]
      reg <- subsample$subsample[,-1]
      prb <- subsample$prob
    } else{       
      subsample <- sketch(fullsample, m, method = method)
      ys <- subsample[,1]
      reg <- subsample[,-1]
    }
    
    if (method == "leverage"){
      submodel_weight <- 1/(m*prb)
      submodel <- lm(ys ~ reg - 1, weights = submodel_weight) 
    } else{
      submodel <- lm(ys ~ reg - 1) 
    }

  # print(summary(submodel))
  
  ztest_tilde <- lmtest::coeftest(submodel, df = Inf)
  # use White's estimator
  ztest_robust_tilde <- lmtest::coeftest(submodel, df = Inf, vcov = sandwich::vcovHC, type = "HC0")

  size_tilde <- abs(ztest_tilde[(d+1),1] - true_beta)/ztest_tilde[(d+1),2]
  size_robust_tilde <- abs(ztest_robust_tilde[(d+1),1] - true_beta)/ztest_robust_tilde[(d+1),2]
  
  power_tilde <- abs(ztest_tilde[(d+1),1] - false_beta)/ztest_tilde[(d+1),2]
  power_robust_tilde <- abs(ztest_robust_tilde[(d+1),1] - false_beta)/ztest_robust_tilde[(d+1),2]  
  
  results[i,,met] <- c(size_tilde, size_robust_tilde, power_tilde, power_robust_tilde)
  
}

  if (i%%5 == 0){
    
    print(c(hetero,i))
    
  }  
  
  if (i%%20 == 0){
    
    test_tmp <- (results[1:i,,] > 1.96)
    if (length(methods) == 1){
      test_summary_tmp <- t(apply(test_tmp, 2, mean))
    } else{
      test_summary_tmp <- t(apply(test_tmp, c(2,3), mean))
    }
    
    print(test_summary_tmp)
    
  }  
  
}

test <- (results > 1.96)
if (length(methods) == 1){
  test_summary <- t(apply(test, 2, mean))
} else{
  test_summary <- t(apply(test, c(2,3), mean))
}

rownames(test_summary) <- methods
colnames(test_summary) <- c("size:non-robust","size:robust",
                            "power:non-robust","power:robust")

save(test_summary, results, file = save_file_name)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")
test_summary_xtable = xtable::xtable(test_summary, digits=3)

sink(file = output_text_file_name, append = TRUE)
options(digits=3)

print("Simulation Results")
print(test_summary_xtable)

time_taken <- proc.time() - start_time 
print("Time Taken (seconds)")
print(time_taken)
sink()

}

closeAllConnections()


