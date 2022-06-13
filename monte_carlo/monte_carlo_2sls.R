# Replication file: Monte Carlo experiments for OLS
# Lee and Ng (2022, ICML)
# Date: June 13, 2022
# Note: The number of replication is set as "nrep <- 5" for checking the code,
#       but it needs to be changed to "nrep <- 1000"  for replication purposes.

rm(list = ls())

library("sketching")
setwd("~/Dropbox/Github/replication-LeeNg-2022-ICML/monte_carlo")

seed <- 53426  
set.seed(seed)  

# Setup

n <- 1e+6
m <- 500
p <- 5
q <- 20
nrep <- 5 # This needs to be changed to "nrep <- 1000"  for replication.
alpha <- 0.05
cv <- qnorm(1 - alpha/2)

mu <- rep(0,q)
rho <- 0.5
Sigma_first_row <- rho^(0:(q-1)) 
Sigma <- stats::toeplitz(Sigma_first_row)
beta <- matrix(rep(1,p), nrow = p, ncol = 1)

true_beta <- beta[p]

heteroset <- c("hetero","homo")

for (i_het in 1:length(heteroset)){ 
  
start_time <- proc.time()    
  
hetero <- heteroset[i_het]  

if (hetero == "hetero"){
  false_beta <- true_beta + 0.10
} else{
  false_beta <- true_beta + 0.05
}   

output_text_file_name <- paste("mc_2sls",hetero,"summary","June2022.txt", sep = "_")
save_file_name <- paste("mc_2sls",hetero,"output","June2022.Rdata", sep = "_")

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

methods <- c("bernoulli","unif","countsketch","srht","fft")

test_summary <- {}

results <- array(NA, dim = c(nrep,4,length(methods)))

for (i in 1:nrep){
  
  # generate data
  Z <- MASS::mvrnorm(n, mu, Sigma)
  V <- matrix(rnorm(n), nrow = n, ncol = 1)
  X_exg <- Z[,1:(p-1)]
  
  if (hetero == "hetero"){
    het <- exp(apply(5*abs(Z),1,mean))/100
  } else{
    het <- 1
  }
  
  X1 <- 0.1*apply(X_exg,1,sum) + 0.5*apply(Z[,p:q],1,sum) + V
  X <- cbind(X_exg, X1)
  
  eps <- matrix(rnorm(n), nrow = n, ncol = 1)
  Y <- X %*% beta + het*(V + eps)

  intercept <- matrix(rep(1,n), nrow = n, ncol = 1)
  fullsample <- cbind(Y,intercept,X,intercept,Z)
  
  for (met in 1:length(methods)){
    
    method <- methods[met]
    
    # generate a sketch
    
    subsample <- sketch(fullsample, m, method = method)
      ys <- subsample[,1]
      reg <- subsample[,2:(p+2)]
      inst <- subsample[,(p+3):ncol(subsample)] 
      
    # estimate 2sls
      
    submodel <- ivreg::ivreg(ys ~ reg - 1 | inst - 1) 

   # print(summary(submodel))
  
  ztest_tilde <- lmtest::coeftest(submodel, df = Inf)
  # use HC-robust estimator
  ztest_robust_tilde <- lmtest::coeftest(submodel, df = Inf, vcov = sandwich::vcovHC, type = "HC0")

  size_tilde <- abs(ztest_tilde[(p+1),1] - true_beta)/ztest_tilde[(p+1),2]
  size_robust_tilde <- abs(ztest_robust_tilde[(p+1),1] - true_beta)/ztest_robust_tilde[(p+1),2]
  
  power_tilde <- abs(ztest_tilde[(p+1),1] - false_beta)/ztest_tilde[(p+1),2]
  power_robust_tilde <- abs(ztest_robust_tilde[(p+1),1] - false_beta)/ztest_robust_tilde[(p+1),2]  
  
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

test <- (results > cv)
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


