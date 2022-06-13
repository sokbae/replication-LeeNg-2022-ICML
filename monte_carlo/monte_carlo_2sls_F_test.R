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

testset <- c("size","power")
heteroset <- c("hetero","homo")


for (i_test in 1:length(testset)){ 

for (i_het in 1:length(heteroset)){ 
  
start_time <- proc.time()   

testtype <- testset[i_test] 
  
hetero <- heteroset[i_het]  

output_text_file_name <- paste("mc_2sls_F_test",testtype,hetero,"summary","June2022.txt", sep = "_")
save_file_name <- paste("mc_2sls_F_test",testtype,hetero,"output","June2022.Rdata", sep = "_")

sink(file = output_text_file_name, append = FALSE)
options(digits=3)

print("Number of Simulations")
print(nrep)

if (testtype == "size"){
  print("size")
} else{
  print("power")
}   

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

results <- array(NA, dim = c(nrep,2,length(methods)))

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
  
  if (testtype == "size"){
    X1 <- 0.1*apply(X_exg,1,sum) + het*V
  } else{
    X1 <- 0.1*apply(X_exg,1,sum) + 0.05*apply(Z[,p:q],1,sum) + het*V
  }   
  
  X <- cbind(X_exg, X1)
  intercept <- matrix(rep(1,n), nrow = n, ncol = 1)
  fullsample <- cbind(intercept,X,intercept,Z)
  
  for (met in 1:length(methods)){
    
    method <- methods[met]
    
    # generate a sketch
    
      subsample <- sketch(fullsample, m, method = method)
      reg <- subsample[,1:(p+1)]
      inst <- subsample[,(p+2):ncol(subsample)] 
      endog_reg <- reg[,(p+1)]
      exog_reg <- inst[,1:p]
      excluded_inst <- inst[,(p+1):ncol(inst)]
    
    # first-stage F test
      
    first_stage <- lm(endog_reg ~ exog_reg + excluded_inst - 1)  
    first_stage_restricted <- lm(endog_reg ~ exog_reg - 1) 
    F_test <- lmtest::waldtest(first_stage_restricted, first_stage, 
                             test = "F")
    F_test_hc <- lmtest::waldtest(first_stage_restricted, first_stage, 
                               test = "F", vcov = sandwich::vcovHC)
    
    results[i,,met] <- c(F_test$`Pr(>F)`[2], F_test_hc$`Pr(>F)`[2])
  
}

  if (i%%5 == 0){
    
    print(c(testtype,hetero,i))
    
  }  
  
  if (i%%20 == 0){
    
    test_tmp <- (results[1:i,,] < alpha)
    if (length(methods) == 1){
      test_summary_tmp <- t(apply(test_tmp, 2, mean))
    } else{
      test_summary_tmp <- t(apply(test_tmp, c(2,3), mean))
    }
    
    print(test_summary_tmp)
    
  }  
  
}

test <- (results < alpha)
if (length(methods) == 1){
  test_summary <- t(apply(test, 2, mean))
} else{
  test_summary <- t(apply(test, c(2,3), mean))
}

rownames(test_summary) <- methods
colnames(test_summary) <- c("non-robust","robust")

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
  
}  

closeAllConnections()


