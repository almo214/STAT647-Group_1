
## Question 1

## b) - Part i
## Note: ALL simulations will use the following R exponential (damped) covariance function parameters as
## such:
## - Model 1 Cov: phi=gamma^2=4, given that gamma=2; alpha=1/range=1/rho=1/.75, given that rho=.75
## - Model 2 Cov: (Haven't looked into this yet, but will update once I do. Task: Define/map exp damped 
## function cov param's to sigma^2=1 and rho=.75) 
## This is just my preliminary edits and annotations to understand the code better. Please double check 
## these facts, and correct me if I'm wrong on anything!
## Look for updated code and annotations with the keyword, "UPDATE: "

#########################################################################################################

## Load required libraries
library(fields)
library(MASS)
library(Matrix)
library(geoR) 
library(geostatsp)
# library(RandomFields) # Library for Exponential damping no longer supported

# Define all Functions ####

# Define Exponential Damped covariance 
exp_Damped_cov <- function(h, rho, phi) {
  (phi^2) * exp(-h / rho^2) * cos(h)
}

# Luis's function, removed nugget
damped_covariance <- function(D, phi, range) {
  exp(-D / range) * (1 + D / (phi * range)) * exp(-D / phi)
}



## Define GMLE function
likelihood.exponential_GMLE <- function(cov.pars) {
  gamma <- cov.pars[1] 
  rho <- cov.pars[2] 
  mu <- cov.pars[3] # Corresponds to Y(s)_hat (estimated mean)
  z = z - mu # Subtracted the estimated mean from z so that the parameter estimates calculations are 
  # accurate when z is centered around close to 0
  
  # Calculate the exponential covariance
  cov <- (gamma^2) * exp(-D / rho^2)
  
  # Calculate Cholesky decomposition
  temp <- chol(cov)
  
  # Calculate log likelihood components
  logpart <- 2 * sum(log(diag(temp)))
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  
  # Negative log-likelihood
  temp4 <- logpart + exponentpart
  
  return(temp4 / 2)

}



## Define GMLE function for exponential damped cov
likelihood.exp_damp_GMLE <- function(cov.pars) {
  gamma <- cov.pars[1] 
  rho <- cov.pars[2] 
  mu <- cov.pars[3] # Corresponds to Y(s)_hat (estimated mean)
  z = z - mu # Subtracted the estimated mean from z so that the parameter estimates calculations are 
  # accurate when z is centered around close to 0
  
  # Calculate the exponential covariance
  cov <- exp_Damped_cov(D, rho, phi)
   # cov <- damped_covariance(D, phi, rho)
  
  # Calculate Cholesky decomposition
  temp <- chol(cov)
  
  # Calculate log likelihood components
  logpart <- 2 * sum(log(diag(temp)))
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  
  # Negative log-likelihood
  temp4 <- logpart + exponentpart
  
  return(temp4 / 2)
  
}




# Simulations for Model 1 (exponential cov) ####
## Generate data
set.seed(42)
n  <- 1000
size = 2 # over range [0,size]

# 1-D
coords <- runif(n,0,size) # Over range [0, size]

## For 2-D 
#coords <- cbind(runif(n,0,size), runif(n,0,size))


D <- rdist(coords) # distance matrix

rho <- .75
alpha <- 1 / rho
nu <- 0.5
phi <- 4.0


## Method 1: GMLE
sims1 <- matrix(NA, nrow=1000, ncol=3) 

## Obtain the estimates and populate sims1
for(i in c(1:1000)){ 
  set.seed(i)
  
  Sigma_grid <- Matern(D, alpha = alpha, nu = nu, phi = phi) # alpha=1/rho=1/.75; phi=gamma^2=2^2=4
  
  z <- rmvnorm(n = 1, mu = rep(5, n), Sigma = Sigma_grid) # Mean: 5; Sigma: cov(ep_1(s_1),ep_1(s_2)); 
  # n = 1 Y(s) output per simulation rep
 
  cov.pars<-c(1,1,1) # initial parameter values
  
  # Use tryCatch to skip any iterations which the likelihood function throws an error. Those will be recorded
  # as NA in the sims1 matrix.
  tryCatch ({
  out1 <- nlm(likelihood.exponential_GMLE, cov.pars, stepmax=2, print.level=0, gradtol=10^(-6))

  out2 <- (out1$est)**2

  sims1[i,] <- out2
  
  ## UPDATE: By squaring all the estimates, the mean gets squared too, which would not be an accurate      
  ## reflection of the estimated mean [ideally, it should equal to approximately 5]. So we go back and 
  ## square root all the mean estimates in col 3 of sims1
  sims1[i,3]=sqrt(sims1[i,3])
  print(i)
  }, error = function(e) {
    sims1[i,] <- c(NA, NA, NA)
  })
}

par(mfrow=c(3,1))
## Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sims1[,1], main ="Histograms GMLE [0, 2], Exponential, n =1000", xlab ='sigma') # Distribution of gamma^2
hist(sims1[,2], main='', xlab ='rho') # Distribution of rho
hist(sims1[,3], main='', xlab ='mean') # UPDATE: Distribution of Y(s)
mean(sims1[,1], na.rm =TRUE) - 4 # Empirical bias of gamma^2: E(gamma^2_hat)-gamma^2 
mean(sims1[,2], na.rm =TRUE) - .75 # Empirical bias of rho: E(rho_hat)-rho
mean(sims1[,3], na.rm =TRUE) - 5 # Empirical bias of Y(s): E(Y(s)_hat)-Y(s)
sd(sims1[,1], na.rm =TRUE) # Standard deviation of gamma^2
sd(sims1[,2], na.rm =TRUE) # Standard deviation of rho
sd(sims1[,3], na.rm =TRUE) # Standard deviation of Y(s)



## Method 2: REML ####


sims2 <- matrix(NA, nrow=1000, ncol=2)

for(i in c(1:1000)){ 
  set.seed(i)
  Sigma_grid <-  Matern(D, alpha = alpha, nu = nu, phi = phi) # Generate the Matern covariance matrix, 
  # Sigma_grid
  
  z <- rmvnorm(n = n, mean = rep(5, n), sigma = Sigma_grid) # Using Exponential cov. function
  # Generate the observed data 'z' using multivariate normal distribution
  
  # Use tryCatch to skip any iterations which the likelihood function throws an error. Those will be recorded
  # as NA in the sims1 matrix.
  tryCatch ({
  # Fit separate Matern models for each column of 'z' using REML
  variogram_fits <- lapply(1:ncol(z), function(col) {
    likfit(
      geodata = list(coords = coords, data = z[, col]), # Use each column of 'z'
      ini.cov.pars = c(1, .5), 
      fix.nugget = TRUE, # Assuming a fixed nugget
      lik.method = "REML", # Use REML estimation
      messages =  FALSE
    )
  })
  
  # Calculate the average values of the model parameters
  avg_params <- sapply(variogram_fits, function(fit) {
    cov.pars <- fit$cov.pars
  })
  
  # Store the average values of the model parameters
  out2 <- rowMeans(avg_params) # Average Model Parameters for Exponential Cov
  print(i)
  sims2[i,] <- out2
  }, error = function(e) {
    sims2[i,] <- c(NA, NA)
  })
}

par(mfrow=c(2,1))
## 5. Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sqrt(sims2[,1]), main ="Histograms REML [0, 2], Exponential, n =100", xlab ='sigma')
hist(sims2[,2], main='', xlab ='rho')
mean(sqrt(sims2[,1]), na.rm =TRUE) - 4
mean(sims2[,2], na.rm =TRUE) - .75
sd(sims2[,1], na.rm =TRUE)
sd(sims2[,2], na.rm =TRUE)

## UPDATE: Code for 1D sims --> HW#1 (but GMLE and REML functions remain the same; both 1D and 2D use the
## same likelihood functions)

#########################################################################################################

### Simulations for Model 2 (exponential damped cov) ####

#########################################################################################################



## Generate data
set.seed(42)
n  <- 100
size = 2 # over range [0,size]

# 1-D
# coords <- runif(n,0,size) # Over range [0, size]

## For 2-D 
coords = cbind(runif(n,0,size), runif(n,0,size))



D <- rdist(coords) # distance matrix



rho <- .75
alpha <- 1 / rho
nu <- 0.5
phi <- 4.0


## Method 1: GMLE
sims1 <- matrix(NA, nrow=1000, ncol=3) # UPDATE: Added a third column in the sims matrix that stores the 
# estimates for the Y(s)_hat

## Obtain the estimates and populate sims1
for(i in c(1:1000)){ 
  set.seed(i)
 
   # signma_grd <- damped_covariance(D, phi, rho)

  Sigma_grid <- exp_Damped_cov(D,  rho, phi)
  
  z = rmvnorm(n = 1, mu = rep(5,n), Sigma = Sigma_grid) # Model 1 --> Mean: 5; Sigma: 
  # cov(ep_1i(s_1),ep_1i(s_2)), n = 1 Y(s) output per simulation rep
 
  cov.pars<-c(1,1,1) # initial parameter values
  
  # Use tryCatch to skip any iterations which the likelihood function throws an error. Those will be recorded
  # as NA in the sims1 matrix.
  tryCatch ({
    out1 <- nlm(likelihood.exp_damp_GMLE, cov.pars, stepmax=2, print.level=0, gradtol=10^(-6))
    
    out2 <- (out1$est)**2
    
    sims1[i,] <- out2
    
    ## UPDATE: By squaring all the estimates, the mean gets squared too, which would not be an accurate      
    ## reflection of the estimated mean [ideally, it should equal to approximately 5]. So we go back and 
    ## square root all the mean estimates in col 3 of sims1
    sims1[i,3]=sqrt(sims1[i,3])
    print(i)
  }, error = function(e) {
    sims1[i,] <- c(NA, NA, NA)
  })
}


## Obtain distribution, mean, and sd for each param estimate of the above simulation

par(mfrow=c(3,1))
## Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sims1[,1], main ="Histograms GMLE [0, 2]^2, Exponential Damped, n =100", xlab ='sigma') # Distribution of gamma^2
hist(sims1[,2], main='', xlab ='rho') # Distribution of rho
hist(sims1[,3], main='', xlab ='mean') # UPDATE: Distribution of Y(s)
mean(sims1[,1], na.rm =TRUE) - 4 # Empirical bias of gamma^2: E(gamma^2_hat)-gamma^2 
mean(sims1[,2], na.rm =TRUE) - .75 # Empirical bias of rho: E(rho_hat)-rho
mean(sims1[,3], na.rm =TRUE) - 5 # Empirical bias of Y(s): E(Y(s)_hat)-Y(s)
sd(sims1[,1], na.rm =TRUE) # Standard deviation of gamma^2
sd(sims1[,2], na.rm =TRUE) # Standard deviation of rho
sd(sims1[,3], na.rm =TRUE) # Standard deviation of Y(s)


## Method 2: REML
sims2 <- matrix(NA, nrow=5, ncol=2)

for(i in c(1:5)){ 
  
  # signma_grd <- damped_covariance(D, phi, rho)
  
  Sigma_grid <- exp_Damped_cov(D,  rho, phi)
  
  z <- rmvnorm(n = n, mean = rep(5, n), sigma = Sigma_grid) # Using Exponential cov. function
  # Generate the observed data 'z' using multivariate normal distribution
  # We are setting the mean = 0 since REML circumvents estimating the mean
  
  
  # Use tryCatch to skip any iterations which the likelihood function throws an error. Those will be recorded
  # as NA in the sims1 matrix.
  tryCatch ({
    # Fit separate Matern models for each column of 'z' using REML
    variogram_fits <- lapply(1:ncol(z), function(col) {
      likfit(
        geodata = list(coords = coords, data = z[, col]), # Use each column of 'z'
        ini.cov.pars = c(1, .5), 
        fix.nugget = TRUE, # Assuming a fixed nugget
        lik.method = "REML", # Use REML estimation
        messages =  FALSE
      )
    })
    
    # Calculate the average values of the model parameters
    avg_params <- sapply(variogram_fits, function(fit) {
      cov.pars <- fit$cov.pars
    })
    
    # Store the average values of the model parameters
    out2 <- rowMeans(avg_params) # Average Model Parameters for Exponential Cov
    print(i)
    sims2[i,] <- out2
  }, error = function(e) {
    sims2[i,] <- c(NA, NA)
  })
}

## 5. Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sqrt(sims2[,1]), main ="Histograms REML [0, 2], Exponential Damped, n =100", xlab ='sigma')
hist(sims2[,2], main='', xlab ='rho')
mean(sqrt(sims2[,1]), na.rm =TRUE) - 4
mean(sims2[,2], na.rm =TRUE) - .75
sd(sqrt(sims2[,1]), na.rm =TRUE)
sd(sims2[,2], na.rm =TRUE)