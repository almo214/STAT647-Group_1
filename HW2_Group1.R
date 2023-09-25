## Question 1

## b) - Part i
## Note: ALL simulations will use the following R exponential (damped) covariance function parameters as
## such:
## - Model 1 Cov: phi=gamma^2=4, given that gamma=2; alpha=1/range=1/rho=1/.75, given that rho=.75
## - Model 2 Cov: (Haven't looked into this yet, but will update once I do. Task: Define/map exp damped 
## function cov param's to sigma^2=1 and rho=.75) 
## This is just my preliminary edits and annotations to understand the code better. Please double check 
## these facts, and correct me if I'm wrong on anything!

#########################################################################################################

## Load required libraries
library(fields)
library(MASS)
library(Matrix)
# library(RandomFields) # Library for Exponential damping no longer supported

# source("REML_Group1.R") # Note: I commented this out on my Rmd file since the code won't run if I kept 
# it here

### Simulations for Model 1 (exponential cov)
## Generate data
n  <- 100

coords <- cbind(runif(n,0,2),runif(n,0,2)) # Over range [0, 2]
D <- rdist(coords) # distance matrix

## Define GMLE function
likelihood.exponential_GMLE <- function(cov.pars) {
  sigma <- cov.pars[1] ## Note: This value corresponds to gamma for Model 1
  rho <- cov.pars[2]
  
  # Calculate the exponential covariance
  cov <- (sigma**2) * exp(-D / rho^2) # UPDATE: The GMLE simulation study code originally had rho^2, and 
  # the outputs made more sense with rho^2 as opposed to rho (not sure why though)
  
  # Calculate Cholesky decomposition
  # temp <- chol(cov) # Original code from GMLE simulation study example - does not work for phi > 1
  temp <- chol(nearPD(cov)$mat)

  
  # Calculate log likelihood components
  logpart <- 2 * sum(log(diag(temp)))
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  
  # Negative log-likelihood
  temp4 <- logpart + exponentpart
  
  return(temp4 / 2)
}

## Define REML function
likelihood.exponential_REML <- function(cov.pars) {
  sigma <- cov.pars[1]
  rho <- cov.pars[2]
  
  # Calculate the exponential covariance
  cov <- (sigma^2) * exp(-D / rho^2)
  
  # Calculate Cholesky decomposition
  temp <- chol(nearPD(cov)$mat)
  
  # Calculate log likelihood components
  logpart <- 2 * sum(log(diag(temp))) + 2 * (sum(log(diag(t(coords) %*% temp %*% coords))))
  # Appears the only difference between GMLE and REML is the second portion of logpart above.
  
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  
  # Negative log-likelihood
  temp4 <- logpart + exponentpart
  
  return(temp4 / 2)
}

## Create a table for the parameter estimates per simulation. Set the number of simulations, M=50 reps
## (Note: Will eventually change M to 1000 once we get the code working properly for small sim size.)
## Method 1: GMLE
sims1 <- matrix(rep(0,50),ncol=2)

## Obtain the estimates and populate the sims1 table
for(i in c(1:50)){ # nrow(sims1)=M=50, 1<=i<=50 reps
  
  Sigma_grid = Matern(D, alpha=1/.75,nu=0.5,phi=4.0) # Is this the correct thing to use for HW2? 
  # A: Not quite! It is correct that you want to use alpha=1/range=1/.75, but phi should be gamma^2=2^2=4
  # since it represents the *marginal variance* of the exponential process
  # Info source: https://rdrr.io/cran/fields/src/R/Matern.R
  
  z = rmvnorm(n = 1, mu = rep(5,n), Sigma = Sigma_grid) # Model 1 --> Mean: 5; Sigma: 
  # cov(ep_1i(s_1),ep_1i(s_2)), n = 1 Y(s) output per simulation rep
  # UPDATE: I was able to get the empirical bias close to 0 *only* when mu=rep(0,n). I believe this is 
  # because a location shift in the mean of Y(s) causes the bias and standard deviation to increase
 
  cov.pars<-c(1,1) # initial parameter values
  
  out1 <- nlm(likelihood.exponential_GMLE, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))

  out2 <- (out1$est)**2

  sims1[i,] <- out2
  
}

## Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sims1[,1]) # Distribution of gamma^2
hist(sims1[,2]) # Distribution of rho
mean(sims1[,1]) - 4 # Empirical bias of gamma^2: E(gamma^2_hat)-gamma^2 
mean(sims1[,2]) - .75 # Empirical bias of rho: E(rho_hat)-rho
sd(sims1[,1]) # Standard deviation of gamma^2
sd(sims1[,2]) # Standard deviation of rho

## Method 2: REML
sims2 <- matrix(rep(0,50),ncol=2)

for(i in c(1:50)){ 
  
  Sigma_grid = fields::Matern(D, alpha=1/.75,nu=0.5,phi=4.0) 
  
  z = rmvnorm(n = 1, mu = rep(5,n), Sigma = Sigma_grid)
 
  cov.pars<-c(1,1)
  
  out1 <- nlm(likelihood.exponential_REML, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))

  out2 <- (out1$est)**2

  sims2[i,] <- out2

}

hist(sims2[,1])
hist(sims2[,2])
mean(sims2[,1]) - 4
mean(sims2[,2]) - .75
sd(sims2[,1])
sd(sims2[,2])

## UPDATE: The values look much different for empirical bias' between GMLE and REML methods, but it seems
## like the standard errors decreased drastically in REML estimates compared to those of
## GMLE, indicating that REML is a better estimator than GMLE. 

#########################################################################################################

### Simulations for Model 2 (exponential damped cov)
## (To be continued...)
