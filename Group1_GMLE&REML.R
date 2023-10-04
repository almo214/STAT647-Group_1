## Title: STAT 647 Group 1 Code for HW#2


## Code Points of Contact:
# 1. Allison Moore amoore214@tamu.edu
# 2. Luis Quilarque lgquilarque@tamu.edu
# 3. Tiffany Chang tiffchang@tamu.edu


## Required libraries
library(fields)
library(MASS)
library(Matrix)
library(geoR) 


## GMLE Function
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

#########################################################################################################

## Here is an example implementation for exponential cov:

## 1. Set a seed for consistent random generation
set.seed(42)

## 2. Generate the data
n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates

## 3. Define the Exponential Matern Cov parameters and distance matrix
rho <- .75
alpha <- 1 / rho
nu <- 0.5
phi <- 4.0
D <- rdist(coords) # Distance matrix

## 4. Create a table for the parameter and mean estimates per simulation. Set the number of simulations, 
## M=50 reps
sims1 <- matrix(NA, nrow=50, ncol=3) # Add a third column in the sims matrix that stores the estimates 
# for Y(s)_hat

## Obtain the estimates and populate the sims1 table
for(i in c(1:50)){ # nrow(sims1)=M=50, 1<=i<=50 reps
  
  Sigma_grid <- Matern(D, alpha = alpha, nu = nu, phi = phi) # alpha=1/rho=1/.75; phi=gamma^2=2^2=4
  
  z <- rmvnorm(n = 1, mu = rep(5, n), Sigma = Sigma_grid) # Mean: 5; Sigma: cov(ep_1(s_1),ep_1(s_2)); 
  # n = 1 Y(s) output per simulation rep
  
  cov.pars <- c(1, 1, 1) # Initial parameter values
  
  out1 <- nlm(likelihood.exponential_GMLE, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))
  
  out2 <- (out1$est)**2
  
  sims1[i,] <- out2
  
  ## By squaring all the estimates, the mean gets squared too, which would not be an accurate reflection 
  ## of the estimated mean [ideally, it should equal to approximately 5]. So we go back and square root 
  ## all the mean estimates in col 3 of sims1
  sims1[i,3] <- sqrt(sims1[i,3])
  
}

## 5. Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sims1[,1]) # Distribution of gamma^2
hist(sims1[,2]) # Distribution of rho
hist(sims1[,3]) # Distribution of Y(s)
mean(sims1[,1]) - 4 # Empirical bias of gamma^2: E(gamma^2_hat)-gamma^2 
mean(sims1[,2]) - .75 # Empirical bias of rho: E(rho_hat)-rho
mean(sims1[,3]) - 5 # Empirical bias of Y(s): E(Y(s)_hat)-Y(s)
sd(sims1[,1]) # Standard deviation of gamma^2
sd(sims1[,2]) # Standard deviation of rho
sd(sims1[,3]) # Standard deviation of Y(s)


## REML Function
# It is not custom defined but here is a layout of what inputs to use
# Note: It fits separate Matern models for each column of 'z' using REML
# Note: We seem to be getting much smaller estimates for range, this is something we will go back and modify!
# variogram_fits <- lapply(1:ncol(z), function(col) {
#   likfit(
#     geodata = list(coords = coords, data = z[, col]), # Use each column of 'z'
#     ini.cov.pars = c(1, 1), 
#     fix.nugget = TRUE, # Assuming a fixed nugget
#     lik.method = "REML" # Use REML estimation
#   )
# })

#########################################################################################################

## Here is an example implementation for exponential cov:

## 1. Set a seed for consistent random generation
set.seed(42)

## 2. Generate the data
n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates

## 3. Define the Exponential Matern Cov parameters and distance matrix
rho <- .75
alpha <- 1 / rho
nu <- 0.5
phi <- 4.0
D <- rdist(coords) # Distance matrix

## 4. Run simulations using REML
## Note: I am only using 5 sims since the computer crashes with M=50. I would *highly* recommend using
## cluster to run larger number of sims.
sims2 <- matrix(NA, nrow=5, ncol=2)

for(i in c(1:5)){ 
  
  Sigma_grid <- Matern(D, alpha = alpha, nu = nu, phi = phi) # Generate the Matern covariance matrix, 
  # Sigma_grid
  
  z <- rmvnorm(n = n, mean = rep(0, n), sigma = Sigma_grid) # Using Exponential cov. function
  # Generate the observed data 'z' using multivariate normal distribution
  # We are setting the mean = 0 since REML circumvents estimating the mean
  
  # Fit separate Matern models for each column of 'z' using REML
  variogram_fits <- lapply(1:ncol(z), function(col) {
    likfit(
      geodata = list(coords = coords, data = z[, col]), # Use each column of 'z'
      ini.cov.pars = c(1, 1), 
      fix.nugget = TRUE, # Assuming a fixed nugget
      lik.method = "REML" # Use REML estimation
    )
  })
  
  # Calculate the average values of the model parameters
  avg_params <- sapply(variogram_fits, function(fit) {
    cov.pars <- fit$cov.pars
  })
  
  # Store the average values of the model parameters
  out2 <- rowMeans(avg_params) # Average Model Parameters for Exponential Cov
  
  sims2[i,] <- out2
  
}

## 5. Obtain distribution, mean, and sd for each param estimate of the above simulation
hist(sims2[,1])
hist(sims2[,2])
mean(sims2[,1]) - 4
mean(sims2[,2]) - .75
sd(sims2[,1])
sd(sims2[,2])

## NOTE: Need to fix REML rho parameter estimation. Would recommend doing research on likfit() function 
## in order to figure how to estimate it correctly.

########################################################################################################
# Load the geostatsp library for damped covariance
library(geostatsp)

# 1. Set a seed for consistent random generation
set.seed(42)

# 2. Generate the data
n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates

# 3. Define the Damped Covariance parameters and distance matrix
rho <- 0.75

D <- rdist(coords) # Distance matrix

# Define a custom damped cosine covariance function
damped_cosine_covariance <- function(D, rho) {
  gamma_sq <- 2^2  
  exp(-D / rho) * gamma_sq * cos(D)
}

# 4. Run simulations using REML
sims_damped <- matrix(NA, nrow = 5, ncol = 2)  # Store results in a matrix

for (i in 1:5) {
  
  # Generate the Damped cosine covariance matrix
  Sigma_grid <- damped_cosine_covariance(D, rho = rho)
  
  z <- rmvnorm(n = n, mean = rep(0, n), sigma = Sigma_grid) # Generate the observed data 'z'
  
  # Fit separate Damped models for each column of 'z' using REML
  variogram_fits <- lapply(1:ncol(z), function(col) {
    likfit(
      geodata = list(coords = coords, data = z[, col]),
      ini.cov.pars = c(1, rho),  # Provide initial values for sigmasq and rho
      fix.nugget = TRUE,
      lik.method = "REML"
    )
  })
  
  # Calculate the average values of the model parameters
  avg_params <- sapply(variogram_fits, function(fit) fit$cov.pars[2])  # Extract rho values
  
  # Store the average values of rho and gamma_sq
  out_damped <- c(mean(avg_params), sqrt(mean(sapply(variogram_fits, function(fit) fit$cov.pars[1]))))
  
  sims_damped[i,] <- out_damped
}

# 5. Obtain distribution, mean, and sd for the parameter estimates of the above simulation
hist(sims_damped[, 1], main = "Damped Cosine Rho")
hist(sims_damped[, 2], main = "Damped Cosine Gamma")

mean(sims_damped[, 1])
mean(sims_damped[, 2])
sd(sims_damped[, 1])
sd(sims_damped[, 2])
