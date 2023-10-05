## Instructions on how to load the RandomFields library (archived version):
## 1. Download the tar.gz files for RandomFields and its dependency, RandomFieldsUtils: 
## --> RandomFieldsUtils_1.2.3.tar.gz (https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/)
## --> RandomFields_3.3.14.tar.gz (https://cran.r-project.org/src/contrib/Archive/RandomFields/)
## 2. Move both files to Documents folder (or just somewhere that is easy to remember)
## 3. Install the packages (RandomFieldsUtils and RandomFields, in that order), using the following 
## commands (right-click on the tar.gz file and select 'Copy as path' for the first argument of 
## install.packages):
## --> install.packages('C:/Users/nioni/Documents/RandomFieldsUtils_1.2.3.tar.gz',repos=NULL,type='source')
## --> install.packages('C:/Users/nioni/Documents/RandomFields_3.3.14.tar.gz',repos=NULL,type='source')
## 4. Load the RandomFields library:

library(RandomFields)

# Define Exponential Damped covariance 
exp_Damped_cov <- function(h, rho, gamma) {
  (gamma^2) * exp(-h / rho^2) * cos(h)
}

## Define GMLE function for exponential damped cov
likelihood.exp_damp_GMLE <- function(cov.pars) {
  
  gamma <- cov.pars[1] 
  rho <- cov.pars[2] 
  mu <- cov.pars[3] # Corresponds to Y(s)_hat (estimated mean)
  z = z - mu # Subtracted the estimated mean from z so that the parameter estimates calculations are 
  # accurate when z is centered around close to 0
  
  # Calculate the exponential covariance
  cov <- exp_Damped_cov(D, rho, gamma)
  
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

## Here is an example implementation for exponential damped cov:

## 1. Set a seed for consistent random generation
set.seed(42)

## 2. Generate the data
n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates

## 3. Define the Exponential Matern Cov parameters and distance matrix
rho <- .75
lambda <- 1 / rho
var <- 4.0 
D <- rdist(coords) # Distance matrix

## 4. Create a table for the parameter and mean estimates per simulation. Set the number of simulations, 
## M=50 reps
sims1 <- matrix(NA, nrow=50, ncol=3) # Add a third column in the sims matrix that stores the estimates 
# for Y(s)_hat

## Obtain the estimates and populate the sims1 table
for(i in c(1:50)){ # nrow(sims1)=M=50, 1<=i<=50 reps
  
  Sigma_grid <- RMdampedcos(lambda=lambda, var=var) # lambda=1/rho; var=gamma^2=2^2=4
  
  z <- rmvnorm(n = 1, mu = rep(5, n), Sigma = Sigma_grid) # Mean: 5; Sigma: cov(ep_1(s_1),ep_1(s_2)); 
  # n = 1 Y(s) output per simulation rep
  
  cov.pars <- c(1, 1, 1) # Initial parameter values
  
  out1 <- nlm(likelihood.exp_damp_GMLE, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))
  
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

## Note: Not sure how to get this to work, but here is the output of RMdampedcos function and its plot
Sigma_grid <- RMdampedcos(lambda=lambda, var=var) # Returns an object of class RMmodel
plot(Sigma_grid)
