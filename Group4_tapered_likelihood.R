library("fields")
library(MASS)

# Generate Data
# n  = 1000
# num_simulations <- 100

# coords = cbind(runif(n,0,4),runif(n,0,4))
# D = rdist(coords)

# Define the Wendland covariance function
wendland_cov <- function(d, support, smoothness) {
  r <- d / support
  cov <- (1 - r) ^ (smoothness + 1) * (1 + (smoothness + 1) * r) * (r < 1)
  return(cov)
}

# Generate the tapering matrix using the Wendland function
#Be careful changing tarpering will change the quality of simulation

# taper_matrix <- wendland_cov(D, support = 1, smoothness = 2)

# Generate the original covariance matrix using the Matern function
# Sigma_original = Matern(D, alpha=3, nu= 0.5, phi=1.0)

# Apply the tapering to the original covariance matrix
# Sigma_tapered = Sigma_original * taper_matrix


# Generate multivariate normal data using the tapered covariance matrix
# z = rmvnorm(1, mu = rep(0, n), Sigma = Sigma_tapered)
# Define the likelihood function
likelihood.matern <- function(cov.pars){ 
  sigma = cov.pars[1]
  Rg = cov.pars[2]   
  mu <- cov.pars[3] # Corresponds to Y(s)_hat (estimated mean)
  z = z - mu # Subtracted the estimated mean from z so that the parameter estimates calculations are 
  # accurate when z is centered around close to 0
  cov = (sigma^2) * exp(-D/(Rg^2)) 
  tapered_cov = cov * taper_matrix
  temp <- chol(tapered_cov)
  logpart <- 2 * sum(log(diag(temp)))
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  return((logpart + exponentpart)/2)
}
# Perform the optimization
# cov.pars = c(1,1) 
# nlm(likelihood.matern, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))

#####################################
#now repeat the process for multiple times, and record the mean of data generated and the parameter estimated each time
# n  = 1000
# num_simulations <- 10

# Create a vector to store the mean of z for each simulation
# mean_z <- numeric(num_simulations)

# Create a matrix to store the estimated parameters for each simulation
# estimated_params <- matrix(NA, nrow=num_simulations, ncol=2) # Assuming 2 parameters are being estimated

# Initializing cov.pars
# cov.pars = c(1,1)

# for (i in 1:num_simulations) {
#   coords = cbind(runif(n, 0, 4), runif(n, 0, 4))
#   D = rdist(coords)
  
#   taper_matrix <- wendland_cov(D, support = 1, smoothness = 2)
#   Sigma_original = Matern(D, alpha=3, nu= 0.5, phi=1.0)
#   Sigma_tapered = Sigma_original * taper_matrix
  
#   z = rmvnorm(1, mu = rep(0, n), Sigma = Sigma_tapered)
#   mean_z[i] <- mean(z)
  
#   likelihood.matern <- function(cov.pars){ 
#     sigma = cov.pars[1]
#     Rg = cov.pars[2]    
#     cov = (sigma^2) * exp(-D/(Rg^2))
#     tapered_cov = cov * taper_matrix
#     temp <- chol(tapered_cov)
#     logpart <- 2 * sum(log(diag(temp)))
#     step1 <- forwardsolve(t(temp), t(z))
#     step2 <- backsolve(temp, step1)
#     exponentpart <- z %*% step2
#     return((logpart + exponentpart)/2)
#   }
  
#   result <- nlm(likelihood.matern, cov.pars, stepmax=5, print.level=0, gradtol=10^(-10))
#   estimated_params[i,] <- result$estimate
# }

# View the results
# list(mean_z = mean_z, estimated_params = estimated_params)
