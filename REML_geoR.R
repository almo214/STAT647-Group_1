
library(geoR)
library(fields)
set.seed(42) 


n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates

# Generate the Matern covariance matrix Sigma_grid
alpha <- 1 / 0.75
nu <- 0.5
phi <- 4.0
D <- rdist(coords) # Distance matrix
Sigma_grid <- Matern(D, alpha = alpha, nu = nu, phi = phi)

# Generate the observed data 'z' using multivariate normal distribution
z <- rmvnorm(n = n, mean = rep(0, n), sigma = Sigma_grid)

# Fit separate Matern models for each column of 'z' using REML
# Not sure if this even remotely something that should be done.
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
  names(cov.pars) <- c("Sigma_sq", "Range")
  cov.pars
})

# Print the average values of the model parameters
cat("Average Model Parameters:\n")
print(rowMeans(avg_params))
