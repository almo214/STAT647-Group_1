# Load the GeoR package
library(geoR)
set.seed(42) 


# Assuming you have the following data
n <- 100
coords <- cbind(runif(n, 0, 2), runif(n, 0, 2)) # Random 2-D coordinates


# Generate the Matern covariance matrix Sigma_grid
alpha <- 1 / 0.75
nu <- 0.5
phi <- 4.0
D <- rdist(coords) # Distance matrix
Sigma_grid <- Matern(D, alpha = alpha, nu = nu, phi = phi)

# Generate the observed data 'z' using multivariate normal distribution

z <- rmvnorm(n = n, mean = 0, sigma = Sigma_grid)


# Fit the Matern model using REML
variogram_fit <- likfit(
  geodata = list(coords = coords, data = z[,1]), # Can't get this to run with more than one column for data...
  ini.cov.pars = c(1, 1),
  fix.nugget = TRUE, # Assuming a fixed nugget
  lik.method = "REML" # Use REML estimation
)

# Print a summary of the REML estimates
summary(variogram_fit)
