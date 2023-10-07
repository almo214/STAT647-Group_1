##### Functions sourced from Group 2's code

#The 0c(0;θ) function represents the covariance at lag 0 same for both exponential and damped exponential
c0_function  <- function(gam) {
  return(gam^2)
}

# Exponentially Covariance Function (subbing r= 1/Rg)
exp_cov  <- function(D, gam, rho) {
  covariance <- (gam^2) * exp(-D / rho)
  return(covariance)
}

# Exponentially Damped Covariance Function
# The exponentially damped covariance function does not have  as a parameter; instead, it uses γ(or gamma)(the rate of decay) and ρ (or rho)(damping parameter)
exp_damped_cov <- function(D, gam, rho) {
  # if(rho > tan(pi / (2 * d))){
  #   stop("rho must be less than or equal to tan(pi / (2 * d))")
  # }
  covariance <- gam^2 * exp(- (D / rho)) * cos(D)
  return(covariance)
}

# Marginal likelihood function (used to evaluate the first term in the summation)
marglik = function(y, gam, sigma, mu){
  return(sum(log(c0_function(gam) + sigma) - (y - mu)^2 / (c0_function(gam) + sigma)))
}

## Exponentially Covariance Function
### Given the parameters μ=5, σ^2 =1 and estimating phi :
vecchia_likelihood <- function(params, y, mu, gamma, sigma, D, damped = F) {
  # sigma = 1 #initialed equal to 1
  rho = params
  # UPDATE: I think y means Y(s), which is z <- rmvnorm(n = 1, mu = rep(5, n), Sigma = Sigma_grid) in our simulations.
  # Please confirm if the above interpretation is correct.
  # mu = 5 #  μ=5
  # gamma = 2
  if(damped){
    covfun = exp_damped_cov
  } else {
    covfun = exp_cov
  }
  
  QV = 0
  for (i in 1:n) {
    if (i == 1) {
      QV = QV + marglik(y[1], gamma, sigma, mu)
      next
    } else if (i <= 20) {
      neighbors = 1:(i - 1)
    } else {
      neighbors = (i - 20):(i - 1)
    }
    
    D_neighbors = matrix(D[neighbors, neighbors], nrow=length(neighbors), ncol=length(neighbors))
    Sigma = matrix(covfun(D_neighbors, gamma, phi), nrow=length(neighbors), ncol=length(neighbors)) + diag(sigma, nrow = length(neighbors))
    
    # Add a small nugget to the diagonal for regularization
    # Sigma = Sigma + diag(1e-6, nrow=length(neighbors))
    #if(any(is.na(Sigma)) || any(is.infinite(Sigma))) {
    #  return(1e+6)  # Large positive value since we're minimizing
    #}
    
    Sigma_inv = solve(Sigma) 
    
    r = matrix(covfun(D[i, neighbors], gamma, phi), ncol=1)
    
    rinv = t(r) %*% Sigma_inv %*% r
    
    term1 = log(c0_function(gamma) + sigma - rinv)
    
    term2_numerator = y[i] - mu - t(r) %*% Sigma_inv %*% matrix(y[neighbors] - mu, ncol=1)
    term2_denominator = c0_function(gamma) + sigma - rinv
    term2 = term2_numerator^2 / term2_denominator
    
    QV = QV + term1 + term2
  }
  return(QV)
}
