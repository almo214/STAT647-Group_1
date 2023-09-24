library( fields)
# library(RandomFields) # Library for Exponential damping no longer supported

source("REML_Group1.R")



# Generate data
n  <- 100

# Part 1 ####
## a ####

# Over range [0, 2]
coords <- cbind(runif(n,0,2),runif(n,0,2))
D <- rdist(coords) # distance matrix



likelihood.exponential_GMLE <- function(cov.pars) {
  sigma <- cov.pars[1]
  rho <- cov.pars[2]
  
  # Calculate the exponential covariance
  cov <- (sigma^2) * exp(-D / rho)
  
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


likelihood.exponential_REML <- function(cov.pars) {
  sigma <- cov.pars[1]
  rho <- cov.pars[2]
  
  # Calculate the exponential covariance
  cov <- (sigma^2) * exp(-D / rho)
  
  # Calculate Cholesky decomposition
  temp <- chol(cov)
  
  # Calculate log likelihood components
  logpart <- 2 * sum(log(diag(temp))) + 2 * (sum(log(diag(t(coords) %*% temp %*% coords ))))
  # Appears the only difference between GMLE and REML is the second portion of logpart above.
  
  step1 <- forwardsolve(t(temp), t(z))
  step2 <- backsolve(temp, step1)
  exponentpart <- z %*% step2
  
  # Negative log-likelihood
  temp4 <- logpart + exponentpart
  
  return(temp4 / 2)
}

sims1 <- matrix(rep(0,50),ncol=2)


for(i in c(1:50)){
  
  Sigma_grid = fields::Matern(D, alpha=1/.75,nu= 0.5,
                              phi=1.0)   # Is this the correct thing to use for HW2?
  z = rmvnorm(1, mu = rep(5,10), Sigma = Sigma_grid)
 
  cov.pars<-c(1,1) ## initial parameter values
  
  out1 <- nlm(likelihood.exponential_REML, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))

  out2 <- (out1$est)**2

    sims1[i,] <- out2

  
}


hist(sims1[,1])
hist(sims1[,2])
mean(sims1[,1]) - 1
mean(sims1[,2]) - 1/3
sd(sims1[,1])
sd(sims1[,2])


hist(sims2[,1])
hist(sims2[,2])
mean(sims2[,1]) - 1
mean(sims2[,2]) - 1/3
sd(sims2[,1])
sd(sims2[,2])


