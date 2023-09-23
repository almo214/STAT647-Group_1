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


likelihood.matern <- function(cov.pars){ 

  sigma <- cov.pars[1]
  Rg <- cov.pars[2]    
  cov<-(sigma**2)*exp(-D/Rg) # Exponential Cov.

  # chol calculates the Cholesky decomposition of a positive definite matrix
  
  temp <- chol(cov)
  
  logpart <- 2*sum(log(diag(temp)))
  
  step1 <- forwardsolve(t(temp),t(z))
  step2 <- backsolve(temp, step1)
  
  exponentpart <- z%*%step2
  
  # negative log-likelihood
  
  temp4 <- logpart+exponentpart

  return(temp4/2)
  
}





sims1 <- matrix(rep(0,50),ncol=2)

for(i in c(1:50)){
  
  Sigma_grid = fields::Matern(D, alpha=1/.75,nu= 0.5,
                              phi=1.0)   # Is this the correct thing to use for HW2?
  z = rmvnorm(1, mu = rep(5,10), Sigma = Sigma_grid)
  
  cov.pars<-c(1,1)
  
  nlm(likelihood.matern, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))
  
  cov.pars<-c(1,1) ## initial parameter values
  
  out1 <- nlm(likelihood.matern, cov.pars, stepmax=5, print.level=2, gradtol=10^(-10))
  
  out2 <- (out1$est)**2
  
  sims1[i,] <- out2
  
}


hist(sims1[,1])
hist(sims1[,2])
mean(sims1[,1]) - 1
mean(sims1[,2]) - 1/3
sd(sims1[,1])
sd(sims1[,2])





