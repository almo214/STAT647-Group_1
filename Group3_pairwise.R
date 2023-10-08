library(fields)
library(MASS)

#######
## suppose x and y coordinates are saved in x and y, and the data are saved in z 
##(x,y,z are Nx1 vectors, when N is the number of locations)
#######

# Generate data with Matern covariance 
#Y=beta'X+e+epsilon
#E(e)=0, cov(e)~Matern(D, alpha=theta0, nu= 0.5, phi=1.0)
# epsilon ~N(0,s^2)


# set.seed(100)

# n=100
# p=3   # number of beta parameters

# coords = cbind(runif(n,0,10),runif(n,0,10))
# D = rdist(coords) 


#### function for pairwise likelihood ####

likelihood.pairwise <- function(param, n){ 
  
  beta = param[1:p]
  s2 = param[p+1]
  theta = param[p+2]
  #nu = param[p+3]
  #phi = param[p+4]
  logpart=0 
  exponentpart=0 
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      EY = c(beta%*%X[,i],beta%*%X[,j])
      pairD = matrix(c(D[i,i],D[i,j],D[j,i],D[j,j]), nrow = 2)
      cov = Matern(pairD, alpha=theta, nu= 0.5, phi=1)
      #cov = Matern(pairD, range=theta, nu= nu0, phi=phi0)
      Mij=matrix(c(cov[1,1]+s2, cov[1,2], cov[2,1], cov[2,2]+s2), nrow = 2) 
      
      logpart=logpart+log(det(Mij))
      exponentpart=exponentpart+ t(c(Y[i]-EY[1],Y[j]-EY[2]))%*%solve(Mij)%*%c(Y[i]-EY[1],Y[j]-EY[2]) 
      
    }
  }
  
  Q  <- logpart+exponentpart
  return(Q)
}

########################################################

### True parameter and corresponding data ###
# beta0 = t(c(1,4,3))
# s20= 1
# theta0= 4
# #nu0 = 3
# #phi0 = 2

# S= Matern(D, alpha=theta0, nu= 0.5, phi=1)
# X=matrix(1, nrow = p, ncol = n) + matrix(rnorm(p*n,0,1), nrow = p, ncol = n)
# Y = t(beta0%*%X+ rmvnorm(1, mu = 0, Sigma = S)+rnorm(n,0,s20))

### Initialization of estimates ###
# beta_initial = t(rep(1,p))
# s2_initial=1
# theta_initial=1
# #nu_initial = 1
# #phi_initial = 1
# param=c(beta_initial, s2_initial ,theta_initial)
# #param=c(beta_initial, s2_initial ,theta_initial, nu_initial, phi_initial)

### Estimation ###
# fit=nlm(likelihood.pairwise, param, stepmax=5, print.level=2, gradtol=10^(-10))
# beta_fit= fit$estimate[1:p]
# s2_fit=fit$estimate[p+1] 
# theta_fit=fit$estimate[p+2] 



######################## Simulation ########################

# simulation = matrix(rep(0,250),ncol=5)

# for(i in c(1:50)){
  
#   ### True parameter and simulated data ###
#   beta0 = t(c(1,4,3))
#   s20= 1
#   theta0= 4
#   #nu0 = 3
#   #phi0 = 2
  
#   S= Matern(D, alpha=theta0, nu= 0.5, phi=1)
#   X=matrix(1, nrow = p, ncol = n) + matrix(rnorm(p* n,0,1), nrow = p, ncol = n)
#   Y = t(beta0%*%X+ rmvnorm(1, mu = 0, Sigma = S)+rnorm(n,0,s20))
  
#   ### Initialization of estimates ###
#   beta_initial = t(rep(1,p))
#   s2_initial=1
#   theta_initial=1
#   #nu_initial = 1
#   #phi_initial = 1
#   param=c(beta_initial, s2_initial, theta_initial)
#   #param=c(beta_initial, s2_initial ,theta_initial, nu_initial, phi_initial)
  
#   ### Estimation ###
#   fit=nlm(likelihood.pairwise, param, stepmax=5, print.level=2, gradtol=10^(-10))
#   #beta_fit= fit$estimate[1:p]
#   #s2_fit=(fit$estimate[p+1]) 
#   #theta_fit=fit$estimate[p+2] 
  
#   est = fit$estimate
  
#   simulation[i,] = est
  
# }

# hist(simulation[,1])
# hist(simulation[,2])
# hist(simulation[,3])
# hist(simulation[,4])
# hist(simulation[,5])

# mean(simulation[,1])
# mean(simulation[,2])
# mean(simulation[,3])
# mean(simulation[,4])
# mean(simulation[,5])

# sd(simulation[,1])
# sd(simulation[,2])
# sd(simulation[,3])
# sd(simulation[,4])
# sd(simulation[,5])

