# Group 1 REML Code adapted from:
# Source -> https://rh8liuqy.github.io/Example_Linear_mixed_model.html


# Using lme4 library ####
# install.packages('lme4')
# install.packages('optimization')
library(lme4)
library(optimization)



## Define the data ####
df.bond <- data.frame(ingot = rep(1:7,each = 3),
                      metal = rep(c("n","i","c"),7),
                      pres = c(67,71.9,72.2,
                               67.5,68.8,66.4,
                               76,82.6,74.5,
                               72.7,78.1,67.3,
                               73.1,74.2,73.2,
                               65.8,70.8,68.7,
                               75.6,84.9,69))
df.bond$ingot <- factor(df.bond$ingot)
df.bond




## Linear mixed model by MLE ####
lmm.bond <- lmer( pres ~ metal + (1|ingot),REML = FALSE,data = df.bond)
summary(lmm.bond)
anova(lmm.bond)

## Linear mixed model by REML  ####
relmm.bond <- lmer( pres ~ metal + (1|ingot),REML = TRUE,data = df.bond)
summary(relmm.bond)
anova(relmm.bond)

# Compare AIC/BIC. Note REML provides lowest AIC and BIC
(AIC_mle <- AIC(lmm.bond))
(BIC_mle <- BIC(lmm.bond))
(AIC_reml <- AIC(relmm.bond))
(BIC_reml <- BIC(relmm.bond))

## Compare covariance estimates ####
# Define X and Z design matrices from model Y = X*Beta + Z*theta + e_iid_noise
matrix.x <- model.matrix(lmm.bond)
matrix.z <- matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1),
                   21, 7)


# Create functions to calculate estimates using MLE / REML
loglikef <- function(x) {
  vector.Y <- as.vector(df.bond$pres)
  matrix.G <- x[1] * diag(1,nrow = 7)
  matrix.R <- x[2] * diag(1,nrow = 21)
  matrix.V <- matrix.z %*% matrix.G %*% t(matrix.z) + matrix.R
  vector.Beta <- solve(t(matrix.x) %*% solve(matrix.V) %*% matrix.x) %*% t(matrix.x) %*% solve(matrix.V) %*% vector.Y
  loglike <- log(det(matrix.V)) + t(vector.Y - matrix.x %*% vector.Beta) %*% solve(matrix.V) %*% (vector.Y - matrix.x %*% vector.Beta)
  return(c(loglike))
}


reloglikef <- function(x) {
  vector.Y <- as.vector(df.bond$pres)
  matrix.G <- x[1] * diag(1,nrow = 7)
  matrix.R <- x[2] * diag(1,nrow = 21)
  matrix.V <- matrix.z %*% matrix.G %*% t(matrix.z) + matrix.R
  vector.Beta <- solve(t(matrix.x) %*% solve(matrix.V) %*% matrix.x) %*% t(matrix.x) %*% solve(matrix.V) %*% vector.Y
  loglike <- log(det(matrix.V)) + log(det(t(matrix.x) %*% solve(matrix.V) %*% matrix.x)) + t(vector.Y - matrix.x %*% vector.Beta) %*% solve(matrix.V) %*% (vector.Y - matrix.x %*% vector.Beta)
  return(c(loglike))
}

# Calculate MLE / REML using our functions defined above and compare to lmm / relmm models.
# Note that both methods produce estimates very close to the lmm / relmm model values.

print(noquote(c("MLE model values: "))) 
print(summary(lmm.bond)$varcor)
print(noquote(c("Maximum Likelihood Estimates: ",sqrt(optim_nm(fun = loglikef, k = 2, start = c(1,1),maximum= FALSE, tol = 0.0000000000001)$par))))

print(noquote(c("REML model values: "))) 
print(summary(relmm.bond)$varcor)
print(noquote(c("Restricted / Residual Maximum Likelihood Estimates: ",sqrt(optim_nm(fun = reloglikef, k = 2, start = c(5,6),maximum= FALSE, tol = 0.0000000000001)$par))))




#  Using nlme library ###############################################################

# install.packages('nlme')
library(nlme)


df.bond <- data.frame(ingot = rep(1:7, each = 3),
                      metal = rep(c("n", "i", "c"), 7),
                      pres = c(67, 71.9, 72.2,
                               67.5, 68.8, 66.4,
                               76, 82.6, 74.5,
                               72.7, 78.1, 67.3,
                               73.1, 74.2, 73.2,
                               65.8, 70.8, 68.7,
                               75.6, 84.9, 69))
df.bond$ingot <- factor(df.bond$ingot)

# Linear mixed model by MLE
lmm.bond_mle <- lme(pres ~ metal, random = ~1 | ingot, data = df.bond, method = "ML")

# Linear mixed model by REML
lmm.bond_reml <- lme(pres ~ metal, random = ~1 | ingot, data = df.bond, method = "REML")

# Compare AIC/BIC
AIC_mle <- AIC(lmm.bond_mle)
BIC_mle <- BIC(lmm.bond_mle)
AIC_reml <- AIC(lmm.bond_reml)
BIC_reml <- BIC(lmm.bond_reml)

# Print AIC and BIC values
cat("MLE AIC:", AIC_mle, "\n")
cat("MLE BIC:", BIC_mle, "\n")
cat("REML AIC:", AIC_reml, "\n")
cat("REML BIC:", BIC_reml, "\n")

# Summary of MLE model
summary(lmm.bond_mle)

# Summary of REML model
summary(lmm.bond_reml)
