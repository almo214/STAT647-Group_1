# Group 1 REML Code adapted from:
# Source -> https://rh8liuqy.github.io/Example_Linear_mixed_model.html


# Using lme4 library ####
# install.packages('lme4')
# install.packages('optimization')
library(lme4)
library(optimization)

# [TODO] Build Function  

## Define the data ####
df <- data.frame(ingot = rep(1:7,each = 3),
                      metal = rep(c("n","i","c"),7),
                      pres = c(67,71.9,72.2,
                               67.5,68.8,66.4,
                               76,82.6,74.5,
                               72.7,78.1,67.3,
                               73.1,74.2,73.2,
                               65.8,70.8,68.7,
                               75.6,84.9,69))

Y <- 'pres'
rfact <- 'ingot'

# Scroll down for example function run

RMLE_Funct <- function(df, Y, rfact, rfact_int = 1, sum_print = FALSE){
  # df = dataframe with all data
  # Y = name of the target variable column
  # rfact = name of the random factor variable
  # rfact_int = number of random factor intercepts

  generate_matrix_z <- function(df, rfact) {
    unique_levels <- unique(df[[rfact]])  # Get unique levels of the random factor
    num_observations <- nrow(df)
    num_levels <- length(unique_levels)
  
    # Initialize an empty matrix.z
    matrix.z <- matrix(0, nrow = num_observations, ncol = num_levels)
  
    # Populate matrix.z based on the grouping structure
    for (i in 1:num_levels) {
      matrix.z[df[[rfact]] == unique_levels[i], i] <- 1
    }
  
    return(matrix.z)
  }
  
  # Given a dynamic input data frame, convert any non-numeric columns to a factor
  to_factor <- names(df)[!sapply(df,is.numeric)]
  df[to_factor] <- lapply(df[to_factor], factor)
  
  # Develop the equation for the linear mixed models
  # List variables that are not the target variable or the random factor variable
  independent_vars <-  colnames(df)[!colnames(df) %in% c(Y, rfact)]
  
  # Create the linear mixed model formula to dynamically accept any dataframe
  formula_str <- paste(Y, "~", paste(independent_vars, collapse = " + "), "+ (1|", rfact, ")")

  # Fit the linear mixed model
  lmm.mle <- lmer(as.formula(formula_str), REML = FALSE, data = df)
  lmm.reml <- lmer(as.formula(formula_str), REML = TRUE, data = df)
  
  # Compare AIC/BIC:  Note REML provides lowest AIC and BIC
  AIC_mle <- AIC(lmm.mle)
  BIC_mle <- BIC(lmm.mle)
  AIC_reml <- AIC(lmm.reml)
  BIC_reml <- BIC(lmm.reml)
  
  # Print AIC and BIC values
  cat("MLE AIC:", AIC_mle, "\n")
  cat("REML AIC:", AIC_reml, "\n")
  cat("MLE BIC:", BIC_mle, "\n")
  cat("REML BIC:", BIC_reml, "\n")
  cat("\n")
  # Print model summaries, if user opts to print
  if (sum_print == TRUE){
    print(summary(lmm.mle))
    print(summary(lmm.reml))
  }
  
  ## Compare covariance estimates ####
  # Define X and Z design matrices from model Y = X*Beta + Z*theta + e_iid_noise
  matrix.x <- model.matrix(lmm.mle)
  
  ###### [TODO] Update with dynamic creation of matrix.z. Work in progress below (outside function)
  # matrix.z <- matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  #                    0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  #                    0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
  #                    0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,
  #                    0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,
  #                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
  #                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1),
  #                  21, 7)
  
  # Create functions to calculate estimates using MLE / REML
  ###### [TODO] Update funcitons with dynamic values based on input
  loglikef <- function(x) {
    vector.Y <- as.vector(df$pres)
    matrix.G <- x[1] * diag(1, nrow = ncol(matrix.z))
    matrix.R <- x[2] * diag(1, nrow = nrow(matrix.z))
    matrix.V <- matrix.z %*% matrix.G %*% t(matrix.z) + matrix.R
    vector.Beta <- solve(t(matrix.x) %*% solve(matrix.V) %*% matrix.x) %*% t(matrix.x) %*% solve(matrix.V) %*% vector.Y
    loglike <- log(det(matrix.V)) + t(vector.Y - matrix.x %*% vector.Beta) %*% solve(matrix.V) %*% (vector.Y - matrix.x %*% vector.Beta)
    return(c(loglike))
  }
  
  
  reloglikef <- function(x) {
    vector.Y <- as.vector(df$pres)
    matrix.G <- x[1] * diag(1, nrow = ncol(matrix.z))
    matrix.R <- x[2] * diag(1, nrow = nrow(matrix.z))
    matrix.V <- matrix.z %*% matrix.G %*% t(matrix.z) + matrix.R
    vector.Beta <- solve(t(matrix.x) %*% solve(matrix.V) %*% matrix.x) %*% t(matrix.x) %*% solve(matrix.V) %*% vector.Y
    loglike <- log(det(matrix.V)) + log(det(t(matrix.x) %*% solve(matrix.V) %*% matrix.x)) + t(vector.Y - matrix.x %*% vector.Beta) %*% solve(matrix.V) %*% (vector.Y - matrix.x %*% vector.Beta)
    return(c(loglike))
  }

  
  # Calculate MLE / REML using our functions defined above and compare to lmm / relmm models.
  # Note that both methods produce estimates very close to the lmm / relmm model values.
  
  print(noquote(c("MLE model values: "))) 
  print(summary(lmm.mle)$varcor)
  
  # Note that the value x applied to the function loglikef is initiliazed using start = in the optim_nm function.
  print(noquote(c("Maximum Likelihood Estimates: ",sqrt(optim_nm(fun = loglikef, k = 2, start = c(1,1),maximum= FALSE, tol = 0.0000000000001)$par))))
  
  
  print(noquote(c("REML model values: "))) 
  print(summary(lmm.reml)$varcor)
  print(noquote(c("Restricted / Residual Maximum Likelihood Estimates: ",sqrt(optim_nm(fun = reloglikef, k = 2, start = c(5,6),maximum= FALSE, tol = 0.0000000000001)$par))))
  

}




# [TODO] Automate the creation of matrix.z based on dynamic input, work in progress
identity_mat <- diag(10)
duplicated_matrix.z <- matrix(0, nrow = nrow(identity_mat), ncol = 3 * ncol(identity_mat))
for (columns in ncol(identity_mat)) {
  start_col <- (i - 1) * ncol(identity_mat) + 1
  end_col <- i * ncol(identity_mat)
  duplicated_matrix.z[, start_col:end_col] <- identity_mat
}



# Run the funciton to produce AICs, BICs, and estimator comparisons ####
RMLE_Funct(df, Y, rfact, sum_print =  FALSE)




#  Using nlme library ###############################################################

# install.packages('nlme')
library(nlme)


df <- data.frame(ingot = rep(1:7, each = 3),
                      metal = rep(c("n", "i", "c"), 7),
                      pres = c(67, 71.9, 72.2,
                               67.5, 68.8, 66.4,
                               76, 82.6, 74.5,
                               72.7, 78.1, 67.3,
                               73.1, 74.2, 73.2,
                               65.8, 70.8, 68.7,
                               75.6, 84.9, 69))
df$ingot <- factor(df$ingot)

# Linear mixed model by MLE
lmm.mle_mle <- lme(pres ~ metal, random = ~1 | ingot, data = df, method = "ML")

# Linear mixed model by REML
lmm.mle_reml <- lme(pres ~ metal, random = ~1 | ingot, data = df, method = "REML")

# Compare AIC/BIC
AIC_mle <- AIC(lmm.mle_mle)
BIC_mle <- BIC(lmm.mle_mle)
AIC_reml <- AIC(lmm.mle_reml)
BIC_reml <- BIC(lmm.mle_reml)

# Print AIC and BIC values
cat("MLE AIC:", AIC_mle, "\n")
cat("MLE BIC:", BIC_mle, "\n")
cat("REML AIC:", AIC_reml, "\n")
cat("REML BIC:", BIC_reml, "\n")

# Summary of MLE model
summary(lmm.mle_mle)

# Summary of REML model
summary(lmm.mle_reml)
