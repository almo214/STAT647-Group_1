# Group 1 REML Code adapted from:
# Source -> https://rh8liuqy.github.io/Example_Linear_mixed_model.html


# Code Points of Contact:
# Allison Moore amoore214@tamu.edu
# Luis Quilarque lgquilarque@tamu.edu
# Tiffany Chang tiffchang@tamu.edu

# Using lme4 library ####
# install.packages('lme4')
# install.packages('optimization')
library(lme4)
library(optimization)

# EXAMPLE IMPLEMENTATION
# RMLE_Funct(df, Y, rfact, sum_print =  FALSE)


## Define the data ####
# Example df used, please replace with your own
df <- data.frame(ingot = rep(1:7,each = 3),
                      metal = rep(c("n","i","c"),7),
                      pres = c(67,71.9,72.2,
                               67.5,68.8,66.4,
                               76,82.6,74.5,
                               72.7,78.1,67.3,
                               73.1,74.2,73.2,
                               65.8,70.8,68.7,
                               75.6,84.9,69))

# You must pass Y, your target variable into our RMLE_Funct as a string.
Y <- 'pres'

# You must pass rfact, the random factor variable(s) as a string or list of strings.
rfact <- 'ingot'


# Function ####



# Function to compare MLE and REML
# df = dataframe with all data
# Y = name of the target variable column
# rfact = name of the random factor variable
# rfact_int = number of random factor intercepts

RMLE_Funct <- function(df, Y, rfact, rfact_int = 1, sum_print = FALSE){


  # User input validation
  if (!is.character(Y) | length(Y) >1) {
    # It's neither a string nor a list of strings
    print("Please submit your Y, target variable, as a single string.")
  }
  
  #User input validation
  if (!is.character(rfact) && !is.list(rfact)) {
    # It's neither a string nor a list of strings
    print("Please submit your rfact, random variable(s), as a string or list of strings.")
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
  
  # Function to generate the z matrix
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
  
  matrix.z <- generate_matrix_z(df, rfact)

  
  # Create functions to calculate estimates using MLE / REML
  loglikef <- function(x) {
    vector.Y <- as.vector(df[[Y]])
    matrix.G <- x[1] * diag(1, nrow = ncol(matrix.z))
    matrix.R <- x[2] * diag(1, nrow = nrow(matrix.z))
    matrix.V <- matrix.z %*% matrix.G %*% t(matrix.z) + matrix.R
    vector.Beta <- solve(t(matrix.x) %*% solve(matrix.V) %*% matrix.x) %*% t(matrix.x) %*% solve(matrix.V) %*% vector.Y
    loglike <- log(det(matrix.V)) + t(vector.Y - matrix.x %*% vector.Beta) %*% solve(matrix.V) %*% (vector.Y - matrix.x %*% vector.Beta)
    return(c(loglike))
  }
  
  
  reloglikef <- function(x) {
    vector.Y <- as.vector(df[[Y]])
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

# Run the function to produce AICs, BICs, and estimator comparisons ####
RMLE_Funct(df, Y, rfact, sum_print =  FALSE)

