##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to process the residuals from ARIMA models and construct a logarithmic transformation
# of the outcome variable for regression analysis. This involves creating a matrix of residuals,
# performing matrix multiplication to obtain interaction terms, and then applying logarithmic
# transformation on eigenvalues.

# Input:
# - resi.acf.list: A list containing residuals for each region of interest.
# - length: The number of time points in the residuals data.

# Output:
# - A matrix of log-transformed outcomes for LS regression.

LS.construct.logY <- function(resi.acf.list, length) {
  # Create a matrix to hold residuals for all regions across time points
  y_v_t <- matrix(nrow = 100, ncol = length)
  
  # Populate the matrix with residuals
  for (i in 1:100) {
    y_v_t[i, ] <- resi.acf.list[[i]] # y_v_t is 100 by T
  }
  
  # Construct a 3D array to hold the product of residuals across time
  vectorbycol <- vector("numeric")
  for (i in 1:length) {
    single.vectorbycol <- as.vector(y_v_t[, i] %*% t(y_v_t[, i]))
    vectorbycol <- c(vectorbycol, single.vectorbycol)
  }
  array.residuals <- array(vectorbycol, dim = c(100, 100, length)) # 100 by 100 by T
  
  # Logarithmically transform the first eigenvalue of each time point
  log.version.vector <- vector("numeric")
  for (i in 1:length) {
    eigen_info <- eigen(array.residuals[,,i])
    # there is only one positive eigenvalue, set the first to logarithm, and the second to end = 0
    eigen_info$values[1] <- log(eigen_info$values[1])
    eigen_info$values[-1] <- 0
    # contruct the 
    log.version.vector <- c(log.version.vector, as.vector(eigen_info$vectors %*% diag(eigen_info$values) %*% t(eigen_info$vectors)))
  }
  log.Y.arrary <- array(log.version.vector, dim = c(100, 100, length))
  
  # Flatten the array to prepare for regression analysis
  log.Y.outcome.vector <- vector("numeric")
  for (i in 1:100) {
    for (j in 1:i) {
      single.vector <- vector()
      for (t in 1:length) {
        scalar <- log.Y.arrary[i,j,t]
        single.vector <- c(single.vector, scalar)
      }
      log.Y.outcome.vector <- c(log.Y.outcome.vector, single.vector)
    }
  }
  # construct T by 5050 outcome matrix
  log.Y.forLS <- array(log.Y.outcome.vector, dim = c(length,5050))
  
  return(log.Y.forLS)
}

# Example of usage:
# residuals_list <- list(replicate(671, rnorm(100)))  # Example list of residuals
# log_Y_matrix <- LS.construct.logY(residuals_list, 671)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to perform Least Squares (LS) regression using robust estimation techniques
# on the log-transformed outcomes. The function fits models for each outcome and returns regression
# coefficients along with robust covariance matrices.

# Input:
# - log.Y.forLS: A matrix of log-transformed outcomes (671 by 5050).
# - totalcovariates.scale: Scaled covariates for inclusion in the regression model.

# Output:
# - A list containing regression results for all models.

LS.robust.estimation <- function(log.Y.forLS, totalcovariates.scale, length) {
  
  # Initialize an empty list to store regression results
  robust.regression.all <- list()
  
  # Perform OLS regression with robust standard errors for each model
  for (i in 1:ncol(log.Y.forLS)) {
    ols.model <- lm(log.Y.forLS[, i] ~ ., data = totalcovariates.scale)
    # robust coef
    robust.ols <- coeftest(ols.model, vcov = vcovHC(ols.model, type = "HC3"))
    # covariance matrix
    cov.blink.fix <- vcovHC(ols.model, type = "HC3")
    # store each coefs in the list
    robust.regression.all[[i]] <- list(robust.ols, cov.blink.fix)
  }
  
  return(robust.regression.all)
}

# Example of usage:
# robust_results <- LS.robust.estiomation(log_y_data, scaled_covariates, 671)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to calculate the z-statistic for differences between eye-blink and eye-fixation
# coefficients across all pairs of regions, and construct a lower triangular matrix
# of the results.

# Input:
# - robust.regression.all: A list of regression results, where each element contains
#   robust coefficients and covariance matrices for each region.

# Output:
# - A lower triangular matrix storing t-statistics for the comparisons between eye-blink
#   and eye-fixation across all regions.

Z.statistic.compare <- function(robust.regression.all) {
  # Initialize a vector to store the t-statistics for the comparison between eye-blink and eye-fixation
  t.stat.blink.fixation <- vector()
  
  # Loop through each regression result to calculate the t-statistics
  for (i in 1:length(robust.regression.all)) {
    # Extract coefficients and covariance information
    robust.coef <- robust.regression.all[[i]][[1]]
    cov.1.2 <- robust.regression.all[[i]][[2]][2,3]
    
    # Calculate the t-statistic for the difference between eye-blink and eye-fixation coefficients
    t.stat <- (robust.coef[3,1] - robust.coef[2,1]) / sqrt(robust.regression.all[[i]][[2]][3,3] +
                                                             robust.regression.all[[i]][[2]][2,2] - 2 * cov.1.2)
    t.stat.blink.fixation <- c(t.stat.blink.fixation, t.stat)
  }
  
  # Construct a lower triangular matrix to store the t-statistics
  diff.blink.fixation <- matrix(0, nrow = 100, ncol = 100)
  index <- 1  # Initialize an index to track the position in the vector
  
  # Populate the lower triangular part of the matrix with t-statistics
  for (i in 1:100) {
    for (j in 1:i) {
      diff.blink.fixation[i, j] <- t.stat.blink.fixation[index]
      index <- index + 1
    }
  }
  
  return(diff.blink.fixation)
}

# Example of usage:
# Assuming 'robust_regression_results' is your list of robust regression outcomes:
# comparison_matrix <- Z.statistic.compare(robust_regression_results)
