#####################################################################################
# Function to conduct model diagnostics for a randomly selected OLS model.
# This function assesses the residuals of the model to check assumptions like normality,
# autocorrelation, etc.

# Input:
# - log.Y.array: Matrix of log-transformed outcomes from function 5.
# - totalcovariates.scale: Scaled design matrix from function 3.2.
# - length: The number of observations in each model.

# Output:
# - A list of diagnostic plots including ACF, PACF, Q-Q plot, and histogram of residuals.

Residual.Diag_OLS <- function(log.Y.array, totalcovariates.scale, length) {
  # Randomly select indices for a single OLS model from the lower triangle of the matrix
  i <- sample(1:100, 1)
  j <- sample(1:i, 1)
  
  # Calculate the column index based on selected i and j
  index.col = ((i - 1) * i / 2) + j
  
  # Construct the design matrix including an intercept
  design.matrix <- as.matrix(cbind(rep(1, length), totalcovariates.scale))
  
  # Calculate the expected Y using the projection matrix formula
  E_Y_t_tY_t <- design.matrix %*% solve(t(design.matrix) %*% design.matrix) %*% t(design.matrix) %*% log.Y.array[, index.col]
  
  # Compute residuals
  residuals_At <- log.Y.array[, index.col] - E_Y_t_tY_t
  
  # Generate diagnostic plots
  acfplot <- acf(residuals_At, main = paste("ACF for Region Pair", i, "-", j))
  pacfplot <- pacf(residuals_At, main = paste("PACF for Region Pair", i, "-", j))
  qqplot <- qqnorm(residuals_At, main = paste("Q-Q plot for Region Pair", i, "-", j))
  hist <- hist(residuals_At, main = paste("Histogram of Residuals for Region Pair", i, "-", j), breaks = 30)
  
  return(list(ACF = acfplot, PACF = pacfplot, QQPlot = qqplot, Histogram = hist))
}

# Example of usage:
# Assuming 'log_y_data' is your log-transformed outcome matrix and 'scaled_covariates' is your scaled design matrix:
# diagnostic_plots <- Residual.Diag_OLS(log_y_data, scaled_covariates, 671)
