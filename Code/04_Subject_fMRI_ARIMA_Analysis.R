#####################################################################################
# Function to perform Polynomial ARIMA regression on fMRI time series data. The function
# returns coefficients and residuals for eyeblink and eyefixation, along with covariance matrices
# for each region of interest.

# Inputs:
# - xii_pmean: fMRI time series data.
# - design.matrix: A pre-scaled design matrix used as regressors in the ARIMA model.
# - num: The number of regions of interest (excluding 19 subcortex regions).

# Output:
# - A list containing:
#   - Coefficients for eyeblink across all regions.
#   - Coefficients for eyefixation across all regions.
#   - Residuals from the ARIMA model.
#   - Covariance matrices for each region.

ARIMAmodel <- function(xii_pmean, design.matrix, num = 100) {
  # Initialize matrices to store the ARIMA model's output for blink and fixation coefficients
  estimate.region.blink <- matrix(data = NA, nrow = num, ncol = 4)
  estimate.region.fixation <- matrix(data = NA, nrow = num, ncol = 4)
  
  # Initialize lists to store residuals and covariance matrices for each region
  resi.acf <- vector("list", num)
  cov.matrix <- vector("list", num)
  
  # Loop through each region to fit the ARIMA model and store the outputs
  for (i in 1:num) {
    # Fit the ARIMA model using the specified order and external regressors
    arima.region <- arima(xii_pmean[i,], order = c(3,0,0), xreg = design.matrix)
    
    # Extract coefficients for eyeblink and eyefixation using their specific positions in the output
    estimate.region.blink[i, ] <- coeftest(arima.region)["blink.ses1.ses2",]
    estimate.region.fixation[i, ] <- coeftest(arima.region)["fixation.ses1.ses2",]
    
    # Store residuals and covariance matrix of the fitted model
    resi.acf[[i]] <- residuals(arima.region)
    cov.matrix[[i]] <- vcov(arima.region)
  }
  
  # Return a list containing all the results
  return(list(estimate.region.blink, estimate.region.fixation, resi.acf, cov.matrix))
}

# Example of usage:
# Assuming 'fMRI_time_series' is your fMRI data matrix and 'scaled_design_matrix' is your design matrix prepared earlier
# results <- ARIMAmodel(xii_pmean = fMRI_time_series, design.matrix = scaled_design_matrix, num = 100)
