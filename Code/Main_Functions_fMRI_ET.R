#########################################################################################################
#########################################################################################################
################################# Functions supporting the project ######################################
#########################################################################################################
#########################################################################################################

# Define a function to compute the mean fMRI time series for specified brain regions.
# The function calculates the mean across vertices for each time point.
# Inputs:
# - dtseries_data: the raw .dtseries.nii MRI data to be processed.
# - brain_parcellation: the label of the parcellation scheme, e.g., 'Schaefer_100', 'Schaefer_400', etc.
# - region_count: the total number of cortical brain regions in the parcellation 
#   (this does not include subcortex regions, which are assumed to be 19).

# Output:
# - A matrix with brain regions as rows and time points as columns (region by time).

fMRI_xii_pmean <- function(dtseries_data, brain_parcellation = "Schaefer_100", region_count = 100) {
  # Load the brain parcellation data
  parc <- load_parc(brain_parcellation)
  
  # Convert parcellation data into a vector, indicating the associated brain region
  parc_vec <- c(as.matrix(parc))
  
  # Adjust for non-cortical vertices to align with the parcellation scheme
  adjusted_data <- move_from_mwall(dtseries_data, NA)
  
  # Convert the adjusted data into a matrix format (voxel by time)
  xii_mat <- as.matrix(adjusted_data)
  
  # Retrieve and process the labels for subcortex regions
  subcortex_labels <- adjusted_data$meta$subcort$labels
  
  # Combine cortical and subcortex labels
  combined_labels <- c(as.character(parc_vec), as.character(subcortex_labels))
  
  # Generate labels for brain regions including cortex and subcortex
  labels_vector <- c(as.character(1:region_count), levels(subcortex_labels)[3:21])
  
  # Initialize a matrix to store the computed mean fMRI signal values
  regional_means <- matrix(nrow = region_count + 19, ncol = ncol(xii_mat))
  
  # Calculate mean signals for each parcellation region across time points
  for (label in labels_vector) {
    region_data <- xii_mat[combined_labels == label, ]
    regional_means[which(labels_vector == label), ] <- colMeans(region_data, na.rm = TRUE)
  }
  
  return(regional_means)
}

#########################################################################################################
#########################################################################################################
# Define a function to process eye-tracking data from an .ASC file and extract
# either blink or fixation data, creating a timestamped dataset for further analysis.

# Inputs:
# - file_path: The path to the eye-tracking .ASC file.
# - blink: Logical flag to indicate if blink data should be extracted.
# - fixation: Logical flag to indicate if fixation data should be extracted.

# Output:
# - A list containing the total time of the eye-tracking task, onsets, durations,
#   and sampling rate accuracy.

ETascDataProcess <- function(file_path, blink = FALSE, fixation = FALSE) {
  # Read the raw .asc data
  ETasc <- read.asc(fname = file_path, samples = TRUE, events = TRUE)
  # Initialize an empty container for the extracted type
  ETtypein <- NULL
  
  # Extract blink or fixation data and calculate real time
  if (blink) {
    ETtypein <- ETasc$blinks %>% 
      dplyr::mutate(ts = (stime - min(ETasc$raw$time)) / 1e3,
                    te = (etime - min(ETasc$raw$time)) / 1e3,
                    duration = dur / 1e3)
  } else if (fixation) {
    ETtypein <- ETasc$fix %>% 
      dplyr::mutate(ts = (stime - min(ETasc$raw$time)) / 1e3,
                    te = (etime - min(ETasc$raw$time)) / 1e3,
                    duration = dur / 1e3)
  } else {
    stop("Please specify blink or fixation data to process.")
  }
  
  # Convert the results to a data.frame
  ETtype <- as.data.frame(ETtypein)
  # Calculate the total time, onsets, and durations for further analysis
  totaltime <- (max(ETasc$raw$time) - min(ETasc$raw$time)) / 1e3
  onsets <- ETtype$ts
  durations <- ETtype$duration
  sampling_rate <- 0.002
  
  return(list(totaltime, onsets, durations, sampling_rate))
}

# Example of usage:
# et_data <- ETascDataProcess(file_path = "path/to/data.asc", blink = TRUE, fixation = FALSE)
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Define a function to perform the convolution of eye-tracking time series data
# with a hemodynamic response function (HRF).

# Inputs:
# - totaltime: Total duration of the stimulus experiment in seconds.
# - onsets: A vector containing the start times for each stimulus event.
# - durations: A vector containing the duration of each stimulus event in seconds.
# - sampling_rate: The sampling rate for the experiment (e.g., 0.002 sec for 500 Hz).

# Output:
# - A list containing the convolution result and the corresponding time vector.

Convolution_function <- function(totaltime, onsets, durations, sampling_rate) {
  # Generate a stimulus boxcar function based on the inputs
  stimulus <- stimfunction(totaltime = totaltime, onsets = onsets, durations = durations, accuracy = sampling_rate)
  # Calculate the time vector for the length of the experiment
  time_vector <- seq(0, totaltime, by = sampling_rate)
  
  # Perform convolution with the stimulus and HRF
  # (Ensure 'stimulus' and 'newhrf' have the same length for later convolution)
  stimulus_padded <- c(stimulus, rep(0, length(time_vector) - length(stimulus)))
  hrf_kernel <- canonicalHRF(time_vector, verbose = FALSE) # Double-gamma Function
  hrf_padded <- c(hrf_kernel, rep(0, length(stimulus)))
  
  # Compute the frequency domain representation via FFT
  fft_stimulus <- fft(stimulus_padded)
  fft_hrf <- fft(hrf_padded)
  # Multiply element-wise in the frequency domain and perform inverse FFT
  convolved_signal <- Re(fft(fft_stimulus * fft_hrf, inverse = TRUE) / length(stimulus))
  # Trim the convolution result to the original time vector length
  convolved_signal <- convolved_signal[1:length(time_vector)]
  
  return(list(convolved_signal, time_vector))
}

# Example of usage:
# convolution_result <- Convolution_function(totaltime = 780, onsets = c(1, 100, 200), durations = c(1, 1, 1),  accuracy = 0.002)
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Define a function to extract time series data from convolved eye-tracking (ET) data
# and align it with fMRI data time points.

# Inputs:
# - conv_data: Convolution result from the ET data processed with HRF.
# - fmri_data: Pre-processed fMRI time series data obtained from earlier analysis.
# - tr: Repetition time or the time interval of fMRI data acquisition (in seconds).

# Output:
# - A vector of the standardized extracted ET time series aligned with fMRI data.

Extraction_ETtime <- function(conv_data, fmri_data, tr = 1.127) {
  # Calculate the indices to extract from the convolved data based on fMRI TR
  indices <- round(seq(from = 1, to = length(conv_data[[1]]), by = tr / 0.002), digits = 0)
  
  # Extract the convolved data at the specified indices
  extracted_data <- conv_data[[1]][indices]
  
  # Ensure the extracted ET time series length matches the fMRI time series
  if (ncol(fmri_data) > length(extracted_data)) {
    # Pad the end with zeros if the fMRI data is longer
    extracted_data <- c(extracted_data, rep(0, ncol(fmri_data) - length(extracted_data)))
  } else {
    # Truncate the extracted data if it is longer than the fMRI data
    extracted_data <- extracted_data[1:ncol(fmri_data)]
  }
  
  # Standardize the extracted time series
  standardized_data <- extracted_data / max(extracted_data)
  
  return(standardized_data)
}

# Example of usage:
# Assuming 'convolved_hrf_data' is your convolved ET data and 'fmri_time_series' is the fMRI time series
# standardized_et_series <- Extraction_ETtime(conv_data = convolved_hrf_data, fmri_data = fmri_time_series, tr = 1.127)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to process confounder covariates and generate the design matrix used in the analysis.
# This function reads the confounder data from the specified path and computes quadratic terms for selected variables.

# Input:
# - path: Path to the file containing head motion data.
# - xii.mean: Reference matrix to align the number of rows of confounder data.

# Output:
# - Data frame with six head motion confounders and their quadratic terms, sized to match xii.mean.

HeadMotionConfounder_process <- function(path, xii.mean) {
  # Read the confounder data
  scrub <- read.table(file = path, header = TRUE)
  # Select the required six confounder variables
  six.confounder <- scrub[, c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z")]
  # Create quadratic terms for each confounder
  six.confounder$trans_x_squared <- six.confounder$trans_x^2
  six.confounder$trans_y_squared <- six.confounder$trans_y^2
  six.confounder$trans_z_squared <- six.confounder$trans_z^2
  six.confounder$rot_x_squared <- six.confounder$rot_x^2
  six.confounder$rot_y_squared <- six.confounder$rot_y^2
  six.confounder$rot_z_squared <- six.confounder$rot_z^2
  
  # Align the confounder data to the size of the fMRI data matrix
  six.confounder <- six.confounder[1:ncol(xii.mean),]
  
  return(six.confounder)
}

# Example usage:
# adjusted_confounders <- HeadMotionConfounder_process(path = "path/to/data.tsv", xii.mean = fmri_matrix)

# NOTE!!!!!This part can be changed based on the specific input file.

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# Function to construct the design matrix for statistical analysis by integrating eye-tracking
# covariates with head motion data, and adding session indicators and interaction terms.

# Input:
# - ET.covariates: Time series data from eye-tracking. NOTE! The input is T by 2, i.e., ET blink and ET fixation.
# - head.motion.covariate.ses1: Head motion data for session 1.
# - head.motion.covariate.ses2: Head motion data for session 2.

# Output:
# - Scaled and combined design matrix with covariates and interaction terms for analysis.

DesignMatrix_process <- function(ET.covariates, head.motion.covariate.ses1, head.motion.covariate.ses2) {
  # Create session indicators and time index
  indicator <- c(rep(0, nrow(head.motion.covariate.ses1)), rep(1, nrow(head.motion.covariate.ses2)))
  indextime <- 1:(nrow(head.motion.covariate.ses1) + nrow(head.motion.covariate.ses2))
  
  # Combine the eye-tracking and head motion data with additional covariates
  covariate.data <- cbind(ET.covariates, rbind(head.motion.covariate.ses1, head.motion.covariate.ses2), indicator, indextime)
  # Add interaction terms for time and session indicators 
  covariate.data = covariate.data %>% mutate(
    tx_indicator = trans_x * indicator,
    ty_indicator = trans_y * indicator,
    tz_indicator = trans_z * indicator,
    rx_indicator = rot_x * indicator,
    ry_indicator = rot_y * indicator,
    rz_indicator = rot_z * indicator,
    tx_square_indicator = trans_x_squared * indicator,
    ty_square_indicator = trans_y_squared * indicator,
    tz_square_indicator = trans_z_squared * indicator,
    rx_square_indicator = rot_x_squared * indicator,
    ry_square_indicator = rot_y_squared * indicator,
    rz_square_indicator = rot_z_squared * indicator
  )
  
  # Scale the covariates within each task session
  totalcovariates.scale <- covariate.data %>% group_by(indicator) %>%
    dplyr::mutate(across(trans_x:rz_square_indicator, ~scale(., center = mean(., na.rm = TRUE), scale = sd(., na.rm = TRUE))))
  # Convert to data.frame and replace all NA values resulting from scaling with zeros
  totalcovariates.scale <- as.data.frame(totalcovariates.scale)
  totalcovariates.scale[is.na(totalcovariates.scale)] <- 0
  
  return(totalcovariates.scale)
}

# Example usage:
# design_matrix <- DesignMatrix_process(ET_covariates, head_motion_data_session1, head_motion_data_session2)

##########################################################################################################################################################################
##########################################################################################################################################################################
# for only have one ses data
DesignMatrix_process.ses1 <- function(ET.covariates = ET.covariates, head.motion.covariate = head.motion.covariate) {
  # create two more columns, indicator for movie session (based on fMRI) and time index
  indextime <- c(1:nrow(head.motion.covariate))
  # just use cbind to merge the Et data and 12 head motion confounders together
  covariate.data <- cbind(ET.covariates, head.motion.covariate, indextime)
  # create the other interaction covariate terms
  covariate.data$time_square <- covariate.data$indextime^2
  
  ########### scaling
  # scale all the covariates within each task session
  totalcovariates.scale <- covariate.data %>% 
    mutate(across(trans_x:time_square, ~scale(., center = mean(., na.rm = TRUE), scale = sd(., na.rm = TRUE))))
  # convert to data.frame
  totalcovariates.scale <- as.data.frame(totalcovariates.scale)
  # Replace all NAs with 0 in the entire dataframe, these 0 is due to indicator for session1 = 0, scale these o gives NAs
  totalcovariates.scale[is.na(totalcovariates.scale)] <- 0
  return(totalcovariates.scale)
}
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

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

ARIMAmodel.ONEses <- function(xii_pmean, design.matrix, num = 100) {
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
    estimate.region.blink[i, ] <- coeftest(arima.region)["sub.convolution.timeseries.blink.ses1",]
    estimate.region.fixation[i, ] <- coeftest(arima.region)["sub.convolution.timeseries.fixation.ses1",]
    
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

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
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
      for (t in 1:length) {
        log.Y.outcome.vector <- c(log.Y.outcome.vector, log.Y.arrary[i, j, t])
      }
    }
  }
  # construct T by 5050 outcome matrix
  log.Y.arrary <- matrix(log.Y.outcome.vector, ncol = 5050, byrow = TRUE)
  
  return(log.Y.arrary)
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
# - log.Y.array: A matrix of log-transformed outcomes (671 by 5050).
# - totalcovariates.scale: Scaled covariates for inclusion in the regression model.

# Output:
# - A list containing regression results for all models.

LS.robust.estiomation <- function(log.Y.array, totalcovariates.scale, length) {
  # Add intercept column to the design matrix
  design.matrix <- cbind(rep(1, length), totalcovariates.scale)
  
  # Initialize an empty list to store regression results
  robust.regression.all <- vector("list", ncol(log.Y.array))
  
  # Perform OLS regression with robust standard errors for each model
  for (i in 1:ncol(log.Y.array)) {
    ols.model <- lm(log.Y.array[, i] ~ ., data = as.data.frame(design.matrix))
    robust.ols <- coeftest(ols.model, vcov = vcovHC(ols.model, type = "HC3"))
    cov.blink.fix <- vcovHC(ols.model, type = "HC3")
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
      diff.blink.fixation[i, j] := t.stat.blink.fixation[index]
      index <- index + 1
    }
  }
  
  return(diff.blink.fixation)
}

# Example of usage:
# Assuming 'robust_regression_results' is your list of robust regression outcomes:
# comparison_matrix <- Z.statistic.compare(robust_regression_results)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

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


# Define a function to map statistical data (e.g., mean signals, coefficients, or t-statistics)
# onto a brain graph. This is often used in neuroimaging to visualize data on a brain model.

# Inputs:
# - dtseries_data: Raw .dtseries.nii MRI data.
# - mapping_data: The statistical data to be mapped onto the brain graph.
#                 Must be of the same length as the number of brain regions.
# - brain_parcellation: The parcellation scheme used for brain region delineation.
# - region_count: The number of brain regions (excluding the 19 subcortex regions).

# Output:
# - A brain graph object with the mapped statistical data.

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################


fMRIBrain_Mapping <- function(dtseries_data, mapping_data, brain_parcellation = "Schaefer_100", region_count = 100) {
  # Load the parcellation data
  parc <- load_parc(brain_parcellation)
  
  # Convert the parcellation to a vector to identify the brain region of each voxel
  parc_vec <- c(as.matrix(parc))
  
  # Adjust subcortex labels to match cortical labeling
  subcortex_labels_adjusted <- as.numeric(dtseries_data$meta$subcort$labels) - 2
  subcortex_labels_adjusted <- region_count + subcortex_labels_adjusted
  
  # Create a combined brain vector of cortical and adjusted subcortex labels
  brain_vec <- c(parc_vec, subcortex_labels_adjusted)
  
  # Adjust the fMRI data for non-cortical vertices to align with parcellation
  adjusted_data <- move_from_mwall(dtseries_data, NA)
  
  # Initialize fMRI data with empty values for subsequent mapping
  prepared_data <- select_xifti(adjusted_data * 0, 1)
  
  # If the data represents p-values, perform -log10 transformation for visualization
  # Uncomment the following line if necessary
  # transformed_data <- -log10(mapping_data)
  
  # Map the statistical data onto the brain vector
  brain_graph <- newdata_xifti(prepared_data, c(NA, mapping_data)[brain_vec + 1])
  
  # Return the brain graph object for visualization
  return(brain_graph)
}

# Example of usage:
# brain_graph_result <- fMRIBrain_Mapping(dtseries_data, statistical_data, 'Schaefer_100', 100)
