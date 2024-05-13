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
  indices <- round(seq(from = 1, to = length(conv_data), by = tr / 0.002), digits = 0)
  
  # Extract the convolved data at the specified indices
  extracted_data <- conv_data[indices]
  
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
    mutate(across(trans_x:rz_square_indicator, ~scale(., center = mean(., na.rm = TRUE), scale = sd(., na.rm = TRUE))))
  # Convert to data.frame and replace all NA values resulting from scaling with zeros
  totalcovariates.scale <- as.data.frame(totalcovariates.scale)
  totalcovariates.scale[is.na(totalcovariates.scale)] <- 0
  
  return(totalcovariates.scale)
}

# Example usage:
# design_matrix <- DesignMatrix_process(ET_covariates, head_motion_data_session1, head_motion_data_session2)


##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
