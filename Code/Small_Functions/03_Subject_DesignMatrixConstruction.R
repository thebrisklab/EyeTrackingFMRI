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
