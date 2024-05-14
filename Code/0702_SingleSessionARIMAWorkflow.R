##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
# Note: Some subject only has one available task data
#####################################################################################
# Function to execute the entire ARIMA modeling workflow for subjects with only one recorded task session.
# This function handles fMRI and Eye Tracking data processing, constructs a design matrix, and performs ARIMA modeling.

# Input:
# - path_fMRI_ses1: Path to the fMRI data file (can be from any session, not limited to session 1).
# - path_HC_ses1: Path to the Head Motion confounder data file.
# - path_ET_ses1: Path to the Eye Tracking data file.

# Output:
# - A list containing ARIMA model results and the scaled design matrix used in the analysis.

Everything_to_ARIMAoutput.ONEses <- function(path_fMRI_ses1, path_HC_ses1, path_ET_ses1) {
  # Re-read the fMRI data for input into subsequent analysis steps
  xii1 <- read_xifti(path_fMRI_ses1, brainstructures = "all")
  
  # Process fMRI data to get the mean time series
  pmean.ses1 <- fMRI_xii_pmean(dtseries_data = xii1)
  
  # Process Eye Tracking data and perform convolution for blink and fixation data
  input.conv.blink.ses1 <- ETascDataProcess(file_path = path_ET_ses1, blink = TRUE)
  convolution.timeseries.blink.ses1 <- Convolution_function(totaltime = input.conv.blink.ses1[[1]],
                                                            onsets = input.conv.blink.ses1[[2]], 
                                                            duration = input.conv.blink.ses1[[3]], 
                                                            accuracy = input.conv.blink.ses1[[4]])
  sub.convolution.timeseries.blink.ses1 <- Extraction_ETtime(conv.data = convolution.timeseries.blink.ses1, xii.mean = pmean.ses1, interval = 1.127)
  
  input.conv.fixation.ses1 <- ETascDataProcess(path = path_ET_ses1, fixation = TRUE)
  convolution.timeseries.fixation.ses1 <- Convolution_function(totaltime = input.conv.fixation.ses1[[1]],
                                                               onsets = input.conv.fixation.ses1[[2]], 
                                                               duration = input.conv.fixation.ses1[[3]], 
                                                               accuracy = input.conv.fixation.ses1[[4]])
  sub.convolution.timeseries.fixation.ses1 <- Extraction_ETtime(conv.data = convolution.timeseries.fixation.ses1, xii.mean = pmean.ses1, interval = 1.127)
  
  # Process head motion confounders
  head.confounder.ses1 <- HeadMotionConfounder_process(path = path_HC_ses1, xii.mean = pmean.ses1)
  
  # Construct the final design matrix for the session
  Scale.design.matrix.ses1 <- DesignMatrix_process(ET.covariates = cbind(sub.convolution.timeseries.blink.ses1, sub.convolution.timeseries.fixation.ses1), 
                                                   head.confounder.ses1)
  
  # Perform ARIMA modeling
  arima.ses1 <- ARIMAmodel(xii_pmean = pmean.ses1, design.matrix = Scale.design.matrix.ses1)
  
  return(list(arima.ses1, Scale.design.matrix.ses1))
}

# Example of usage:
# result <- Everything_to_ARIMAoutput.ONEses(path_fMRI_ses1 = "path_to_fMRI_data.nii", path_HC_ses1 = "path_to_head_confounder_data.csv", path_ET_ses1 = "path_to_eye_tracking_data.csv")
