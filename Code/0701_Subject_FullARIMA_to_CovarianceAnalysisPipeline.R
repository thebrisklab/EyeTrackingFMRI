##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
#######################################################################################.  Brain Activation!!!!!!!#################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################

# Comprehensive function to do the entire workflow from reading fMRI and Eye Tracking data,
# processing confounders, constructing a design matrix, and finally executing ARIMA modeling.
# This function integrates various preprocessing and modeling steps into a seamless pipeline.

# Input:
# - path_fMRI_ses1, path_fMRI_ses2: Paths to fMRI data files for two sessions.
# - path_ET_ses1, path_ET_ses2: Paths to Eye Tracking data files for two sessions.
# - path_HC_ses1, path_HC_ses2: Paths to Head Motion confounder data files for two sessions.

# Output:
# - A list containing ARIMA model results and the scaled design matrix used in the analysis.

Everything_to_ARIMAoutput <- function(path_fMRI_ses1, path_fMRI_ses2, 
                                      path_ET_ses1, path_ET_ses2, 
                                      path_HC_ses1, path_HC_ses2) {
  #################################################################################################################
  # Process fMRI data to get mean time series for each session
  xii.ses1 <- read_xifti(path_fMRI_ses1, brainstructures = "all")
  xii.ses2 <- read_xifti(path_fMRI_ses2, brainstructures = "all")
  xii_pmean.ses1 <- fMRI_xii_pmean(dtseries_data = xii.ses1)
  xii_pmean.ses2 <- fMRI_xii_pmean(dtseries_data = xii.ses2)
  
  #################################################################################################################
  # Process Eye Tracking data and compute convolution signals for blinks and fixations
  input.conv.blink.ses1 <- ETascDataProcess(file_path = path_ET_ses1, blink = TRUE)
  input.conv.fixation.ses1 <- ETascDataProcess(file_path = path_ET_ses1, fixation = TRUE)
  input.conv.blink.ses2 <- ETascDataProcess(file_path = path_ET_ses2, blink = TRUE)
  input.conv.fixation.ses2 <- ETascDataProcess(file_path = path_ET_ses2, fixation = TRUE)
  
  # Compute and extract convolution time series for each session
  convolution.timeseries.blink.ses1 <- Convolution_function(totaltime = input.conv.blink.ses1[[1]],
                                                            onsets = input.conv.blink.ses1[[2]], 
                                                            duration = input.conv.blink.ses1[[3]], 
                                                            sampling_rate = input.conv.blink.ses1[[4]])
  sub.convolution.timeseries.blink.ses1 <- Extraction_ETtime(conv.data = convolution.timeseries.blink.ses1, 
                                                             xii.mean = xii_pmean.ses1, interval = 1.127)
  
  convolution.timeseries.fixation.ses1 <- Convolution_function(totaltime = input.conv.fixation.ses1[[1]],
                                                               onsets = input.conv.fixation.ses1[[2]], 
                                                               duration = input.conv.fixation.ses1[[3]], 
                                                               sampling_rate = input.conv.fixation.ses1[[4]])
  sub.convolution.timeseries.fixation.ses1 <- Extraction_ETtime(conv.data = convolution.timeseries.fixation.ses1, 
                                                                xii.mean = xii_pmean.ses1, interval = 1.127)
  
  convolution.timeseries.blink.ses2 <- Convolution_function(totaltime = input.conv.blink.ses2[[1]],
                                                            onsets = input.conv.blink.ses2[[2]], 
                                                            duration = input.conv.blink.ses2[[3]], 
                                                            sampling_rate = input.conv.blink.ses2[[4]])
  sub.convolution.timeseries.blink.ses2 <- Extraction_ETtime(conv.data = convolution.timeseries.blink.ses2, 
                                                             xii.mean = xii_pmean.ses2, interval = 1.127)
  
  convolution.timeseries.fixation.ses2 <- Convolution_function(totaltime = input.conv.fixation.ses2[[1]],
                                                               onsets = input.conv.fixation.ses2[[2]], 
                                                               duration = input.conv.fixation.ses2[[3]], 
                                                               sampling_rate = input.conv.fixation.ses2[[4]])
  sub.convolution.timeseries.fixation.ses2 <- Extraction_ETtime(conv_data = convolution.timeseries.fixation.ses2, 
                                                                fmri_data = xii_pmean.ses2, tr = 1.127)
  
  #################################################################################################################
  # Construct and scale the final design matrix with confounders from both sessions
  head.confounder.ses1 <- HeadMotionConfounder_process(path = path_HC_ses1, xii.mean = xii_pmean.ses1)
  head.confounder.ses2 <- HeadMotionConfounder_process(path = path_HC_ses2, xii.mean = xii_pmean.ses2)
  blink.ses1.ses2 <- c(sub.convolution.timeseries.blink.ses1, sub.convolution.timeseries.blink.ses2)
  fixation.ses1.ses2 <- c(sub.convolution.timeseries.fixation.ses1, sub.convolution.timeseries.fixation.ses2)
  Scale.design.matrix <- DesignMatrix_process(ET.covariates = cbind(blink.ses1.ses2, fixation.ses1.ses2), 
                                              head.confounder.ses1, head.confounder.ses2)
  
  #################################################################################################################
  # Combine session fMRI time series and apply ARIMA model
  xii_pmean.ses1.ses2 <- cbind(xii_pmean.ses1, xii_pmean.ses2)
  ARIMAmodel.results <- ARIMAmodel(xii_pmean = xii_pmean.ses1.ses2, design.matrix = Scale.design.matrix)
  
  return(list(ARIMAmodel.results, Scale.design.matrix))
}

# Example of usage:
# final_results <- Everything_to_ARIMAoutput(path_fMRI_ses1, path_fMRI_ses2, path_ET_ses1, path_ET_ses2, path_HC_ses1, path_HC_ses2)

##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
############################################################################################################################ Covariance Regression !!!############################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################

#####################################################################################
# Function to execute the OLS model fitting and compare t-statistics across outcomes.
# This function integrates the process of transforming ARIMA model residuals to a logarithmic
# scale, fitting robust OLS models, and comparing the resulting coefficients statistically.

# Input:
# - ARIMAoutput: A list containing ARIMA model results and a scaled design matrix.

# Output:
# - A matrix of t-statistics comparing different regression outcomes.

Function_get_Z_stat_covariance <- function(ARIMAoutput) {
  # Unpack ARIMA model results and scaled design matrix from the input list
  ARIMAmodel.results <- ARIMAoutput[[1]]
  Scale.design.matrix <- ARIMAoutput[[2]]
  length <- nrow(Scale.design.matrix)  # Determine the number of observations
  
  ###################################################################################################################
  # Generate the log-transformed outcome matrix from ARIMA residuals for OLS regression
  # This step converts the residuals into a log Y matrix of shape 100x100x671, then reshapes it to 671x5050 for OLS modeling
  log.Y.LS <- LS.construct.logY(resi.acf.list = ARIMAmodel.results[[3]], length = length)
  
  #######################################################################################################################
  # Perform robust estimation using the log-transformed outcomes and the scaled design matrix
  robust.estimation <- LS.robust.estiomation(log.Y.array = log.Y.LS, totalcovariates.scale = Scale.design.matrix, length = length)
  
  #######################################################################################################################
  # Calculate the matrix of t-statistics for comparisons between different model outcomes
  z.statistic.matrix <- Z.statistic.compare(robust.regression.all = robust.estimation)
  
  return(z.statistic.matrix)
}

# Example of usage:
# Assuming 'ARIMAmodel_results_1887002' is your input list containing ARIMA results and the scaled design matrix:
# final_comparison_matrix <- Function_get_Z_stat_covariance(ARIMAmodel_results_1887002)
