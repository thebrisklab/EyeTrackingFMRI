# Call function to get the ARIMA output for each participant
# Required R packages
# Function to check and install required packages
install_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      message(sprintf("Package '%s' is already installed.", pkg))
    }
  }
}

# Define the packages needed
required_packages <- c(
  "eyelinker", "dplyr", "tidyr", "ggplot2", "intervals", "stringr", "httr", "eyeQuality", "tidyverse",
  "neuRosim", "lmtest", "ciftiTools", "sandwich", "readxl", "RCurl", "git2r"
)

# Run the function to install and load the packages
install_and_load_packages(required_packages)
# Setting Workbench path for ciftiTools
ciftiTools.setOption("wb_path", "/Users/fredhuang/Downloads/workbench")

# Source the R scripts with functions needed
source("Code/Main_Functions_fMRI_ET.R")
source("Code/0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R")
source("Code/0702_SingleSessionARIMAWorkflow.R")

############################################################################################################################################################
# Single subject level model - Aim 1: ARIMA model
# Call Everything_to_ARIMAoutput function to get the arima model output
# Note that read_cifti does not support direct reading from URLs
############################################################################################################################################################


############################################################################################################################################################
################### Data download: ET data; fMRI data; head motion data ###################
options(timeout = 6000) # large data file needs more time 

# Define the Dropbox URL
url <- "https://www.dropbox.com/scl/fo/yvxr7s8ivy9vei4plozgr/ADvFl4wxAiFx367d7vZzUWo?rlkey=50c1mrlvjftkdl9ifcp5c7thw&st=t9s5mx4h&dl=1"
# Define the destination file path
zipfile <- "Data/Download_data/EyeTrackingFMRI_Data.zip"
# Download the zip file from Dropbox
download.file(url, zipfile, mode = "wb")

# Define the destination directory to unzip the files
unzip_dir <- "Data/Download_data/EyeTrackingFMRI_Data"
# Unzip the downloaded file
unzip(zipfile, exdir = unzip_dir)
############################################################################################################################################################

# Call function to do single subject ARIMA analyses

# 1917203 TD ET
ARIMAmodel.results_1917203 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/sub-1917203_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/sub-1917203_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/sub-1917203_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/sub-1917203_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/1917201.asc",
  path_ET_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1917203/1917202.asc"
)

saveRDS(ARIMAmodel.results_1917203, "Data/Processed_data/ARIMAmodel.results_1917203.rds")

############################################################################################################################################################

# 1879403 TD ET
ARIMAmodel.results_1879403 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/sub-1879403_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/sub-1879403_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/sub-1879403_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/sub-1879403_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/18794_03_01.asc",
  path_ET_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1879403/18794_03_02.asc"
)

saveRDS(ARIMAmodel.results_1879403, "Data/Processed_data/ARIMAmodel.results_1879403.rds")

############################################################################################################################################################

# 1908702 ASD ET okay, but movie-1 has some jumpy around the edges of the screen 
ARIMAmodel.results_1908702 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/sub-1908702_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/sub-1908702_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/sub-1908702_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/sub-1908702_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/19087_02_01.asc",
  path_ET_ses2 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1908702/19087_02_02.asc"
)

saveRDS(ARIMAmodel.results_1908702, "Data/Processed_data/ARIMAmodel.results_1908702.rds")

############################################################################################################################################################

# 1878002 ASD, usability sheet sys et-2 sleep/crash
ARIMAmodel.results_1878002 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1878002/sub-1878002_ses-02_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1878002/sub-1878002_ses-02_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1878002/18780_02_01.asc"
)

saveRDS(ARIMAmodel.results_1878002, "Data/Processed_data/ARIMAmodel.results_1878002.rds")

############################################################################################################################################################

# Sub-1876402 ASD. no usable ET-2 data. For movie-1, can only analyze the current 246 timepoints.
ARIMAmodel.results_1876402 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1876402/sub-1876402_ses-02_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 =  "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1876402/sub-1876402_ses-02_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 =  "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1876402/1876401.asc"
)

saveRDS(ARIMAmodel.results_1876402, "Data/Processed_data/ARIMAmodel.results_1876402.rds")

############################################################################################################################################################

# 1873503 ASD
ARIMAmodel.results_1873503 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1873503/sub-1873503_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1873503/sub-1873503_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "Data/Download_data/EyeTrackingFMRI_Data/fMRI_data/1873503/1873501.asc"
)

saveRDS(ARIMAmodel.results_1873503, "Data/Processed_data/ARIMAmodel.results_1873503.rds")

############################################################################################################################################################

# 1908004 TD

