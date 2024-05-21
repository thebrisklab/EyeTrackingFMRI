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
  "eyelinker", "dplyr", "tidyr", "ggplot2", "intervals", "stringr",
  "neuRosim", "lmtest", "ciftiTools", "sandwich", "readxl", "RCurl", "git2r"
)

# Run the function to install and load the packages
install_and_load_packages(required_packages)
# Setting Workbench path for ciftiTools
ciftiTools.setOption("wb_path", "/Users/fredhuang/Downloads/workbench") # your real path

# Clone the repository
# repo <- git2r::clone(url = "https://github.com/thebrisklab/EyeTrackingFMRI.git", 
#                      local_path = "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis/Thesis_fMRI/For Github/Clone",
#                      credentials = git2r::cred_token()) # maybe due to private

# Source the R scripts with functions needed
setwd("Code")
source("Main_Functions_fMRI_ET.R")
source("0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R")
source("0702_SingleSessionARIMAWorkflow.R")

############################################################################################################################################################
############################################################################################################################################################
# Single subject level model - Aim 1: ARIMA model
# Call Everything_to_ARIMAoutput function to get the arima model output
# 1917203 TD ET
ARIMAmodel.results_1917203 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-01/func/sub-1917203_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-01/func/sub-1917203_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-02/func/sub-1917203_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-02/func/sub-1917203_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19172-03/1917201.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19172-03/1917202.asc"
)

saveRDS(ARIMAmodel.results_1917203, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1917203.rds")

# 1879403 TD ET
ARIMAmodel.results_1879403 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-01/func/sub-1879403_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-01/func/sub-1879403_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-02/func/sub-1879403_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-02/func/sub-1879403_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18794_03/18794_03_01.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18794_03/18794_03_02.asc"
)

saveRDS(ARIMAmodel.results_1879403, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1879403.rds")

# 1908702 ASD ET okay, but movie-1 has some jumpy around the edges of the screen 
ARIMAmodel.results_1908702 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-01/func/sub-1908702_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-01/func/sub-1908702_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-02/func/sub-1908702_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-02/func/sub-1908702_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19087-02/19087_02_01.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19087-02/19087_02_02.asc"
)

saveRDS(ARIMAmodel.results_1908702, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1908702.rds")

# 1878002 ASD, usability sheet sys et-2 sleep/crash
ARIMAmodel.results_1878002 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1878002_Ses4/ses-02/func/sub-1878002_ses-02_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1878002_Ses4/ses-02/func/sub-1878002_ses-02_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18780_02/18780_02_01.asc"
)

saveRDS(ARIMAmodel.results_1878002, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1878002.rds")

# Sub-1876402 ASD. no usable ET-2 data. For movie-1, can only analyze the current 246 timepoints.
ARIMAmodel.results_1876402 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1876402_Ses4/ses-02/func/sub-1876402_ses-02_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 =  "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1876402_Ses4/ses-02/func/sub-1876402_ses-02_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 =  "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18764-02/1876401.asc"
)

saveRDS(ARIMAmodel.results_1876402, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1876402.rds")

# 1873503 ASD
ARIMAmodel.results_1873503 <- Everything_to_ARIMAoutput.ONEses(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1873503_OnlySes1/ses-01/func/sub-1873503_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1873503_OnlySes1/ses-01/func/sub-1873503_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18735-03/1873501.asc"
)

saveRDS(ARIMAmodel.results_1873503, "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis Output/ARIMA_results/ARIMAmodel.results_1873503.rds")
