########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
# Code for reproduce the poster "Analysis of Simultaneous Eye-tracking and fMRI data Collected in Children with ASD".
# AIM 1 - GLM model to examine the relationship between eye movement events and brain activation.
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
ciftiTools.setOption("wb_path", "/Users/fredhuang/Downloads/workbench")

# Clone the repository
# repo <- git2r::clone(url = "https://github.com/thebrisklab/EyeTrackingFMRI.git", 
#                      local_path = "/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis/Thesis_fMRI/For Github/Clone",
#                      credentials = git2r::cred_token()) # maybe due to private

# Source the R scripts with functions needed
setwd("/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis/Thesis_fMRI/For Github")
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

########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
# To do the population level analysis - Aim 1: Model 1/2/3

# Function to convert each subject's ARIMA output as the data.frame
########## Function input is ARIMAmodel.results_1873503, output is dataframe for this subject ###########
data.frame.for.lm <- function(ARIMAmodel.results.input = ARIMAmodel.results_1873503) {
  
  ARIMAmodel.results <- ARIMAmodel.results.input[[1]]
  variable.names <- deparse(substitute(ARIMAmodel.results.input)) # for extracting the subject id later
  
  # below is covariance for the subject
  ARIMA.cov.blink.fix <- sapply(1:100, function(i) ARIMAmodel.results[[4]][[i]][5, 6]) 
  # eyeblink & eyefixation merge and convert to data.frame, change column names
  dat.lm <- as.data.frame(cbind(ARIMAmodel.results[[1]], ARIMAmodel.results[[2]], ARIMA.cov.blink.fix))
  colnames(dat.lm) <- c("blinkcoef","blinkstd","blinkt","blinkpval","fixcoef","fixstd","fixt","fixpval", "cov")
  
  # extract the subject ID 
  # Use gsub to remove everything except the last numbers
  last_numbers <- gsub("^.*_([0-9]+)$", "\\1", variable.names)
  # Convert the result to numeric if needed
  last_numbers_numeric <- as.numeric(last_numbers)
  
  # add factor columns
  dat.lm = dat.lm %>% dplyr::mutate(Region = c(1:100), subjID = rep(last_numbers_numeric, 100))
  return(dat.lm)
}

# for each subject, call function to get the output 
# Remove variables that start with the pattern "dat.lm", just for clean up
to_remove <- grep("^dat\\.lm", ls(), value = TRUE)
rm(list = to_remove)

dat.lm.1917203 <- data.frame.for.lm(ARIMAmodel.results_1917203) # TD, 
dat.lm.1879403 <- data.frame.for.lm(ARIMAmodel.results_1879403) # TD, 
dat.lm.1908702 <- data.frame.for.lm(ARIMAmodel.results_1908702) # ASD, et-1 a bit jumpy
dat.lm.1878002 <- data.frame.for.lm(ARIMAmodel.results_1878002) # ASD, movie-1 only, et-2 sleep
dat.lm.1876402 <- data.frame.for.lm(ARIMAmodel.results_1876402) # ASD, movie-1 only 246 timepoints
dat.lm.1873503 <- data.frame.for.lm(ARIMAmodel.results_1873503) # ASD, movie-1 only, et-1 medium

####### rbind together
# Create a list of all data frame names following the pattern dat.lm.xxxx
data_frame_names <- ls(pattern = "^dat\\.lm\\.[0-9]+")
# Use mget to get the data frames in a list
data_frames_list <- mget(data_frame_names)
# Combine all data frames into one
combined_data_frame <- do.call(rbind, data_frames_list)

#### read ASD, TD file for merging together to add the demographic data
asd.td.info <- read.csv("/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/ASDStatisticalMethod-ASDAndNonASDParticip_DATA_2024-02-27_1144.csv", header = T)
# remove the "-" in the first column
asd.td.info$src_subject_id <- as.numeric(gsub("-", "", asd.td.info$src_subject_id))
# merge with combined_data_frame
combined_data_frame.asdtd <- left_join(combined_data_frame, asd.td.info, by = c("subjID" = "src_subject_id"))
# add one more column to do lm easily
combined_data_frame.asdtd$ASD <- ifelse(combined_data_frame.asdtd$phenotype == "Autism Spectrum Disorder", 1, 0)


########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
# LM model / t-test between ASD & TD
# Model 1: blinkcoef_iv1 ~ \alpha_v1 + \gamma_v1 * ASD + \nu_v1
# Model 2: fixcoef_iv2 ~ \alpha_v2 + \gamma_v2 * ASD + \nu_v2
# Model 3: (fixcoef_iv2 - blinkcoef_iv1) ~ (\alpha_v2 - \alpha_v1) + (\gamma_v2 - \gamma_v1) + (\nu_v2 - \nu_v1)

# source the function
setwd("/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/Data Analysis/Thesis_fMRI/For Github")
source("ThreeLMs_BrainMapping.R")

# read the xii for brain mapping template
xii <- read_xifti("/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-01/func/sub-1879403_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii", 
                  brainstructures = "all")

#######################
####### Model 1 #######
#######################

# call function to get the LM results
# NOTE, the output is a list with 1 indicates population effects (non ASD); 2 indicates ASD modification on effects
LMresults.blink <- LM.pop.model(combined_data_frame.asdtd, outcome = "blinkcoef")
lm1.nonASD <- LMresults.blink[[1]]
lm1.ASD <- LMresults.blink[[2]]

########################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1: non.ASD 
brain.plot.nonASD.coef.lm1 <- Brainmap.coefs(lm1.nonASD[,1], xii = xii)
title.name.nonASD.lm1 = "nonASD: Effects of Eyeblink events on fMRI Signals"
plot(brain.plot.nonASD.coef.lm1, zlim = c(-2,2), title = title.name.nonASD.lm1, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.nonASD.fdr.lm1 <- Brainmap.pval(lm1.nonASD[,4], fdr.threshold = 0.05, xii = xii)
title.name.nonASD.fdr.lm1 = "nonASD: Effects of Eyeblink events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.05 Threshold"
plot(brain.plot.nonASD.fdr.lm1, zlim = c(0,3), title = title.name.nonASD.fdr.lm1, borders = "black")

#######################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1 - ASD vs. non.ASD
brain.plot.ASD.coef.lm1 <- Brainmap.coefs(lm1.ASD[,1], xii = xii)
title.name.ASD.lm1 = "nonASD: Effects of Eyeblink events on fMRI Signals"
plot(brain.plot.ASD.coef.lm1, zlim = c(-2,2), title = title.name.ASD.lm1, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.ASD.fdr.lm1 <- Brainmap.pval(lm1.ASD[,4], fdr.threshold = 0.2, xii = xii)
title.name.ASD.fdr,lm1 = "nonASD: Effects of Eyeblink events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.2 Threshold"
plot(brain.plot.ASD.fdr.lm1, zlim = c(0,3), title = title.name.ASD.fdr, borders = "black")

#######################
####### Model 2 #######
#######################

# call function to get the LM results
# NOTE, the output is a list with 1 indicates population effects (non ASD); 2 indicates ASD modification on effects
LMresults.fixation <- LM.pop.model(combined_data_frame.asdtd, outcome = "fixcoef")
lm2.nonASD <- LMresults.fixation[[1]]
lm2.ASD <- LMresults.fixation[[2]]

########################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1: non.ASD 
brain.plot.nonASD.coef.lm2 <- Brainmap.coefs(lm2.nonASD[,1], xii = xii)
title.name.nonASD.lm2 = "nonASD: Effects of Eyefixation events on fMRI Signals"
plot(brain.plot.nonASD.coef.lm2, zlim = c(-2,2), title = title.name.nonASD.lm2, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.nonASD.fdr.lm2 <- Brainmap.pval(lm2.nonASD[,4], fdr.threshold = 0.05, xii = xii)
title.name.nonASD.fdr.lm2 = "nonASD: Effects of Eyefixation events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.05 Threshold"
plot(brain.plot.nonASD.fdr.lm2, zlim = c(0,3), title = title.name.nonASD.fdr.lm2, borders = "black")

#######################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1 - ASD vs. non.ASD
brain.plot.ASD.coef.lm2 <- Brainmap.coefs(lm2.ASD[,1], xii = xii)
title.name.ASD.lm2 = "nonASD: Effects of Eyefixation events on fMRI Signals"
plot(brain.plot.ASD.coef.lm2, zlim = c(-2,2), title = title.name.ASD.lm2, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.ASD.fdr.lm2 <- Brainmap.pval(lm2.ASD[,4], fdr.threshold = 0.2, xii = xii)
title.name.ASD.fdr.lm2 = "nonASD: Effects of Eyefixation events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.2 Threshold"
plot(brain.plot.ASD.fdr.lm2, zlim = c(0,3), title = title.name.ASD.fdr.lm2, borders = "black")

#######################
####### Model 3 #######
#######################

# call function to get the LM results
# NOTE, the output is a list with 1 indicates population effects (non ASD); 2 indicates ASD modification on effects
LMresults.coef_blink <- LM.pop.model(combined_data_frame.asdtd, model3 = TRUE)
lm3.nonASD <- LMresults.coef_blink[[1]]
lm3.ASD <- LMresults.coef_blink[[2]]

########################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1: non.ASD 
brain.plot.nonASD.coef.lm3 <- Brainmap.coefs(lm3.nonASD[,1], xii = xii)
title.name.nonASD.lm3 = "nonASD: Effects of Eyeblink events on fMRI Signals"
plot(brain.plot.nonASD.coef.lm3, zlim = c(-2,2), title = title.name.nonASD.lm3, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.nonASD.fdr.lm3 <- Brainmap.pval(lm3.nonASD[,4], fdr.threshold = 0.05, xii = xii)
title.name.nonASD.fdr.lm3 = "nonASD: Effects of Eyeblink events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.05 Threshold"
plot(brain.plot.nonASD.fdr.lm3, zlim = c(0,3), title = title.name.nonASD.fdr.lm3, borders = "black")

#######################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots, i.e, Results Aim1 - ASD vs. non.ASD
brain.plot.ASD.coef.lm3 <- Brainmap.coefs(lm3.ASD[,1], xii = xii)
title.name.ASD.lm3 = "nonASD: Effects of Eyeblink events on fMRI Signals"
plot(brain.plot.ASD.coef.lm3, zlim = c(-2,2), title = title.name.ASD.lm3, borders = "black")

#### FDR adjusted P-value and set those region greater than 0.2 as NAs
brain.plot.ASD.fdr.lm3 <- Brainmap.pval(lm3.ASD[,4], fdr.threshold = 0.2, xii = xii)
title.name.ASD.fdr.lm3 = "nonASD: Effects of Eyeblink events  on fMRI Signals: FDR-Adjusted -log10(P-Values) with 0.2 Threshold"
plot(brain.plot.ASD.fdr.lm3, zlim = c(0,3), title = title.name.ASD.fdr.lm3, borders = "black")
