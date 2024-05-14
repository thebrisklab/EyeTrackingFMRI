####################################################################################
####################################################################################
#### This R script is used to call the function to reproduce the project results ###
####################################################################################
####################################################################################

# source the R scripts with functions needed
source("Main_Functions_fMRI_ET.R")
source("0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R")
source("0702_SingleSessionARIMAWorkflow.R")

# Call Everything_to_ARIMAoutput funtion to get the arima model output
# 1917203 TD ET
ARIMAmodel.results_1917203 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-01/func/sub-1917203_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-01/func/sub-1917203_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-02/func/sub-1917203_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1917203/ses-02/func/sub-1917203_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19172-03/1917201.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19172-03/1917202.asc"
)

# 1879403 TD ET
ARIMAmodel.results_1879403 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-01/func/sub-1879403_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-01/func/sub-1879403_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-02/func/sub-1879403_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1879403/ses-02/func/sub-1879403_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18794_03/18794_03_01.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/18794_03/18794_03_02.asc"
)

# 1908702 ASD ET okay, but movie-1 has some jumpy around the edges of the screen 
ARIMAmodel.results_1908702 <- Everything_to_ARIMAoutput(
  path_fMRI_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-01/func/sub-1908702_ses-01_task-movie1_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-01/func/sub-1908702_ses-01_task-movie1_run-01_desc-confounds_timeseries.tsv",
  path_fMRI_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-02/func/sub-1908702_ses-02_task-movie2_run-01_space-fsLR_den-91k_bold.dtseries.nii",
  path_HC_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ImproveData/sub-1908702/ses-02/func/sub-1908702_ses-02_task-movie2_run-01_desc-confounds_timeseries.tsv",
  path_ET_ses1 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19087-02/19087_02_01.asc",
  path_ET_ses2 = "/Users/fredhuang/Library/CloudStorage/OneDrive-EmoryUniversity/ImproveFConn_FredXuchengHuang/ET Data/19087-02/19087_02_02.asc"
)

####################################################################################
####################################################################################
# Data processing outputs from ARIMA model. The goal is to construct a data.frame and do t-test / lm
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

####################################################################################
####################################################################################
# Remove variables that start with the pattern "dat.lm"
to_remove <- grep("^dat\\.lm", ls(), value = TRUE)
rm(list = to_remove)

# for each subject, call function to get the output 
dat.lm.1917203 <- data.frame.for.lm(ARIMAmodel.results_1917203) # TD, 
dat.lm.1879403 <- data.frame.for.lm(ARIMAmodel.results_1879403) # TD, 
dat.lm.1908702 <- data.frame.for.lm(ARIMAmodel.results_1908702) # ASD, et-1 a bit jumpy

####### rbind together
# Create a list of all data frame names following the pattern dat.lm.xxxx
data_frame_names <- ls(pattern = "^dat\\.lm\\.[0-9]+")
# Use mget to get the data frames in a list
data_frames_list <- mget(data_frame_names)
# Combine all data frames into one
combined_data_frame <- do.call(rbind, data_frames_list)

#### read ASD, TD file for merging together.
asd.td.info <- read.csv("/Users/fredhuang/Desktop/Research_Project/Dr. Risk/thesis/ASDStatisticalMethod-ASDAndNonASDParticip_DATA_2024-02-27_1144.csv", header = T)
# remove the "-" in the first column
asd.td.info$src_subject_id <- as.numeric(gsub("-", "", asd.td.info$src_subject_id))
# merge with combined_data_frame
combined_data_frame.asdtd <- left_join(combined_data_frame, asd.td.info, by = c("subjID" = "src_subject_id"))
# add one more column to do lm easily
combined_data_frame.asdtd$ASD <- ifelse(combined_data_frame.asdtd$phenotype == "Autism Spectrum Disorder", 1, 0)

####################################################################################
####################################################################################
####  for loop to do t-test for each region
ASD.TD.lm <- data.frame()
for (i in 1:100) {
  test <- lm( (fixcoef) ~ ASD, data = subset(combined_data_frame.asdtd, Region == i))
  pvalue <- summary(test)$coefficient[2,] ##!!!!!### 1 is beta0, i.e., TD only. 2 is beta1, the difference between TD and ASD 
  ASD.TD.lm <- rbind(ASD.TD.lm, pvalue)
}
                                
