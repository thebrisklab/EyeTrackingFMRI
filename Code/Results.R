####################################################################################
####################################################################################
#### This R script is used to call the function to reproduce the project results ###
####################################################################################
####################################################################################

# source the R scripts with functions needed
source("Main_Functions_fMRI_ET.R")
source("0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R")
source()

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
