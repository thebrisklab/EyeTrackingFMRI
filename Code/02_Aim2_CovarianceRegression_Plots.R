########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
# Code for reproduce the poster "Analysis of Simultaneous Eye-tracking and fMRI data Collected in Children with ASD".
# AIM 2 - Covariance regression model to examine the relationship between eye movement events and brain functional connectivity. (Chord Plots, Heatmap Plots)
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
  "neuRosim", "lmtest", "ciftiTools", "sandwich", "readxl", "RCurl", "git2r", "ggraph", "igraph", "tidyverse", "RColorBrewer"
)

# Run the function to install and load the packages
install_and_load_packages(required_packages)
# Setting Workbench path for ciftiTools
ciftiTools.setOption("wb_path", "/Users/fredhuang/Downloads/workbench")

# source the function
# Source the R scripts with functions needed

source("Code/Main_Functions_fMRI_ET.R")
source("Code/0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R")
source("Code/Heatmap_HierarchicalCirclePlot.R") 

# Load the ARIMA output, which is stored in a 'list' variable that includes: 1. ARIMA output; 2. Scale design matrix. 
# The ARIMA output itself is also a 'list' variable containing: 
# 1. a V by 4 matrix of eyeblink coefficients; 2. a V by 4 matrix of eyefixation coefficients; 3. T length residuals for V regions ; 4. V covariance matrices.

# read the data from ARIMA output
# can skip this section when these variables exist!
ARIMAmodel.results_1917203 <- readRDS("Data/Processed_data/ARIMAmodel.results_1917203.rds")
ARIMAmodel.results_1873503 <- readRDS("Data/Processed_data/ARIMAmodel.results_1873503.rds")
ARIMAmodel.results_1876402 <- readRDS("Data/Processed_data/ARIMAmodel.results_1876402.rds")
ARIMAmodel.results_1878002 <- readRDS("Data/Processed_data/ARIMAmodel.results_1878002.rds")
ARIMAmodel.results_1879403 <- readRDS("Data/Processed_data/ARIMAmodel.results_1879403.rds")
ARIMAmodel.results_1908702 <- readRDS("Data/Processed_data/ARIMAmodel.results_1908702.rds")

# call function to conduct 1. Log transformation; 2. LS regression; 3. Get the Z matrices
# NOTE: TIME CONSUMING!!!
zstat.mat.1917203 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1917203) # TD, 
zstat.mat.1873503 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1873503) # ASD, movie-1 only, et-1 medium 
zstat.mat.1876402 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1876402) # ASD, movie-1 only 246 timepoints
zstat.mat.1878002 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1878002) # ASD, movie-1 only, et-2 sleep
zstat.mat.1879403 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1879403) # TD, 
zstat.mat.1908702 <- Function_get_Z_stat_covariance(ARIMAmodel.results_1908702) # ASD, et-1 a bit jumpy 

########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
###### Heatmap plot #######
###########################

# call function to get the heatmap plot with customized color range
Heatmapplot <- Heatmap.plot(mat = zstat.mat.1917203, lower.bound = -2, upper.bound = +2)
print(Heatmapplot)

########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
###### Hierarchical Circle Plot #######
#######################################
# call function to get the circle plot with customized threshold value
# threshold.value is the quantile value, the function calculates the data value based on quantile

Circleplot <- Hierarchical.Circle.plot(z_statistic_matrix = zstat.mat.1917203, threshold.value = 0.999)
print(Circleplot)


########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
###### Population level, i.e., mean #####
# Subject to be revised.
