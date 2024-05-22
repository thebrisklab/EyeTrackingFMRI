# EyeTrackingFMRI - repository
Analysis of the simultaneous eye-tracking and movie-watching fMRI data. Includes 1) task activation modeling using the general linear model and 2) covariance regression; and the necessary data processing procedures for eye-tracking (ET) data. The repository contains the code for the replication of "Simultaneous Analysis of Eye-tracking and fMRI Data Collected in Children with ASD".

# R Code Usage Guidelines

## Main functions
- `Main_Functions_fMRI_ET.R` contains the level 2 main functions used in this project. Individual sub-functions files can be found in the folder `Small_Functions`.
- `0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R` contains the level 1 functions that complete the workflow from 1. reading fMRI and eye-tracking data to the generalized regression output (Aim 1), and 2. covariance regression output (Aim 2).
- `0702_SingleSessionARIMAWorkflow.R` contains functions to analyze subjects with only one useful task session data as a complement to the **0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R**.
- `ThreeLMs_BrainMapping.R` contains functions to do population-level analyses in Aim 1.
- `Heatmap_HierarchicalCirclePlot.R` contains functions to create the circle and heatmap plots in Aim 2. 

## Working functions
- `00_Aim1_GLMmodel.R` calls the functions to do subject-level ARIMA regression (Aim 1)
- `01_Aim1_GLMmodel_to_Population_Results.R` to do the population-level analyses from ARIMA regression and generate the brain mapping plots (Aim 1).
- `02_Aim2_CovarianceRegression_Plots.R` to do the covariance regression and generate the circle plots (Aim 2).

**When reproducing the project's results, open and run the three working functions in order.**

# Overall analytical workflow (sub-functional description)

### A. Data Processing:

#### 1. fMRI data processing function
1.1 `fMRI_xii_pmean`: Extract the fMRI time series for each brain region by calculating mean values across voxels by time points.
  - __Input:__ (cifti object, parcellation, the number of brain region V) $\rightarrow$ __Return:__ fMRI Time Series: V by T matrix

1.2 `fMRIBrain_Mapping`: Map data to brain regions to visualize data such as mean, p-values, or t-statistics on a brain map.
  - __Input:__ (cifti object, mapping data (V length vector), parcellation, the number of brain region V) $\rightarrow$ __Return:__ brain graph object with mapped data

#### 2. Eye Tracking (ET) Data Processing
2.1 `ETascDataProcess`: Process the eye-tracking data to obtain input parameters needed for convolution. Used `detectBlinks` function from `eyeQuality` R package to detect the eye blink events.
  - __Input:__ (path for eye-tracking data, blink or fixation)
  - __Return:__ list object includes:
    1. Total time of the ET task
    2. Onsets of the ET events
    3. Durations of the ET events
    4. ET data sampling rate

2.2 `Convolution_function`: Performs convolution between eyeblink & eyefixation events and the Double-gamma Hemodynamic Response Function (HRF) using FFT.
  - __Input:__ (1. total time for the task, 2. onsets of the events, 3. durations of the events, 4. sampling interval)
  - __Return:__ 1. ET convolution time series vector; 2. real-time vector

2.3 `Extraction_ETtime`: Align the time points from eye-tracking convolution data with fMRI time points (T length), standardize the extracted convolution time series to ensure the max value is equal to 1.
   - __Input:__ (1. convolution time series from 2.2, 2. fMRI time series from 1.1, 3. sampling interval of fMRI data, 1,127 in this project)
   - __Output:__ T length convolution time series of ET events.

#### 3. Confounder Processing

3.1 **confounder process function**: To process the nuisance covariates of head motion (In this case: 6 head motion control, 6 quadratics of 6 head motion control.)

   - Function Signature: {Function input:} (path of the file, fMRI mean time series for 1.1)

3.2 **design matrix function**: To construct the design matrix [Note, need to construct the eye blink, fixation, and nuisance covariates (Need to be scaled).]

   - Function Signature: {Function input:} (ET convolution from 2.3, nuisance covariates form 3.1) $\rightarrow$ {Returns:} T by J Design Matrix

### B.1 General Linear Regression:

4. **Brain Activation**: To execute autoregressive ARIMA regression to examine the relationship between ET data and brain region activation.

   -  Function Signature: {Function input:} (fMRI mean time series form 1.1, design matrix from 3.2 ) $\rightarrow$ {Returns:} (list of 1. ET-blink coefs, 2. ET-fixation coefs, 3. Residuals 4. Covariance matrix)

### B.2. Covariance Regression:

5. **OLSoutcome process function**: To construct the logarithm of \(Y\) as the outcome variable for regression analysis.

   - Function Signature: {Function input:} (residuals from 4) $\rightarrow$ {Returns:} (logarithm of T by $\frac{V(V+1)}{2}$ matrix)

6.1 **OLSfitting function**: To perform LS regression for all $\frac{V(V+1)}{2}$ regions with robust estimation techniques. [HC 3]
   - Function Signature: {Function input:} (Log.Y from 5, Design matrix from 3.2) $\rightarrow$ {Returns:} (list of 1. Coefficients, 2. Robust covariance matrix)

6.2 **Comparison statistic function**: To compare the z-statistic between the slope of Eye-blink and Eye-fixation estimated coefficients from 6.1
   - Function Signature: {Function input:} (list from 6.1) $\rightarrow$ {Returns:} (V by V Z-statistic matrix)

7. **Model Diagnostics for OLS function**: To randomly pick the OLS model to do the model diagnostics on residuals
   - Function Signature: {Function input:} (1. Log.Y from 5, 2. Design matrix form 3.2) $\rightarrow$ {Returns:} (list of 1. ACF plot, 2.PACF plot, 3. Q-Q plot, 4.Histogram)
  
#### Note: A "big" function was also constructed to call the above functions all at once in real analysis. However, I recommend making such separate subfunctions because 1. it is easier to debug, and 2. some subjects have different available data points to analyze, and such subfunctions can be used to change the parameters for a single participant.

### C. Population effects models

### D. Visualisation: Heatmap & Chord Diagram

# The required R packages
- `eyelinker` - For processing eye-tracking data from Eyelink 1000 Plus devices.
- `stringr` - For manipulation of string objects.
- `sandwich` - For robust estimation of covariance matrices.
- `neuRosim` - For simulation of fMRI data in neuroimaging.
- `dplyr` - For data manipulation and transformation within data frames.
- `tidyr` - For tidying data, making it suitable for analysis.
- `ggplot2` - For creating sophisticated data visualizations.
- `intervals` - For working with interval data such as confidence intervals.
- `lmtest` - For diagnostic tests in linear regression models.
- `ciftiTools` - For working with CIFTI files used in neuroimaging data.
- `readxl` - For reading Excel files into R.
- `RCurl` - For network (HTTP/FTP/...) client interface.
- `ggraph` - For creating network graphs using the grammar of graphics (Chord Diagram).
- `igraph` - For network analysis and graph theory operations (Chord Diagram).
  
Additional commands:

- `ciftiTools.setOption("wb_path", "your local path")` - Sets the Workbench command-line applications path in `ciftiTools`.


# Simple Example
