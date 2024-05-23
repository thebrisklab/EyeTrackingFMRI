# EyeTrackingFMRI - repository
Analysis of the simultaneous eye-tracking and movie-watching fMRI data. Includes 1) task activation modeling using the general linear model and 2) covariance regression; and the necessary data processing procedures for eye-tracking (ET) data. The repository contains the code for the replication of "Simultaneous Analysis of Eye-tracking and fMRI Data Collected in Children with ASD".

# R Code usage guidelines

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

2.3 `Extraction_ETtime`: Align the time points from eye-tracking convolution data with fMRI time points (T length), and standardize the extracted convolution time series to ensure the max value is equal to 1.
   - __Input:__ (1. convolution time series from **2.2**, 2. fMRI time series from **1.1**, 3. sampling interval of fMRI data, 1.127 in this project)
   - __Output:__ T length convolution time series of ET events.

#### 3. Confounder processing

3.1 `HeadMotionConfounder_process`: To process the nuisance covariates of head motion (In this case: 6 head motion control, 6 quadratics of 6 head motion control.)
   - __Input:__ (location of the head motion file; fMRI mean time series from **1.1**.)
   - __Output:__ T by 12 data.frame (6 "trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z", and 6 quadatic terms).

3.2 `DesignMatrix_process`: To construct the scaled design matrix. The default is two useful task sessions' (ses1, ses2) data.
   - __Input:__ (1. ET covariates from **2.3**, 2. Head motion covariates for ses1 from **3.1**, 3. Head motion covariates for ses2 from **3.1**.)
   - __Output:__ T by J data.frame, ready to be used in regression models (no intercept column!)

3.3 `DesignMatrix_process.ses1`: To construct the scaled design matrix for those who only have one useful ses data
   - __Input:__  (1. ET covariates from **2.3**, 2. Head motion covariates for useful ses from **3.1**)
   - __Output:__ T by J` data.frame, no interaction terms across different sessions like **3.2**.

### B.1 General linear regression:

4. `ARIMAmodel`: To execute autoregressive AR(3) regression for each of V brain regions to examine the relationship between ET data and brain activation.
   - __Input:__ (1. fMRI time series from **1.1**, 2. design matrix from **3**, 3.number of brain region V)
   - __Output:__ A list object includes:
      1. Estimated ET blink coefficients matrix (V by 4)
      2. Estimated ET fixation coefficients matrix (V by 4)
      3. Residuals for each of the V regressions (V by T)
      4. A list object storing covariance matrics for each of the V regressions (V elements, each being a (J+4) by (J+4) matrix, where the number 4 corresponds to the inclusion of ar1, ar2, ar3, intercept)

### B.2 Covariance regression:

5. `LS.construct.logY`: To process the residuals from ARIMA models and construct a logarithmic transformation of the outcome variable for Least Squares (LS) regression analysis.
   - __Input:__ (1.residuals from **4.3**, 2. the length of time series T)
   - __Output:__ T by $\frac{V(V+1)}{2}$ matrix, each column being dependent variables for LS regression below.

6.1 `LS.robust.estimation`: To perform LS regression using robust estimation (HC3) on T by $\frac{V(V+1)}{2}$ = 5050 regions.
   - __Input:__ (1. T by 5050 matrix from **5**, 2. scaled design matrix from **3**, 3. time series length T)
   - __Output:__ A list object with 5050 elements, each storing:
      1. Robust estimation coefficients matrix (J+1 by 4)
      2. Robust covariance matrix (J+1 by J+1)

6.2 `Z.statistic.compare`: To compare the z-statistics between the estimated coefficients for Eye-blink and Eye-fixation from **6.1**.
   - __Input:__ The list object from **6.1**.
   - __Output:__ V by V matrix of z-statistics.

6.3. `Residual.Diag_OLS`: To randomly pick the LS models to do the model diagnostics on residuals.
   - __Input:__ (1. Log Y from **5**, 2. scaled design matrix from **3**, 3. time series length T)
   - __Output:__ A list object includes ACF, PACF, Q-Q norm, and Histogram plots for model diagnostics. 
  
#### In practice
A level 1 function `0701_Subject_FullARIMA_to_CovarianceAnalysisPipeline.R` was constructed to call the above functions in sequence in real analysis. However, I recommend making such separate subfunctions because 1. it is easier to debug, and 2. some subjects have different available data points to analyze, and such subfunctions can be used to change the parameters for a single participant.

### C. General linear regression: Population effects models

#### LM models between ASD & non-ASD:
 1. Model 1: $blinkcoef_{iv1} = \alpha_{v1} + \gamma_{v1} * ASD + \nu_{v1}$
 2. Model 2: $fixcoef_{iv2} = \alpha_{v2} + \gamma_{v2} * ASD + \nu_{v2}$
 3. Model 3: $(fixcoef_{iv2} - blinkcoef_{iv1}) = (\alpha_{v2} - \alpha_{v1}) + (\gamma_{v2} - \gamma_{v1}) + (\nu_{v2} - \nu_{v1})$

`ThreeLMs_BrainMapping.R` includes the following three functions for conducting three LM models to quantify the population effects:

`LM.pop.model`: To conduct the three LM models; 
`Brainmap.coefs`: Returns the mapping data for coefficients;
`Brainmap.pval`: Returns the mapping data for log10 FDR P-values.

### D. Covariance regression visualization: Heatmap & chord diagram

`Heatmap_HierarchicalCirclePlot.R` includes the following two functions for creating the heatmap and circle plots for covariance regression:

`Heatmap.plot`: Returns graph object storing heatmap figure;
   - __Input:__ 1. V by V matrix, 2. lower bound of color range, 3. upper bound of color range
     
`Hierarchical.Circle.plot`: Returns graph object storing circle plots;
   - __Input:__ 1. V by V matrix, 2. quantile value, 3. subject = TRUE
   - **Note:** When subject == TRUE, the quantile value is used based on _qnorm()_, when subject != TRUE, quantile value is based on the data quantile.

# The required R packages

## Data processing and manipulation
  - `dplyr`,        # For data manipulation and transformation within data frames.
  - `tidyr`,        # For tidying data, making it suitable for analysis.
  - `tidyverse`,    # For an assortment of data science packages.
  - `stringr`,      # For manipulation of string objects.
  - `readxl`,       # For reading Excel files into R.

## Statistical analysis and modeling
  - `neuRosim`,     # For simulation of fMRI data in neuroimaging.
  - `lmtest`,       # For diagnostic tests in linear regression models.
  - `sandwich`,     # For robust estimation of covariance matrices.
  - `intervals`,    # For working with interval data such as confidence intervals.

## Visualization
  - `ggplot2`,      # For creating sophisticated data visualizations.
  - `ggraph`,       # For creating network graphs using the grammar of graphics (Chord Diagram).
  - `igraph`,       # For network analysis and graph theory operations (Chord Diagram).
  - `RColorBrewer`  # For color palettes for visualizations.

## Neuroimaging
  - `ciftiTools`,   # For working with CIFTI files used in neuroimaging data.

## Networking and version control
  - `httr`,         # For tools to work with URLs and HTTP.
  - `git2r`,        # For working with Git repositories.
  - `RCurl`,        # For network (HTTP/FTP/...) client interface.

## Specialized packages
  - `eyelinker`,    # For processing eye-tracking data from Eyelink 1000 Plus devices.
  - `eyeQuality`,   # For assessing the quality of eye-tracking data.
  
Additional commands:

- `ciftiTools.setOption("wb_path", "your local path")` - Sets the Workbench command-line applications path in `ciftiTools`.


# Simple Example
