# EyeTrackingFMRI - repository
Analysis of the simultaneous eye-tracking and movie-watching fMRI data. Includes 1) task activation modeling using the general linear model and 2) covariance regression; and the necessary data processing procedures for eye-tracking (ET) data. The repository contains the code for the replication of "Simultaneous Analysis of Eye-tracking and fMRI Data Collected in Children with ASD".

# The workflow for the overall process in this repository

### A. Data Processing:

#### 1. fMRI Data Processing
- **fMRI Process Function**: Extracts the fMRI time series for each brain region by calculating mean values across voxels by time points.
  - Function Signature: {Function input}(cifiti object, parcellation) $\rightarrow$ {Function output: }Mean fMRI Time Series: V by T matrix

- **fMRI Mapping Plot Function**: Maps brain region-specific data such as mean, p-value, coefficients, or t-statistics to a brain graph.
  - Returns: cifiti object for plotting

#### 2. Eye Tracking (ET) Data Processing
- **ET Data Process Function**: Processes the eye tracking data to obtain parameters needed for convolution.
  - Returns:
    1. Total time of the ET task
    2. Onset of the ET events
    3. Duration of the ET events
    4. ET data sampling rate

- **Convolution Function**: Performs convolution between eyeblink & eyefixation events and the Double-gamma Hemodynamic Response Function (HRF) using FFT.
  - Function Signature: `$\text{convolution}(eye-blink, eye-fixation, HRF)$`
  - Returns: list object with 1. ET convolution time series vector; 2. real-time vector

2.3 **time series extraction function**: To align the time points from eye tracking data with fMRI time points.

   - **Equation**: $\text{timeseries subtraction}(eye\_tracking\_data, fMRI\_time\_points)$

3.1 **confounder process function**: To process the confounder covariates of head motion 

3.2 **design matrix function**: To construct the design matrix [Note, need to construct the eye blink, fixation, and other covariates (Need to be scaled).]

   - **Equation**: $\text{confounder process}(covariates) \rightarrow \text{Design Matrix}$

#### B. Regression:

4. **PolynomialARIMA function**: To execute Polynomial ARIMA regression, returning a list of eyeblink & eyefixation coefficients and residuals.

   - **Equation**: $\text{PolynomialARIMA}(data) \rightarrow \{\text{coefficients}, \text{residuals}\}$

5. **OLSoutcome process function**: To construct the logarithm of \(Y\) as the outcome variable for regression analysis.

   - **Equation**: $\log(Y_{vivj}) \text{ as outcome}$

6.1 **OLSfitting function**: To perform OLS regression for all "5050" regions with robust estimation techniques. [Robust Estimation]
   - **Equation**: $\text{OLSfitting}(data, \text{regions}=5050) \rightarrow \text{OLS Outcome}$

6.2 **Comparison statistic function**: To compare the t-statistic between Eyeblink and Eyefixation estimated coefficients
    -**Input**: robust regression coefficients for all 5050 OLS models
    **Output**: lower-left triangle matrix storing 5050 comparisons 

7. **Model Diagnostics for OLS function**: To randomly pick the OLS model to do the model diagnostics on residuals
   - **Input**: is the same as function 6.1, i.e., use log. Y and design matrix to calculate the residuals and make plots

8. **Got the z-statistics for Eye-fixation - Eye-blink.**
