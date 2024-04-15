# EyeTrackingFMRI - repository
Analysis of the simultaneous eye-tracking and movie-watching fMRI data. Includes 1) task activation modeling using the general linear model and 2) covariance regression; and the necessary data processing procedures for eye-tracking (ET) data. The repository contains the code for the replication of "Simultaneous Analysis of Eye-tracking and fMRI Data Collected in Children with ASD".

# The workflow for the overall process in this repository

### A. Data Processing:

#### 1. fMRI Data Processing
1.1 **fMRI Process Function**: Extracts the fMRI time series for each brain region by calculating mean values across voxels by time points.
  - Function Signature: {Function input:} (cifti object, parcellation) $\rightarrow$ {Returns:} Mean fMRI Time Series: V by T matrix

1.2 **fMRI Mapping Plot Function**: Maps brain region-specific data such as mean, p-value, coefficients, or t-statistics to a brain graph.
  - Returns: cifti object for plotting

#### 2. Eye Tracking (ET) Data Processing
2.1 **ET Data Process Function**: Processes the eye-tracking data to obtain input parameters needed for convolution.
  - Returns:
    1. Total time of the ET task
    2. Onset of the ET events
    3. Duration of the ET events
    4. ET data sampling rate

2.2 **Convolution Function**: Performs convolution between eyeblink & eyefixation events and the Double-gamma Hemodynamic Response Function (HRF) using FFT.
  - Function Signature: {Function input} (Returns from the 2.1 ET Data Process Function )
  - Returns: 1. ET convolution time series vector; 2. real-time vector

2.3 **Time series extraction function**: To align the time points from eye-tracking convolution data with fMRI time points.

   - Function Signature: {Function input: }(ET convolution time series from 2.2, fMRI mean time series from 1.1, fMRI sampling rate (1.127))

#### 3. Confounder Processing

3.1 **confounder process function**: To process the nuisance covariates of head motion (In this case: 6 head motion control, 6 quadratics of 6 head motion control.)

   - Function Signature: {Function input:} (path of the file, fMRI mean time series for 1.1)

3.2 **design matrix function**: To construct the design matrix [Note, need to construct the eye blink, fixation, and nuisance covariates (Need to be scaled).]

   - Function Signature: {Function input:} (ET convolution from 2.3, nuisance covariates form 3.1) $\rightarrow$ {Returns:} T by J Design Matrix

#### B.1 General Linear Regression:

4. **Brain Activation**: To execute autoregressive ARIMA regression to examine the relationship between ET data and brain region activation.

   -  Function Signature: {Function input:} (fMRI mean time series form 1.1, design matrix from 3.2 ) $\rightarrow$ {Returns:} (list of 1. ET-blink coefs, 2. ET-fixation coefs, 3. Residuals 4. Covariance matrix)

#### B.2. Covariance Regression:

5. **OLSoutcome process function**: To construct the logarithm of \(Y\) as the outcome variable for regression analysis.

   - Function Signature: {Function input:} (residuals from 4) $\rightarrow$ {Returns:} (logarithm of T by $\frac{V(V+1)}{2}$ matrix)

6.1 **OLSfitting function**: To perform LS regression for all $\frac{V(V+1)}{2}$ regions with robust estimation techniques. [HC 3]
   - Function Signature: {Function input:} (Log.Y from 5, Design matrix from 3.2) $\rightarrow$ {Returns:} (list of 1. Coefficients, 2. Robust covariance matrix)

6.2 **Comparison statistic function**: To compare the z-statistic between the slope of Eye-blink and Eye-fixation estimated coefficients from 6.1
   - Function Signature: {Function input:} (list from 6.1) $\rightarrow$ {Returns:} (V by V Z-statistic matrix)

7. **Model Diagnostics for OLS function**: To randomly pick the OLS model to do the model diagnostics on residuals
   - Function Signature: {Function input:} (1. Log.Y from 5, 2. Design matrix form 3.2) $\rightarrow$ {Returns:} (list of 1. ACF plot, 2.PACF plot, 3. Q-Q plot, 4.Histogram)

