# EyeTrackingFMRI
Analysis of the simultaneous eye-tracking and movie-watching fMRI data. Includes 1) task activation modeling using the general linear model and 2) covariance regression.

#### A. Data Processing:

1.1 **fMRI process function**: To extract the fMRI time series for each brain region, calculating mean values across voxel by time points.

   - **Equation**: $\text{fMRIprocess}(vertex, time\_points) \rightarrow \text{Mean fMRI Time Series}$

1.2 **fMRIMapping Plot function**: To map the xii_pmean/coefficients/t-statistics to our brain graph 

   - **Return** return xii_seed_long to plot the data in -log(10) scale
 
2.1 **ET Data Process function**: To process the eye tracking data and get the input parameter for *function 1*.

2.2 **Convolution function**: To perform the convolution between eyeblink & eyefixation and the HRF (Hemodynamic Response Function).

    - **Equation**: $\text{convolution}(eye-blink, eye-fixation, HRF)$

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
