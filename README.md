# Quantifying Under-Reporting in Longitudinal Health Data: Genital Warts Example

This repository reproduces the results from the study "Quantifying the Under-Reporting of Uncorrelated Longitudinal Data: The Genital Warts Example" by Moriña et al. (2021). The aim is to estimate the true incidence of genital warts (GW) cases in Catalunya by addressing the under-reporting in health data. The study uses a mixture model to quantify under-reporting and adjusts for missing data, particularly for the period 2009-2016. I did this project during my mathematical modeling workshop in UAB under DR. David Moriña.

## Project Overview

This project investigates the use of mixture models for quantifying under-reporting in longitudinal health data. It focuses on genital warts incidence in Catalunya, where the disease is often under-reported due to its asymptomatic nature. By employing a mixture model approach, we estimate the true burden of genital warts in the population.

The analysis uses data from the Catalan public health system (SIDIAP) for the period 2009-2016. The goal is to reconstruct the actual incidence of genital warts by modeling the under-reporting process as a mixture of two distributions.


## Contents

### 1. **Data Preprocessing and Model Fitting**
   - **Population**: All residents of Catalunya assigned to an ICS primary care center.
   - **Model**: A mixture model with two components—one for the observed incidence and one for under-reported incidence. The under-reporting frequency and intensity are modeled as parameters that vary over time.

### 2. **Mixture Model and Expectation-Maximization (EM) Algorithm**
   - **Mixture Model**: A weighted sum of two normal distributions to account for under-reported cases. The mixing proportion (ωt) and intensity (q) of under-reporting are estimated.
   - **Parameter Estimation**: Using Maximum Likelihood Estimation (MLE) via the EM algorithm to estimate the mixing proportion and the parameters of the distributions.

### 3. **Model Validation**
   - **Goodness of Fit**: Akaike’s Information Criterion (AIC) is used to compare the mixture model against a single normal distribution model.
   - **Residuals Analysis**: Residuals are checked to ensure the model adequately captures the data's structure and that no significant autocorrelations remain.

### 4. **Results**
   - The study estimates that only around 80% of the actual genital warts incidence was recorded in the SIDIAP database.
   - The incidence of under-reporting was found to be more significant among women over 30 years old.
   - **Estimated vs. Registered Cases**: The reconstruction showed a 23% underestimation of GW cases.

### 5. **Impact of Under-Reporting**
   - The study found that under-reporting had a substantial impact on the true burden of genital warts in Catalunya, with an estimated under-reporting of almost 10,000 cases (23% of the registered cases).
   - The total annual cost for genital warts treatment in Catalunya was underestimated by about 10 million Euros.

### 6. **Code and Data**
   - All code used for fitting the models and obtaining results and figures is available in the repository.

---

## Methods

### 1. **Mixture Model for Under-Reporting**
   - The registered data (Yt) is modeled as a mixture of two normal distributions: one for the observed incidence and one for under-reported cases.
   - The probability of under-reporting (ωt) and the intensity of under-reporting (q) are estimated as part of the mixture model.

### 2. **Expectation-Maximization (EM) Algorithm**
   - The initial parameters for the EM algorithm are obtained using an Expectation-Maximization approach.
   - The algorithm iteratively estimates the component membership and updates parameters to maximize the log-likelihood function.

### 3. **Maximum Likelihood Estimation (MLE)**
   - The parameters are estimated by maximizing the log-likelihood of the observed data.
   - The log-likelihood is computed using the normal distribution for each observation in the dataset.

---

## Results and Analysis

- **Under-Reporting**: On average, 20% of genital warts cases were under-reported, with significant variation based on age and sex.
- **Estimation Accuracy**: The reconstructed incidence showed an increase in the reported cases, especially for women over 30 years old.
- **Impact on Healthcare Costs**: The under-reporting led to an underestimation of healthcare costs by around 10 million Euros annually.

---

## Code and Implementation

- The models and estimation methods are implemented in **R**.
- The key R package used for parameter estimation is `mixtools`, which is employed for fitting the mixture model and obtaining the maximum likelihood estimates.
## References
Moriña, D., Fernández-Fontelo, A., Cabaña, A., Puig, P., Monfil, L., Brotons, M., & Diaz, M. (2020). Quantifying the Under-Reporting of Uncorrelated Longitudinal Data: The Genital Warts Example. BMC Medical Research Methodology.

## Contact

For any questions or inquiries, please contact:
- Vishal Nair
- Email: [v1292002@gmail.com](mailto:v1292002@gmail.com)
