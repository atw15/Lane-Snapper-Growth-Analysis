# Lane Snapper Growth Analysis

This repository contains the R code and data analysis workflow for studying the growth patterns of Lane Snapper. The analysis involves data cleaning, exploratory data analysis, model fitting using von Bertalanffy growth functions, and Bayesian modeling to compare different models and estimate growth parameters.

## Table of Contents

- [Introduction](#introduction)
- [Data Description](#data-description)
- [Data Preparation and Loading Libraries](#data-preparation-and-loading-libraries)
- [Data Exploration and Cleaning](#data-exploration-and-cleaning)
- [Exploratory Data Analysis](#exploratory-data-analysis)
- [Modeling and Statistical Analysis](#modeling-and-statistical-analysis)
- [Predictions and Plotting](#predictions-and-plotting)
- [Bayesian Modeling](#bayesian-modeling)
- [Residual Analysis](#residual-analysis)
- [Model Comparison](#model-comparison)
- [Summarizing Data](#summarizing-data)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)
- [Repository Structure](#repository-structure)
- [Running the Script](#running-the-script)

## Introduction

This project aims to analyze the growth patterns of Lane Snapper using a dataset of fish lengths and ages. The analysis includes data cleaning, exploratory data analysis, fitting growth models, and comparing the fits of different models. Bayesian methods are used to obtain credible intervals for the growth parameters.

## Data Description

The dataset `Age Master Final_5_12.csv` includes the following columns:
- **Final.Age**: The final age of the fish in years.
- **FFL.mm**: Fork length of the fish in millimeters.
- **Fishery**: The type of fishery (e.g., fishery-independent (FI) or fishery-dependent).
- **Macro.sex**: The sex of the fish (e.g., F for female, M for male, U for unknown).

## Data Preparation and Loading Libraries

- **Clearing Environment**: The script starts by removing all existing R objects to ensure a clean workspace.
- **Loading Data**: The dataset `Age Master Final_5_12.csv` is loaded, and summary statistics for the `Final.Age` column are displayed.
- **Installing and Loading Libraries**: Required packages (`devtools`, `magrittr`, `FSA`, `dplyr`, `car`, `nlstools`, `patchwork`, and `FSAmisc` from GitHub) are installed and loaded.

## Data Exploration and Cleaning

- **Examine the Data**: The `Final.Age` column is inspected to check its values.
- **Data Cleaning**: Rows with missing values are removed using `na.omit`.
- **Summary Statistics**: Summary statistics for the `Final.Age` column are displayed after cleaning.

## Exploratory Data Analysis

- **Visualizing Distributions and Relationships**:
  - Histograms and scatter plots are created to visualize the distribution of ages and the relationship between fork length (`FFL.mm`) and age (`Final.Age`).

## Modeling and Statistical Analysis

- **Model Fitting and Comparison**:
  - Functions for the von Bertalanffy growth model with truncated normal and regular normal distributions are defined.

## Predictions and Plotting

- **Predictions and Plotting**:
  - Predictions are made using the fitted models.
  - The data and model predictions are plotted to visualize the fit.

## Bayesian Modeling

To examine the potential bias associated with the minimum size limit on the Lane Snapper fishery, Bayesian models were used to re-estimate growth parameters, focusing on the models preferred in the preliminary model evaluation (Burton et al., 2019; McGarvey et al., 2002). A truncated normal probability density function was used for the residual variability in length, assuming a zero probability of capture below the minimum size limit of the fishery-dependent samples.

To incorporate the 203.2 mm (8 inch) TL size limit into the truncated normal probability function, we converted the TL size limit to FL using a linear regression model fitted to the TL and FL data from the samples so that the lower limit of the truncated normal would be consistent with the data. The Deviance Information Criterion (DIC) was used for model selection (Lunn et al., 2012). The full, untruncated normal likelihood was applied to all fishery-independent samples, which were not subject to the minimum size limit.

The models were fit using a Gibbs Sampler implemented in JAGS, and the Gelman-Rubin diagnostic and effective sample size were used to evaluate convergence. R version 4.2.2 was used for statistical analysis and the R-package ‘FSA’ was used for among-group comparisons, model fitting, and selection for the non-Bayesian analysis (Ogle et al., 2023; R Core Team, 2022). JAGS version 4.3.1 and the R-packages “rjags” and “R2jags” were used for the analysis of the Bayesian models (Plummer).

### Truncated Normal Model with Sex Dependency

- Models sex-specific growth parameters for fishery-dependent and fishery-independent samples.
- Outputs: Residuals, predicted values, and model parameters.

### Combined Sexed Model Truncated

- Models combined sexed samples using truncated normal distribution.
- Outputs: Residuals, predicted values, and model parameters.

### Combined Sexed Model Normal

- Models combined sexed samples using regular normal distribution.
- Outputs: Residuals, predicted values, and model parameters.

### Sex-Specific Model Normal

- Models sex-specific growth parameters without truncation.
- Outputs: Residuals, predicted values, and model parameters.

## Residual Analysis

- **Residual Analysis**:
  - Residuals from the models are analyzed to check for goodness-of-fit.
  - Q-Q plots are created to check the normality of residuals.

## Model Comparison

- **Model Comparison**:
  - DIC values are calculated for each model to compare model fits.
  - The model with the smallest DIC is considered the best-fitting model.

## Summarizing Data

- **Summarizing Data**:
  - Data are summarized by fishery type, including the number of samples, mean, standard deviation, standard error, and confidence intervals.

## Citation

This work includes sections inspired by Derek H. Ogle's book:
- Ogle, Derek H. Introductory Fisheries Analyses with R. 1st ed. vol. 32. Milton: CRC Press, 2016. Web.

## Acknowledgments

This project was developed under the guidance of Dr. Elizabeth Babcock. Her expertise and support were invaluable in completing this analysis.

---

### Repository Structure

- `README.md`: This file.
- `Annotated Age and Growth Analysis.R`: The R script containing the annotated code.
- `Annotated Truncated Growth Models.R`: The R script containing the truncated normal and regular normal growth models.
