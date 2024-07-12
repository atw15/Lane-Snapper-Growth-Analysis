# Lane Snapper Growth Analysis

This repository contains the R code and data analysis workflow for studying the growth patterns of Lane Snapper. The analysis involves data cleaning, exploratory data analysis, model fitting using von Bertalanffy growth functions, and Bayesian modeling to compare different models and estimate growth parameters.

## Table of Contents

- [Introduction](#introduction)
- [Data Preparation and Loading Libraries](#data-preparation-and-loading-libraries)
- [Data Exploration and Cleaning](#data-exploration-and-cleaning)
- [Exploratory Data Analysis](#exploratory-data-analysis)
- [Modeling and Statistical Analysis](#modeling-and-statistical-analysis)
- [Predictions and Plotting](#predictions-and-plotting)
- [Bayesian Modeling](#bayesian-modeling)
- [Residual Analysis](#residual-analysis)
- [Model Comparison](#model-comparison)
- [Summarizing Data](#summarizing-data)

## Introduction

This project aims to analyze the growth patterns of Lane Snapper using a dataset of fish lengths and ages. The analysis includes data cleaning, exploratory data analysis, fitting growth models, and comparing the fits of different models. Bayesian methods are used to obtain credible intervals for the growth parameters.

## Data Preparation and Loading Libraries

- **Clearing Environment:** The script starts by removing all existing R objects to ensure a clean workspace.
- **Loading Data:** The dataset `Age Master Final_5_12.csv` is loaded, and summary statistics for the `Final.Age` column are displayed.
- **Installing and Loading Libraries:** Required packages (`devtools`, `magrittr`, `FSA`, `dplyr`, `car`, `nlstools`, `patchwork`, and `FSAmisc` from GitHub) are installed and loaded.

## Data Exploration and Cleaning

- **Examine the Data:** The `Final.Age` column is inspected to check its values.
- **Data Cleaning:** Rows with missing values are removed using `na.omit`.
- **Summary Statistics:** Summary statistics for the `Final.Age` column are displayed after cleaning.

## Exploratory Data Analysis

- **Visualizing Distributions and Relationships:**
  - Histograms and scatter plots are created to visualize the distribution of ages and the relationship between fork length (`FFL.mm`) and age (`Final.Age`).

## Modeling and Statistical Analysis

- **Model Fitting and Comparison:**
  - Functions for the von Bertalanffy growth model with truncated normal and regular normal distributions are defined.
  - Models are fitted using the `nlminb` function, and parameter estimates and log-likelihood values are compared.

## Predictions and Plotting

- **Predictions and Plotting:**
  - Predictions are made using the fitted models.
  - The data and model predictions are plotted to visualize the fit.

## Bayesian Modeling

- **Bayesian Modeling:**
  - JAGS is used for Bayesian modeling to obtain credible intervals for the growth parameters.
  - Fishery-independent and fishery-dependent data are separated to account for size limits.

## Residual Analysis

- **Residual Analysis:**
  - Residuals from the models are analyzed to check for goodness-of-fit.
  - Q-Q plots are created to check the normality of residuals.

## Model Comparison

- **Model Comparison:**
  - Models are compared using Deviance Information Criterion (DIC) to determine the best-fitting model.

## Summarizing Data

- **Summarizing Data:**
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
- `Age Master Final_5_12.csv`: The dataset used in the analysis.
