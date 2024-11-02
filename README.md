# Smooth Transition Regression Model Analysis

## Overview
This repository contains an implementation of Smooth Transition Regression (STR) models in MATLAB. The focus is on simulating data, estimating model parameters, and comparing different regression techniques, including logistic and exponential smoothing functions.

## Tasks

### Task 4.1: Simulate Smooth Transition Model and Compute Loss Score
- **Transition Variable**
- **Parameters**: Set lambda = 3 and c = 0.3.
- **Computation**: The function is derived.
- **Matrix Definition**: Constructed matrix X_t with two columns.
- **Model Simulation**: Generated data.
- **Loss Function**: Implemented a function to compute the loss score for the least squares estimator of the exponential STR model.

### Task 4.2: Prepare a Grid Search of Initial Parameters
- **Data Loading**: Loaded the dataset 'Report4 2'.
- **Parameter Grids**:
  - **Lambda Grid**: Logarithmically spaced values between 0.01 and 100.
  - **C Grid**: Linearly spaced values from the 10th to the 90th percentile of the transition variable \( s \).
- **Estimation Process**:
  - For each pair of lambda and c, estimated beta_1 and beta_2^*.
  - Calculated beta_2 = beta_2^* - beta_1 and computed the loss score.
- **Optimization**: Selected parameters lambda, c that minimize the loss score and re-estimated the corresponding beta values.

### Task 4.3: Estimating Parameters of Smooth Transition Regression
- **Data Loading**: Loaded the dataset 'Report4 3'.
- **Initial Parameters**: Based on data, initialized the parameters vector containing beta_1, beta_2, lambda, c.
- **Parameter Estimation**: Used `fminunc` and `fmincon` for optimization under specified constraints.
- **Comparison Table**: Prepared a comparison table for estimated parameters and loss scores.

### Task 4.4: Additional Model Comparison
- **Model Definition**: Evaluated the model for GDP growth as a function of unemployment, inflation, and renewable energy sources.
- **Data Preparation**: 
  - Downloaded data from [World Bank](https://data.worldbank.org).
  - Ensured completeness of the dataset, excluding missing values for each country specified by 'countrycode'.
- **Model Estimation**:
  - Implemented logistic and exponential STR models with grid searches for parameters.
  - Estimated parameters using `fmincon` under constraints.
- **GDP Estimation**: Calculated GDP estimates using both STR models and a baseline linear regression model.
- **Visualization**: Plotted actual GDP against estimated values from all models, including legends and axis labels.
- **BIC Comparison**: Calculated the Bayesian Information Criterion (BIC) for each model to facilitate comparison.

## Programming Language
This project was implemented in **MATLAB**.
