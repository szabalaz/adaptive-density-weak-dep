# Adaptive Density Estimation for Weakly Dependent Data

![Language](https://img.shields.io/badge/Language-R-blue)
![License](https://img.shields.io/badge/License-MIT-green)

This repository contains the **R** implementation of non-parametric kernel density estimation methods, developed as part of my academic research.

The main focus is on **adaptive bandwidth selection** using the **Goldenshluger-Lepski (GL) method**, specifically tailored for contexts where data is not independent and identically distributed (i.i.d.), but exhibits **weak dependence** (e.g., time series, stochastic processes).

## 📋 Overview

Standard methods for bandwidth selection (like Cross-Validation) often fail or perform poorly when the independence assumption is violated. This project implements theoretical estimators that ensure optimal convergence rates even under weak dependence conditions, without prior knowledge of the density's regularity.

## ✨ Key Features

*   **Kernel Density Estimation (KDE):** Implementation of standard kernel estimators.
*   **Adaptive Bandwidth Selection:** Full implementation of the Goldenshluger-Lepski (GL) algorithm.
*   **Weak Dependence Handling:** Adjustments and simulations for time series data.
*   **Visualization:** High-quality plots using `ggplot2` to illustrate:
    *   The Bias-Variance tradeoff.
    *   The "Anchor Point" effect in histograms vs. Kernels.
    *   Convergence of adaptive estimators.

## 🛠️ Prerequisites

To run the scripts, you will need **R** installed along with the following packages:

```r
install.packages(c("ggplot2", "KernSmooth", "dplyr"))
```
## 📂 Repository Structure
/scripts: R source codes for the estimators and GL method.
/plots: Generated figures used in the thesis/paper.
/simulations: Scripts to generate synthetic weakly dependent data (e.g., AR(1), MA(1) processes) for testing purposes.
## 📄 License
This project is licensed under the MIT License - see the LICENSE file for details.
