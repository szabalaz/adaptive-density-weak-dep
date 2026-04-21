# Adaptive Density Estimation for Weakly Dependent Data

![Language](https://img.shields.io/badge/Language-R-blue)
![License](https://img.shields.io/badge/License-MIT-green)

This repository contains **R** and **Python** implementations of nonparametric kernel density and regression estimation methods, developed as part of academic research.

The main focus is on **adaptive density and regression estimation** using the **Goldenshluger-Lepski (GL) method**, specifically for data with **weak dependence** (e.g., time series or stochastic processes) where standard i.i.d. assumptions do not hold.

## 📋 Overview

Classical bandwidth selection methods assume independence and that the smoothness of the density (regularity parameter $B$) is known, usually $B = 2$. In practice, we usually do not know the density function, and therefore we certainly do not know its derivatives. When these assumptions fail, classical methods should not be expected to perform well.

This project contrasts classical density estimation for independent data with adaptive density and regression estimation for weakly dependent data, using the GL method. The approach constructs kernel estimators that belong to more general function classes, such as Hölder classes, where the optimal bandwidth does not depend on the unknown regularity $B$ and the estimator converges at nearly optimal rates.

## ✨ Key Features

* **Kernel Density Estimation (KDE):** Standard kernel-based density estimators.
* **Adaptive Bandwidth Selection:** Implementations and experiments with the Goldenshluger-Lepski algorithm.
* **Weak Dependence:** Simulations and examples for dependent time series models.
* **Visualization:** Plots comparing estimators, bandwidth choices, and kernel effects.

## 🛠️ Prerequisites

### R
Install R and the required packages:

```r
install.packages(c(
  "ggplot2", "KernSmooth", "dplyr", "np", "forecast", "tseries", "latex2exp"
))
```

### Python
If you use the Python examples, create and activate a virtual environment:

```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## 📂 Repository Structure

* `scripts/adaptive_estimation/`: Core adaptive estimation implementations and validation scripts.
* `scripts/examples/`: Jupyter/R notebooks and reproducible examples.
* `scripts/functions/`: Utility functions used across scripts.
* `scripts/simulations/`: Synthetic data generation for AR, MA, and white noise models.
* `scripts/plots/`: Generated figures and visualizations.

## ▶️ Usage

Run a script from the repository root, for example:

```powershell
Rscript scripts/adaptive_estimation/validacion_completa_estimacion_adaptativa_densidad.R
```

Or open one of the notebook examples:

```powershell
jupyter notebook scripts/examples/histograma_estimador_densidad_python.ipynb
```

## 📄 License
This repository is licensed under the MIT License. See the `LICENSE` file for full terms.
