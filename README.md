# Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data

This repository contains the R code used in:

**Quan, Z., Iosif, A-M., & Chen, S. (2025).  
_Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data_.**

It includes:

- **Simulation study** evaluating ALRD under multiple data-generating mechanisms  
- **Real-data analysis** applying ALRD to MIA monkey MRI + cytokine data  

---

## File Structure
```
/
├── real_data_analysis.R # ALRD pipeline for real MIA dataset
├── simulation_scenario.R # SIM() + data generation for simulation study
└── README.md

---
```

## Requirements

Install required R packages:

```r
install.packages(c(
  "tidyverse", "readxl", "glmnet", "glmmLasso", "lme4",
  "mgcv", "gratia", "caret", "MASS", "nlme", "Matrix",
  "ggpubr", "grid"
))
```
## Real-Data Analysis

Place your dataset at:

```
DATA/hyp_all.xlsx
```

Run the ALRD analysis:
```r
source("real_data_analysis.R")

fit <- ALRD(
  subregion      = "B_frontal",   # or "B_prefrontal"
  response       = "GM",          # or "WM"
  DATA           = DATA,
  B              = 2,
  lambda         = seq(1,100,1),
  lambda.min_aug = 89
)

fit$coef_table_orig   # coefficients on original cm³ scale
fit$f_1               # partial ITE plot 1
fit$f_2               # partial ITE plot 2
```

## Simulation Study

Run a simulation scenario:
```r

source("simulation_scenario.R")

result_list <- SIM(
  ICC_type    = "b",
  main_effect = "big",
  data_type   = "lin",
  n           = 100,
  K           = 5,
  N           = 100
)

```
Simulation output includes:

- Accuracy of individualized treatment rule

- Spearman correlation with true ITEs

- Average prediction error

- Optimal λ for LASSO and augmented LASSO


## Contact

Zhikuan Quan
PhD Candidate, Biostatistics, UC Davis
zkquan@ucdavis.edu
