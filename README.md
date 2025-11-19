# Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data

This repository contains the R code used in:

**Quan, Z., Iosif, A-M., & Chen, S. (2025).  
_Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data_.**

It includes:

- **Simulation study** evaluating ALRD under multiple data-generating mechanisms  

---

## File Structure
```
/
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
