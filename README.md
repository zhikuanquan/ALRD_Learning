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

