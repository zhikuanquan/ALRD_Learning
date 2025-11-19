\section*{README: Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data}

This repository contains the R code used in:

\textbf{Quan, Z., Iosif, A.-M., \& Chen, S. (2025).  
\textit{Doubly Robust Estimation of Individualized Treatment Effects in Longitudinal Data}.}

It includes:
\begin{itemize}
    \item \textbf{Simulation study} evaluating ALRD under multiple data-generating mechanisms.
    \item \textbf{Real-data analysis} applying ALRD to MIA monkey MRI and cytokine data.
\end{itemize}

\subsection*{File Structure}
\begin{verbatim}
/ 
├── real_data_analysis.R       # ALRD pipeline for real MIA dataset
├── simulation_scenario.R      # SIM() + data generation for simulation study
└── README.tex
\end{verbatim}

\subsection*{Required R Packages}
\begin{verbatim}
install.packages(c(
  "tidyverse", "readxl", "glmnet", "glmmLasso", "lme4",
  "mgcv", "gratia", "caret", "MASS", "nlme", "Matrix",
  "ggpubr", "grid"
))
\end{verbatim}

\subsection*{Real-Data Analysis}
Place the dataset at:
\begin{verbatim}
DATA/hyp_all.xlsx
\end{verbatim}

Run the ALRD analysis:
\begin{verbatim}
source("real_data_analysis.R")

fit <- ALRD(
  subregion      = "B_frontal",   # or "B_prefrontal"
  response       = "GM",          # or "WM"
  DATA           = DATA,
  B              = 2,
  lambda         = seq(1,100,1),
  lambda.min_aug = 89
)

fit$coef_table_orig   # coefficients on original cm^3 scale
fit$f_1               # partial ITE plot 1
fit$f_2               # partial ITE plot 2
\end{verbatim}

\subsection*{Simulation Study}
Run one simulation scenario:
\begin{verbatim}
source("simulation_scenario.R")

result_list <- SIM(
  ICC_type    = "b",
  main_effect = "big",
  data_type   = "lin",
  n           = 100,
  K           = 5,
  N           = 100
)
\end{verbatim}

Simulation outputs include:
\begin{itemize}
    \item Accuracy of individualized treatment rule
    \item Spearman correlation with true ITEs
    \item Average prediction error
    \item Optimal lambda for LASSO and augmented LASSO
\end{itemize}

\subsection*{Citation}
\begin{verbatim}
Quan Z., Iosif A.-M., Chen S. (2025).
Doubly Robust Estimation of Individualized Treatment Effects 
in Longitudinal Data.
\end{verbatim}

\subsection*{Contact}
\textbf{Zhikuan (Zach) Quan} \\
PhD Candidate, Biostatistics, UC Davis \\
Email: \texttt{zkquan@ucdavis.edu}
