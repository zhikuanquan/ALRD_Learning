#################################
###     LOADING PACKAGE       ###
#################################
library(readxl)
library(tidyverse)
library(glmnet)
library(caret)
library(glmmLasso)
library(lme4)
library(ggpubr)
library(grid)  
set.seed(100)
#################################
###       INPUT DATA          ###
#################################
DATA <- read_excel("Desktop/DATA/hyp_all.xlsx")
DATA <- na.omit(DATA)
#################################
###           EDA             ###
#################################

DATA_BASE <- DATA %>%
  group_by(Subject_ID) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup
DATA_BASE <- as.data.frame(DATA_BASE)
DATA_BASE <- as.data.frame(DATA_BASE[-6,])
for (i in 1:dim(DATA)[1]){
  DATA[i,63] <- ifelse(DATA$time[i]==1,1,0)
  DATA[i,64] <- ifelse(DATA$time[i]==2,1,0)
  DATA[i,65] <- ifelse(DATA$time[i]==3,1,0)
  DATA[i,66] <- ifelse(DATA$time[i]==4,1,0)
  DATA[i,67] <- ifelse(DATA$time[i]==5,1,0)
}
colnames(DATA)[63:67]  <- paste('t',unique(DATA$time),sep='')
####################################################
###        MODEL: Response ~ T*(X)+T*time        ###
####################################################
ALRD_Lasso <- function(B=2,lambda = c(seq(1,100,1)),DATA_ALL){
  ###### Find the best lambda for Lasso by 2-fold CV with the smallest MSE
  # B = 2 ## 2-folds
  # lambda = c(seq(1,100,1))
  # DATA_ALL: Fitting Data with Augmentation
  MSE_aug <- matrix(NA,length(lambda),B)
  set.seed(101)
  for(i in 1:length(lambda)){
    try({
      flds <- createFolds(unique(DATA_ALL$Subject_ID), k = B, list = TRUE, returnTrain = FALSE)
      for(b in 1:B)
      {
        sel_ind <- which(DATA_ALL$Subject_ID%in%unique(DATA_ALL$Subject_ID)[flds[[b]]])
        data_in <- DATA_ALL[sel_ind,]
        data_in$Subject_ID <- factor(data_in$Subject_ID) 
        data_out <- DATA_ALL[-sel_ind,]
        data_out$Subject_ID <- factor(data_out$Subject_ID) 
        try({
          fit_aug <- glmmLasso(as.formula(form_aug),rnd = list(Subject_ID=~1),lambda=lambda[i],data=data_out,control = list(center=FALSE))
          MSE_aug[i,b] = mean((predict(fit_aug,newdata=data_in)-data_in$N_Y)^2)
        },silent=TRUE)
      }
    },silent=TRUE)
    print(i)
  }
  
  MSE_aug1=apply(MSE_aug,1,function(x){mean(x,na.rm=T)})
  MSE_aug1[apply(is.na(MSE_aug),1,mean)>0.8]=NA
  
  lambda.min_aug = min(lambda[MSE_aug1==min(MSE_aug1,na.rm=T)],na.rm=T)
  resultsss <- data.frame(lambda=lambda,MSE = MSE_aug1)
  return(list(lambda.min_aug,resultsss))
}

ALRD <- function(subregion, response, DATA, 
                 B = 2, lambda = c(seq(1,100,1)), lambda.min_aug = 0) {
  
  # --- 0) Prep / filtering ----------------------------------------------------
  DATA_ALL <- as.data.frame(DATA[DATA$Subject_ID != 6, ])
  colnames(DATA_ALL)[colnames(DATA_ALL) == response] <- "Y"
  DATA_ALL <- DATA_ALL[DATA_ALL$hemi_region == subregion, ]
  DATA_ALL <- na.omit(DATA_ALL)
  
  # Encode treatment as +/-1 if needed
  for (i in 1:nrow(DATA_ALL)) {
    DATA_ALL[i,58] <- as.numeric(ifelse(DATA_ALL[i,58] == "PolyICLC", 1, -1))
  }
  DATA_ALL$group <- as.numeric(DATA_ALL$group)
  
  # Keep a copy (you were using this earlier)
  DATA_o <- DATA_ALL
  
  # --- 1) Put Y on your GM/1000 scale *before* standardization ---------------
  DATA_ALL$Y <- DATA_ALL$Y / 1000
  DATA_o$Y   <- DATA_o$Y   / 1000
  DATA_ALL$base_IL_8 <- DATA_ALL$base_IL_8/1000
  DATA_o$base_IL_8 <- DATA_o$base_IL_8/1000
  
  # Indices / names
  idx_cont   <- c(2:3, 5:26)     # continuous baseline X to standardize
  idx_other  <- c(63:67)         # time dummies (not standardized)
  Xcont_names  <- colnames(DATA_ALL)[idx_cont]
  Xother_names <- colnames(DATA_ALL)[idx_other]
  
  # --- 2) Compute and store scaling info -------------------------------------
  meanY <- mean(DATA_ALL$Y)
  sdY   <- sd(DATA_ALL$Y)
  
  meanX <- sapply(DATA_ALL[, idx_cont, drop = FALSE], mean)
  sdX   <- sapply(DATA_ALL[, idx_cont, drop = FALSE], sd)
  
  scale_info <- list(
    meanY = meanY, sdY = sdY,
    meanX = meanX, sdX = sdX,
    Xcont_names = Xcont_names,
    Xother_names = Xother_names
  )
  
  # --- 3) Standardize Y and continuous X (leave dummies untouched) -----------
  DATA_STD <- DATA_ALL
  DATA_STD$Y <- (DATA_STD$Y - meanY) / sdY
  for (nm in Xcont_names) {
    DATA_STD[[nm]] <- (DATA_STD[[nm]] - meanX[[nm]]) / sdX[[nm]]
  }
  # dummies in Xother_names remain as-is
  
  # --- 4) Build formulas ------------------------------------------------------
  X.name0  <- c(Xcont_names, Xother_names)             # covariates allowed in interactions
  AX.name  <- paste(paste("group", X.name0, sep=":"), collapse="+")
  form_non <- paste("Y ~ -1 + group +", AX.name)       # lasso (non-aug) form
  form_aug <- paste("N_Y ~ -1 + group +", AX.name)     # augmented lasso form
  
  # --- 5) Augmentation (lasso on standardized scale) -------------------------
  # Fit separate lasso models within treatment groups on standardized data
  data_1 <- DATA_STD[DATA_STD$group ==  1, ]
  data_0 <- DATA_STD[DATA_STD$group == -1, ]
  
  cv_1 <- cv.glmnet(x=as.matrix(cbind(1,data_1[,c(2:3,5:26)])), y=data_1$Y, alpha = 1)
  m1 <- glmnet(x=as.matrix(cbind(1,data_1[,c(2:3,5:26)])), y=(data_1$Y), alpha = 1,lambda = cv_1$lambda.min)
  Y_1 <- 0.5*(as.matrix(cbind(1,DATA_STD[,c(2:3,5:26)]))%*%as.matrix(m1$beta,ncol=1))
  
  cv_0 <- cv.glmnet(x=as.matrix(cbind(1,data_0[,c(2:3,5:26)])), y=data_0$Y, alpha = 1)
  m0 <- glmnet(x=as.matrix(cbind(1,data_0[,c(2:3,5:26)])), y=data_0$Y, alpha = 1,lambda = cv_0$lambda.min)
  Y_0 <- 0.5*(as.matrix(cbind(1,DATA_STD[,c(2:3,5:26)]))%*%as.matrix(m0$beta,ncol=1))
  
  DATA_STD$N_Y <- as.numeric(DATA_STD$Y - (Y_1 + Y_0))
  
  # --- 6) Variable selection with glmmLasso on standardized scale ------------
  DATA_STD$Subject_ID <- factor(DATA_STD$Subject_ID)
  fit_dam_aug <- glmmLasso(as.formula(form_aug),
                           rnd = list(Subject_ID = ~1),
                           lambda = lambda.min_aug,
                           data = DATA_STD,
                           control = list(center = FALSE))
  
  # Keep only nonzero interactions
  kept_cov <- X.name0[fit_dam_aug$coefficients[-1] != 0]
  AX.keep  <- paste(paste("group", kept_cov, sep=":"), collapse="+")
  form_final <- paste("Y ~ -1 + group +", AX.keep, " + (1|Subject_ID)")
  
  # --- 7) Final lmer on standardized scale -----------------------------------
  model_std <- lmer(as.formula(form_final), data = DATA_STD)
  ssum <- summary(model_std)
  
  # Fixed effects (standardized)
  fe_std <- as.data.frame(ssum$coefficients)
  fe_std$term <- rownames(fe_std)
  rownames(fe_std) <- NULL
  names(fe_std)[1:2] <- c("estimate_std", "se_std")  # keep the essentials
  
  # --- 8) Rescale fixed effects back to original Y (/1000) scale -------------
  rescale_factor <- function(term) {
    if (term == "group") {
      # ITE intercept: only Y was standardized
      return(sdY)
    }
    if (startsWith(term, "group:")) {
      var <- substring(term, 7)  # after "group:"
      if (var %in% Xcont_names)  return(sdY / sdX[[var]])
      if (var %in% Xother_names) return(sdY)        # unstandardized dummies
    }
    # default: no rescaling (shouldn't hit here)
    return(1)
  }
  
  fe_std$scale_factor <- vapply(fe_std$term, rescale_factor, numeric(1))
  fe_std$estimate_orig <- fe_std$estimate_std * fe_std$scale_factor
  fe_std$se_orig       <- fe_std$se_std       * fe_std$scale_factor
  
  # Arrange nice output
  coef_table_orig <- fe_std[, c("term","estimate_orig","se_orig","estimate_std","se_std","scale_factor")]
  
  ## GAMM
  ID_gamm <- DATA_STD$Subject_ID
  group_gamm <- DATA_STD$group
  if(subregion == "B_frontal"&response  == "GM"){
    formu = paste("Y ~-1+ group_gamm+s(base_IL_17a,k=3,by=group_gamm)+s(base_IL_8,k=3,by=group_gamm)+s(ID_gamm, bs = 're')")
    y_label <- "Partial ITE in Frontal GM (cm³)"
  }
  if(subregion == "B_prefrontal"&response  == "GM"){
    formu = paste("Y ~-1+ group_gamm+s(base_IL_17a,k=3,by=group_gamm)+s(base_IL_8,k=3,by=group_gamm)+s(base_MIP_1a,k=3,by=group_gamm)+s(ID_gamm, bs = 're')")
    y_label <- "Partial ITE in Prefrontal GM (cm³)"
  }
  if(subregion == "B_frontal"&response  == "WM"){
    formu = paste("Y ~-1+ group_gamm+base_G_CSF:group_gamm+Dam_GD40_weight__kg_:group_gamm+s(base_IL_17a,k=3,by=group_gamm)+s(base_IL_8,k=2,by=group_gamm)+s(ID_gamm, bs = 're')")
    y_label <- "Partial ITE in Frontal WM (cm³)"
  }
  if(subregion == "B_prefrontal"&response  == "WM"){
    formu = paste("Y ~-1+ group_gamm+s(base_IL_17a,k=3,by=group_gamm)+s(base_IL_8,k=2,by=group_gamm)+s(Dam_Age_at_conception__yr_,by=group_gamm)+s(Dam_GD40_weight__kg_,by=group_gamm)+s(ID_gamm, bs = 're')")
    y_label <- "Partial ITE in Prefrontal WM (cm³)"
  }
  
  fit_gamm <- gam(formula=as.formula(formu), data = cbind(DATA_STD,ID_gamm,group_gamm))
  
  sm1 <- smooth_estimates(fit_gamm, smooth = "s(base_IL_17a):group_gamm")
  f_1 <- sm1 %>%
    add_confint() %>%
    ggplot(aes(y = .estimate*scale_info$sdY, x = base_IL_17a*scale_info$sdX["base_IL_17a"]+scale_info$meanX["base_IL_17a"])) +
    geom_ribbon(aes(ymin = .lower_ci*scale_info$sdY, ymax = .upper_ci*scale_info$sdY),
                alpha = 0.2, fill = "grey") +
    geom_line(colour = "black", size = 1.5) +
    labs(y = y_label,
         x = 'IL-17A (pg/mL)')+
    theme_light()
  
  sm2 <- smooth_estimates(fit_gamm, smooth = "s(base_IL_8):group_gamm")
  f_2 <- sm2 %>%
    add_confint() %>%
    ggplot(aes(y = .estimate*scale_info$sdY, x = base_IL_8*scale_info$sdX["base_IL_8"]+scale_info$meanX["base_IL_8"])) +
    geom_ribbon(aes(ymin = .lower_ci*scale_info$sdY, ymax = .upper_ci*scale_info$sdY),
                alpha = 0.2, fill = "grey") +
    geom_line(colour = "black", size = 1.5) +
    labs(y = y_label,
         x = 'IL-8 (10³ pg/mL)')+
    theme_light()
  
  
  list(
    model_std         = model_std,          # lmer fit on standardized scale
    coef_table_orig   = coef_table_orig,    # rescaled to original (/1000) Y units
    scale_info        = scale_info,         # means/SDs used for back-transforms
    kept_covariates   = kept_cov,
    formula_final_std = form_final,
    f_1 = f_1,
    f_2 = f_2
  )
}

## B_frontal GM: "base_IL_8"   "base_IL_17a"
fit <- ALRD(
  subregion = "B_frontal",
  response  = "GM",
  DATA = DATA,
  B = 2,
  lambda = c(seq(1,100,1)),
  lambda.min_aug = 89
)

# Coefficients and SEs on the original (/1000) scale:
fit$coef_table_orig
coef_table_orig_fmt <- fit$coef_table_orig
coef_table_orig_fmt[,2:5] <- lapply(coef_table_orig_fmt[,2:5],
                                    function(x) formatC(x, format = "f", digits = 4))

coef_table_orig_fmt

# visualization 
f_gm_1 <- fit$f_1
f_gm_2 <- fit$f_2

## B_prefrontal GM: "base_IL_8"   "base_IL_17a" "base_MIP_1a"
fit <- ALRD(
  subregion = 'B_prefrontal',
  response  = "GM",
  DATA = DATA,
  B = 2,
  lambda = c(seq(1,100,1)),
  lambda.min_aug = 91
)

# Coefficients and SEs on the original (/1000) scale:
fit$coef_table_orig
coef_table_orig_fmt <- fit$coef_table_orig
coef_table_orig_fmt[,2:5] <- lapply(coef_table_orig_fmt[,2:5],
                                    function(x) formatC(x, format = "f", digits = 4))

coef_table_orig_fmt

# visualization 
pf_gm_1 <- fit$f_1
pf_gm_2 <- fit$f_2

## B_frontal WM: "Dam_GD40_weight__kg_" "base_G_CSF" "base_IL_8" "base_IL_17a"    
fit <- ALRD(
  subregion = 'B_frontal',
  response  = "WM",
  DATA = DATA,
  B = 2,
  lambda = c(seq(1,100,1)),
  lambda.min_aug = 20
)

# Coefficients and SEs on the original (/1000) scale:
fit$coef_table_orig
coef_table_orig_fmt <- fit$coef_table_orig
coef_table_orig_fmt[,2:5] <- lapply(coef_table_orig_fmt[,2:5],
                                    function(x) formatC(x, format = "f", digits = 4))

coef_table_orig_fmt

# visualization 
f_wm_1 <- fit$f_1
f_wm_2 <- fit$f_2

## B_prefrontal WM: "Dam_Age_at_conception__yr_" "Dam_GD40_weight__kg_" "base_G_CSF" "base_IL_8" "base_IL_17a" 
fit <- ALRD(
  subregion = 'B_prefrontal',
  response  = "WM",
  DATA = DATA,
  B = 2,
  lambda = c(seq(1,100,1)),
  lambda.min_aug = 17.3
)

# Coefficients and SEs on the original (/1000) scale:
fit$coef_table_orig
coef_table_orig_fmt <- fit$coef_table_orig
coef_table_orig_fmt[,2:5] <- lapply(coef_table_orig_fmt[,2:5],
                                    function(x) formatC(x, format = "f", digits = 4))

coef_table_orig_fmt

# visualization 
pf_wm_1 <- fit$f_1 
pf_wm_2 <- fit$f_2

# Helper function for consistent formatting
format_panel <- function(p, xlab_text) {
  p +
    labs(y = "Estimated Effect on ITE", x = xlab_text, title = NULL) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 13),
      axis.text  = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    )
}

# X labels
x_il17 <- "IL-17A (pg/mL)"
x_il8  <- expression(paste("IL-8 (10"^3, " pg/mL)"))

# Apply formatting
f_gm_1  <- format_panel(f_gm_1,  x_il17)
pf_gm_1 <- format_panel(pf_gm_1, x_il17)
f_wm_1  <- format_panel(f_wm_1,  x_il17)
pf_wm_1 <- format_panel(pf_wm_1, x_il17)

f_gm_2  <- format_panel(f_gm_2,  x_il8)
pf_gm_2 <- format_panel(pf_gm_2, x_il8)
f_wm_2  <- format_panel(f_wm_2,  x_il8)
pf_wm_2 <- format_panel(pf_wm_2, x_il8)

f_gm_1 <- f_gm_1 + scale_x_continuous(limits = c(0,7))
pf_gm_1 <- pf_gm_1 + scale_x_continuous(limits = c(0,7))
f_wm_1 <- f_wm_1 + scale_x_continuous(limits = c(0,7))
pf_wm_1 <- pf_wm_1 + scale_x_continuous(limits = c(0,7))
f_gm_2 <- f_gm_2 + scale_x_continuous(limits = c(0,7))
pf_gm_2 <- pf_gm_2 + scale_x_continuous(limits = c(0,7))
f_wm_2 <- f_wm_2 + scale_x_continuous(limits = c(0,7))
pf_wm_2 <- pf_wm_2 + scale_x_continuous(limits = c(0,7))

# --- Build the 2x4 panel grid (collect legend at bottom) ---
panel_grid <- ggarrange(
  f_gm_1,  pf_gm_1,  f_wm_1,  pf_wm_1,
  f_gm_2,  pf_gm_2,  f_wm_2,  pf_wm_2,
  nrow = 2, ncol = 4,
  common.legend = TRUE, legend = "bottom"
)

# --- Create a top row with ONE title per column (bold, centered) ---
title_row <- ggarrange(
  as_ggplot(textGrob("Frontal GM (cm\u00B3)",    gp = gpar(fontface = "bold", cex = 1.1))),
  as_ggplot(textGrob("Prefrontal GM (cm\u00B3)", gp = gpar(fontface = "bold", cex = 1.1))),
  as_ggplot(textGrob("Frontal WM (cm\u00B3)",    gp = gpar(fontface = "bold", cex = 1.1))),
  as_ggplot(textGrob("Prefrontal WM (cm\u00B3)", gp = gpar(fontface = "bold", cex = 1.1))),
  ncol = 4
)

# --- Combine titles + grid ---
final_fig <- ggarrange(
  title_row,
  panel_grid,
  ncol = 1,
  heights = c(0.08, 1)  # adjust spacing for your layout
)

print(final_fig)