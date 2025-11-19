# Data Generation: Y_it = t*\beta_t + \beta*X + INTERACTIONS + \epsilon_it + alpha_i
# X = (1, X_1,...,X_30)
# beta = (beta_0, beta_1, ... , beta_30)
# INTERACTIONS: gamma*(X*X_T)
# gamma = (\gamma_T, \gamma_1,...,\gamma_30)
# Time: t = 1, 2, ..., K;
# i = 1,2,...,n
# n: number of individuals
# K: number of repeated measurements of each individual
# X_T: treatment ~ Bernoulli(0.5)

# p_con: the number of continuous covariates
# p_bin: the number of binary covariates

n <- 100
K <- 5
p_x <- 0.5
rho_x <- 0.6
sigma_e <- sqrt(0.5)
sigma_a <- sqrt(1.5)
p_con <- 15
p_bin <- 15

beta <- c(0.3, c(0.5,0.4,0.6,-0.3,-0.6,0.3,0.1,-0.2,-0.1,0.2)*4,rep(0,p_con+p_con-10)) # big main effect
beta_t <- 0.2
beta_tg <- 2
gamma <- c(1.4,c(c(0,0,0,0,0),0.4,-0.6,0.8,-0.8,0)*10,rep(0,p_con+p_con-10))

# lambda: possible values of lambda in lasso
# N: number of simulation
lambda <- c(seq(1,500,1),600,800)
N <- 100

# ICC_type = c('b','s')
# main_effect = c('big','small')
# data_type = c('lin','non')
ICC_type = 'b'
main_effect = 'big'
data_type = 'lin'
cor_type = 'cor'

##############################################################################
##############################################################################
###################           Simulation Study             ###################
##############################################################################
##############################################################################
SIM <- function(ICC_type, main_effect, data_type,
                n, K, beta_t, beta, gamma, p_x, rho_x, sigma_e, sigma_a,
                lambda,beta_tg,
                N,p_con,p_bin){
  
  ## Loading Packages 
  library(MASS)
  library(nlme)
  library(glmmLasso)
  library(ggplot2)
  library(caret)
  library(mgcv)
  library(dplyr)
  library(Matrix)
  library(glmnet)
  library(lme4)
  
  ##############################################################################
  ###################          DATA Simulation              ####################
  ##############################################################################
  
  #### Pre-defined function
  # AR(1) Covariance Structure
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    return(rho^exponent)
  }
  
  # Average Prediction Error
  ape <- function(estimate,true){
    pe <- c()
    for (i in 1:dim(true)[1]) {
      pe[i] <- norm(as.matrix(estimate[i,]-true[i,]),type="2")
    }
    mean_pe <- mean(pe)
    return(mean_pe)
  }
  
  # return correlation over time
  cor_time <- function(score_T,true_score_T,K){
    result <- c()
    for (i in 1:K){
      result[i] <- cor(true_score_T[,i],score_T[,i],method = 'spearman')
    }
    return(result)
  }
  
  DATA_GEN <- function(n=200, K=5, beta_t, beta, beta_tg, gamma, p_x=0.5, rho_x = 0, sigma_e, sigma_a, p_con, p_bin){
    alpha_i <- rnorm(n,0,sigma_a)
    time <- 1:K
    DATA <- matrix(0, nrow = 1, ncol = (4+p_con+p_bin))
    for (i in 1:n) {
      epsilon_it <- rnorm(K,0,sigma_e)
      # binary covariates
      X_bin <- rbinom(p_bin,1,p_x) 
      # continuous covariates
      X_con <- mvrnorm(n = 1, mu = rep(0,p_con), Sigma = ar1_cor(p_con, rho_x)) 
      # treatment group
      X_T <- sample(c(1,-1),size =1)
      # baseline covariates for i-th individual
      X <- c(1, X_bin[1:5],X_con[1:5],X_bin[6:p_bin],X_con[6:p_con])
      X_inter <- X*X_T
      if (data_type == 'non'){
        Y <- rep(beta%*%X +beta[2:5]%*%X[c(7,8,9,10)]^2+((X_inter%*%gamma-1*X_T*(X[c(7)])^2-1*X_T*(X[c(8)])^2))/2,K)+beta_t*time+epsilon_it+alpha_i[i]+beta_tg*X_T*time/2
      }
      if (data_type == 'lin'){
        Y <- rep(beta%*%X + ((X_inter%*%gamma-1*X_T*(X[c(7)])^2-1*X_T*(X[c(8)])^2))/2,K)+beta_t*time+epsilon_it+alpha_i[i]+beta_tg*X_T*time/2
      }
      DATA_i <- cbind(Y,X_T,i,time,matrix(X[-1],nrow=K,ncol=p_con+p_bin,byrow = TRUE))
      DATA <- rbind(DATA,DATA_i)
    }
    DATA <- DATA[-1,]
    colnames(DATA)[c(1:4)]<-c('Y','group','ID','time')
    X_name <- c()
    for (j in 1:(p_con+p_bin)){
      X_name[j] <- paste("X",as.character(j),sep = '')
    }
    colnames(DATA)[5:(p_con+p_bin+4)]<-X_name
    DATA <- as.data.frame(DATA)
    return(DATA)
  }
  
  ##############################################################################
  ###################            Best Lambda                ####################
  ##############################################################################
  
  ### Create formula for model
  data <- DATA_GEN(n, K, beta_t, beta, beta_tg, gamma, p_x, rho_x, sigma_e, sigma_a, p_con, p_bin)
  X.name0 <- colnames(data)[-c(1:4)] # select all covariates
  X.name = paste(X.name0,collapse = "+")
  AX.name = paste(paste("I(group/2)",X.name0,sep=":"),collapse = "+") # interaction between treatment and covariates
  
  form <- paste("Y~I(group/2)+time+time:I(group/2)+",X.name,"+",AX.name,sep="") ## full linear regression
  form_int <- paste("N_Y~-1+I(group/2)+time:I(group/2)+",AX.name,sep="") ## augmentation-type new method
  
  ### Data Augmentation
  data$ID <- as.factor(data$ID)
  data_1 <- data[data$group==1,]
  cv_model_1 <- cv.glmnet(x=as.matrix(cbind(1,data_1[,4:(p_con+p_bin+4)])), y=data_1$Y, alpha = 1)
  model_1 <- glmnet(x=as.matrix(cbind(1,data_1[,4:(p_con+p_bin+4)])), y=data_1$Y, alpha = 1, lambda = cv_model_1$lambda.min)
  Y1 <- 0.5*(as.matrix(cbind(1,data[,4:(p_con+p_bin+4)]))%*%as.matrix(model_1$beta,ncol=1))
  
  data_0 <- data[data$group==-1,]
  cv_model_0 <- cv.glmnet(x=as.matrix(cbind(1,data_0[,4:(p_con+p_bin+4)])), y=data_0$Y, alpha = 1)
  model_0 <- glmnet(x=as.matrix(cbind(1,data_0[,4:(p_con+p_bin+4)])), y=data_0$Y, alpha = 1, lambda = cv_model_0$lambda.min)
  Y0 <- 0.5*(as.matrix(cbind(1,data[,4:(p_con+p_bin+4)]))%*%as.matrix(model_0$beta,ncol=1))
  
  data <- cbind(data,Y1,Y0)
  colnames(data)[(p_con+p_bin+4+1):(p_con+p_bin+4+2)] <- c('Y1','Y0')
  
  data_int <- data %>%
    group_by(time) %>%
    mutate(C_Y = Y-mean(Y), N_Y = Y-mean(Y1+Y0))
  data_int <- as.data.frame(data_int)
  
  ### Creat 5-folds
  B = 5 # 5-folds
  flds <- createFolds(unique(data$ID), k = B, list = TRUE, returnTrain = FALSE)
  
  MSE <- matrix(NA,length(lambda),B)
  MSE_int <- matrix(NA,length(lambda),B)
  
  for(i in 1:length(lambda)){
    for(b in 1:B)
    {
      sel_ind <- which(data$ID%in%flds[[b]])
      data_in <- data_int[sel_ind,]
      data_in$ID <- droplevels(data_in$ID) 
      data_out <- data_int[-sel_ind,]
      data_out$ID <- droplevels(data_out$ID) 
      try({
        fit_all <- glmmLasso(as.formula(form),rnd = list(ID=~1),lambda=lambda[i],data=data_out)
        fit_int <- glmmLasso(as.formula(form_int),rnd = list(ID=~1),lambda=lambda[i],data=data_out,control = list(center=FALSE))
        MSE[i,b] = mean((predict(fit_all,newdata=data_in)-data_in$Y)^2)
        MSE_int[i,b] = mean((predict(fit_int,newdata=data_in)-data_in$N_Y)^2)
      },silent=TRUE)
    }
  }
  MSE1 <- apply(MSE,1,function(x){mean(x,na.rm=T)})
  MSE1[apply(is.na(MSE),1,mean)>0.8] <- NA	# for those with too many failed runs, does not include the MSE
  MSE_int1=apply(MSE_int,1,function(x){mean(x,na.rm=T)})
  MSE_int1[apply(is.na(MSE_int),1,mean)>0.8]=NA
  
  lambda.min = min(lambda[MSE1==min(MSE1,na.rm=T)],na.rm=T)
  lambda.min_int = min(lambda[MSE_int1==min(MSE_int1,na.rm=T)],na.rm=T)
  
  ##############################################################################
  ###################          MODEL COMPARISON             ####################
  ##############################################################################
  
  ### True Test Group
  data_test <- DATA_GEN(n=1000, K, beta_t, beta, beta_tg, gamma, p_x, rho_x, sigma_e, sigma_a, p_con, p_bin)
  x_base_test <- data_test %>%
    group_by(ID) %>%
    filter(row_number()==1)
  time <- 1:K
  base_group <- as.numeric(x_base_test$group)
  true_score <- apply((as.matrix(cbind(1,x_base_test[,c(5:(p_con+p_bin+4))]))%*%as.matrix(gamma,1,(p_con+p_bin+1))-1*(x_base_test[,c(10)])^2-1*(x_base_test[,c(11)])^2),1,mean)+beta_tg*mean(time)
  true_score_T <- matrix(apply((as.matrix(cbind(1,x_base_test[,c(5:(p_con+p_bin+4))]))%*%as.matrix(gamma,1,(p_con+p_bin+1))-1*(x_base_test[,c(10)])^2-1*(x_base_test[,c(11)])^2),1,mean),ncol=K,nrow=1000,byrow = F) + 
    beta_tg*matrix(c(1:K),ncol=K,nrow=1000,byrow = T)
  true_group <- sign(true_score>0)
  true_group_T <- sign(true_score_T>0)
  bs_true_T <- colMeans(x_base_test$Y*ifelse(sign(true_score_T>0)==true_group_T,1,-1))
  
  acc_group_all <- c()
  acc_group_int <- c()
  acc_group_gamm <- c()
  acc_true <- c()
  
  ape_group_int <- c()
  ape_group_all <- c()
  ape_group_gamm <- c()
  
  scc_group_all <- c()
  scc_group_int <- c()
  scc_group_gamm <- c()
  
  acc_group_all_T <- matrix(0, nrow = 1, ncol = K)
  acc_group_int_T <- matrix(0, nrow = 1, ncol = K)
  acc_group_gamm_T <- matrix(0, nrow = 1, ncol = K)
  
  scc_group_all_T <- matrix(0, nrow = 1, ncol = K)
  scc_group_int_T <- matrix(0, nrow = 1, ncol = K)
  scc_group_gamm_T <- matrix(0, nrow = 1, ncol = K)
  
  bs_all_T <- matrix(0, nrow = 1, ncol = K)
  bs_int_T <- matrix(0, nrow = 1, ncol = K)
  bs_gamm_T <- matrix(0, nrow = 1, ncol = K)
  
  for (i in 1:N) {
    data <- DATA_GEN(n, K, beta_t, beta, beta_tg, gamma, p_x, rho_x, sigma_e, sigma_a, p_con, p_bin)
    ## create design matrix for each formula
    data$ID <- as.factor(data$ID)
    # long
    x_base <- data %>%
      group_by(ID) %>%
      filter(row_number()==1)
    x.list <-  rep(list(as.matrix(x_base[,5:(4+p_con+p_bin)])), K)
    y.list <- split(data[,1],as.factor(data$time))
    idx.list <- rep(list(1:n), K)
    trt <- x_base$group
    # others
    data$ID <- as.factor(data$ID)
    data_1 <- data[data$group==1,]
    cv_model_1 <- cv.glmnet(x=as.matrix(cbind(1,data_1[,4:(4+p_con+p_bin)])), y=data_1$Y, alpha = 1)
    model_1 <- glmnet(x=as.matrix(cbind(1,data_1[,4:(4+p_con+p_bin)])), y=data_1$Y, alpha = 1, lambda = cv_model_1$lambda.min)
    Y1 <- (as.matrix(cbind(1,data[,4:(4+p_con+p_bin)]))%*%as.matrix(model_1$beta,ncol=1))
    
    data_0 <- data[data$group==-1,]
    cv_model_0 <- cv.glmnet(x=as.matrix(cbind(1,data_0[,4:(4+p_con+p_bin)])), y=data_0$Y, alpha = 1)
    model_0 <- glmnet(x=as.matrix(cbind(1,data_0[,4:(4+p_con+p_bin)])), y=data_0$Y, alpha = 1, lambda = cv_model_0$lambda.min)
    Y0 <- (as.matrix(cbind(1,data[,4:(4+p_con+p_bin)]))%*%as.matrix(model_0$beta,ncol=1))
    data <- cbind(data,Y1,Y0)
    colnames(data)[(4+p_con+p_bin+1):(4+p_con+p_bin+2)] <- c('Y1','Y0')
    
    data_int <- data %>%
      group_by(time) %>%
      mutate(C_Y = Y-mean(Y), N_Y = Y-mean(Y1+Y0))
    data_int <- as.data.frame(data_int)
    
    
    ## int
    fit_int_best <- glmmLasso(as.formula(form_int),rnd = list(ID=~1),lambda=lambda.min_int,data=data_int,control = list(center=FALSE))
    
    ## all
    fit_best <- glmmLasso(as.formula(form),rnd = list(ID=~1),lambda=lambda.min,data=data_int)
    
    ## GAMM
    # create formula used by gamm
    #var_id <- c(6:10,21:30)
    #sel_var <- var_id[fit_int_best$coefficients[c(8:12,23:32)]!=0]
    var_id <- c(6:9)
    sel_var <- var_id[fit_int_best$coefficients[c(8:11)]!=0]
    
    if(length(sel_var)==0){formu="N_Y~-1+I(group_gamm/2)+time:I(group_gamm/2)+s(ID_gamm, bs = 're')"
    }else if (length(sel_var)==1){
      formu=paste("N_Y~-1+I(group_gamm/2)+time:I(group_gamm/2)+s(ID_gamm, bs = 're')+","s(X",sel_var[1],",by=I(group_gamm/2))",sep="")
    }else if(length(sel_var)>1){
      formu=paste("N_Y~-1+I(group_gamm/2)+time:I(group_gamm/2)+s(ID_gamm, bs = 're')+","s(X",sel_var[1],",by=I(group_gamm/2))",paste("+s(X",sel_var[-1],",by=I(group_gamm/2))",sep="",collapse=""),sep="")
    }
    formu<- paste(paste("N_Y~-1+I(group_gamm/2)",sep=""),'+time:I(group_gamm/2)+I(group_gamm/2):X6+I(group_gamm/2):X7+I(group_gamm/2):X8+I(group_gamm/2):X9',"+s(ID_gamm, bs = 're')+s(X6,by=I(group_gamm/2))+s(X7,by=I(group_gamm/2))",sep="")
    ID_gamm <- rep(c(1:n),each=K)
    group_gamm <- data_int$group
    fit_gamm <- gam(formula=as.formula(formu), data = cbind(data_int,ID_gamm,group_gamm))
    
    
    ## Scores
    
    score_all_T <- matrix(as.matrix(x_base_test[,5:(4+p_con+p_bin)])%*%fit_best$coefficients[(5+p_con+p_bin):(4+p_con+p_bin+p_con+p_bin)]+fit_best$coefficients[2],ncol=K,nrow=1000,byrow = F) +
      fit_best$coefficients[(4+p_con+p_bin)]*matrix(c(1:K),ncol=K, nrow=1000, byrow = T)
    score_all <- apply(score_all_T,1,mean)
    
    score_int_T <- matrix(as.matrix(x_base_test[,5:(4+p_con+p_bin)])%*%fit_int_best$coefficients[3:(2+p_con+p_bin)]+fit_int_best$coefficients[1],ncol=K,nrow=1000,byrow = F) +
      fit_int_best$coefficients[2]*matrix(c(1:K), ncol=K, nrow=1000, byrow = T)
    score_int <- apply(score_int_T,1,mean)
    
    # create test data for gamm
    data_test_gamm <- cbind(data_test,ID_gamm = 0,group_gamm=1)
    score_gamm_T <- 2*matrix(as.numeric(predict(fit_gamm,newdata=data_test_gamm)),ncol=K,nrow=1000,byrow=TRUE)
    score_gamm <- apply(score_gamm_T,1,mean)
    
    #accuracy
    acc_true[i] <- mean(true_group)
    acc_group_all[i] <- mean(sign(score_all>0)==true_group) 
    acc_group_int[i] <- mean(sign(score_int>0)==true_group) 
    acc_group_gamm[i] <- mean(sign(score_gamm>0)==true_group) 
    
    # Spearman cor
    scc_group_all[i] <- cor(true_score,score_all,method = 'spearman')
    scc_group_int[i] <- cor(true_score,score_int,method = 'spearman')
    scc_group_gamm[i] <- cor(true_score,score_gamm,method = 'spearman')
    
    # Average Prediction Error
    ape_group_int[i] <- ape(score_int_T,true_score_T)
    ape_group_all[i] <- ape(score_all_T,true_score_T)
    ape_group_gamm[i] <- ape(score_gamm_T,true_score_T)
    
  }
  
  result <- rbind(cbind(acc = acc_group_all,scc = scc_group_all,ape = ape_group_all,model_type='all'),
                  cbind(acc = acc_group_int,scc = scc_group_int,ape = ape_group_int,model_type='int'),
                  cbind(acc = acc_group_gamm,scc = scc_group_gamm,ape = ape_group_gamm,type='gamm')
  )
  result <- as.data.frame(cbind(result,data_type = data_type,ICC_type=ICC_type,main_effect=main_effect))
  
  ##############################################################################
  #####################          SAVE RESULTS             ######################
  ##############################################################################
  # result_save
  true_result <- data.frame(lambda.min = lambda.min,lambda.min_int=lambda.min_int, mean_true_group = mean(true_group))
  
  result_list <- list(result,true_result)
  return(result_list)
}
result_list <- SIM(ICC_type, main_effect, data_type,n, K, beta_t, beta, gamma, p_x, rho_x, sigma_e, sigma_a,lambda,beta_tg,N,p_con,p_bin)

result_save <- paste('/home/zkquan/icc_',ICC_type,'_',main_effect,'_',data_type,'_',cor_type,'.csv',sep='')
true_save <- paste('/home/zkquan/icc_',ICC_type,'_',main_effect,'_',data_type,'_',cor_type,'_true_result.csv',sep='')
write.csv(as.data.frame(result_list[1]),result_save, row.names = FALSE)
write.csv(as.data.frame(result_list[2]),true_save, row.names = FALSE)