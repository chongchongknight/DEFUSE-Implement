source("./Source/computing_function.R")
source("./Source/main_function_single_continuous.R")
source("./Source/DEFUSE.R")
source("./Source/data_generation.R")
library(matrixStats) # efficient computation of row/column statistics (e.g., column variances)
library(stepPlr)     # penalized logistic regression with stepwise variable selection
library(evd)         # extreme value distributions and related functions
library(methods)     # formal methods and S4 object system (base R)
library(MASS)        # statistical functions (e.g., mvrnorm, lda, glm.nb)
library(dplyr)       # data manipulation and cleaning
library(caret)       # unified interface for machine learning models and model evaluation
library(splines)     # spline basis functions for non-linear modeling
library(ggplot2)     # data visualization and plotting
library(glmnet)      # regularized generalized linear models (lasso, ridge, elastic net)
library(nloptr)      # nonlinear optimization (e.g., for calibration parameter tuning)
library(kernlab)     # kernel-based machine learning methods (e.g., KRR, SVM)
library(tensor)      # tensor algebra and operations for multi-dimensional optimization
library(xgboost)     # gradient boosting and density ratio estimation
## sample size, degress of density shift, nonlinearity
parameter = list(c(500, 1000, 5000, 0.1, 0.5),
                 c(500, 1000, 5000, 0.2, 0.5),
                 c(500, 1000, 5000, 0.3, 0.5),
                 c(500, 1000, 5000, 0.2, 0.25),
                 c(500, 1000, 5000, 0.2, 0.75))
for (taskid in 1:5) {
  n_site = 3 # number of sites
  n_label = parameter[[taskid]][1:3] # labeled data in each site, f
  n_total = parameter[[taskid]][1:3] # sample size in each site
  dim = 5 # dimension of X
  M = vector("list", length = n_site)
  M[[1]] = NA # missing index
  M[[2]] = c(4,5)
  ## first M is label complete, unlabel data is not included 
  mu = c(0,0,0,0,0) # multivariate normal distribution parameters
  rho = 0 # correlation
  cor_1 = matrix(rho, nrow = dim, ncol = dim)
  for (i in 1:dim) {
    cor_1[i,i] = 1
  }
  sigma = cor_1
  gamma.bar = t(matrix(rep(c(-1, -1, 1, 1, 1), 3), ncol = 3))
  y_correct = 2
  imput_linear = 1
  n_coe = c(4, 5)
  n_lc = n_total[1] # nrow(data_lc)
  n_lic = n_total[2] # nrow(data_lic)
  n_uc = n_total[3] # nrow(data_uc)
  source_index = list(
    c(1:n_lc), 
    c((n_lc + 1):(n_lc + n_lic)),
    c((n_lc + n_lic + 1): (n_lc + n_lic + n_uc))
  )
  miss_source = list(
    NULL,
    n_coe,
    NULL
  )
  nrep = 100
  res_table1 = matrix(0, nrow = nrep, ncol = 5)
  res_table2 = matrix(0, nrow = nrep, ncol = 5)
  res_table3 = matrix(0, nrow = nrep, ncol = 5)
  reduces_res_table2 = matrix(0, nrow = nrep, ncol = 5)
  res_check = matrix(0, nrow = nrep, ncol = 5)
  var_table_lm1 = matrix(0, nrow = nrep, ncol = 5)
  var_table_lm2 = matrix(0, nrow = nrep, ncol = 5)
  var_table_uc1 = matrix(0, nrow = nrep, ncol = 5)
  var_table_uc2 = matrix(0, nrow = nrep, ncol = 5)
  reduced_var_table_lm1 = matrix(0, nrow = nrep, ncol = 5)
  reduced_var_table_lm2 = matrix(0, nrow = nrep, ncol = 5)
  ## population level estimate
  n_total2 = c(500, 1000, 300000)
  dat = mix_gaussian(n_site = n_site, n_total = n_total2, n_label = n_total2, mu = mu, 
                     sigma = sigma, imput_linear = imput_linear, y_correct = y_correct, 
                     M = M, gamma.bar = gamma.bar, seednum = 123, 
                     parameter[[taskid]][4], parameter[[taskid]][5])
  real = lm(dat$Y[[3]] ~ dat$X[[3]]- 1)$coefficients
  ## replication
  for (i in 1:nrep) {
    set.seed(i)
    #i=1
    dat = mix_gaussian(n_site = n_site, n_total = n_total, n_label = n_total, mu = mu, 
                       sigma = sigma, imput_linear = imput_linear, y_correct = y_correct, 
                       M = M, gamma.bar = gamma.bar, seednum = i, 
                       parameter[[taskid]][4], parameter[[taskid]][5])
    data_lc = cbind(dat$Y[[1]], dat$X[[1]]) ### first column is Y
    data_lic = cbind(dat$Y[[2]], dat$X[[2]][, 1:3], NA, NA)
    data_uc = cbind(NA, dat$X[[3]])
    data_lc = data.frame(data_lc)
    data_lic = data.frame(data_lic)
    data_uc = data.frame(data_uc)
    tryCatch({
      Combined_X = as.matrix(rbind(data_lc[, -1], data_lic[, -1] , data_uc[, -1]))
      Combined_Y = c(data_lc[, 1], data_lic[, 1], data_uc[, 1])
      start = Sys.time()
      lm_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "lm", i)
      end = Sys.time()
      reduce_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "reduced", i)
      res_table1[i, ] = lm_res$gamma_tilt # LC-only
      res_table2[i, ] = lm_res$gamma_1_ad # DEFUSE_LM CM
      res_table3[i, ] = lm_res$gamma_2_ad # DEFUSE 
      reduces_res_table2[i, ] = reduce_res$gamma_1_ad # DEFUSE_LM RM
      var_table_lm1[i, ] = lm_res$var_tilt # LC-only theoretical variance
      var_table_lm2[i, ] = lm_res$var_ad1 # DEFUSE_LM CM theoretical variance
      var_table_uc1[i, ] = lm_res$T_part_var0 # DEFUSE_LM CM theoretical variance
      var_table_uc2[i, ] = lm_res$T_part_var # DEFUSE_LM theoretical variance
      reduced_var_table_lm1[i, ] = reduce_res$var_tilt # LC-only theoretical variance
      reduced_var_table_lm2[i, ] = reduce_res$var_ad1 # DEFUSE_LM RM theoretical variance
    })
  }
}