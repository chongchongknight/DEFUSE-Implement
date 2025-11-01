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
## number of X_Omegar and number of X_Gammar
parameter = list(c(2000, 2000, 5000, 2, 22),
                 c(2000, 2000, 5000, 6, 22),
                 c(2000, 2000, 5000, 4, 12),
                 c(2000, 2000, 5000, 4, 22),
                 c(2000, 2000, 5000, 4, 32))
set.seed(123)
XOmgea_mimusA_coef = matrix(round(rnorm(18, 0, 0.5), 1), ncol = 6)
XGamma_coef = matrix(round(rnorm(64, 0, 0.5), 1), ncol = 32)
for (taskid in 1:5) {
  n_site = 3 # number of sites
  n_total = parameter[[taskid]][1:3] # sample size in each site
  A = c(1:3)
  XOmgea_mimusA = c(4:(3 + parameter[[taskid]][4]))
  XOmgea_nonsense = c((4 + parameter[[taskid]][4]):(13 + parameter[[taskid]][4]))
  XGamma_use = c((14 + parameter[[taskid]][4]):(13 + parameter[[taskid]][4] + parameter[[taskid]][5]))
  XGamma_nonsense = c((14 + parameter[[taskid]][4] + parameter[[taskid]][5]):(23 + parameter[[taskid]][4] + parameter[[taskid]][5]))
  XOmega = c(1:(13 + parameter[[taskid]][4]))
  XGamma = c(1:(23 + parameter[[taskid]][4] + parameter[[taskid]][5]))
  Xunobs = c((24 + parameter[[taskid]][4] + parameter[[taskid]][5]):(33 + parameter[[taskid]][4] + parameter[[taskid]][5]))
  dim = (33 + parameter[[taskid]][4] + parameter[[taskid]][5])
  gamma.bar = t(matrix(rep(c(-1, -3, -1, 1, 1, -1, 1, -2, 2), 3), ncol = 3))
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
    Xunobs,
    NULL
  )
  nrep = 100
  res_table1 = matrix(0, nrow = nrep, ncol = length(A))
  res_table2 = matrix(0, nrow = nrep, ncol = length(A))
  res_table3 = matrix(0, nrow = nrep, ncol = length(A))
  reduces_res_table2 = matrix(0, nrow = nrep, ncol = length(A))
  res_check = matrix(0, nrow = nrep, ncol = length(A))
  var_table_lm1 = matrix(0, nrow = nrep, ncol = length(A))
  var_table_lm2 = matrix(0, nrow = nrep, ncol = length(A))
  var_table_uc1 = matrix(0, nrow = nrep, ncol = length(A))
  var_table_uc2 = matrix(0, nrow = nrep, ncol = length(A))
  reduced_var_table_lm1 = matrix(0, nrow = nrep, ncol = length(A))
  reduced_var_table_lm2 = matrix(0, nrow = nrep, ncol = length(A))
  ## population level estimate
  n_total2 = c(500, 1000, 300000)
  dat = mix_gaussian(n_site, n_total, A, XOmgea_mimusA, XOmgea_nonsense, XGamma_use,
                     XGamma_use, XGamma_nonsense, Xunobs, 123, 0.4, dim, 
                     XOmgea_mimusA_coef, XGamma_coef)
  real = lm(dat$Y[[3]] ~ dat$X[[3]][, A] - 1)$coefficients
  ## replication
  for (i in 1:nrep) {
    seednum = i + nrep * (taskid - 1)
    dat = mix_gaussian(n_site, n_total, A, XOmgea_mimusA, XOmgea_nonsense, XGamma_use,
                       XGamma_use, XGamma_nonsense, Xunobs, seednum, 0.4, dim,
                       XOmgea_mimusA_coef, XGamma_coef)
    data_lc = cbind(dat$Y[[1]], dat$X[[1]]) ### first column is Y
    data_lic = cbind(dat$Y[[2]], dat$X[[2]][, -Xunobs], NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    data_uc = cbind(NA, dat$X[[3]])
    data_lc = data.frame(data_lc)
    data_lic = data.frame(data_lic)
    data_uc = data.frame(data_uc)
    tryCatch({
      Combined_X = as.matrix(rbind(data_lc[, -1], data_lic[, -1] , data_uc[, -1]))
      Combined_Y = c(data_lc[, 1], data_lic[, 1], data_uc[, 1])
      start = Sys.time()
      lm_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, 
                            A, XOmega, XGamma, "lm", seednum)
      end = Sys.time()
      #print(end - start)
      res_table1[i, ] = lm_res$gamma_tilt # LC-only
      res_table2[i, ] = lm_res$gamma_1_ad # DEFUSE_LM CM
      res_table3[i, ] = lm_res$gamma_2_ad # DEFUSE 
      var_table_lm1[i, ] = lm_res$var_tilt # LC-only theoretical variance
      var_table_lm2[i, ] = lm_res$var_ad1 # DEFUSE_LM CM theoretical variance
      var_table_uc1[i, ] = lm_res$T_part_var0 # DEFUSE_LM CM theoretical variance
      var_table_uc2[i, ] = lm_res$T_part_var # DEFUSE_LM theoretical variance
    })
  }
}
