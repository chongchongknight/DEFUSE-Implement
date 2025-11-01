source("./Source/DEFUSE.R")
source("./Source/computing_function.R")
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
library(pROC)        # receiver operating characteristic (ROC) analysis and AUC computation
set.seed(100)
n_site = 6 
n_sub = 1000
n_total = c(n_sub, n_sub, n_sub, n_sub, n_sub, 5000) # sample size in each site
n_label = c(n_sub, n_sub, n_sub, n_sub, n_sub, NA) # labeled data in each site, f
dim = 5 # dimension of X
Miss = c(4, 5)
n_lc = n_total[1] # nrow(data_lc)
n_lic1 = n_total[2] # nrow(data_lic)
n_lic2 = n_total[3] # nrow(data_lic)
n_lic3 = n_total[4] # nrow(data_lic)
n_lic4 = n_total[5] # nrow(data_lic)
n_uc = n_total[6] # nrow(data_uc)
source_index = list(
  c(1:n_lc), 
  c((n_lc + 1):(n_lc + n_lic1)),
  c((n_lc + n_lic1 + 1):(n_lc+ n_lic1 + n_lic2)),
  c((n_lc + n_lic1 + n_lic2 + 1):(n_lc+ n_lic1 + n_lic2 + n_lic3)),
  c((n_lc + n_lic1 + n_lic2 + n_lic3 + 1):(n_lc+ n_lic1 + n_lic2 + n_lic3 + n_lic4)),
  c((n_lc + n_lic1 + n_lic2 + n_lic3 + n_lic4 + 1): (n_lc+ n_lic1 + n_lic2 + n_lic3 + n_lic4 + n_uc))
)
gamma.bar = t(matrix(rep(c(-1, -1, 1, 1, 1), n_site), ncol = n_site))
## coefficient shift, nonlinearity 
parameter_list = list(
  c(0.6, 0.25),
  c(0.6, 0.5),
  c(0.6, 0.75),
  c(0.6, 1),
  c(0.3, 0.5),
  c(0.6, 0.5),
  c(0.9, 0.5),
  c(1.2, 0.5)
)
n_cc = 4
AUC_F1_table = matrix(0, nrow = n_cc, ncol = 8)
for (cc in 1:n_cc) {
  nrep = 100
  res_table1 = matrix(0, nrow = nrep, ncol = 4)
  res_table2 = matrix(0, nrow = nrep, ncol = 4)
  res_table3 = matrix(0, nrow = nrep, ncol = 4)
  res_table4 = matrix(0, nrow = nrep, ncol = 4)
  AUC_F1_subtable = matrix(0, nrow = nrep, ncol = 8)
  para1 = parameter_list[[cc]][1]
  para2 = parameter_list[[cc]][2]
  for (i in 1:nrep) {
    seednum = i + (cc - 1) * nrep
    dat = mix_gaussian(n_site = n_site, n_total = n_total, dim = dim, Miss = Miss, 
                       gamma.bar = gamma.bar, seednum = seednum, shift = para1, nonlinear = para2)
    data_lc   = data.frame(cbind(dat$Y[[1]], dat$X[[1]]))
    data_lic1 = data.frame(cbind(dat$Y[[2]], dat$X[[2]][, 1:3], NA, dat$X[[2]][, 5]))
    data_lic2 = data.frame(cbind(dat$Y[[3]], dat$X[[3]][, 1:3], NA, dat$X[[3]][, 5]))
    data_lic3 = data.frame(cbind(dat$Y[[4]], dat$X[[4]][, 1:4], NA))
    data_lic4 = data.frame(cbind(dat$Y[[5]], dat$X[[5]][, 1:4], NA))
    data_uc   = data.frame(cbind(NA,         dat$X[[6]]))
    Combined_X = as.matrix(rbind(data_lc[, -1], data_lic1[, -1], data_lic2[, -1],
                                 data_lic3[, -1], data_lic4[, -1], data_uc[, -1]))
    Combined_Y = c(data_lc[, 1], data_lic1[, 1], data_lic2[, 1], 
                   data_lic3[, 1], data_lic4[, 1], data_uc[, 1])
    start = Sys.time()
    res = violation(Combined_X, Combined_Y, source_index, seednum)
    end = Sys.time()
    res_table1[i, ] = res$keep_id1
    res_table2[i, ] = res$keep_id2
    res_table3[i, ] = res$keep_id3
    res_table4[i, ] = res$keep_id4
    suppressMessages({
      AUC_F1_subtable[i, 1] = auc(roc(c(0, 1, 0, 1), as.numeric(res$keep_id1), levels = c(0, 1), direction = "<"))
      AUC_F1_subtable[i, 2] = auc(roc(c(0, 1, 0, 1), as.numeric(res$keep_id2), levels = c(0, 1), direction = "<"))
      AUC_F1_subtable[i, 3] = auc(roc(c(0, 1, 0, 1), as.numeric(res$keep_id3), levels = c(0, 1), direction = "<"))
      AUC_F1_subtable[i, 4] = auc(roc(c(0, 1, 0, 1), as.numeric(res$keep_id4), levels = c(0, 1), direction = "<"))
    })
    AUC_F1_subtable[i, 5] = f1_score(c(0, 1, 0, 1), as.numeric(res$keep_id1))
    AUC_F1_subtable[i, 6] = f1_score(c(0, 1, 0, 1), as.numeric(res$keep_id2))
    AUC_F1_subtable[i, 7] = f1_score(c(0, 1, 0, 1), as.numeric(res$keep_id3))
    AUC_F1_subtable[i, 8] = f1_score(c(0, 1, 0, 1), as.numeric(res$keep_id4))
  }
}

