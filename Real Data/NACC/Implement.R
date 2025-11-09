source("./Source/computing_function.R")
source("./Source/main_function_multiple_binary.R")
source("./Source/DEFUSE.R")
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
## sample size, degrees of density shift, nonlinearity

#set.seed(1) the result is not based on seed, I forget to do so
data = read.csv("data/score_memory_language_function.csv") %>% 
  select(NACCID, NACCMMSE, SEX, RACE, NACCAGE, WEIGHT, HEIGHT, CDRLANG) %>%
  mutate(NACCMMSE = case_when(
    is.na(NACCMMSE) ~ NA,
    NACCMMSE >= 25 ~ 0,
    NACCMMSE < 25 ~ 1
  )) %>% 
  filter(RACE <= 10 & SEX <= 10 & WEIGHT <= 500 & HEIGHT<= 500) %>% 
  mutate(NACCAGE = scale(NACCAGE),
         WEIGHT = scale(WEIGHT),
         HEIGHT = scale(HEIGHT))
data2 = read.csv("data/added_demo.csv") %>% 
  select(-X)

## selected center ids
id_center = read.csv("id_center.csv")[, -1]
id_center = id_center[!duplicated(id_center$NACCID), ]
adc_select = read.csv("data/adc_select.csv")[, -1]

data = left_join(data, data2, by = "NACCID") %>% select(NACCMMSE, SEX, RACE, NACCAGE, WEIGHT, HEIGHT, NACCLIVS, ALCOHOL, TOBAC100, CDRLANG)
## select data
set.seed(3)
data_lc = data %>% 
  filter(!is.na(NACCMMSE) & !is.na(CDRLANG)) 
data_lc = data_lc[sample(1:nrow(data_lc), 1000, replace = F), ]
data_lic = data %>% 
  filter(!is.na(NACCMMSE) & is.na(CDRLANG)) %>% 
  filter(NACCADC %in% adc_select) %>% 
  select(-NACCADC)
data_uc = data %>% 
  filter(is.na(NACCMMSE) & !is.na(CDRLANG))
data_uc = data_uc[sample(1:nrow(data_uc), 10000, replace = F), ]
colnames = colnames(data_uc)
## implement code
for (taskid in 1:1) {
  n_lc  = nrow(data_lc)
  n_lic = nrow(data_lic)
  n_uc  = nrow(data_uc)
  source_index = list(
    c(1:n_lc), 
    c((n_lc + 1):(n_lc + n_lic)),
    c((n_lc + n_lic + 1): (n_lc + n_lic + n_uc))
  )
  miss_source = list(
    NULL,
    9,
    NULL
  )
  nrep = 1
  res_table1 = matrix(0, nrow = nrep, ncol = 9)
  res_table2 = matrix(0, nrow = nrep, ncol = 9)
  res_table3 = matrix(0, nrow = nrep, ncol = 9)
  constant_res_table2 = matrix(0, nrow = nrep, ncol = 9)
  constant_res_table3 = matrix(0, nrow = nrep, ncol = 9)
  res_check = matrix(0, nrow = nrep, ncol = 9)
  var_table_lm1 = matrix(0, nrow = nrep, ncol = 9)
  var_table_lm2 = matrix(0, nrow = nrep, ncol = 9)
  var_table_uc1 = matrix(0, nrow = nrep, ncol = 9)
  var_table_uc2 = matrix(0, nrow = nrep, ncol = 9)
  constant_var_table_lm1 = matrix(0, nrow = nrep, ncol = 9)
  constant_var_table_lm2 = matrix(0, nrow = nrep, ncol = 9)
  constant_var_table_uc1 = matrix(0, nrow = nrep, ncol = 9)
  constant_var_table_uc2 = matrix(0, nrow = nrep, ncol = 9)
  for (i in 1:nrep) {
    tryCatch({
      Combined_X = as.matrix(rbind(data_lc[, -1], data_lic[, -1] , data_uc[, -1]))
      Combined_Y = c(data_lc[, 1], data_lic[, 1], data_uc[, 1])
      lm_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "lm", i)
      constant_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "constant", i)
      res_table1[i, ] = lm_res$gamma_tilt
      res_table2[i, ] = lm_res$gamma_1_ad
      res_table3[i, ] = lm_res$gamma_2_ad
      constant_res_table2[i, ] = constant_res$gamma_1_ad
      constant_res_table3[i, ] = constant_res$gamma_2_ad
      var_table_lm1[i, ] = lm_res$var_tilt
      var_table_lm2[i, ] = lm_res$var_ad1
      var_table_uc1[i, ] = lm_res$T_part_var0
      var_table_uc2[i, ] = lm_res$T_part_var
      reduced_var_table_lm1[i, ] = reduce_res$var_tilt
      reduced_var_table_lm2[i, ] = reduce_res$var_ad1
      reduced_var_table_uc1[i, ] = reduce_res$T_part_var0
      reduced_var_table_uc2[i, ] = reduce_res$T_part_var
    })
  }
}