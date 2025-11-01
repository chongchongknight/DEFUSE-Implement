source("./misalign/computing_function.R")
source("./misalign/violation.R")
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
data = read.csv("data/score_memory_language_function.csv") %>%
  select(NACCID, NACCMMSE, SEX, RACE, NACCAGE, WEIGHT, HEIGHT, CDRLANG) %>%
  mutate(
    NACCMMSE = case_when(
      is.na(NACCMMSE) ~ NA,
      NACCMMSE >= 25 ~ 0,
      NACCMMSE < 25 ~ 1
    )
  ) %>%
  filter(RACE <= 10 & SEX <= 10 & WEIGHT <= 500 & HEIGHT <= 500) %>%
  mutate(
    NACCAGE = scale(NACCAGE),
    WEIGHT = scale(WEIGHT),
    HEIGHT = scale(HEIGHT)
  )
data2 = read.csv("data/added_demo.csv") %>% select(-X)
data = left_join(data, data2, by = "NACCID") %>%
  select(NACCMMSE, SEX, RACE, NACCAGE, WEIGHT, HEIGHT, NACCLIVS, ALCOHOL, TOBAC100, CDRLANG)
set.seed(3)
data_lc = data %>%
  filter(!is.na(NACCMMSE) & !is.na(CDRLANG)) %>%
  slice_sample(n = 1000)
data_lic = data %>%
  filter(!is.na(NACCMMSE) & is.na(CDRLANG)) %>%
  slice_sample(n = 1000)
data_uc = data %>%
  filter(is.na(NACCMMSE) & !is.na(CDRLANG)) %>%
  slice_sample(n = 10000)
colnames = colnames(data_uc)
for (taskid in 1:1) {
  n_lc = nrow(data_lc)
  n_lic = nrow(data_lic)
  n_uc = nrow(data_uc)
  source_index = list(
    1:n_lc,
    (n_lc + 1):(n_lc + n_lic),
    (n_lc + n_lic + 1):(n_lc + n_lic + n_uc)
  )
  miss_source = list(NULL, 9, NULL)
  nrep = 1
  for (i in 1:nrep) {
    seednum = i + (taskid - 1) * nrep
    tryCatch({
      Combined_X = as.matrix(rbind(data_lc[, -1], data_lic[, -1], data_uc[, -1]))
      Combined_Y = c(data_lc[, 1], data_lic[, 1], data_uc[, 1])
      res = violation(Combined_X, Combined_Y, source_index, seednum)
    })
  }
}
