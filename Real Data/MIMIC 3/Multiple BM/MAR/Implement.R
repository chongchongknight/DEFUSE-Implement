source("./Source/computing_function.R")
source("./Source/main_function_single_binary.R")
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

## ---- Load Data ----
psychiatric_labeled_freq   <- read_csv("data/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("data/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq   <- read_csv("data/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("data/heart.unlabeled.non.freq.csv")

psychiatric_labeled   <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq)
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq, psychiatric_unlabeled_nonfreq)

auxi <- read.csv("data/main_phenotype_all_each_icd_cpt.csv")[, c(3, 12:959)]

## ---- Select and Merge Basic Variables ----
psychiatric_labeled <- psychiatric_labeled %>%
  select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease) %>%
  mutate(weight = NA)

psychiatric_unlabeled <- psychiatric_unlabeled %>%
  select(subject_id, recent_age, HDL, diabetes, l50993, l51233, weight) %>%
  mutate(advanced_heart_disease = NA) %>%
  select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease, weight)

sum_dat <- rbind(psychiatric_labeled, psychiatric_unlabeled) %>%
  mutate(
    recent_age = as.numeric(scale(as.numeric(recent_age))),
    l50993 = as.numeric(scale(l50993)),
    l51233 = as.numeric(scale(l51233))
  )

psychiatric_labeled   <- sum_dat %>% filter(is.na(weight))
psychiatric_unlabeled <- sum_dat %>% filter(!is.na(weight))

psychiatric_labeled   <- left_join(psychiatric_labeled, auxi, by = "subject_id")
psychiatric_unlabeled <- left_join(psychiatric_unlabeled, auxi, by = "subject_id")

## ---- Construct Labeled/Unlabeled Blocks ----
X_lic1 <- psychiatric_labeled %>%
  select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>%
  filter(!is.na(l50993), is.na(l51233))
X_lic2 <- psychiatric_labeled %>%
  select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>%
  filter(is.na(l50993), !is.na(l51233))
X_lic3 <- psychiatric_labeled %>%
  select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>%
  filter(is.na(l50993), is.na(l51233))
X_lc <- psychiatric_labeled %>%
  select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>%
  filter(!is.na(l50993) & !is.na(l51233))
X_uc <- psychiatric_unlabeled %>%
  select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_unlabeled)) %>%
  filter(!is.na(l50993), !is.na(l51233))

## Binary transform
for (df in list(X_lic1, X_lic2, X_lic3, X_lc, X_uc)) {
  df$HDL <- ifelse(df$HDL > 0, 1, 0)
  df$diabetes <- ifelse(df$diabetes > 0, 1, 0)
}

## Outcome & Weights
weight <- pull(psychiatric_unlabeled %>% filter(!is.na(l50993), !is.na(l51233)), weight)
weight <- weight / mean(weight)

Y_lc   <- pull(psychiatric_labeled %>% filter(!is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic1 <- pull(psychiatric_labeled %>% filter(!is.na(l50993), is.na(l51233)), advanced_heart_disease)
Y_lic2 <- pull(psychiatric_labeled %>% filter(is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic3 <- pull(psychiatric_labeled %>% filter(is.na(l50993), is.na(l51233)), advanced_heart_disease)

## ---- Pre-Screen Informative Covariates ----
pre_lm1 <- lm(l50993 ~ . - 1, data = X_uc %>% select(-l51233))
p_values1 <- summary(pre_lm1)$coefficients[, "Pr(>|t|)"]
significant_covariates1 <- names(p_values1)[p_values1 < 0.01][-1]

pre_lm2 <- lm(l51233 ~ . - 1, data = X_uc %>% select(-l50993))
p_values2 <- summary(pre_lm2)$coefficients[, "Pr(>|t|)"]
significant_covariates2 <- names(p_values2)[p_values2 < 0.01][-1]

significant_covariates <- union(significant_covariates1, significant_covariates2)
select_covariates <- c("recent_age", "HDL", "diabetes", "l50993", "l51233", significant_covariates)

## Select and add intercept
add_intercept <- function(df) {
  df %>%
    select(all_of(select_covariates)) %>%
    mutate(intercept = 1) %>%
    select(intercept, everything())
}
X_lc   <- add_intercept(X_lc)
X_lic1 <- add_intercept(X_lic1)
X_lic2 <- add_intercept(X_lic2)
X_lic3 <- add_intercept(X_lic3)
X_uc   <- add_intercept(X_uc)

## ---- Run DEFUSE ----
for (taskid in 1:1) {
  n_lc   <- nrow(X_lc)
  n_lic1 <- nrow(X_lic1)
  n_lic2 <- nrow(X_lic2)
  n_lic3 <- nrow(X_lic3)
  n_uc   <- nrow(X_uc)
  
  source_index <- list(
    c(1:n_lc),
    c((n_lc + 1):(n_lc + n_lic1)),
    c((n_lc + n_lic1 + 1):(n_lc + n_lic1 + n_lic2)),
    c((n_lc + n_lic1 + n_lic2 + 1):(n_lc + n_lic1 + n_lic2 + n_lic3)),
    c((n_lc + n_lic1 + n_lic2 + n_lic3 + 1):(n_lc + n_lic1 + n_lic2 + n_lic3 + n_uc))
  )
  
  miss_source <- list(NULL, 6, 5, c(5, 6), NULL)
  nrep <- 1
  A <- 1:6
  
  ## Storage matrices
  res_table1 <- res_table2 <- res_table3 <- matrix(0, nrep, 6)
  reduces_res_table2 <- reduces_res_table3 <- matrix(0, nrep, 6)
  var_table_lm1 <- var_table_lm2 <- var_table_uc1 <- var_table_uc2 <- matrix(0, nrep, 6)
  reduced_var_table_lm1 <- reduced_var_table_lm2 <- reduced_var_table_uc1 <- reduced_var_table_uc2 <- matrix(0, nrep, 6)
  
  for (i in 1:nrep) {
    tryCatch({
      Combined_X <- as.matrix(rbind(X_lc, X_lic1, X_lic2, X_lic3, X_uc))
      Combined_Y <- c(Y_lc, Y_lic1, Y_lic2, Y_lic3, NA)
      
      lm_res     <- DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "lm", i, A)
      reduce_res <- DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "constant", i, A)
      
      res_table1[i, ] <- lm_res$gamma_tilt
      res_table2[i, ] <- lm_res$gamma_1_ad
      res_table3[i, ] <- lm_res$gamma_2_ad
      reduces_res_table2[i, ] <- reduce_res$gamma_1_ad
      reduces_res_table3[i, ] <- reduce_res$gamma_2_ad
      
      var_table_lm1[i, ] <- lm_res$var_tilt
      var_table_lm2[i, ] <- lm_res$var_ad1
      var_table_uc1[i, ] <- lm_res$T_part_var0
      var_table_uc2[i, ] <- lm_res$T_part_var
      reduced_var_table_lm1[i, ] <- reduce_res$var_tilt
      reduced_var_table_lm2[i, ] <- reduce_res$var_ad1
      reduced_var_table_uc1[i, ] <- reduce_res$T_part_var0
      reduced_var_table_uc2[i, ] <- reduce_res$T_part_var
    })
  }
}

