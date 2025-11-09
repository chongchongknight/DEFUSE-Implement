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

## load data
psychiatric_labeled_freq <- read_csv("data/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("data/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq <- read_csv("data/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("data/heart.unlabeled.non.freq.csv") 
most_lab = read_csv("data/most_lab.csv")[, -1]

psychiatric_labeled <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq) 
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq,psychiatric_unlabeled_nonfreq)
psychiatric_labeled = psychiatric_labeled[, -c(11:12, 14, 16)]
psychiatric_unlabeled = psychiatric_unlabeled[, -c(12:13, 15:16)]
psychiatric_labeled = left_join(psychiatric_labeled, most_lab, by = "subject_id") %>% 
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, advanced_heart_disease) %>% 
  mutate(weight = NA)

psychiatric_unlabeled = left_join(psychiatric_unlabeled, most_lab, by = "subject_id") %>% 
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, weight) %>% 
  mutate(advanced_heart_disease = NA) %>% 
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, advanced_heart_disease, weight)

sum_dat = rbind(psychiatric_labeled, psychiatric_unlabeled) %>% 
  mutate(recent_age = as.numeric(recent_age)) %>% 
  mutate(l51252 = as.numeric(scale(l51252)), 
         l50820 = as.numeric(scale(l50820)),
         recent_age = as.numeric(scale(recent_age)))
psychiatric_labeled = sum_dat %>% 
  filter(is.na(weight))
psychiatric_unlabeled = sum_dat %>% 
  filter(!is.na(weight))
psychiatric_labeled_freq <- read_csv("data/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("data/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq <- read_csv("data/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("data/heart.unlabeled.non.freq.csv") 
psychiatric_labeled <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq) 
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq,psychiatric_unlabeled_nonfreq)
## load auxiliary
auxi = read.csv("./data/main_phenotype_all_each_icd_cpt.csv")[, c(3, 12:959)]
psychiatric_labeled = psychiatric_labeled %>% 
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease) %>% 
  mutate(weight = NA)
psychiatric_unlabeled = psychiatric_unlabeled %>% 
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, weight) %>% 
  mutate(advanced_heart_disease = NA) %>% 
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease, weight)
sum_dat = rbind(psychiatric_labeled, psychiatric_unlabeled) %>% 
  mutate(recent_age = as.numeric(recent_age)) %>% 
  mutate(l50993 = as.numeric(scale(l50993)), 
         l51233 = as.numeric(scale(l51233)),
         recent_age = as.numeric(scale(recent_age)))
psychiatric_labeled = sum_dat %>% 
  filter(is.na(weight))
psychiatric_unlabeled = sum_dat %>% 
  filter(!is.na(weight))
psychiatric_labeled = left_join(psychiatric_labeled, auxi, by = "subject_id")
psychiatric_unlabeled = left_join(psychiatric_unlabeled, auxi, by = "subject_id")
## data modefication to our required shape
X_lic1 <- psychiatric_labeled %>% 
  dplyr::select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>% 
  filter(!is.na(l50993), is.na(l51233)) 
X_lic2 <- psychiatric_labeled %>% 
  dplyr::select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>% 
  filter(is.na(l50993), !is.na(l51233)) 
X_lic3 <- psychiatric_labeled %>% 
  dplyr::select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>% 
  filter(is.na(l50993), is.na(l51233)) 
X_lic = rbind(X_lic1, X_lic2, X_lic3)
X_lic$HDL<- ifelse(X_lic$HDL> 0, 1, 0)
X_lic$diabetes<- ifelse(X_lic$diabetes> 0, 1, 0)
X_lc <- psychiatric_labeled %>% 
  dplyr::select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_labeled)) %>% 
  filter(!is.na(l50993) & !is.na(l51233)) 
X_lc$HDL<- ifelse(X_lc$HDL> 0, 1, 0)
X_lc$diabetes<- ifelse(X_lc$diabetes> 0, 1, 0)
X_uc <- psychiatric_unlabeled %>% 
  dplyr::select(recent_age, HDL, diabetes, l50993, l51233, 9:ncol(psychiatric_unlabeled)) %>% 
  filter(!is.na(l50993), !is.na(l51233)) 
X_uc$HDL<- ifelse(X_uc$HDL> 0, 1, 0)
X_uc$diabetes<- ifelse(X_uc$diabetes> 0, 1, 0)
weight = pull(psychiatric_unlabeled %>% filter(!is.na(l50993), !is.na(l51233)),weight)
weight = weight/mean(weight)
Y_lc <- pull(psychiatric_labeled %>% filter(!is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic1 <- pull(psychiatric_labeled %>% filter(!is.na(l50993), is.na(l51233)), advanced_heart_disease)
Y_lic2 <- pull(psychiatric_labeled %>% filter(is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic3 <- pull(psychiatric_labeled %>% filter(is.na(l50993), is.na(l51233)), advanced_heart_disease)
## pre-screening
pre_lm1 = lm(l50993 ~ .-1, data = X_uc %>% select(-l51233))
summary_model1 <- summary(pre_lm1)
p_values1 <- summary_model1$coefficients[, "Pr(>|t|)"]
significant_covariates1 <- names(p_values1)[p_values1 < 0.01]
significant_covariates1 = significant_covariates1[-1]
pre_lm2 = lm(l51233 ~ .-1, data = X_uc %>% select(-l50993))
summary_model2 <- summary(pre_lm2)
p_values2 <- summary_model2$coefficients[, "Pr(>|t|)"]
significant_covariates2 <- names(p_values2)[p_values2 < 0.01]
significant_covariates2 = significant_covariates2[-1]
significant_covariates = union(significant_covariates1, significant_covariates2)
select_covariates = c(c("recent_age", "HDL", "diabetes", "l50993", "l51233"), significant_covariates)
## final data
X_lc = X_lc %>% 
  select(select_covariates)
X_lic1 = X_lic1 %>% 
  select(select_covariates)
X_lic2 = X_lic2 %>% 
  select(select_covariates)
X_lic3 = X_lic3 %>% 
  select(select_covariates)
X_uc = X_uc %>% 
  select(select_covariates)
for (taskid in 1:1) {
  n_lc  = nrow(X_lc)
  n_lic1 = nrow(X_lic1)
  n_lic2 = nrow(X_lic2)
  n_lic3 = nrow(X_lic3)
  n_uc  = nrow(X_uc)
  source_index = list(
    c(1:n_lc), 
    c((n_lc + 1):(n_lc + n_lic1)),
    c((n_lc + n_lic1 + 1): (n_lc + n_lic2 + n_lic2)),
    c((n_lc + n_lic1 + n_lic2 + 1): (n_lc + n_lic1 + n_lic2 + n_lic3)),
    c((n_lc + n_lic1 + n_lic2 + n_lic3 + 1): (n_lc + n_lic1 + n_lic2 + n_lic3 + n_uc))
  )
  miss_source = list(
    NULL,
    5,
    4,
    c(4, 5),
    NULL
  )
  nrep = 1
  for (i in 1:nrep) {
    seednum = i + (taskid - 1) * nrep
    tryCatch({
      Combined_X = as.matrix(rbind(X_lc, X_lic1, X_lic2, X_lic3, X_uc))
      Combined_Y = c(Y_lc, Y_lic1, Y_lic2, Y_lic3, NA)
      res = violation(Combined_X, Combined_Y, source_index, seednum) ## MAR case
      res2 = violation2(Combined_X, Combined_Y, source_index, seednum) ## MCAR case
    })
  }
}
