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

## ---- Load Data (Heart) ----
psychiatric_labeled_freq <- read_csv("data/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("data/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq <- read_csv("data/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("data/heart.unlabeled.non.freq.csv")

psychiatric_labeled <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq)
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq, psychiatric_unlabeled_nonfreq)
most_lab <- read_csv("data/most_lab.csv")[, -1]

psychiatric_labeled <- psychiatric_labeled[, -c(11:12, 14, 16)]
psychiatric_unlabeled <- psychiatric_unlabeled[, -c(12:13, 15:16)]

psychiatric_labeled <- left_join(psychiatric_labeled, most_lab, by = "subject_id") %>%
  mutate(
    recent_age = as.numeric(recent_age),
    na = l51279 + l50820 + l50971 + l51252,
    l51279 = ifelse(!is.na(na), l51279, NA),
    l50820 = ifelse(!is.na(na), l50820, NA),
    l50971 = ifelse(!is.na(na), l50971, NA),
    l51252 = ifelse(!is.na(na), l51252, NA)
  )

psychiatric_unlabeled <- left_join(psychiatric_unlabeled, most_lab, by = "subject_id") %>%
  mutate(
    recent_age = as.numeric(recent_age),
    na = l51279 + l50820 + l50971 + l51252,
    l51279 = ifelse(!is.na(na), l51279, NA),
    l50820 = ifelse(!is.na(na), l50820, NA),
    l50971 = ifelse(!is.na(na), l50971, NA),
    l51252 = ifelse(!is.na(na), l51252, NA)
  )

## ---- Load Lab Scores ----
main_labeled_freq <- read_csv("data/advanced_heart_disease_l50993_data.labeled.freq.csv") %>%
  select(subject_id, adheart_count)
main_unlabeled_freq <- read_csv("data/advanced_heart_disease_l50993_data.unlabeled.freq.csv") %>%
  select(subject_id, adheart_count)
main_labeled_nonfreq <- read_csv("data/advanced_heart_disease_l50993_data.labeled.non.freq.csv") %>%
  select(subject_id, adheart_count)
main_unlabeled_nonfreq <- read_csv("data/advanced_heart_disease_l50993_data.unlabeled.non.freq.csv") %>%
  select(subject_id, adheart_count)

main_labeled <- rbind(main_labeled_freq, main_labeled_nonfreq)
main_unlabeled <- rbind(main_unlabeled_freq, main_unlabeled_nonfreq)

psychiatric_labeled <- left_join(psychiatric_labeled, main_labeled, by = "subject_id")
psychiatric_unlabeled <- left_join(psychiatric_unlabeled, main_unlabeled, by = "subject_id")

## ---- Score Model ----
sub1 <- psychiatric_labeled %>%
  filter(!is.na(na)) %>%
  select(subject_id, l51252, l51279, l50820, l50971, adheart_count)
sub2 <- psychiatric_unlabeled %>%
  filter(!is.na(na)) %>%
  select(subject_id, l51252, l51279, l50820, l50971, adheart_count)

data <- rbind(sub1, sub2)
glmfit <- glm(adheart_count ~ l51252 + l51279 + l50820 + l50971, data = data, family = "binomial")
p <- glmfit$fitted.values
data2 <- data %>% mutate(score = p) %>% select(subject_id, score)

psychiatric_labeled <- left_join(psychiatric_labeled, data2, by = "subject_id")
psychiatric_unlabeled <- left_join(psychiatric_unlabeled, data2, by = "subject_id")

## ---- Basic Variables ----
dat1 <- psychiatric_labeled %>% select(subject_id, recent_age, HDL, diabetes, na, score)
dat2 <- psychiatric_unlabeled %>% select(subject_id, recent_age, HDL, diabetes, na, score)

sum <- rbind(dat1, dat2) %>%
  select(subject_id, recent_age, HDL, diabetes) %>%
  mutate(recent_age = as.numeric(scale(recent_age)))

dat_label1 <- left_join(psychiatric_labeled %>% select(subject_id, advanced_heart_disease), sum)
dat_unlabel1 <- left_join(psychiatric_unlabeled %>% select(subject_id, weight), sum)

sum1 <- rbind(dat1, dat2) %>%
  select(subject_id, na, score) %>%
  filter(!is.na(na)) %>%
  mutate(score = as.numeric(scale(score)))

## ---- Combine Auxillary Data ----
dat_label <- left_join(dat_label1, sum1)
dat_unlabel <- left_join(dat_unlabel1, sum1)

auxi <- read.csv("data/main_phenotype_all_each_icd_cpt.csv")[, c(3, 12:959)]
dat_label <- left_join(dat_label, auxi, by = "subject_id")
dat_unlabel <- left_join(dat_unlabel, auxi, by = "subject_id")

## ---- Construct Feature Matrices ----
X_lic <- dat_label %>%
  filter(is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_label))
X_lc <- dat_label %>%
  filter(!is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_label))
X_uc <- dat_unlabel %>%
  filter(!is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_unlabel))

X_lic$HDL <- ifelse(X_lic$HDL > 0, 1, 0)
X_lic$diabetes <- ifelse(X_lic$diabetes > 0, 1, 0)
X_lc$HDL <- ifelse(X_lc$HDL > 0, 1, 0)
X_lc$diabetes <- ifelse(X_lc$diabetes > 0, 1, 0)
X_uc$HDL <- ifelse(X_uc$HDL > 0, 1, 0)
X_uc$diabetes <- ifelse(X_uc$diabetes > 0, 1, 0)

weight <- pull(dat_unlabel %>% filter(!is.na(na)), weight)
weight <- weight / mean(weight)

Y_lc <- pull(dat_label %>% filter(!is.na(na)), advanced_heart_disease)
Y_lic <- pull(dat_label %>% filter(is.na(na)), advanced_heart_disease)

## ---- Pre-Screen Covariates ----
pre_lm <- lm(score ~ . - 1, data = X_uc)
summary_model <- summary(pre_lm)
p_values <- summary_model$coefficients[, "Pr(>|t|)"]
significant_covariates <- names(p_values)[p_values < 0.01]
significant_covariates <- significant_covariates[-(1:2)]
select_covariates <- c("recent_age", "HDL", "diabetes", "score", significant_covariates)

X_lc <- X_lc %>%
  select(all_of(select_covariates)) %>%
  mutate(intercept = 1) %>%
  select(intercept, everything())
X_lic <- X_lic %>%
  select(all_of(select_covariates)) %>%
  mutate(intercept = 1) %>%
  select(intercept, everything())
X_uc <- X_uc %>%
  select(all_of(select_covariates)) %>%
  mutate(intercept = 1) %>%
  select(intercept, everything())

## ---- Run Violation Analysis ----
for (taskid in 1:1) {
  n_lc <- nrow(X_lc)
  n_lic <- nrow(X_lic)
  n_uc <- nrow(X_uc)
  
  source_index <- list(
    c(1:n_lc),
    c((n_lc + 1):(n_lc + n_lic)),
    c((n_lc + n_lic + 1):(n_lc + n_lic + n_uc))
  )
  
  nrep <- 1
  for (i in 1:nrep) {
    seednum <- i + (taskid - 1) * nrep
    tryCatch({
      Combined_X <- as.matrix(rbind(X_lc, X_lic, X_uc))
      Combined_Y <- c(Y_lc, Y_lic, NA)
      res <- violation(Combined_X, Combined_Y, source_index, seednum)   # MAR case
      res2 <- violation2(Combined_X, Combined_Y, source_index, seednum) # MCAR case
    })
  }
}
