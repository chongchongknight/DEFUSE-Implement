source("./Source/computing_function.R")
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
library(tidyverse)   # data loading and cleaning


# ---- Load Heart Datasets ----
psychiatric_labeled_freq <- read_csv("data/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("data/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq <- read_csv("data/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("data/heart.unlabeled.non.freq.csv")

psychiatric_labeled <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq)
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq, psychiatric_unlabeled_nonfreq)
most_lab <- read_csv("data/most_lab.csv")[, -1]

psychiatric_labeled <- psychiatric_labeled[, -c(11:12, 14, 16)]
psychiatric_unlabeled <- psychiatric_unlabeled[, -c(12:13, 15:16)]

# ---- Merge lab results ----
psychiatric_labeled <- left_join(psychiatric_labeled, most_lab, by = "subject_id") %>%
  mutate(
    recent_age = as.numeric(recent_age),
    na = l51279 + l50820 + l50971 + l51252,
    l51279 = case_when(!is.na(na) ~ l51279, TRUE ~ NA),
    l50820 = case_when(!is.na(na) ~ l50820, TRUE ~ NA),
    l50971 = case_when(!is.na(na) ~ l50971, TRUE ~ NA),
    l51252 = case_when(!is.na(na) ~ l51252, TRUE ~ NA)
  )

psychiatric_unlabeled <- left_join(psychiatric_unlabeled, most_lab, by = "subject_id") %>%
  mutate(
    recent_age = as.numeric(recent_age),
    na = l51279 + l50820 + l50971 + l51252,
    l51279 = case_when(!is.na(na) ~ l51279, TRUE ~ NA),
    l50820 = case_when(!is.na(na) ~ l50820, TRUE ~ NA),
    l50971 = case_when(!is.na(na) ~ l50971, TRUE ~ NA),
    l51252 = case_when(!is.na(na) ~ l51252, TRUE ~ NA)
  )

# ---- Score Construction ----
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

psychiatric_unlabeled <- left_join(psychiatric_unlabeled, main_unlabeled, by = "subject_id")
psychiatric_labeled <- left_join(psychiatric_labeled, main_labeled, by = "subject_id")

# ---- Build GLM score model ----
sub1 <- psychiatric_labeled %>%
  filter(!is.na(na)) %>%
  select(subject_id, l51252, l51279, l50820, l50971, adheart_count)
sub2 <- psychiatric_unlabeled %>%
  filter(!is.na(na)) %>%
  select(subject_id, l51252, l51279, l50820, l50971, adheart_count)

data <- rbind(sub1, sub2)
glmfit <- glm(adheart_count ~ l51252 + l51279 + l50820 + l50971, data = data, family = "binomial")
p <- glmfit$fitted.values
data2 <- data %>%
  mutate(score = p) %>%
  select(subject_id, score)

psychiatric_unlabeled <- left_join(psychiatric_unlabeled, data2, by = "subject_id")
psychiatric_labeled <- left_join(psychiatric_labeled, data2, by = "subject_id")

# ---- Select features ----
dat1 <- psychiatric_labeled %>%
  select(subject_id, recent_age, HDL, diabetes, na, score)
dat2 <- psychiatric_unlabeled %>%
  select(subject_id, recent_age, HDL, diabetes, na, score)

sum <- rbind(dat1, dat2) %>%
  select(subject_id, recent_age, HDL, diabetes) %>%
  mutate(recent_age = as.numeric(scale(recent_age)))

dat_label1 <- left_join(psychiatric_labeled %>% select(subject_id, advanced_heart_disease), sum)
dat_unlabel1 <- left_join(psychiatric_unlabeled %>% select(subject_id, weight), sum)

sum1 <- rbind(dat1, dat2) %>%
  select(subject_id, na, score) %>%
  filter(!is.na(na)) %>%
  mutate(score = as.numeric(scale(score)))

# ---- Combine auxiliary data ----
dat_label <- left_join(dat_label1, sum1)
dat_unlabel <- left_join(dat_unlabel1, sum1)
auxi <- read.csv("data/main_phenotype_all_each_icd_cpt.csv")[, c(3, 12:959)]

dat_label <- left_join(dat_label, auxi, by = "subject_id")
dat_unlabel <- left_join(dat_unlabel, auxi, by = "subject_id")

# ---- Split Labeled/Unlabeled subsets ----
X_lic <- dat_label %>%
  filter(is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_label))
X_lic$HDL <- ifelse(X_lic$HDL > 0, 1, 0)
X_lic$diabetes <- ifelse(X_lic$diabetes > 0, 1, 0)

X_lc <- dat_label %>%
  filter(!is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_label))
X_lc$HDL <- ifelse(X_lc$HDL > 0, 1, 0)
X_lc$diabetes <- ifelse(X_lc$diabetes > 0, 1, 0)

X_uc <- dat_unlabel %>%
  filter(!is.na(na)) %>%
  select(recent_age, HDL, diabetes, score, 8:ncol(dat_unlabel))
X_uc$HDL <- ifelse(X_uc$HDL > 0, 1, 0)
X_uc$diabetes <- ifelse(X_uc$diabetes > 0, 1, 0)

weight <- pull(dat_unlabel %>% filter(!is.na(na)), weight)
weight <- weight / mean(weight)
Y_lc <- pull(dat_label %>% filter(!is.na(na)), advanced_heart_disease)
Y_lic <- pull(dat_label %>% filter(is.na(na)), advanced_heart_disease)

# ---- Pre-screen informative variables ----
pre_lm <- lm(score ~ . - 1, data = X_uc)
summary_model <- summary(pre_lm)
p_values <- summary_model$coefficients[, "Pr(>|t|)"]
significant_covariates <- names(p_values)[p_values < 0.01]
significant_covariates <- significant_covariates[-(1:2)]
select_covariates <- c("recent_age", "HDL", "diabetes", "score", significant_covariates)

set.seed(1)
folds <- sample(cut(seq(1, n_lc), breaks = 5, labels = FALSE))
mod <- vector("list", 5)
ctrl <- trainControl(method = "cv", number = 5)

for (i in 1:5) {
  mod[[i]][[1]] <- train(
    score ~ ., data = X_lc[folds != i, ],
    method = "svmRadial",
    trControl = ctrl,
    tuneGrid = expand.grid(C = 2.297397, sigma = 0.01991501)
  )
  mod[[i]][[2]] <- mod[[i]][[1]]
}

# ---- Step 2: Main loop for selection and estimation ----
n <- 1
for (i in 1:n) {
  tryCatch({
    set.seed(i)
    
    # For no-bootstrap theoretical setting
    sam_c <- 1:n_lc
    sam_ic <- 1:n_lic
    
    # Split data into train/test
    Y_lc_train <- Y_lc[sam_c]
    Y_lic_train <- Y_lic[sam_ic]
    X_lc_train <- X_lc[sam_c, 1:4]
    X_lic_train <- X_lic[sam_ic, 1:4]
    X_uc_train <- X_uc[, 1:4]
    
    X_lc_train_withauxi <- X_lc[sam_c, ]
    X_lic_train_withauxi <- X_lic[sam_ic, ]
    X_uc_train_withauxi <- X_uc
    
    X_lc_train_withauxi_select <- X_lc_train_withauxi %>%
      dplyr::select(select_covariates)
    X_lic_train_withauxi_select <- X_lic_train_withauxi %>%
      dplyr::select(select_covariates)
    X_uc_train_withauxi_select <- X_uc_train_withauxi %>%
      dplyr::select(select_covariates)
    
    # Define test sets (if needed later)
    Y_lc_test <- Y_lc[-sam_c]
    Y_lic_test <- Y_lic[-sam_ic]
    X_lc_test <- X_lc[-sam_c, 1:4]
    X_lic_test <- X_lic[-sam_ic, 1:4]
    
    n_lc_train <- length(sam_c)
    n_lic_train <- length(sam_ic)
    n_uc <- nrow(X_uc)
    M <- 4
    
    # ---- Step 3: Get Estimators ----
    ## (a) Benchmark methods
    mice_res <- mice_get()
    ssl_res <- ssl_get()
    
    ## (b) Proposed methods
    sip_dop_res <- sip_dop_get(seed = i, cali_method = "intercept")
    sip_dop_res2 <- sip_dop_get(seed = i, cali_method = "linear")
    
  }, error = function(e) {
    cat("Error occurred for seed number", i, ":", conditionMessage(e), "\n")
  })
}
