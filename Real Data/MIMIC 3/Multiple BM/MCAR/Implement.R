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

## load data
psychiatric_labeled_freq   <- read_csv("other_phe/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("other_phe/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq   <- read_csv("other_phe/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("other_phe/heart.unlabeled.non.freq.csv")
most_lab <- read_csv("other_phe/most_lab.csv")[, -1]

psychiatric_labeled   <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq)
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq, psychiatric_unlabeled_nonfreq)

psychiatric_labeled   <- psychiatric_labeled[, -c(11:12, 14, 16)]
psychiatric_unlabeled <- psychiatric_unlabeled[, -c(12:13, 15:16)]

psychiatric_labeled <- left_join(psychiatric_labeled, most_lab, by = "subject_id") %>%
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, advanced_heart_disease) %>%
  mutate(weight = NA)
psychiatric_unlabeled <- left_join(psychiatric_unlabeled, most_lab, by = "subject_id") %>%
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, weight) %>%
  mutate(advanced_heart_disease = NA) %>%
  dplyr::select(recent_age, HDL, HBP, diabetes, l51252, l50820, advanced_heart_disease, weight)

sum_dat <- rbind(psychiatric_labeled, psychiatric_unlabeled) %>%
  mutate(
    recent_age = as.numeric(scale(recent_age)),
    l51252 = as.numeric(scale(l51252)),
    l50820 = as.numeric(scale(l50820))
  )

psychiatric_labeled   <- sum_dat %>% filter(is.na(weight))
psychiatric_unlabeled <- sum_dat %>% filter(!is.na(weight))

psychiatric_labeled_freq   <- read_csv("multi_block/heart.labeled.freq.csv")
psychiatric_unlabeled_freq <- read_csv("multi_block/heart.unlabeled.freq.csv")
psychiatric_labeled_nonfreq   <- read_csv("multi_block/heart.labeled.non.freq.csv")
psychiatric_unlabeled_nonfreq <- read_csv("multi_block/heart.unlabeled.non.freq.csv")

psychiatric_labeled   <- rbind(psychiatric_labeled_freq, psychiatric_labeled_nonfreq)
psychiatric_unlabeled <- rbind(psychiatric_unlabeled_freq, psychiatric_unlabeled_nonfreq)

auxi <- read.csv("./auxillary_main/main_phenotype_all_each_icd_cpt.csv")[, c(3, 12:959)]

psychiatric_labeled <- psychiatric_labeled %>%
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease) %>%
  mutate(weight = NA)
psychiatric_unlabeled <- psychiatric_unlabeled %>%
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, weight) %>%
  mutate(advanced_heart_disease = NA) %>%
  dplyr::select(subject_id, recent_age, HDL, diabetes, l50993, l51233, advanced_heart_disease, weight)

sum_dat <- rbind(psychiatric_labeled, psychiatric_unlabeled) %>%
  mutate(
    recent_age = as.numeric(scale(recent_age)),
    l50993 = as.numeric(scale(l50993)),
    l51233 = as.numeric(scale(l51233))
  )

psychiatric_labeled   <- sum_dat %>% filter(is.na(weight))
psychiatric_unlabeled <- sum_dat %>% filter(!is.na(weight))

psychiatric_labeled   <- left_join(psychiatric_labeled, auxi, by = "subject_id")
psychiatric_unlabeled <- left_join(psychiatric_unlabeled, auxi, by = "subject_id")
## data construction
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

for (dat in list(X_lc, X_lic1, X_lic2, X_lic3, X_uc)) {
  dat$HDL      <- ifelse(dat$HDL > 0, 1, 0)
  dat$diabetes <- ifelse(dat$diabetes > 0, 1, 0)
}

weight <- pull(psychiatric_unlabeled %>% filter(!is.na(l50993), !is.na(l51233)), weight)
weight <- weight / mean(weight)

Y_lc   <- pull(psychiatric_labeled %>% filter(!is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic1 <- pull(psychiatric_labeled %>% filter(!is.na(l50993), is.na(l51233)), advanced_heart_disease)
Y_lic2 <- pull(psychiatric_labeled %>% filter(is.na(l50993), !is.na(l51233)), advanced_heart_disease)
Y_lic3 <- pull(psychiatric_labeled %>% filter(is.na(l50993), is.na(l51233)), advanced_heart_disease)


set.seed(510)
n_lc <- length(Y_lc)
folds <- sample(cut(seq(1, n_lc), breaks = 5, labels = FALSE))
ctrl  <- trainControl(method = "cv", number = 5)

mod4 <- mod5 <- mod45_4 <- mod45_5 <- vector("list", 5)

for (i in 1:5) {
  mod5[[i]][[1]] <- train(l51233 ~ ., data = X_lc[folds != i, ], method = "svmRadial", trControl = ctrl)
  mod4[[i]][[1]] <- train(l50993 ~ ., data = X_lc[folds != i, ], method = "svmRadial", trControl = ctrl)
  mod45_4[[i]][[1]] <- train(l50993 ~ ., data = X_lc[folds != i, -5], method = "svmRadial", trControl = ctrl)
  mod45_5[[i]][[1]] <- train(l51233 ~ ., data = X_lc[folds != i, -4], method = "svmRadial", trControl = ctrl)
  
  mod5[[i]][[2]]   <- mod5[[i]][[1]]
  mod4[[i]][[2]]   <- mod4[[i]][[1]]
  mod45_4[[i]][[2]] <- mod45_4[[i]][[1]]
  mod45_5[[i]][[2]] <- mod45_5[[i]][[1]]
}


n <- 1

gamma_tilt_df   <- gamma_tilt1_df <- gamma_tilt2_df <- matrix(0, nrow = n, ncol = 6)
gamma_tilt1_df2 <- gamma_tilt2_df2 <- gamma_mice_df <- matrix(0, nrow = n, ncol = 6)
se_tilt_df   <- se_tilt1_df <- se_tilt2_df <- matrix(0, nrow = n, ncol = 6)
se_tilt1_df2 <- se_tilt2_df2 <- se_mice_df <- matrix(0, nrow = n, ncol = 6)

nlm <- c(length(Y_lic1), length(Y_lic2), length(Y_lic3))

for (i in 1:n) {
  tryCatch({
    set.seed(i)
    sam_c   <- 1:n_lc
    sam_lm1 <- 1:nlm
    sam_lm2 <- 1:nlm
    sam_lm3 <- 1:nlm
    
    ```
    Y_lc_train  <- Y_lc[sam_c]
    Y_lm1_train <- Y_lic1[sam_lm1]
    Y_lm2_train <- Y_lic2[sam_lm2]
    Y_lm3_train <- Y_lic3[sam_lm3]
    
    X_lc_train  <- X_lc[sam_c, 1:5]
    X_lm1_train <- X_lic1[sam_lm1, 1:5]
    X_lm2_train <- X_lic2[sam_lm2, 1:5]
    X_lm3_train <- X_lic3[sam_lm3, 1:5]
    X_uc_train  <- X_uc[, 1:5]
    
    mice_res     <- mice_get()
    ssl_res      <- ssl_get(length(sam_c), X_lc_train, Y_lc_train, X_uc_train)
    sip_dop_res  <- sip_dop_get(seed = i)
    sip_dop_res2 <- sip_dop_get2(seed = i)
    
    gamma_tilt_df[i, ]   <- sip_dop_res$gamma
    gamma_tilt1_df[i, ]  <- sip_dop_res$gamma1
    gamma_tilt2_df[i, ]  <- sip_dop_res$gamma2
    gamma_tilt1_df2[i, ] <- sip_dop_res2$gamma1
    gamma_tilt2_df2[i, ] <- sip_dop_res2$gamma2
    gamma_mice_df[i, ]   <- mice_res$gamma
    
    se_tilt_df[i, ]   <- sip_dop_res$se
    se_tilt1_df[i, ]  <- sip_dop_res$se1
    se_tilt2_df[i, ]  <- sip_dop_res$se2
    se_tilt1_df2[i, ] <- sip_dop_res2$se1
    se_tilt2_df2[i, ] <- sip_dop_res2$se2
    se_mice_df[i, ]   <- mice_res$se
    ```
    
  }, error = function(e) {
    cat("Error occurred for seed number", i, ":", conditionMessage(e), "\n")
  })
}
