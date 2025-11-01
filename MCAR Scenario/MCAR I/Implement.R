## use the simulation folder
source("./Source/computing_function.R")
source("./Source/main_function_single_continuous.R")
source("./Source/DEFUSE.R")
source("./Source/data_generation.R")

library(matrixStats)      # efficient computation of row/column statistics for matrices
library(stepPlr)          # stepwise penalized logistic regression
library(evd)              # extreme value distributions (for tail analysis, max-stability)
library(methods)          # formal S4 methods and classes (base R dependency)
library(MASS)             # classical statistics functions and multivariate distributions
library(glmnet)           # penalized generalized linear models (LASSO, ridge)
library(dplyr)            # data manipulation and transformation (select, mutate, filter, etc.)
library(caret)            # implementation of machine learning methods (model training, tuning)
library(mice)             # multiple imputation for missing data
library(splines)          # spline basis generation for nonlinear regression
library(ggplot2)          # visualization and plotting
library(glmnet)           # (duplicate load) penalized regression models — LASSO / Elastic Net
library(nloptr)           # nonlinear optimization (for M-estimation and constrained minimization)
library(mvtnorm)          # multivariate normal and t-distribution functions (simulation)


set.seed(100)
sample = c(3000, 500, 2500)
K = 2 # number of sites
N = c(rep(sample[1],K)) # sample size in each site
n = c(sample[2], sample[3])# labeled data in each site
dim = 5 # number of predictors
M = vector("list", length = K)
M[[1]] = NA # missing index
M[[2]] = c(4,5)
Mu = c(0,0,0,0,0) # multivariate normal distribution parameters
rho = 0
cor_1 = matrix(rho, nrow = dim, ncol = dim)
var = 1
cor_1[dim,] = 0
cor_1[,dim] = 0
for (i in 1:dim) {
  cor_1[i,i] = var
}
Sigma = cor_1
# Set the value for gamma.true
gamma.bar = c(-1, -1, 1, 1, 1) # adjust c(-1, -1, 1, 1* w, 1* w) to have plot S1, S2
y_form = 1
imput = 1
## imput = 1： miss X is linear correlated with x obs, 
## y_form = 1： model coefficient is correctly specified
# Simulation result
n_coe = c(4, 5)
n_lc = 500 # nrow(data_lc)
n_lic = 1500 # nrow(data_lic), change this from 1000 to 2500 to have Table 1
n_uc = 2500 # nrow(data_uc)
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
nrep = 1
for (i in 1:nrep) {
  set.seed(i)
  dat = mix_gaussian(n_site = K, n_total = N, n_label = n, mu = Mu, 
                     sigma = Sigma, imput_linear = imput, y_correct = y_form, 
                     M = M, gamma.bar = gamma.bar, seednum = i)
  data_lc = cbind(dat$Y[[1]][1:500], dat$X[[1]][1:500, ]) ### first column is Y
  data_lic = cbind(dat$Y[[2]][1:1500], dat$X[[2]][1:1500, ])
  data_lic[, n_coe + 1] = NA ## missing
  data_uc = cbind(dat$Y[[1]][501:3000], dat$X[[1]][501:3000, ])
  data_lc = data.frame(data_lc)
  data_lic = data.frame(data_lic)
  data_uc = data.frame(data_uc)
  tryCatch({
    Combined_X = rbind(data_lc[, -1], data_lic[, -1] , data_uc[, -1])
    Combined_Y = c(data_lc[, 1], data_lic[, 1], data_uc[, 1])
    lm_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "lm", i)
    svm_res = DEFUSE_apply(Combined_X, Combined_Y, source_index, miss_source, "svmLinear", i)
  })
}
