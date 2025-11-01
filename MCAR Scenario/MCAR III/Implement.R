source("./Source/data_generation.R")
source("./Source/computing_function.R")
source("./Source/main_function_multiple_continuous.R")
source("./Source/DEFUSe.R")

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

sample = vector(mode = "list", 4)
coefficient = vector(mode = "list", 10)
for (i in 1:4) {
  sample[[i]] = c(3000, 500, 500 * (i + 1))
}

for (j in 1:10) {
  coefficient[[j]] = c(-1, -1, 1, 0.1 * j, 0.1 * j)
}
seednum = 100
res = double_projection(1, 1, c(sample[[i]]), c(coefficient[[j]]), seednum)
