## implement the DEFUSE 
DEFUSE_apply = function(Combined_X, Combined_Y, source_index, miss_source, imput_method, seednum) {
  ## Combined_X = Combined_X
  ## Combined_Y = Combined_Y
  ## source_index = source_index
  ## miss_source = miss_source
  ## seednum = 1
  ## imput_method = "svmLinear"
  set.seed(seednum)
  K = length(source_index) # number of sites
  N = rep(0, K) # sample size in each site
  for (i in 1:K) {
    N[i] = length(source_index[[i]])
  } 
  n = N # labeled data in each site
  dim = ncol(Combined_X) # number of predictors
  M = miss_source
  X_complete = as.matrix(Combined_X[source_index[[1]],], nrow = n[1])
  Y_complete = Combined_Y[source_index[[1]]]
  # Step 1: estimate the gamma
  mod = lm(Y_complete ~ X_complete - 1)
  gamma.tilt = mod$coefficients
  naive_sd = mean((mod$residuals) ^ 2)
  # Step 2: Estimation of Score Function
  # get the score function
  H_inv = solve(H(X_complete, gamma.tilt))
  S = t(H_inv %*% t(X_complete)) * as.numeric(Y_complete - g_x(X_complete, gamma.tilt))
  set.seed(100)
  # cross fitting, 5 folds
  folds <- sample(cut(seq(1,n[[1]]),breaks = 5,labels = FALSE))
  # Use the main function to do compute
  res = main_func(Combined_X, Combined_Y, n, K, H_inv, S, gamma.tilt, naive_sd, folds, source_index, miss_source, imput_method, seednum)
  return(res)
}
