DEFUSE_apply = function(Combined_X, Combined_Y, source_index, miss_source, imput_method, seednum, A) {
  ## Combined_X = Combined_X
  ## Combined_Y = Combined_Y
  ## source_index = source_index
  ## miss_source = miss_source
  ## seednum = 1
  ## A = 1:5
  set.seed(seednum)
  
  ## ---- Basic Setup ----
  K  <- length(source_index)                 # number of sites
  N  <- sapply(source_index, length)         # sample size per site
  n  <- N                                    # labeled data counts
  dim <- ncol(Combined_X)                    # number of predictors
  M   <- miss_source                         # missing pattern list
  
  X_complete <- as.matrix(Combined_X[source_index[[1]], A], nrow = n[1])
  Y_complete <- Combined_Y[source_index[[1]]]
  gamma_tilt <- glm(
    Y_complete ~ X_complete - 1,
    family = binomial()
  )$coefficients
  J_tilt_inv_lc <- solve(H(Combined_X[source_index[[1]], A], gamma_tilt, rep(1, n[1])))
  f1_tilt <- X_complete[, A] * as.numeric(Y_complete - g_x(X_complete[, A], matrix(gamma_tilt)))
  S_DR_tilt <- f1_tilt %*% J_tilt_inv_lc
  if (imput_method == "constant") {
    res = main_func2(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum, A)
  } else {
    res = main_func(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum, A)
  }
  
  return(res)
}
