### link function for bianry response
g_x = function(x, beta){
  ratio = rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    xbeta = matrix(x[i, ], nrow = 1) %*% beta
    ratio[i] = exp(xbeta) / (exp(xbeta) + 1)
  }
  return(ratio)
}
### derivative of link function for binary response
g_dot <- function(x,beta){
  ratio = rep(0, nrow(x))
  for (i in 1:nrow(x)) {
    xbeta = matrix(x[i, ], nrow = 1) %*% beta
    ratio[i] = exp(xbeta) / (exp(xbeta) + 1) ^ 2
  }
  return(ratio)
}
### Hessian matrix
H <- function(x, beta, weights){
  h = matrix(0, nrow = ncol(x), ncol = ncol(x))
  deriv = g_dot(x,beta)
  for (i in 1:nrow(x)) {
    h = h + weights[i] * deriv[i] * (x[i,] %*% t(x[i,]))
  }
  return(h/nrow(x))
}
### Hessian array
H_array <- function(x, beta){
  h_array = array(0, dim = c(nrow(x), ncol(x), ncol(x)))
  deriv = g_dot(x, beta)
  for (i in 1:nrow(x)) {
    h_array[i,, ] = deriv[i] * x[i, ] %*% t(x[i, ])
  }
  return(h_array)
}
### inverse of Hessian, also see our definition in main paper Section 2.1
J <- function(X_complete, gamma_tilt, weight_lc, j1_DR_tilt_cond_LC, j1_DR_tilt_cond_UC) {
  j = matrix(0, nrow = ncol(X_complete), ncol = ncol(X_complete))
  deriv = g_dot(X_complete, gamma_tilt)
  for (i in 1:nrow(X_complete)) {
    j = j + weight_lc[i] * (deriv[i] * X_complete[i,] %*% t(X_complete[i,]))
  }
  j2 = j / nrow(X_complete) - apply(weight_lc * j1_DR_tilt_cond_LC, c(2,3), mean) + apply(j1_DR_tilt_cond_UC, c(2,3), mean)
  return(j2)
}
### Density Ratio for simple logistic 
density_ratio = function(target, source) {
  combined_Y = c(rep(0, nrow(source)), rep(1, nrow(target)))
  combined_X = rbind(source, target)
  fit = glm(combined_Y ~ combined_X, family = binomial)
  odds = fit$fitted.values/(1 - fit$fitted.values)
  weight = odds[combined_Y == 0] * nrow(source) / nrow(target)
  return(weight)
}
### Density Ratio for simple logistic with evaluation data and removal of extreme values 
density_ratio2 = function(target, source, evaluate) {
  combined_Y = c(rep(0, nrow(source)), rep(1, nrow(target)))
  combined_X = rbind(source, target)
  combined_X = cbind(combined_X)
  fit = glm(combined_Y ~ combined_X, family = binomial)
  X_test <- cbind(1, evaluate)
  pred_prob <- 1 / (1 + exp(- X_test %*% fit$coefficients))
  odds = pred_prob  / (1 - pred_prob)
  w = as.numeric(odds * nrow(evaluate) / nrow(target))
  q <- quantile(w, probs = c(0.05, 0.95), na.rm = TRUE)
  weight <- pmin(pmax(w, q[1]), q[2])
  return(weight)
}
### Selection of tuning parameter for Gaussian process
tune_gausspr_cv <- function(X_tr, y_tr, w = NULL, K = 2,
                            sigma_grid = NULL,
                            var_grid   = c(1e-3, 3e-3, 1e-2, 3e-2, 1e-1),
                            seed = 1) {
  stopifnot(is.matrix(X_tr))
  n <- nrow(X_tr)
  y_vec <- as.numeric(y_tr)
  if (length(y_vec) != n) stop("length(y_tr) must equal nrow(X_tr).")
  if (is.null(w)) w <- rep(1, n)
  if (length(w) != n) stop("length(w) must equal nrow(X_tr).")
  if (is.null(sigma_grid)) {
    D2 <- as.matrix(dist(X_tr))^2
    med2 <- median(D2[upper.tri(D2)])
    if (!is.finite(med2) || med2 <= 0) med2 <- mean(D2[upper.tri(D2)])
    base <- 1 / max(med2, .Machine$double.eps)
    sigma_grid <- c(0.5, 1, 2, 4) * base
  }
  set.seed(seed)
  K <- max(2, min(K, n))
  fold_in <- sample(rep(seq_len(K), length.out = n))
  best <- list(score = Inf, sigma = sigma_grid[1], var = var_grid[1])
  for (sg in sigma_grid) {
    for (vg in var_grid) {
      se_list <- numeric(K)
      for (k in seq_len(K)) {
        idx_tr <- which(fold_in != k); idx_va <- which(fold_in == k)
        fit_ok <- TRUE
        pred <- tryCatch({
          fit <- gausspr(
            x = X_tr[idx_tr, , drop = FALSE],
            y = y_vec[idx_tr],
            kernel = "rbfdot",
            kpar   = list(sigma = sg),
            var    = vg
          )
          stats::predict(fit, newdata = X_tr[idx_va, , drop = FALSE])
        }, error = function(e) { fit_ok <<- FALSE; rep(NA_real_, length(idx_va)) })
        
        if (!fit_ok || anyNA(pred)) {
          se_list[k] <- Inf
        } else {
          se_list[k] <- sum(w[idx_va] * (pred - y_vec[idx_va])^2) / sum(w[idx_va])
        }
      }
      score <- mean(se_list)
      if (score < best$score) best <- list(score = score, sigma = sg, var = vg)
    }
  }
  best[c("sigma", "var")]
}
### Asymptotic variance calculation for the gamma1
var_cal = function(S_DR_tilt, K1_cond_Z_LC,
                   weight_LC, weight_LM1, weight_LM2, weight_LM3, alpha, rho1, rho2, rho3, 
                   K1_cond_XgammaY1_LC, K1_cond_XgammaY2_LC,K1_cond_XgammaY3_LC,
                   K1_cond_XgammaY1_Z_LC, K1_cond_XgammaY2_Z_LC, K1_cond_XgammaY3_Z_LC,
                   K1_cond_XgammaY1_LM1, K1_cond_XgammaY1_Z_LM1, 
                   K1_cond_XgammaY2_LM2, K1_cond_XgammaY2_Z_LM2,
                   K1_cond_XgammaY3_LM3, K1_cond_XgammaY3_Z_LM3,
                   X_LC, X_LM1, X_LM2, X_LM3, covariate_z, nK2) {
  if (nK2 > 1) {
    res = colVars(weight_LC * (S_DR_tilt - K1_cond_Z_LC) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[1:nK2, k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[1:nK2, k])) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[(nK2 + 1):(2 * nK2), k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[(nK2 + 1):(2 * nK2), k]))- 
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[(2 * nK2 + 1):(3 * nK2), k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[(2 * nK2 + 1):(3 * nK2), k]))) +
      1 / rho1 * colVars(
        weight_LM1 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LM1[, k] * cbind(1, X_LM1[, covariate_z]) %*% alpha[1:nK2, k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1[, k] * cbind(1, X_LM1[, covariate_z]) %*% alpha[1:nK2, k]))) + 
      1 / rho2 * colVars(
        weight_LM2 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_LM2[, k] * cbind(1, X_LM2[, covariate_z]) %*% alpha[(nK2 + 1):(2 * nK2), k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_Z_LM2[, k] * cbind(1, X_LM2[, covariate_z]) %*% alpha[(nK2 + 1):(2 * nK2), k]))) +
      1 / rho3 * colVars(
        weight_LM3 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_LM3[, k] * cbind(1, X_LM3[, covariate_z]) %*% alpha[(2 * nK2 + 1):(3 * nK2), k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_Z_LM3[, k] * cbind(1, X_LM3[, covariate_z]) %*% alpha[(2 * nK2 + 1):(3 * nK2), k])))
  } else {
    res = colVars(as.numeric(weight_LC) * (S_DR_tilt - K1_cond_Z_LC) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * alpha[1:nK2, k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC[, k] * alpha[1:nK2, k])) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_LC[, k] * alpha[(nK2 + 1):(2 * nK2), k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_Z_LC[, k] * alpha[(nK2 + 1):(2 * nK2), k]))- 
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_LC[, k] * alpha[(2 * nK2 + 1):(3 * nK2), k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_Z_LC[, k] * alpha[(2 * nK2 + 1):(3 * nK2), k]))) +
      1 / rho1 * colVars(
        weight_LM1 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LM1[, k] * alpha[1:nK2, k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1[, k] * alpha[1:nK2, k]))) + 
      1 / rho2 * colVars(
        weight_LM2 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_LM2[, k] * alpha[(nK2 + 1):(2 * nK2), k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY2_Z_LM2[, k] * alpha[(nK2 + 1):(2 * nK2), k]))) +
      1 / rho3 * colVars(
        weight_LM3 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_LM3[, k] * alpha[(2 * nK2 + 1):(3 * nK2), k]) - 
                        sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY3_Z_LM3[, k] * alpha[(2 * nK2 + 1):(3 * nK2), k])))
  }
  return(res)
}
### Asymptotic variance calculation for the gamma2
var_cal_2 = function(S_DR_tilt, K1_cond_Z_LC, 
                     S1, T1_cond_Z_LC,
                     weight_LC, zeta, rho4,
                     T1_cond_X_LC, T1_cond_X_Z_LC, 
                     T1_cond_X_UC,
                     X_LC, X_UC, covariate_z, nT2) {
  if (nT2 > 1) {
    res = colVars(weight_LC * (S1 - T1_cond_Z_LC) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% zeta[ , k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% zeta[ , k]))
    ) + 
      1 / rho4 * colVars(sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_UC[, k] * cbind(1, X_UC[, covariate_z]) %*% zeta[ , k]))
  } else {
    res = colVars(weight_LC * (S1 - T1_cond_Z_LC) -
                    weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_LC[, k] * zeta[ , k]) - 
                                   sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC[, k] * zeta[ , k]))
    ) + 
      1 / rho4 * colVars(sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_UC[, k] * zeta[ , k]))
  }
  return(res)
}
### Density ration using Xgboost
density_ratio_xgboost <- function(target_X, source_X, apply_X) {
  target_X <- as.matrix(target_X)
  source_X <- as.matrix(source_X)
  apply_X  <- as.matrix(apply_X)
  n_t <- nrow(target_X); n_s <- nrow(source_X)
  pi_t <- n_t / (n_t + n_s); pi_s <- 1 - pi_t
  prior_ratio <- pi_s / pi_t
  train_X <- rbind(target_X, source_X)
  y <- c(rep(1L, n_t), rep(0L, n_s))
  set.seed(2025)
  idx <- sample(seq_along(y))
  n_valid <- max(1L, floor(0.2 * length(y)))
  valid_id <- idx[seq_len(n_valid)]
  train_id <- idx[-seq_len(n_valid)]
  dtrain <- xgb.DMatrix(train_X[train_id, , drop = FALSE], label = y[train_id])
  dvalid <- xgb.DMatrix(train_X[valid_id, , drop = FALSE], label = y[valid_id])
  grid <- expand.grid(
    eta = c(0.1, 0.2),
    max_depth = c(3, 4),
    nrounds = c(800, 1200)
  )
  best_loss <- Inf
  best_model <- NULL
  for (i in seq_len(nrow(grid))) {
    pars <- grid[i, ]
    bst <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "logloss",
        max_depth = pars$max_depth,
        eta = pars$eta,
        subsample = 0.8,
        colsample_bytree = 1.0,
        verbosity = 0
      ),
      data = dtrain,
      nrounds = pars$nrounds,
      watchlist = list(eval = dvalid),
      verbose = 0
    )
    pred_valid <- predict(bst, dvalid)
    loss <- -mean(y[valid_id] * log(pred_valid + 1e-8) +
                    (1 - y[valid_id]) * log(1 - pred_valid + 1e-8))
    if (loss < best_loss) {
      best_loss <- loss
      best_model <- bst
    }
  }
  dapply <- xgb.DMatrix(apply_X)
  p <- predict(best_model, dapply)
  p <- pmax(1e-3, pmin(1 - 1e-3, p))
  w <- prior_ratio * p / (1 - p)
  w <- w / mean(w)
  return(w)
}
### Density ration using Xgboost linear
density_ratio_xgboost_linear <- function(target_X, source_X, apply_X) {
  target_X <- as.matrix(target_X)
  source_X <- as.matrix(source_X)
  apply_X  <- as.matrix(apply_X)
  feat_names <- paste0("x", seq_len(ncol(apply_X)))
  colnames(target_X) <- colnames(source_X) <- colnames(apply_X) <- feat_names
  n_t <- nrow(target_X); n_s <- nrow(source_X)
  pi_t <- n_t / (n_t + n_s); pi_s <- 1 - pi_t
  prior_ratio <- pi_s / pi_t
  train_X <- rbind(target_X, source_X)
  y <- c(rep(1L, n_t), rep(0L, n_s))
  set.seed(2025)
  idx <- sample(seq_along(y))
  n_valid <- max(1L, floor(0.2 * length(y)))
  valid_id <- idx[seq_len(n_valid)]
  train_id <- idx[-seq_len(n_valid)]
  dtrain <- xgb.DMatrix(train_X[train_id, , drop = FALSE], label = y[train_id])
  dvalid <- xgb.DMatrix(train_X[valid_id, , drop = FALSE], label = y[valid_id])
  grid <- expand.grid(
    eta    = c(0.1, 0.3, 0.5),
    lambda = c(0, 1e-3, 1e-2),
    alpha  = c(0, 1e-3, 1e-2),
    nrounds = c(800, 1200)
  )
  best_loss  <- Inf
  best_model <- NULL
  for (i in seq_len(nrow(grid))) {
    pars <- grid[i, ]
    bst <- xgb.train(
      params = list(
        objective = "binary:logistic",
        booster   = "gblinear",
        eta       = pars$eta,
        lambda    = pars$lambda,
        alpha     = pars$alpha,
        eval_metric = "logloss",
        verbosity   = 0
      ),
      data = dtrain,
      nrounds = pars$nrounds,
      watchlist = list(eval = dvalid),
      early_stopping_rounds = 50,
      verbose = 0
    )
    loss <- bst$best_score
    if (!is.null(loss) && loss < best_loss) {
      best_loss  <- loss
      best_model <- bst
    }
  }
  dapply <- xgb.DMatrix(apply_X)
  p <- predict(best_model, dapply)
  #p <- pmax(1e-3, pmin(1 - 1e-3, p))
  w <- prior_ratio * p / (1 - p)
  qs <- quantile(w, c(0.01, 0.99), na.rm = TRUE)
  w  <- pmin(pmax(w, qs[1]), qs[2])
  w <- w / mean(w)
  return(w)
}
