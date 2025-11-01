### implementation code for DEFUSE
DEFUSE_apply = function(Combined_X, Combined_Y, source_index, miss_source, imput_method, seednum) {
  ## initial setup
  set.seed(seednum)
  K  = length(source_index)      # number of sites
  N  = sapply(source_index, length)
  n  = N                         # labeled data in each site
  dim = ncol(Combined_X)         # number of predictors
  M   = miss_source
  X_complete = as.matrix(Combined_X[source_index[[1]],], nrow = n[1])
  Y_complete = Combined_Y[source_index[[1]]]
  fold_id = sample(rep(0:2, length.out = nrow(X_complete)))
  model_list_LC = vector("list", 2)
  f1_DR_check_cond_LC = matrix(0, nrow = n[1], ncol = dim)
  f1_DR_check_cond_UC = matrix(0, nrow = n[3], ncol = dim)
  j1_DR_check_cond_LC = array(0, dim = c(n[1], dim, dim))
  j1_DR_check_cond_UC = array(0, dim = c(n[3], dim, dim))
  weight_lc = rep(0, nrow(X_complete))
  gamma_check = 0
  ## cross fitting starts, construct density ratio function, conditional score model and conditional hessian model
  for (ff in 1:3) {
    density_id = ff %% 3
    score_id   = (ff + 1) %% 3
    hessian_id = (ff + 2) %% 3
    weight_lc_1 = density_ratio_xgboost_linear(
      Combined_X[source_index[[3]],][, 1:3],                          # UC X
      Combined_X[source_index[[1]],][fold_id == density_id, 1:3],     # LC train (density_id)
      Combined_X[source_index[[1]],][fold_id == score_id,    1:3]     
    )
    weight_lc[fold_id == score_id] = weight_lc_1
    mod_1 = glm(
      Y_complete[fold_id == score_id] ~ X_complete[fold_id == score_id, ] - 1,
      weights = weight_lc_1,
      family  = gaussian()
    )
    gamma_check = gamma_check + mod_1$coefficients
    f1_DR_check_1 = X_complete * as.numeric(Y_complete - g_x(X_complete, mod_1$coefficients))
    j1_DR_check_1 = H_array(X_complete, mod_1$coefficients)
    X_tr_score <- X_complete[fold_id == score_id, 1:3, drop = FALSE]
    for (col in 1:ncol(X_complete)) {
      Y_tr_score <- f1_DR_check_1[fold_id == score_id, col, drop = FALSE]
      best_score_hp <- tune_gausspr_cv(X_tr_score, Y_tr_score)
      invisible(capture.output({
        fit_col <- gausspr(
          x = X_tr_score,
          y = Y_tr_score,
          kernel = "rbfdot",
          kpar   = list(sigma = best_score_hp$sigma),
          var    = best_score_hp$var
        )
      }))
      f1_DR_check_cond_LC[fold_id == hessian_id, col] <- predict(fit_col, newdata = X_complete[fold_id == hessian_id, 1:3, drop = FALSE])
      f1_DR_check_cond_UC[, col] <- f1_DR_check_cond_UC[, col] + predict(fit_col, newdata = Combined_X[source_index[[3]],][, 1:3, drop = FALSE])
    }
    X_tr_hess <- X_complete[fold_id == hessian_id, 1:3, drop = FALSE]
    for (col in 1:ncol(X_complete)) {
      for (row in col:ncol(X_complete)) {
        Y_tr_hess <- j1_DR_check_1[fold_id == hessian_id, row, col]
        best_hess_hp = tune_gausspr_cv(X_tr_hess, Y_tr_hess)
        invisible(capture.output({
          fit_rc <- gausspr(
            x = X_tr_hess,
            y = Y_tr_hess,
            kernel = "rbfdot",
            kpar = list(sigma = best_hess_hp$sigma),
            var  = best_hess_hp$var
          )
        }))
        pred_test <- predict(fit_rc, newdata = X_complete[fold_id == density_id, 1:3, drop = FALSE])
        j1_DR_check_cond_LC[fold_id == density_id, row, col] <- pred_test
        j1_DR_check_cond_LC[fold_id == density_id, col, row] <- pred_test 
        pred_uc <- predict(fit_rc, newdata = Combined_X[source_index[[3]],][, 1:3, drop = FALSE])
        j1_DR_check_cond_UC[, row, col] <- j1_DR_check_cond_UC[, row, col] + pred_uc
        j1_DR_check_cond_UC[, col, row] <- j1_DR_check_cond_UC[, row, col]
      }
    }
  }
  gamma_check <- gamma_check / 3
  f1_DR_check_cond_UC <- f1_DR_check_cond_UC / 3
  j1_DR_check_cond_UC <- j1_DR_check_cond_UC / 3
  ## One step Augmentation
  J_check_inv_lc <- solve(J(X_complete, gamma_check, weight_lc, j1_DR_check_cond_LC, j1_DR_check_cond_UC))
  f1_check  <- X_complete * as.numeric(Y_complete - g_x(X_complete, gamma_check))
  gamma_tilt <- gamma_check + J_check_inv_lc %*% (colMeans(weight_lc * (f1_check - f1_DR_check_cond_LC)) + colMeans(f1_DR_check_cond_UC))
  J_tilt_inv_lc <- solve(H(Combined_X[source_index[[3]],], gamma_tilt, rep(1, n[3])))
  f1_tilt <- X_complete * as.numeric(Y_complete - g_x(X_complete, matrix(gamma_tilt)))
  S_DR_tilt <- f1_tilt %*% J_tilt_inv_lc
  if (imput_method == "reduced") {
    res = main_func2(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum)
  } else {
    res = main_func(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum)
  }
  res = append(res, list(gamma_check = as.numeric(gamma_check)))
  return(res)
}


