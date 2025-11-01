### calculation of gamma1 and gamma2 for different calibration and conditional models
main_func <- function(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum){
  ## initial setup
  set.seed(seednum)
  X_LC = Combined_X[source_index[[1]],]
  X_LM = Combined_X[source_index[[2]],]
  X_UC = Combined_X[source_index[[3]],]
  Y_LC = Combined_Y[source_index[[1]]]
  Y_LM = Combined_Y[source_index[[2]]]
  Y_UC = Combined_Y[source_index[[3]]]
  fold_id_LC = sample(rep(0:3, length.out = n[1])) 
  fold_id_LM = sample(rep(0:1, length.out = n[2])) 
  nK2 = 1 + ncol(X_LC[, 1:3])
  K1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_Z_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_XgammaY_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY_Z_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_XgammaY_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_Z_LM_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  weight_LC = rep(0, n[1])
  weight_LM = rep(0, n[2])
  ## column-wisely construct conditional models
  for (col in 1:ncol(S_DR_tilt)) {
    K1 = S_DR_tilt[, col]
    for (ff in 1:4) {
      density_id = ff %% 4
      score_z_id = (ff + 1) %% 4
      score_xsuby_id = (ff + 2) %% 4
      score_xsuby_z_id = (ff + 3) %% 4
      weight_lc = density_ratio_xgboost_linear(X_UC[ ,1:3], X_LC[fold_id_LC == density_id, 1:3], X_LC[fold_id_LC == score_xsuby_id, 1:3])
      weight_LC[fold_id_LC == score_xsuby_id] = weight_lc
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, 1:3]
      K1_tr_score_z <- K1[fold_id_LC == score_z_id]
      XY_LC_tr_score_xsuby <- cbind(X_LC[fold_id_LC == score_xsuby_id, 1:3], Y_LC[fold_id_LC == score_xsuby_id])
      K1_tr_score_score_xsuby <- K1[fold_id_LC == score_xsuby_id]
      best_hess_hp_score_z <- tune_gausspr_cv(X_LC_tr_score_z, K1_tr_score_z) #tuning parameter selection
      best_hess_hp_score_xsuby <- tune_gausspr_cv(XY_LC_tr_score_xsuby, K1_tr_score_score_xsuby)
      invisible(capture.output({
        mod_K1_Z_LC = gausspr(X_LC_tr_score_z, 
                              K1_tr_score_z, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_z$sigma), 
                              var = best_hess_hp_score_z$var)
        mod_K1_XgammaY_LC = gausspr(XY_LC_tr_score_xsuby, 
                                    K1_tr_score_score_xsuby, 
                                    kernel = "rbfdot", 
                                    kpar = list(sigma = best_hess_hp_score_xsuby$sigma), 
                                    var = best_hess_hp_score_xsuby$var)
      }))
      K1_cond_XgammaY_LC[fold_id_LC == score_xsuby_z_id, col] = predict(mod_K1_XgammaY_LC, newdata = cbind(X_LC[fold_id_LC == score_xsuby_z_id, 1:3], Y_LC[fold_id_LC == score_xsuby_z_id]))
      X_LC_tr_score_xsuby_z <- X_LC[fold_id_LC == score_xsuby_z_id, 1:3]
      K1_tr_score_xsuby_z <- K1_cond_XgammaY_LC[fold_id_LC == score_xsuby_z_id, col]
      best_hess_hp_score_xsuby_z <- tune_gausspr_cv(X_LC_tr_score_xsuby_z, K1_tr_score_xsuby_z)
      invisible(capture.output({
        mod_K1_XgammaY_Z_LC = gausspr(X_LC_tr_score_xsuby_z, 
                                      K1_tr_score_xsuby_z,
                                      kernel = "rbfdot", 
                                      kpar = list(sigma = best_hess_hp_score_xsuby_z$sigma), 
                                      var = best_hess_hp_score_xsuby_z$var)
      }))
      K1_cond_XgammaY_Z_LC[fold_id_LC == score_z_id, col] = predict(mod_K1_XgammaY_Z_LC, newdata = X_LC[fold_id_LC == score_z_id, 1:3])
      K1_cond_Z_LC[fold_id_LC == density_id, col] = predict(mod_K1_Z_LC, newdata = X_LC[fold_id_LC == density_id, 1:3])
      K1_cond_Z_UC[, col] = predict(mod_K1_Z_LC, newdata = X_UC[,1:3]) + K1_cond_Z_UC[, col]
      K1_cond_XgammaY_LM[, col] = predict(mod_K1_XgammaY_LC, newdata = cbind(X_LM[,1:3], Y_LM)) + K1_cond_XgammaY_LM[, col]
      K1_cond_XgammaY_Z_LC_UC[, col] = predict(mod_K1_XgammaY_Z_LC, newdata = X_UC[,1:3]) + K1_cond_XgammaY_Z_LC_UC[, col]
    }
    K1_cond_Z_UC[, col] = K1_cond_Z_UC[, col] / 4
    K1_cond_XgammaY_LM[, col] = K1_cond_XgammaY_LM[, col] / 4
    K1_cond_XgammaY_Z_LC_UC[, col] = K1_cond_XgammaY_Z_LC_UC[, col] / 4
    ## LM source training
    for (ff in 1:2) {
      density_id = ff %% 2
      score_xsuby_z_id = (ff + 1) %% 2
      weight_lm = density_ratio_xgboost_linear(X_UC[ ,1:3], X_LM[fold_id_LM == density_id, 1:3], X_LM[fold_id_LM == score_xsuby_z_id, 1:3])
      weight_LM[fold_id_LM == score_xsuby_z_id] = weight_lm
      X_LM_tr_score_xsuby_z <- X_LM[fold_id_LM == score_xsuby_z_id, 1:3]
      K1_tr_score_xsuby_z <- K1_cond_XgammaY_LM[fold_id_LM == score_xsuby_z_id, col]
      best_hess_hp_score_xsuby_z <- tune_gausspr_cv(X_LM_tr_score_xsuby_z, K1_tr_score_xsuby_z)
      invisible(capture.output({
        mod_K1_XgammaY_Z_LM = gausspr(X_LM_tr_score_xsuby_z, 
                                      K1_tr_score_xsuby_z, 
                                      kernel = "rbfdot", 
                                      kpar = list(sigma = best_hess_hp_score_xsuby_z$sigma), 
                                      var = best_hess_hp_score_xsuby_z$var)
      }))
      K1_cond_XgammaY_Z_LM[fold_id_LM == density_id, col] = predict(mod_K1_XgammaY_Z_LM, newdata = X_LM[fold_id_LM == density_id, 1:3])
      K1_cond_XgammaY_Z_LM_UC[, col] = predict(mod_K1_XgammaY_Z_LM, newdata = X_UC[,1:3]) + K1_cond_XgammaY_Z_LM_UC[, col]
    }
    K1_cond_XgammaY_Z_LM_UC[, col] = K1_cond_XgammaY_Z_LM_UC[, col] / 2
  }
  rho1 = n[2] / n[1]
  ## calibration weight alpha
  alpha = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = cov(weight_LC * (K1_cond_XgammaY_LC[, col] * cbind(1, X_LC[, 1:3]) - K1_cond_XgammaY_Z_LC[, col] * cbind(1, X_LC[, 1:3])))
    var2 = cov(weight_LM * (K1_cond_XgammaY_LM[, col] * cbind(1, X_LM[, 1:3]) - K1_cond_XgammaY_Z_LM[, col] * cbind(1, X_LM[, 1:3])))
    cov1 = cov(weight_LC * (S_DR_tilt[, col] - K1_cond_Z_LC[, col]), 
               weight_LC * (K1_cond_XgammaY_LC[, col] * cbind(1, X_LC[, 1:3]) - K1_cond_XgammaY_Z_LC[, col] * cbind(1, X_LC[, 1:3])))
    alpha[, col] = solve(var1 + 1 / rho1 * var2) %*% as.vector(cov1)
  }
  M_phi1_LC = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_LC[, k] * cbind(1, X_LC[, 1:3]) %*% alpha[ , k]) - 
                             sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_Z_LC[, k] * cbind(1, X_LC[, 1:3]) %*% alpha[ , k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_Z_LC_UC[, k] * cbind(1, X_UC[, 1:3]) %*% alpha[ , k])), n[1]), nrow = 5))
  
  M_phi1_LM = weight_LM * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_LM[, k] * cbind(1, X_LM[, 1:3]) %*% alpha[ , k]) - 
                             sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_Z_LM[, k] * cbind(1, X_LM[, 1:3]) %*% alpha[ , k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_Z_LM_UC[, k] * cbind(1, X_UC[, 1:3]) %*% alpha[ , k])), n[2]), nrow = 5))
  ## gamma1
  gamma_1_ad = gamma_tilt - colMeans(M_phi1_LC) + colMeans(M_phi1_LM)
  ## incorporate unlabel data
  ## initial setup
  S1 = S_DR_tilt - sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY_LC[, k] * cbind(1, X_LC[, 1:3]) %*% alpha[ , k])
  nT2 = 1 + ncol(X_LC[, 1:3])
  T1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  for (col in 1:ncol(S1)) {
    T1 = S1[, col]
    for (ff in 1:4) {
      density_id = ff %% 4
      score_z_id = (ff + 1) %% 4
      score_x_id = (ff + 2) %% 4
      score_x_z_id = (ff + 3) %% 4
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, 1:3]
      T1_tr_score_z <- T1[fold_id_LC == score_z_id]
      X_LC_tr_score_x <- X_LC[fold_id_LC == score_x_id, 1:5] 
      T1_tr_score_x <- T1[fold_id_LC == score_x_id]
      best_hess_hp_score_z <- tune_gausspr_cv(X_LC_tr_score_z, T1_tr_score_z)
      best_hess_hp_score_x <- tune_gausspr_cv(X_LC_tr_score_x, T1_tr_score_x)
      invisible(capture.output({
        mod_T1_Z_LC = gausspr(X_LC_tr_score_z, 
                              T1_tr_score_z, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_z$sigma), 
                              var = best_hess_hp_score_z$var)
        mod_T1_X_LC = gausspr(X_LC_tr_score_x, 
                              T1_tr_score_x, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_x$sigma), 
                              var = best_hess_hp_score_x$var)
      }))
      T1_cond_X_LC[fold_id_LC == score_x_z_id, col] = predict(mod_T1_X_LC, newdata = X_LC[fold_id_LC == score_x_z_id, 1:5])
      X_LC_tr_score_x_z <- X_LC[fold_id_LC == score_x_z_id, 1:3]
      T1_tr_score_x_z <- T1_cond_X_LC[fold_id_LC == score_x_z_id, col]
      best_hess_hp_score_x_z <- tune_gausspr_cv(X_LC_tr_score_x_z, T1_tr_score_x_z)
      invisible(capture.output({
        mod_T1_X_Z_LC = gausspr(X_LC_tr_score_x_z, 
                                T1_tr_score_x_z,
                                kernel = "rbfdot", 
                                kpar = list(sigma = best_hess_hp_score_x_z$sigma), 
                                var = best_hess_hp_score_x_z$var)
      }))
      T1_cond_X_Z_LC[fold_id_LC == score_z_id, col] = predict(mod_T1_X_Z_LC, newdata = X_LC[fold_id_LC == score_z_id, 1:3])
      T1_cond_Z_LC[fold_id_LC == density_id, col] = predict(mod_T1_Z_LC, newdata = X_LC[fold_id_LC == density_id, 1:3])
      T1_cond_X_UC[, col] = predict(mod_T1_X_LC, newdata = X_UC[,1:5]) + T1_cond_X_UC[, col]
      T1_cond_X_Z_LC_UC[, col] = predict(mod_T1_X_Z_LC, newdata = X_UC[,1:3]) + T1_cond_X_Z_LC_UC[, col]
    }
    T1_cond_X_UC[, col] = T1_cond_X_UC[, col] / 4
    T1_cond_X_Z_LC_UC[, col] = T1_cond_X_Z_LC_UC[, col] / 4
  }
  rho2 = n[3] / n[1]
  ## calibration function weight zeta
  zeta = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = cov(weight_LC * (T1_cond_X_LC[, col] * cbind(1, X_LC[, 1:3]) - T1_cond_X_Z_LC[, col] * cbind(1, X_LC[, 1:3])))
    var2 = cov(T1_cond_X_UC[, col] * cbind(1, X_UC[, 1:3]))
    cov1 = cov(weight_LC * (S1[, col] - T1_cond_Z_LC[, col]),
               weight_LC * (T1_cond_X_LC[, col] * cbind(1, X_LC[, 1:3]) - T1_cond_X_Z_LC[, col] * cbind(1, X_LC[, 1:3])))
    zeta[, col] = solve(var1 + 1 / rho2 * var2) %*% as.vector(cov1)
  }
  M_phi2_LC = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_LC[, k] * cbind(1, X_LC[, 1:3]) %*% zeta[ , k]) - 
                             sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC[, k] * cbind(1, X_LC[, 1:3]) %*% zeta[ , k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC_UC[, k] * cbind(1, X_UC[, 1:3]) %*% zeta[ , k])), n[1]), nrow = 5))
  
  phi2_UC = sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_UC[, k] * cbind(1, X_UC[, 1:3]) %*% zeta[ , k]) 
  ## gamma2
  gamma_2_ad = gamma_1_ad - colMeans(M_phi2_LC) + colMeans(phi2_UC)
  ## without calibration or zero weight
  alpha0 = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  zeta0 = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  var_tilt = var_cal(S_DR_tilt, K1_cond_Z_LC,
                     weight_LC, weight_LM, alpha0, rho1,
                     K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, 
                     K1_cond_XgammaY_LM, K1_cond_XgammaY_Z_LM,
                     X_LC, X_LM)
  var_ad1 = var_cal(S_DR_tilt, K1_cond_Z_LC,
                    weight_LC, weight_LM, alpha, rho1,
                    K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, 
                    K1_cond_XgammaY_LM, K1_cond_XgammaY_Z_LM,
                    X_LC, X_LM)
  T_part_var0 = var_cal_2(S_DR_tilt, K1_cond_Z_LC,
                          K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, alpha, 
                          S1, T1_cond_Z_LC,
                          weight_LC, zeta0, rho2,
                          T1_cond_X_LC, T1_cond_X_Z_LC, 
                          T1_cond_X_UC,
                          X_LC, X_UC)
  T_part_var = var_cal_2(S_DR_tilt, K1_cond_Z_LC,
                         K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, alpha, 
                         S1, T1_cond_Z_LC,
                         weight_LC, zeta, rho2,
                         T1_cond_X_LC, T1_cond_X_Z_LC, 
                         T1_cond_X_UC,
                         X_LC, X_UC)
  var_ad2 = var_ad1 - (T_part_var0 - T_part_var)
  res = list(gamma_tilt = gamma_tilt, gamma_1_ad = gamma_1_ad, gamma_2_ad = gamma_2_ad,
             var_tilt = var_tilt, var_ad1 = var_ad1, T_part_var0 = T_part_var0,
             T_part_var = T_part_var)
  return(res)
}
### calculation of gamma1 and gamma2 for different calibration and reduced models
main_func2 <- function(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum){
  ## initial setup
  set.seed(seednum)
  X_LC = Combined_X[source_index[[1]],]
  X_LM = Combined_X[source_index[[2]],]
  X_UC = Combined_X[source_index[[3]],]
  Y_LC = Combined_Y[source_index[[1]]]
  Y_LM = Combined_Y[source_index[[2]]]
  Y_UC = Combined_Y[source_index[[3]]]
  fold_id_LC = sample(rep(0:3, length.out = n[1])) 
  fold_id_LM = sample(rep(0:1, length.out = n[2])) 
  K1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_Z_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_XgammaY_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY_Z_LM = matrix(0, nrow = nrow(X_LM), ncol = ncol(X_LM))
  K1_cond_XgammaY_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY_Z_LM_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  weight_LC = rep(0, n[1])
  weight_LM = rep(0, n[2])
  ## column-wisely train the model
  for (col in 1:ncol(S_DR_tilt)) {
    K1 = S_DR_tilt[, col]
    ## best linear projection for the reduced model, simple parametric without cross-fitting
    mod_K1_XgammaY_LC = lm(Y_LC ~ X_LC[, 1:3] - 1)
    K1_cond_XgammaY_LC = t(t(X_LC) %*% X_LC[, 1:3] %*% solve(t(X_LC[, 1:3]) %*% X_LC[, 1:3]) %*% t(X_LC[, 1:3] * as.numeric(mod_K1_XgammaY_LC$residuals)))
    K1_cond_XgammaY_LM = t(t(X_LC) %*% X_LC[, 1:3] %*% solve(t(X_LC[, 1:3]) %*% X_LC[, 1:3]) %*% t(X_LM[, 1:3] * as.numeric(Y_LM - X_LM[, 1:3] %*% mod_K1_XgammaY_LC$coefficients)))
    for (ff in 1:4) {
      density_id = ff %% 4
      score_z_id = (ff + 1) %% 4
      score_xsuby_id = (ff + 2) %% 4
      score_xsuby_z_id = (ff + 3) %% 4
      weight_lc = density_ratio_xgboost_linear(X_UC[ ,1:3], X_LC[fold_id_LC == density_id, 1:3], X_LC[fold_id_LC == score_xsuby_id, 1:3])
      weight_LC[fold_id_LC == score_xsuby_id] = weight_lc
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, 1:3]
      K1_tr_score_z <- K1[fold_id_LC == score_z_id]
      best_hess_hp_score_z <- tune_gausspr_cv(X_LC_tr_score_z, K1_tr_score_z)
      invisible(capture.output({
        mod_K1_Z_LC = gausspr(X_LC_tr_score_z, 
                              K1_tr_score_z, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_z$sigma), 
                              var = best_hess_hp_score_z$var)
      }))
      X_LC_tr_score_xsuby_z <- X_LC[fold_id_LC == score_xsuby_z_id, 1:3]
      K1_tr_score_xsuby_z <- K1_cond_XgammaY_LC[fold_id_LC == score_xsuby_z_id, col]
      best_hess_hp_score_xsuby_z <- tune_gausspr_cv(X_LC_tr_score_xsuby_z, K1_tr_score_xsuby_z)
      invisible(capture.output({
        mod_K1_XgammaY_Z_LC = gausspr(X_LC_tr_score_xsuby_z, 
                                      K1_tr_score_xsuby_z,
                                      kernel = "rbfdot", 
                                      kpar = list(sigma = best_hess_hp_score_xsuby_z$sigma), 
                                      var = best_hess_hp_score_xsuby_z$var)
      }))
      K1_cond_XgammaY_Z_LC[fold_id_LC == score_z_id, col] = predict(mod_K1_XgammaY_Z_LC, newdata = X_LC[fold_id_LC == score_z_id, 1:3])
      K1_cond_Z_LC[fold_id_LC == density_id, col] = predict(mod_K1_Z_LC, newdata = X_LC[fold_id_LC == density_id, 1:3])
      K1_cond_Z_UC[, col] = predict(mod_K1_Z_LC, newdata = X_UC[,1:3]) + K1_cond_Z_UC[, col]
      K1_cond_XgammaY_Z_LC_UC[, col] = predict(mod_K1_XgammaY_Z_LC, newdata = X_UC[,1:3]) + K1_cond_XgammaY_Z_LC_UC[, col]
    }
    K1_cond_Z_UC[, col] = K1_cond_Z_UC[, col] / 4
    K1_cond_XgammaY_Z_LC_UC[, col] = K1_cond_XgammaY_Z_LC_UC[, col] / 4
    for (ff in 1:2) {
      density_id = ff %% 2
      score_xsuby_z_id = (ff + 1) %% 2
      weight_lm = density_ratio_xgboost_linear(X_UC[ ,1:3], X_LM[fold_id_LM == density_id, 1:3], X_LM[fold_id_LM == score_xsuby_z_id, 1:3])
      weight_LM[fold_id_LM == score_xsuby_z_id] = weight_lm
      X_LM_tr_score_xsuby_z <- X_LM[fold_id_LM == score_xsuby_z_id, 1:3]
      K1_tr_score_xsuby_z <- K1_cond_XgammaY_LM[fold_id_LM == score_xsuby_z_id, col]
      best_hess_hp_score_xsuby_z <- tune_gausspr_cv(X_LM_tr_score_xsuby_z, K1_tr_score_xsuby_z)
      invisible(capture.output({
        mod_K1_XgammaY_Z_LM = gausspr(X_LM_tr_score_xsuby_z, 
                                      K1_tr_score_xsuby_z, 
                                      kernel = "rbfdot", 
                                      kpar = list(sigma = best_hess_hp_score_xsuby_z$sigma), 
                                      var = best_hess_hp_score_xsuby_z$var)
      }))
      K1_cond_XgammaY_Z_LM[fold_id_LM == density_id, col] = predict(mod_K1_XgammaY_Z_LM, newdata = X_LM[fold_id_LM == density_id, 1:3])
      K1_cond_XgammaY_Z_LM_UC[, col] = predict(mod_K1_XgammaY_Z_LM, newdata = X_UC[,1:3]) + K1_cond_XgammaY_Z_LM_UC[, col]
    }
    K1_cond_XgammaY_Z_LM_UC[, col] = K1_cond_XgammaY_Z_LM_UC[, col] / 2
  }
  rho1 = n[2] / n[1]
  nK2 = 1 # constant calibration for reduced model
  alpha = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = var(weight_LC * (K1_cond_XgammaY_LC[, col] - K1_cond_XgammaY_Z_LC[, col]))
    var2 = var(weight_LM * (K1_cond_XgammaY_LM[, col] - K1_cond_XgammaY_Z_LM[, col]))
    cov1 = cov(weight_LC * (S_DR_tilt[, col] - K1_cond_Z_LC[, col]), 
               weight_LC * (K1_cond_XgammaY_LC[, col] - K1_cond_XgammaY_Z_LC[, col]))
    alpha[, col] = solve(var1 + 1 / rho1 * var2) %*% as.vector(cov1)
  }
  M_phi1_LC = weight_LC * (K1_cond_XgammaY_LC * diag(alpha) - K1_cond_XgammaY_Z_LC * diag(alpha)) + 
    t(matrix(rep(colMeans(K1_cond_XgammaY_Z_LC_UC * diag(alpha)), n[1]), nrow = 5))
  
  M_phi1_LM = weight_LM * (K1_cond_XgammaY_LM * diag(alpha) - K1_cond_XgammaY_Z_LM * diag(alpha)) + 
    t(matrix(rep(colMeans(K1_cond_XgammaY_Z_LM_UC * diag(alpha)), n[2]), nrow = 5))
  gamma_1_ad = gamma_tilt - colMeans(M_phi1_LC) + colMeans(M_phi1_LM)
  alpha0 = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  var_tilt = var_cal(S_DR_tilt, K1_cond_Z_LC,
                     weight_LC, weight_LM, alpha0, rho1,
                     K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, 
                     K1_cond_XgammaY_LM, K1_cond_XgammaY_Z_LM,
                     X_LC, X_LM)
  var_ad1 = var_cal(S_DR_tilt, K1_cond_Z_LC,
                    weight_LC, weight_LM, alpha, rho1,
                    K1_cond_XgammaY_LC, K1_cond_XgammaY_Z_LC, 
                    K1_cond_XgammaY_LM, K1_cond_XgammaY_Z_LM,
                    X_LC, X_LM)
  res = list(gamma_tilt = gamma_tilt, gamma_1_ad = gamma_1_ad,
             var_tilt = var_tilt, var_ad1 = var_ad1)
  ## reduced model only have gamma1
  return(res)
}