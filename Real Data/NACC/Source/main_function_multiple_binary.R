### calculation of gamma1 and gamma2 for different calibration and conditional models
main_func <- function(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum){
  ## initial setup
  set.seed(seednum)
  X_LC = Combined_X[source_index[[1]],]
  X_LM1 = Combined_X[source_index[[2]],]
  X_UC = Combined_X[source_index[[3]],]
  Y_LC = Combined_Y[source_index[[1]]]
  Y_LM1 = Combined_Y[source_index[[2]]]
  Y_UC = Combined_Y[source_index[[3]]]
  fold_id_LC = sample(rep(0:2, length.out = n[1])) 
  fold_id_LM1 = sample(rep(0:1, length.out = n[2])) 
  K1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_Z_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY1_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_XgammaY1_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY1_Z_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_XgammaY1_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_Z_LM1_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  weight_LC = rep(0, n[1])
  weight_LM1 = rep(0, n[2])
  covariate_lic1 = which(!is.na(X_LM1[1,]))
  covariate_z = c(1:5)
  ## column-wisely construct conditional models
  for (col in 1:ncol(S_DR_tilt)) {
    K1 = S_DR_tilt[, col]
    for (ff in 1:3) {
      score_xsuby1_id = (ff) %% 3
      density_id = (ff + 1) %% 3
      score_z_id = (ff + 1) %% 3
      score_xsuby1_z_id = (ff + 1) %% 3
      eva_id = (ff + 2) %% 3
      weight_lc = density_ratio_xgboost_linear(X_UC[, covariate_z], 
                                               X_LC[fold_id_LC == density_id, covariate_z], 
                                               X_LC[fold_id_LC == eva_id, covariate_z])
      weight_LC[fold_id_LC == eva_id] = weight_lc
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, covariate_z]
      K1_tr_score_z <- K1[fold_id_LC == score_z_id]
      XY_LC_tr_score_xsuby1 <- cbind(X_LC[fold_id_LC == score_xsuby1_id, covariate_lic1], Y_LC[fold_id_LC == score_xsuby1_id])
      K1_tr_score_xsuby1 <- K1[fold_id_LC == score_xsuby1_id]
      best_hess_hp_score_z <- tune_gausspr_cv(X_LC_tr_score_z, K1_tr_score_z)
      best_hess_hp_score_xsuby1 <- tune_gausspr_cv(XY_LC_tr_score_xsuby1, K1_tr_score_xsuby1)
      invisible(capture.output({
        mod_K1_Z_LC = gausspr(X_LC_tr_score_z, 
                              K1_tr_score_z, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_z$sigma), 
                              var = best_hess_hp_score_z$var)
        mod_K1_XgammaY1_LC = gausspr(XY_LC_tr_score_xsuby1, 
                                     K1_tr_score_xsuby1, 
                                     kernel = "rbfdot", 
                                     kpar = list(sigma = best_hess_hp_score_xsuby1$sigma), 
                                     var = best_hess_hp_score_xsuby1$var)
      }))
      K1_cond_XgammaY1_LC[fold_id_LC == eva_id, col] = predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LC[fold_id_LC == eva_id, covariate_lic1], Y_LC[fold_id_LC == eva_id]))
      X_LC_tr_score_xsuby_z <- X_LC[fold_id_LC == score_xsuby1_z_id, covariate_z]
      K1_tr_score_xsuby1_z <- predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LC[fold_id_LC == score_xsuby1_z_id, covariate_lic1], Y_LC[fold_id_LC == score_xsuby1_z_id]))
      best_hess_hp_score_xsuby1_z <- tune_gausspr_cv(X_LC_tr_score_xsuby_z, K1_tr_score_xsuby1_z)
      invisible(capture.output({
        mod_K1_XgammaY1_Z_LC = gausspr(X_LC_tr_score_xsuby_z, 
                                       K1_tr_score_xsuby1_z,
                                       kernel = "rbfdot", 
                                       kpar = list(sigma = best_hess_hp_score_xsuby1_z$sigma), 
                                       var = best_hess_hp_score_xsuby1_z$var)
        
      }))
      K1_cond_XgammaY1_Z_LC[fold_id_LC == eva_id, col] = predict(mod_K1_XgammaY1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      K1_cond_Z_LC[fold_id_LC == eva_id, col] = predict(mod_K1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      K1_cond_Z_UC[, col] = predict(mod_K1_Z_LC, newdata = X_UC[, covariate_z]) + K1_cond_Z_UC[, col]
      K1_cond_XgammaY1_LM1[, col] = predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LM1[, covariate_lic1], Y_LM1)) + K1_cond_XgammaY1_LM1[, col]
      K1_cond_XgammaY1_Z_LC_UC[, col] = predict(mod_K1_XgammaY1_Z_LC, newdata = X_UC[, covariate_z]) + K1_cond_XgammaY1_Z_LC_UC[, col]
    }
    K1_cond_Z_UC[, col] = K1_cond_Z_UC[, col] / 3
    K1_cond_XgammaY1_LM1[, col] = K1_cond_XgammaY1_LM1[, col] / 3
    K1_cond_XgammaY1_Z_LC_UC[, col] = K1_cond_XgammaY1_Z_LC_UC[, col] / 3
    for (ff in 1:2) {
      density_id = ff %% 2
      score_xsuby1_z_id = (ff) %% 2
      eva_id = (ff + 1) %% 2
      weight_lm1 = density_ratio_xgboost_linear(X_UC[, covariate_z],
                                                X_LM1[fold_id_LM1 == density_id, covariate_z], 
                                                X_LM1[fold_id_LM1 == eva_id, covariate_z])
      weight_LM1[fold_id_LM1 == eva_id] = weight_lm1
      X_LM1_tr_score_xsuby1_z <- X_LM1[fold_id_LM1 == score_xsuby1_z_id, covariate_z]
      K1_tr_score_xsuby1_z <- K1_cond_XgammaY1_LM1[fold_id_LM1 == score_xsuby1_z_id, col]
      best_hess_hp_score_xsuby1_z <- tune_gausspr_cv(X_LM1_tr_score_xsuby1_z, K1_tr_score_xsuby1_z)
      invisible(capture.output({
        mod_K1_XgammaY1_Z_LM1 = gausspr(X_LM1_tr_score_xsuby1_z, 
                                        K1_tr_score_xsuby1_z, 
                                        kernel = "rbfdot", 
                                        kpar = list(sigma = best_hess_hp_score_xsuby1_z$sigma), 
                                        var = best_hess_hp_score_xsuby1_z$var)
        
      }))
      K1_cond_XgammaY1_Z_LM1[fold_id_LM1 == eva_id, col] = predict(mod_K1_XgammaY1_Z_LM1, newdata = X_LM1[fold_id_LM1 == eva_id, covariate_z])
      K1_cond_XgammaY1_Z_LM1_UC[, col] = predict(mod_K1_XgammaY1_Z_LM1, newdata = X_UC[, covariate_z]) + K1_cond_XgammaY1_Z_LM1_UC[, col]
    }
    K1_cond_XgammaY1_Z_LM1_UC[, col] = K1_cond_XgammaY1_Z_LM1_UC[, col] / 2
  }
  rho1 = n[2] / n[1]
  nK2 = 1 + ncol(X_LC[, covariate_z])
  ## calibration weight alpha
  alpha = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = cov(weight_LC * (K1_cond_XgammaY1_LC[, col] * cbind(1, X_LC[, covariate_z]) - K1_cond_XgammaY1_Z_LC[, col] * cbind(1, X_LC[, covariate_z])))
    var2 = cov(weight_LM1 * (K1_cond_XgammaY1_LM1[, col] * cbind(1, X_LM1[, covariate_z]) - K1_cond_XgammaY1_Z_LM1[, col] * cbind(1, X_LM1[, covariate_z])))
    cov1 = cov(weight_LC * (S_DR_tilt[, col] - K1_cond_Z_LC[, col]), 
               weight_LC * (K1_cond_XgammaY1_LC[, col] * cbind(1, X_LC[, covariate_z]) - K1_cond_XgammaY1_Z_LC[, col] * cbind(1, X_LC[, covariate_z])))
    alpha[, col] = solve(var1 + 1 / rho1 * var2) %*% as.vector(cov1)
  }
  M_phi1_LC1 = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[, k]) - 
                              sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[, k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC_UC[, k] * cbind(1, X_UC[, covariate_z]) %*% alpha[, k])), n[1]), nrow = ncol(X_UC)))
  M_phi1_LM1 = weight_LM1 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LM1[, k] * cbind(1, X_LM1[, covariate_z]) %*% alpha[, k]) - 
                               sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1[, k] * cbind(1, X_LM1[, covariate_z]) %*% alpha[, k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1_UC[, k] * cbind(1, X_UC[, covariate_z]) %*% alpha[, k])), n[2]), nrow = ncol(X_UC)))
  gamma_1_ad = gamma_tilt - colMeans(M_phi1_LC1) + colMeans(M_phi1_LM1) 
  ## incorporate unlabel data
  ## initial setup
  S1 = S_DR_tilt - sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% alpha[ , k])
  nT2 = 1 + ncol(X_LC[, covariate_z])
  T1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  for (col in 1:ncol(S1)) {
    T1 = S1[, col]
    for (ff in 1:3) {
      score_x_id = (ff) %% 3
      density_id = (ff + 1) %% 3
      score_z_id = (ff + 1) %% 3
      score_x_z_id = (ff + 1) %% 3
      eva_id = (ff + 2) %% 3
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, covariate_z]
      T1_tr_score_z <- T1[fold_id_LC == score_z_id]
      X_LC_tr_score_x <- X_LC[fold_id_LC == score_x_id, ] 
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
      T1_cond_X_LC[fold_id_LC == score_x_z_id, col] = predict(mod_T1_X_LC, newdata = X_LC[fold_id_LC == score_x_z_id, ])
      X_LC_tr_score_x_z <- X_LC[fold_id_LC == score_x_z_id, covariate_z]
      T1_tr_score_x_z <- T1_cond_X_LC[fold_id_LC == score_x_z_id, col]
      best_hess_hp_score_x_z <- tune_gausspr_cv(X_LC_tr_score_x_z, T1_tr_score_x_z)
      invisible(capture.output({
        mod_T1_X_Z_LC = gausspr(X_LC_tr_score_x_z, 
                                T1_tr_score_x_z,
                                kernel = "rbfdot", 
                                kpar = list(sigma = best_hess_hp_score_x_z$sigma), 
                                var = best_hess_hp_score_x_z$var)
      }))
      T1_cond_X_Z_LC[fold_id_LC == eva_id, col] = predict(mod_T1_X_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      T1_cond_Z_LC[fold_id_LC == eva_id, col] = predict(mod_T1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      T1_cond_X_UC[, col] = predict(mod_T1_X_LC, newdata = X_UC[, ]) + T1_cond_X_UC[, col]
      T1_cond_X_Z_LC_UC[, col] = predict(mod_T1_X_Z_LC, newdata = X_UC[, covariate_z]) + T1_cond_X_Z_LC_UC[, col]
    }
    T1_cond_X_UC[, col] = T1_cond_X_UC[, col] / 3
    T1_cond_X_Z_LC_UC[, col] = T1_cond_X_Z_LC_UC[, col] / 3
  }
  rho2 = n[3] / n[1]
  zeta = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = cov(weight_LC * (T1_cond_X_LC[, col] * cbind(1, X_LC[, covariate_z]) - T1_cond_X_Z_LC[, col] * cbind(1, X_LC[, covariate_z])))
    var2 = cov(T1_cond_X_UC[, col] * cbind(1, X_UC[, covariate_z]))
    cov1 = cov(weight_LC * (S1[, col] - T1_cond_Z_LC[, col]),
               weight_LC * (T1_cond_X_LC[, col] * cbind(1, X_LC[, covariate_z]) - T1_cond_X_Z_LC[, col] * cbind(1, X_LC[, covariate_z])))
    zeta[, col] = solve(var1 + 1 / rho2 * var2) %*% as.vector(cov1)
  }
  M_phi2_LC = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% zeta[ , k]) - 
                             sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC[, k] * cbind(1, X_LC[, covariate_z]) %*% zeta[ , k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC_UC[, k] * cbind(1, X_UC[, covariate_z]) %*% zeta[ , k])), n[1]), nrow = ncol(X_UC)))
  phi2_UC = sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_UC[, k] * cbind(1, X_UC[, covariate_z]) %*% zeta[ , k]) 
  gamma_2_ad = gamma_1_ad - colMeans(M_phi2_LC) + colMeans(phi2_UC)
  alpha0 = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  zeta0 = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  var_tilt = var_cal(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                     weight_LC, weight_LM1, alpha0, rho1,
                     K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, 
                     K1_cond_XgammaY1_LM1, K1_cond_XgammaY1_Z_LM1,
                     X_LC, X_LM1)
  var_ad1 = var_cal(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                    weight_LC, weight_LM1, alpha, rho1,
                    K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, 
                    K1_cond_XgammaY1_LM1, K1_cond_XgammaY1_Z_LM1,
                    X_LC, X_LM1)
  T_part_var0 = var_cal_2(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                          K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, alpha, 
                          S1, T1_cond_Z_LC,
                          weight_LC, zeta0, rho2,
                          T1_cond_X_LC, T1_cond_X_Z_LC, 
                          T1_cond_X_UC,
                          X_LC, X_UC)
  T_part_var = var_cal_2(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                         K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, alpha, 
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
### calculation of gamma1 and gamma2 for different calibration
main_func2 <- function(Combined_X, Combined_Y, n, S_DR_tilt, gamma_tilt, source_index, seednum){
  set.seed(seednum)
  X_LC = Combined_X[source_index[[1]],]
  X_LM1 = Combined_X[source_index[[2]],]
  X_UC = Combined_X[source_index[[3]],]
  Y_LC = Combined_Y[source_index[[1]]]
  Y_LM1 = Combined_Y[source_index[[2]]]
  Y_UC = Combined_Y[source_index[[3]]]
  fold_id_LC = sample(rep(0:2, length.out = n[1])) 
  fold_id_LM1 = sample(rep(0:1, length.out = n[2])) 
  K1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_Z_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY1_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_XgammaY1_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  K1_cond_XgammaY1_Z_LM1 = matrix(0, nrow = nrow(X_LM1), ncol = ncol(X_LM1))
  K1_cond_XgammaY1_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  K1_cond_XgammaY1_Z_LM1_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  weight_LC = rep(0, n[1])
  weight_LM1 = rep(0, n[2])
  covariate_lic1 = which(!is.na(X_LM1[1,]))
  covariate_z = c(1:5)
  for (col in 1:ncol(S_DR_tilt)) {
    K1 = S_DR_tilt[, col]
    for (ff in 1:3) {
      score_xsuby1_id = (ff) %% 3
      density_id = (ff + 1) %% 3
      score_z_id = (ff + 1) %% 3
      score_xsuby1_z_id = (ff + 1) %% 3
      eva_id = (ff + 2) %% 3
      weight_lc = density_ratio_xgboost_linear(X_UC[, covariate_z], 
                                               X_LC[fold_id_LC == density_id, covariate_z], 
                                               X_LC[fold_id_LC == eva_id, covariate_z])
      weight_LC[fold_id_LC == eva_id] = weight_lc
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, covariate_z]
      K1_tr_score_z <- K1[fold_id_LC == score_z_id]
      XY_LC_tr_score_xsuby1 <- cbind(X_LC[fold_id_LC == score_xsuby1_id, covariate_lic1], Y_LC[fold_id_LC == score_xsuby1_id])
      K1_tr_score_xsuby1 <- K1[fold_id_LC == score_xsuby1_id]
      best_hess_hp_score_z <- tune_gausspr_cv(X_LC_tr_score_z, K1_tr_score_z)
      best_hess_hp_score_xsuby1 <- tune_gausspr_cv(XY_LC_tr_score_xsuby1, K1_tr_score_xsuby1)
      invisible(capture.output({
        mod_K1_Z_LC = gausspr(X_LC_tr_score_z, 
                              K1_tr_score_z, 
                              kernel = "rbfdot", 
                              kpar = list(sigma = best_hess_hp_score_z$sigma), 
                              var = best_hess_hp_score_z$var)
        mod_K1_XgammaY1_LC = gausspr(XY_LC_tr_score_xsuby1, 
                                     K1_tr_score_xsuby1, 
                                     kernel = "rbfdot", 
                                     kpar = list(sigma = best_hess_hp_score_xsuby1$sigma), 
                                     var = best_hess_hp_score_xsuby1$var)
      }))
      K1_cond_XgammaY1_LC[fold_id_LC == eva_id, col] = predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LC[fold_id_LC == eva_id, covariate_lic1], Y_LC[fold_id_LC == eva_id]))
      X_LC_tr_score_xsuby_z <- X_LC[fold_id_LC == score_xsuby1_z_id, covariate_z]
      K1_tr_score_xsuby1_z <- predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LC[fold_id_LC == score_xsuby1_z_id, covariate_lic1], Y_LC[fold_id_LC == score_xsuby1_z_id]))
      best_hess_hp_score_xsuby1_z <- tune_gausspr_cv(X_LC_tr_score_xsuby_z, K1_tr_score_xsuby1_z)
      invisible(capture.output({
        mod_K1_XgammaY1_Z_LC = gausspr(X_LC_tr_score_xsuby_z, 
                                       K1_tr_score_xsuby1_z,
                                       kernel = "rbfdot", 
                                       kpar = list(sigma = best_hess_hp_score_xsuby1_z$sigma), 
                                       var = best_hess_hp_score_xsuby1_z$var)
        
      }))
      K1_cond_XgammaY1_Z_LC[fold_id_LC == eva_id, col] = predict(mod_K1_XgammaY1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      K1_cond_Z_LC[fold_id_LC == eva_id, col] = predict(mod_K1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      K1_cond_Z_UC[, col] = predict(mod_K1_Z_LC, newdata = X_UC[, covariate_z]) + K1_cond_Z_UC[, col]
      K1_cond_XgammaY1_LM1[, col] = predict(mod_K1_XgammaY1_LC, newdata = cbind(X_LM1[, covariate_lic1], Y_LM1)) + K1_cond_XgammaY1_LM1[, col]
      K1_cond_XgammaY1_Z_LC_UC[, col] = predict(mod_K1_XgammaY1_Z_LC, newdata = X_UC[, covariate_z]) + K1_cond_XgammaY1_Z_LC_UC[, col]
    }
    K1_cond_Z_UC[, col] = K1_cond_Z_UC[, col] / 3
    K1_cond_XgammaY1_LM1[, col] = K1_cond_XgammaY1_LM1[, col] / 3
    K1_cond_XgammaY1_Z_LC_UC[, col] = K1_cond_XgammaY1_Z_LC_UC[, col] / 3
    for (ff in 1:2) {
      density_id = ff %% 2
      score_xsuby1_z_id = (ff) %% 2
      eva_id = (ff + 1) %% 2
      weight_lm1 = density_ratio_xgboost_linear(X_UC[, covariate_z],
                                                X_LM1[fold_id_LM1 == density_id, covariate_z], 
                                                X_LM1[fold_id_LM1 == eva_id, covariate_z])
      weight_LM1[fold_id_LM1 == eva_id] = weight_lm1
      X_LM1_tr_score_xsuby1_z <- X_LM1[fold_id_LM1 == score_xsuby1_z_id, covariate_z]
      K1_tr_score_xsuby1_z <- K1_cond_XgammaY1_LM1[fold_id_LM1 == score_xsuby1_z_id, col]
      best_hess_hp_score_xsuby1_z <- tune_gausspr_cv(X_LM1_tr_score_xsuby1_z, K1_tr_score_xsuby1_z)
      invisible(capture.output({
        mod_K1_XgammaY1_Z_LM1 = gausspr(X_LM1_tr_score_xsuby1_z, 
                                        K1_tr_score_xsuby1_z, 
                                        kernel = "rbfdot", 
                                        kpar = list(sigma = best_hess_hp_score_xsuby1_z$sigma), 
                                        var = best_hess_hp_score_xsuby1_z$var)
      }))
      K1_cond_XgammaY1_Z_LM1[fold_id_LM1 == eva_id, col] = predict(mod_K1_XgammaY1_Z_LM1, newdata = X_LM1[fold_id_LM1 == eva_id, covariate_z])
      K1_cond_XgammaY1_Z_LM1_UC[, col] = predict(mod_K1_XgammaY1_Z_LM1, newdata = X_UC[, covariate_z]) + K1_cond_XgammaY1_Z_LM1_UC[, col]
    }
    K1_cond_XgammaY1_Z_LM1_UC[, col] = K1_cond_XgammaY1_Z_LM1_UC[, col] / 2
  }
  rho1 = n[2] / n[1]
  nK2 = 1
  alpha = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = var(weight_LC * (K1_cond_XgammaY1_LC[, col] - K1_cond_XgammaY1_Z_LC[, col]))
    var2 = var(weight_LM1 * (K1_cond_XgammaY1_LM1[, col] - K1_cond_XgammaY1_Z_LM1[, col]))
    cov1 = cov(weight_LC * (S_DR_tilt[, col] - K1_cond_Z_LC[, col]), 
               weight_LC * (K1_cond_XgammaY1_LC[, col] - K1_cond_XgammaY1_Z_LC[, col]))
    alpha[, col] = solve(var1 + 1 / rho1 * var2) %*% as.vector(cov1)
  }
  M_phi1_LC1 = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * alpha[, k]) - 
                              sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC[, k] * alpha[, k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LC_UC[, k] * alpha[, k])), n[1]), nrow = ncol(X_UC)))
  M_phi1_LM1 = weight_LM1 * (sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LM1[, k] * alpha[, k]) - 
                               sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1[, k] * alpha[, k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_Z_LM1_UC[, k] * alpha[, k])), n[2]), nrow = ncol(X_UC)))
  gamma_1_ad = gamma_tilt - colMeans(M_phi1_LC1) + colMeans(M_phi1_LM1) 
  S1 = S_DR_tilt - sapply(seq_len(ncol(S_DR_tilt)), function(k) K1_cond_XgammaY1_LC[, k] * alpha[ , k])
  T1_cond_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_X_Z_LC = matrix(0, nrow = nrow(X_LC), ncol = ncol(X_LC))
  T1_cond_Z_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  T1_cond_X_Z_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = ncol(X_UC))
  for (col in 1:ncol(S1)) {
    T1 = S1[, col]
    for (ff in 1:3) {
      score_x_id = (ff) %% 3
      density_id = (ff + 1) %% 3
      score_z_id = (ff + 1) %% 3
      score_x_z_id = (ff + 1) %% 3
      eva_id = (ff + 2) %% 3
      X_LC_tr_score_z <- X_LC[fold_id_LC == score_z_id, covariate_z]
      T1_tr_score_z <- T1[fold_id_LC == score_z_id]
      X_LC_tr_score_x <- X_LC[fold_id_LC == score_x_id, ] 
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
      T1_cond_X_LC[fold_id_LC == score_x_z_id, col] = predict(mod_T1_X_LC, newdata = X_LC[fold_id_LC == score_x_z_id, ])
      X_LC_tr_score_x_z <- X_LC[fold_id_LC == score_x_z_id, covariate_z]
      T1_tr_score_x_z <- T1_cond_X_LC[fold_id_LC == score_x_z_id, col]
      best_hess_hp_score_x_z <- tune_gausspr_cv(X_LC_tr_score_x_z, T1_tr_score_x_z)
      invisible(capture.output({
        mod_T1_X_Z_LC = gausspr(X_LC_tr_score_x_z, 
                                T1_tr_score_x_z,
                                kernel = "rbfdot", 
                                kpar = list(sigma = best_hess_hp_score_x_z$sigma), 
                                var = best_hess_hp_score_x_z$var)
      }))
      T1_cond_X_Z_LC[fold_id_LC == eva_id, col] = predict(mod_T1_X_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      T1_cond_Z_LC[fold_id_LC == eva_id, col] = predict(mod_T1_Z_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
      T1_cond_X_UC[, col] = predict(mod_T1_X_LC, newdata = X_UC[, ]) + T1_cond_X_UC[, col]
      T1_cond_X_Z_LC_UC[, col] = predict(mod_T1_X_Z_LC, newdata = X_UC[, covariate_z]) + T1_cond_X_Z_LC_UC[, col]
    }
    T1_cond_X_UC[, col] = T1_cond_X_UC[, col] / 3
    T1_cond_X_Z_LC_UC[, col] = T1_cond_X_Z_LC_UC[, col] / 3
  }
  rho2 = n[3] / n[1]
  nT2 = 1 
  zeta = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  for (col in 1:ncol(S_DR_tilt)) {
    var1 = var(weight_LC * (T1_cond_X_LC[, col] - T1_cond_X_Z_LC[, col]))
    var2 = var(T1_cond_X_UC[, col])
    cov1 = cov(weight_LC * (S1[, col] - T1_cond_Z_LC[, col]),
               weight_LC * (T1_cond_X_LC[, col] - T1_cond_X_Z_LC[, col]))
    zeta[, col] = solve(var1 + 1 / rho2 * var2) %*% as.vector(cov1)
  }
  M_phi2_LC = weight_LC * (sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_LC[, k] * zeta[ , k]) - 
                             sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC[, k] * zeta[ , k])) + 
    t(matrix(rep(colMeans(sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_Z_LC_UC[, k] * zeta[ , k])), n[1]), nrow = ncol(X_UC)))
  phi2_UC = sapply(seq_len(ncol(S_DR_tilt)), function(k) T1_cond_X_UC[, k] * zeta[ , k]) 
  gamma_2_ad = gamma_1_ad - colMeans(M_phi2_LC) + colMeans(phi2_UC)
  alpha0 = matrix(0, nrow = nK2, ncol = ncol(S_DR_tilt))
  zeta0 = matrix(0, nrow = nT2, ncol = ncol(S_DR_tilt))
  var_tilt = var_cal(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                     weight_LC, weight_LM1, alpha0, rho1,
                     K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, 
                     K1_cond_XgammaY1_LM1, K1_cond_XgammaY1_Z_LM1,
                     X_LC, X_LM1)
  var_ad1 = var_cal(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                    weight_LC, weight_LM1, alpha, rho1,
                    K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, 
                    K1_cond_XgammaY1_LM1, K1_cond_XgammaY1_Z_LM1,
                    X_LC, X_LM1)
  T_part_var0 = var_cal_2(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                          K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, alpha, 
                          S1, T1_cond_Z_LC,
                          weight_LC, zeta0, rho2,
                          T1_cond_X_LC, T1_cond_X_Z_LC, 
                          T1_cond_X_UC,
                          X_LC, X_UC)
  T_part_var = var_cal_2(S_DR_tilt, K1_cond_Z_LC, covariate_z,
                         K1_cond_XgammaY1_LC, K1_cond_XgammaY1_Z_LC, alpha, 
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