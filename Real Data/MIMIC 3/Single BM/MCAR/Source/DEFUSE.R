# Get the proposed estimators
sip_dop_get <- function(seed, cali_method) {
  set.seed(seed)
  model = glm(Y_lc_train ~ ., cbind(Y_lc_train, X_lc_train), family = "binomial")
  gamma.tilt = model$coefficients
  H_inv = solve(H(as.matrix(cbind(1, X_lc_train)), gamma.tilt))
  S = sweep(t(H_inv %*% t(as.matrix(cbind(1, X_lc_train)))),
            MARGIN = 1,
            (Y_lc_train - g_x(as.matrix(cbind(1, X_lc_train)), gamma.tilt)), '*')
  
  mu_hat_C = c()
  mu_hat_IC = matrix(0, nrow = n_lic_train, ncol = 5)
  mu_hat_C_U_y0 = matrix(0, nrow = (n_lc_train + n_uc), ncol = 5)
  mu_hat_C_U_y1 = matrix(0, nrow = (n_lc_train + n_uc), ncol = 5)
  
  for (i in 1:5) {
    mu_hat_C[which((Y_lc_train == 0) & (folds == i))] =
      predict(mod[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)), -M])
    mu_hat_C[which((Y_lc_train == 1) & (folds == i))] =
      predict(mod[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)), -M])
    mu_hat_IC[which(Y_lic_train == 0), i] =
      predict(mod[[i]][[1]], newdata = X_lic_train_withauxi[which(Y_lic_train == 0), -M])
    mu_hat_IC[which(Y_lic_train == 1), i] =
      predict(mod[[i]][[2]], newdata = X_lic_train_withauxi[which(Y_lic_train == 1), -M])
    mu_hat_C_U_y0[, i] =
      predict(mod[[i]][[1]], newdata = rbind(X_lc_train_withauxi, X_uc_train_withauxi)[, -M])
    mu_hat_C_U_y1[, i] =
      predict(mod[[i]][[2]], newdata = rbind(X_lc_train_withauxi, X_uc_train_withauxi)[, -M])
  }
  
  mu_hat_IC = rowMeans(mu_hat_IC)
  mu_hat_C_U_y0 = rowMeans(mu_hat_C_U_y0)
  mu_hat_C_U_y1 = rowMeans(mu_hat_C_U_y1)
  mu_hat_C = colMeans(X_lc[, 4])
  mu_hat_IC = colMeans(X_lc[, 4])
  mu_hat_C_U_y0 = colMeans(X_lc[, 4])
  mu_hat_C_U_y1 = colMeans(X_lc[, 4])
  
  resi = pull(X_lc_train, M) - mu_hat_C
  sigma_hat = c()
  sigma_hat[which(Y_lc_train == 0)] = sqrt(mean(resi[which(Y_lc_train == 0)]^2))
  sigma_hat[which(Y_lc_train == 1)] = sqrt(mean(resi[which(Y_lc_train == 1)]^2))
  sigma_hat_IC = c()
  sigma_hat_IC[which(Y_lic_train == 0)] = sqrt(mean(resi[which(Y_lc_train == 0)]^2))
  sigma_hat_IC[which(Y_lic_train == 1)] = sqrt(mean(resi[which(Y_lc_train == 1)]^2))
  
  monte_sample = 300
  base_sample = rnorm(monte_sample)
  rows = rep(1:n_lc_train, each = monte_sample)
  obs_C = cbind(1, X_lc_train)[rows, ]
  obs_C[, (M + 1)] = rep(base_sample, n_lc_train) * rep(sigma_hat, each = monte_sample) + rep(mu_hat_C, each = monte_sample)
  S_monte_C = sweep(t(H_inv %*% t(as.matrix(obs_C))),
                    MARGIN = 1,
                    rep(Y_lc_train, each = monte_sample) - g_x(as.matrix(obs_C), gamma.tilt), '*')
  phi_1_C = aggregate(S_monte_C, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1]
  
  rowsI = rep(1:n_lic_train, each = monte_sample)
  obs_IC = cbind(1, X_lic_train)[rowsI, ]
  obs_IC[, (M + 1)] = rep(base_sample, n_lic_train) * rep(sigma_hat_IC, each = monte_sample) + rep(mu_hat_IC, each = monte_sample)
  S_monte_IC = sweep(t(H_inv %*% t(as.matrix(obs_IC))),
                     MARGIN = 1,
                     rep(Y_lic_train, each = monte_sample) - g_x(as.matrix(obs_IC), gamma.tilt), '*')
  phi_1_IC = aggregate(S_monte_IC, list(rep(1:n_lic_train, each = monte_sample)), mean)[-1]
  
  allocation_score = allocation(
    cali_method, phi_1_C, phi_1_IC, S, M,
    X_lc_train_withauxi_select, X_lic_train_withauxi_select,
    Y_lc_train, Y_lic_train
  )
  delta = allocation_score$delta
  weight_CIC = allocation_score$weight
  adjusted_phi_1 = allocation_score$adjusted_phi_1
  
  alpha_1 = rep(1, ncol(phi_1_C))
  for (j in 1:ncol(phi_1_C)) {
    alpha_1[j] = lm(S[, j] ~ phi_1_C[, j])$coefficients[2]
  }
  alpha_1 = diag(alpha_1)
  
  semi_gamma_tilt_1 = gamma.tilt - colMeans(phi_1_C) + colMeans(rbind(phi_1_C, phi_1_IC))
  gamma_tilt_1 = gamma.tilt - rowMeans(t(alpha_1) %*% t(phi_1_C)) + rowMeans(t(alpha_1) %*% t(rbind(phi_1_C, phi_1_IC)))
  adjusted_gamma_tilt_1 = gamma.tilt - colMeans(adjusted_phi_1[1:n_lc_train, ]) + colMeans(adjusted_phi_1)
  
  rows2 = rep(1:(n_lc_train + n_uc), each = monte_sample)
  obs_C_U_y0 = cbind(1, rbind(X_lc_train, X_uc_train))[rows2, ]
  obs_C_U_y1 = cbind(1, rbind(X_lc_train, X_uc_train))[rows2, ]
  obs_C_U_y0[, (M + 1)] =
    rep(base_sample, (n_lc_train + n_uc)) * rep(unique(sigma_hat[which(Y_lc_train == 0)]), monte_sample * (n_lc_train + n_uc)) +
    rep(mu_hat_C_U_y0, each = monte_sample)
  obs_C_U_y1[, (M + 1)] =
    rep(base_sample, (n_lc_train + n_uc)) * rep(unique(sigma_hat[which(Y_lc_train == 1)]), monte_sample * (n_lc_train + n_uc)) +
    rep(mu_hat_C_U_y1, each = monte_sample)
  
  S_monte_C_U_y0 = sweep(t(H_inv %*% t(as.matrix(obs_C_U_y0))), MARGIN = 1, rep(0, each = monte_sample) - g_x(as.matrix(obs_C_U_y0), gamma.tilt), '*')
  S_monte_C_U_y1 = sweep(t(H_inv %*% t(as.matrix(obs_C_U_y1))), MARGIN = 1, rep(1, each = monte_sample) - g_x(as.matrix(obs_C_U_y1), gamma.tilt), '*')
  phi_1_C_U_y0 = aggregate(S_monte_C_U_y0, list(rep(1:(n_lc_train + n_uc), each = monte_sample)), mean)[-1]
  phi_1_C_U_y1 = aggregate(S_monte_C_U_y1, list(rep(1:(n_lc_train + n_uc), each = monte_sample)), mean)[-1]
  p_y = g_x(as.matrix(cbind(1, rbind(X_lc_train, X_uc_train))), gamma.tilt)
  phi_2 = sweep(phi_1_C_U_y1, MARGIN = 1, p_y, '*') + sweep(phi_1_C_U_y0, MARGIN = 1, (1 - p_y), '*')
  weight_adjust = c(rep(1, length(Y_lc_train)), weight)
  phi_2 = sweep(phi_2, MARGIN = 1, weight_adjust, "*")
  
  if (cali_method == "without_Y") {
    XY_CU = rbind(as.matrix(cbind(X_lc_train_withauxi_select[, -M], 1)),
                  as.matrix(cbind(X_uc_train_withauxi_select[, -M], 1)))
    colnames(XY_CU) = c()
    weight_CU = as.matrix(XY_CU, ncol = ncol(XY_CU)) %*% delta + 1
    ad_phi_2 = phi_2 * weight_CU
    alpha_2 = rep(1, ncol(ad_phi_2))
    for (j in 1:ncol(ad_phi_2)) {
      alpha_2[j] = solve(t(ad_phi_2[1:n_lc_train, j]) %*% ad_phi_2[1:n_lc_train, j]) %*%
        t(ad_phi_2[1:n_lc_train, j]) %*%
        (S[, j] - n_lic_train / (n_lc_train + n_lic_train) * adjusted_phi_1[1:n_lc_train, j])
    }
    alpha_2 = diag(alpha_2)
    tem_gamma_tilt1 = adjusted_gamma_tilt_1
  } else {
    ad_phi_2 = phi_2
    alpha_2 = rep(1, ncol(ad_phi_2))
    for (j in 1:ncol(ad_phi_2)) {
      alpha_2[j] = solve(t(ad_phi_2[1:n_lc_train, j]) %*% ad_phi_2[1:n_lc_train, j]) %*%
        t(ad_phi_2[1:n_lc_train, j]) %*%
        (S[, j] - n_lic_train / (n_lc_train + n_lic_train) * (as.matrix(phi_1_C) %*% alpha_1)[, j])
    }
    alpha_2 = diag(alpha_2)
    tem_gamma_tilt1 = gamma_tilt_1
  }
  
  gamma_tilt_2 = tem_gamma_tilt1 - rowMeans(t(alpha_2) %*% t(ad_phi_2[1:n_lc_train, ])) + rowMeans(t(alpha_2) %*% t(ad_phi_2))
  a = colVars(S - t(t(phi_1_C)) %*% alpha_1)
  b = colVars(t(t(phi_1_C)) %*% alpha_1)
  var = (a + b) / n_lc_train
  
  semi_phi1 = n_lic_train / (n_lc_train + n_lic_train) * t(t(rbind(phi_1_C, phi_1_IC)))
  semi_S1 = S - semi_phi1[1:n_lc_train, ]
  const_phi1 = n_lic_train / (n_lc_train + n_lic_train) * t(t(rbind(phi_1_C, phi_1_IC))) %*% alpha_1
  const_S1 = S - const_phi1[1:n_lc_train, ]
  linear_phi1 = n_lic_train / (n_lc_train + n_lic_train) * adjusted_phi_1
  linear_S1 = S - linear_phi1[1:n_lc_train, ]
  const_phi2 = n_uc / (n_lc_train + n_uc) * t(t(ad_phi_2)) %*% alpha_2
  
  semi_var1 = colVars(semi_S1) / n_lc_train + colVars(semi_phi1[1:n_lc_train, ]) / n_lic_train
  var1 = colVars(const_S1) / n_lc_train + colVars(const_phi1[1:n_lc_train, ]) / n_lic_train
  adjusted_var1 = colVars(linear_S1) / n_lc_train + colVars(linear_phi1[1:n_lc_train, ]) / n_lic_train
  
  if (cali_method == "without_Y") {
    var2 = colVars(linear_S1 - const_phi2[1:n_lc_train, ]) / n_lc_train +
      colVars(const_phi1[1:n_lc_train, ]) / n_lic_train +
      colVars(const_phi2[1:n_lc_train, ]) / n_uc
  } else {
    var2 = colVars(const_S1 - const_phi2[1:n_lc_train, ]) / n_lc_train +
      colVars(const_phi1[1:n_lc_train, ]) / n_lic_train +
      colVars(const_phi2[1:n_lc_train, ]) / n_uc
  }
  
  return(list(
    gamma = gamma.tilt, se = sqrt(var),
    gamma1 = gamma_tilt_1, se1 = sqrt(var1),
    gamma2 = gamma_tilt_2, se2 = sqrt(var2),
    semi_gamma1 = semi_gamma_tilt_1, semi_se = sqrt(semi_var1),
    ad_gamma1 = adjusted_gamma_tilt_1, ad_se1 = sqrt(adjusted_var1)
  ))
}

