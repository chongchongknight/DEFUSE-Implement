main_func <- function(Combined_X, Combined_Y, n, K, H_inv, S, gamma.tilt, naive_sd, folds, source_index, miss_source, imput_method, seednum){
  set.seed(seednum)
  ctrl <- trainControl(method = "cv", number = 5)
  pppx = ncol(Combined_X)
  colnames = c(colnames(Combined_X)[c(miss_source[[2]], setdiff(1:pppx, miss_source[[2]]))], "Y")
  X_complete = as.matrix(Combined_X[source_index[[1]],], nrow = n[1])
  Y_complete = as.matrix(Combined_Y[source_index[[1]]], ncol = 1)
  phi_1_C_list = vector("list", length = K - 2)
  phi_1_IC_list = vector("list", length = K - 2)
  phi_1_C0_list = vector("list", length = K - 2)
  phi_1_IC0_list = vector("list", length = K - 2)
  phi_1_monte_list = vector("list", length = K - 2)
  phi_1_monte0_list = vector("list", length = K - 2)
  phi_2_list = vector("list", length = K - 2)
  for (i in 2:(K - 1)) {
    n_missing = length(miss_source[[i]])
    mu_hat_C = matrix(0, nrow = n[1], ncol = n_missing) 
    mu_hat_IC = matrix(0, nrow = n[i], ncol = n_missing)
    X_incomplete = as.matrix(Combined_X[source_index[[i]], -miss_source[[i]]], n[i])
    Y_incomplete = as.matrix(Combined_Y[source_index[[i]]], ncol = 1)
    monte_sample = 100
    n_complete_X = n[1] + n[3]
    mu_hat_monte_C = matrix(0, nrow = n_complete_X * monte_sample, ncol = n_missing)
    monte_error = rnorm(monte_sample, sd = naive_sd)
    X_lc_uc = as.matrix(rbind(Combined_X[source_index[[1]], ], Combined_X[source_index[[3]], ]), nrow = n_complete_X)
    monte_X = X_lc_uc[rep(1:n_complete_X, each = monte_sample), ]
    monte_y = X_lc_uc[rep(1:n_complete_X, each = monte_sample), ] %*% gamma.tilt + rep(monte_error, n_complete_X)
    monte_data = data.frame(cbind(monte_X[, c(miss_source[[2]], setdiff(1:pppx, miss_source[[2]]))], monte_y))
    colnames(monte_data) = colnames
    for (j in 1:5) {
      complete_miss = X_complete[which(folds != j), miss_source[[i]]]
      X_complete_obs = X_complete[which(folds != j), -miss_source[[i]]]
      complete_obs = cbind(X_complete_obs, Y_complete[which(folds != j)])
      complete_data = cbind(complete_miss, complete_obs)
      colnames(complete_data) = colnames
      len_miss = length(miss_source[[i]])
      for (n_miss in 1:len_miss) {
        target = colnames[n_miss]
        formula_later = paste(colnames[(len_miss + 1):length(colnames)], collapse = "+")
        formula = as.formula(paste(target, "~ ", formula_later))
        mod = train(formula, data = complete_data, method = imput_method, trControl = ctrl)
        prefict_complete_obs = cbind(matrix(0, nrow = sum(folds == j), ncol = len_miss), cbind(X_complete[which(folds == j), -miss_source[[i]]], Y_complete[which(folds == j)]))
        prefict_incomplete_obs = cbind(matrix(0, nrow = nrow(X_incomplete), ncol = len_miss), cbind(X_incomplete, Y_incomplete))
        colnames(prefict_complete_obs) = colnames
        colnames(prefict_incomplete_obs) = colnames
        mu_hat_C[which(folds == j), n_miss] =  predict(mod, newdata = prefict_complete_obs)
        mu_hat_IC[, n_miss] = mu_hat_IC[, n_miss] + predict(mod, newdata = prefict_incomplete_obs)
        mu_hat_monte_C[, n_miss] = mu_hat_monte_C[, n_miss] + predict(mod, newdata = monte_data)
      }
    }
    mu_hat_IC = mu_hat_IC / 5
    mu_hat_monte_C = mu_hat_monte_C / 5
    resi = X_complete[, miss_source[[i]]] - mu_hat_C
    Cov = t(resi) %*% resi / nrow(X_complete)
    phi_1_C = t(Phi_1(X_complete[, -miss_source[[i]]], Y_complete, 
                      matrix(gamma.tilt[-miss_source[[i]]], ncol = 1), 
                      H_inv, gamma.tilt[miss_source[[i]]], mu_hat_C, Cov, 
                      obs_index = c(1:ncol(S))[-miss_source[[i]]], miss_index = miss_source[[i]]))
    phi_1_IC = t(Phi_1(X_incomplete, Y_incomplete, 
                       matrix(gamma.tilt[-miss_source[[i]]], ncol = 1), 
                       H_inv, gamma.tilt[miss_source[[i]]], mu_hat_IC, Cov, 
                       obs_index = c(1:ncol(S))[-miss_source[[i]]], miss_index = miss_source[[i]]))
    phi_1_C_list[[i - 1]] = phi_1_C
    phi_1_IC_list[[i - 1]] = phi_1_IC
    phi_1_C0_list[[i - 1]] = phi_1_C - t(matrix(rep(colMeans(phi_1_C), n[1]), nrow = ncol(phi_1_C)))
    phi_1_IC0_list[[i - 1]] = phi_1_IC - t(matrix(rep(colMeans(phi_1_IC), n[i]), nrow = ncol(phi_1_IC)))
    phi_1_monte = t(Phi_1(monte_X[, setdiff(1:pppx, miss_source[[2]])], monte_y, 
                          matrix(gamma.tilt[setdiff(1:pppx, miss_source[[2]])], ncol = 1), H_inv, gamma.tilt[miss_source[[2]]], 
                          mu_hat_monte_C, Cov, obs_index = c(1:ncol(S))[-miss_source[[i]]], miss_index = miss_source[[i]]))
    phi_1_monte_list[[i - 1]] = phi_1_monte
    phi_1_monte0_list[[i - 1]] = phi_1_monte - t(matrix(rep(colMeans(phi_1_monte), nrow(phi_1_monte)), nrow = ncol(phi_1_monte)))
  }
  method1 = "intercept"
  method2 = "linear"
  X_C = Combined_X[1:n[1], ]
  X_IC = Combined_X[(n[1] + 1):(n[1] + n[2]), ]
  X_U = Combined_X[(n[1] + n[2] + 1):(n[1] + n[2] + n[3]), ]
  Y_C = Combined_Y[1:n[1]]
  Y_IC = Combined_Y[(n[1] + 1):(n[1] + n[2])]
  miss = miss_source[[2]]
  allocation_intercept = allocation_real(method1, phi_1_C, phi_1_IC, S, cbind(X_C[, -miss], Y_C), cbind(X_IC[, -miss], Y_IC))
  allocation_linear = allocation_real(method2, phi_1_C, phi_1_IC, S, cbind(X_C[, -miss], Y_C), cbind(X_IC[, -miss], Y_IC))
  weight_intercept_C = matrix(0, ncol = ncol(phi_1_C), nrow = nrow(phi_1_C))
  weight_intercept_IC = matrix(0, ncol = ncol(phi_1_IC), nrow = nrow(phi_1_IC))
  weight_intercept_U = matrix(0, ncol = ncol(phi_1_C), nrow = n[3])
  weight_linear_C = matrix(0, ncol = ncol(phi_1_C), nrow = nrow(phi_1_C))
  weight_linear_IC = matrix(0, ncol = ncol(phi_1_IC), nrow = nrow(phi_1_IC))
  weight_linear_U = matrix(0, ncol = ncol(phi_1_C), nrow = n[3])
  for (nn in 1:ncol(S)) {
    weight_intercept_C[, nn] = (rep(allocation_intercept[[1]][[nn]], n[1]) + rep(allocation_intercept[[2]][[nn]], n[1])) / 2 + 1
    weight_intercept_IC[, nn] = (rep(allocation_intercept[[1]][[nn]], n[2]) + rep(allocation_intercept[[2]][[nn]], n[2])) / 2 + 1
    weight_linear_C[, nn] = (as.matrix(cbind(X_C[, -miss], Y_C, 1)) %*% as.matrix(allocation_linear[[1]][[nn]]) + 
                               as.matrix(cbind(X_C[, -miss], Y_C, 1)) %*% as.matrix(allocation_linear[[2]][[nn]])) / 2 + 1
    weight_linear_IC[, nn] = (as.matrix(cbind(X_IC[, -miss], Y_IC, 1)) %*% as.matrix(allocation_linear[[1]][[nn]]) + 
                                as.matrix(cbind(X_IC[, -miss], Y_IC, 1)) %*% as.matrix(allocation_linear[[2]][[nn]])) / 2 + 1
  }
  ad_phi_1_C_intercept = weight_intercept_C * phi_1_C
  ad_phi_1_C_linear = weight_linear_C * phi_1_C
  ad_phi_1_IC_intercept = weight_intercept_IC * phi_1_IC
  ad_phi_1_IC_linear = weight_linear_IC * phi_1_IC
  gamma_1_intercept = gamma.tilt - colMeans(ad_phi_1_C_intercept) + colMeans(rbind(ad_phi_1_C_intercept, ad_phi_1_IC_intercept))
  gamma_1_linear = gamma.tilt - colMeans(ad_phi_1_C_linear) + colMeans(rbind(ad_phi_1_C_linear, ad_phi_1_IC_linear))
  gamma_1_semi = gamma.tilt - colMeans(phi_1_C) + colMeans(rbind(phi_1_C, phi_1_IC))  
  phi_2 = apply(phi_1_monte_list[[1]], 2, function(col) tapply(col, rep(1:n_complete_X, each = nrow(phi_1_monte_list[[1]])/n_complete_X), mean))
  S_CU_intetcept = S - n[2] / (n[1] + n[2]) * ad_phi_1_C_intercept
  S_CU_linear = S - n[2] / (n[1] + n[2]) * ad_phi_1_C_linear
  allocation_intercept_2 = allocation_real(method1, phi_2[1:n[1], ], phi_2[(n[1] + 1):(n[1] + n[3]), ], S_CU_intetcept, X_C, X_U)
  allocation_linear_2 = allocation_real(method2, phi_2[1:n[1], ], phi_2[(n[1] + 1):(n[1] + n[3]), ], S_CU_linear, X_C, X_U)
  for (nn in 1:ncol(S)) {
    weight_intercept_C[, nn] = (rep(allocation_intercept_2[[1]][[nn]], n[1]) + rep(allocation_intercept_2[[2]][[nn]], n[1])) / 2 + 1
    weight_intercept_U[, nn] = (rep(allocation_intercept_2[[1]][[nn]], n[3]) + rep(allocation_intercept_2[[2]][[nn]], n[3])) / 2 + 1
    weight_linear_C[, nn] = (as.matrix(cbind(X_C, 1)) %*% as.matrix(allocation_linear_2[[1]][[nn]]) + 
                               as.matrix(cbind(X_C, 1)) %*% as.matrix(allocation_linear_2[[2]][[nn]])) / 2 + 1
    weight_linear_U[, nn] = (as.matrix(cbind(X_U, 1)) %*% as.matrix(allocation_linear_2[[1]][[nn]]) + 
                               as.matrix(cbind(X_U, 1)) %*% as.matrix(allocation_linear_2[[2]][[nn]])) / 2 + 1
  }
  ad_phi_2_intercept = rbind(weight_intercept_C, weight_intercept_U) * phi_2
  ad_phi_2_linear = rbind(weight_linear_C, weight_linear_U) * phi_2
  gamma_2_intercept = gamma_1_intercept - colMeans(ad_phi_2_intercept[1:n[1], ]) + colMeans(ad_phi_2_intercept)
  gamma_2_linear = gamma_1_linear - colMeans(ad_phi_2_linear[1:n[1], ]) + colMeans(ad_phi_2_linear)
  gamma_2_semi = gamma_1_semi - colMeans(phi_2[1:n[1], ]) + colMeans(phi_2)
  var_semi = var_cal2(S, phi_1_C, phi_1_IC, phi_2[1:n[1], ], phi_2[(n[1] + 1):(n[1] + n[3]), ], correct = "real")
  var_intercept = var_cal2(S, ad_phi_1_C_intercept, ad_phi_1_IC_intercept, ad_phi_2_intercept[1:n[1], ], ad_phi_2_intercept[(n[1] + 1):(n[1] + n[3]), ], correct = "real")
  var_linear = var_cal2(S, ad_phi_1_C_linear, ad_phi_1_IC_linear, ad_phi_2_linear[1:n[1], ], ad_phi_2_linear[(n[1] + 1):(n[1] + n[3]), ], correct = "real")
  res = list(gamma_1_intercept = gamma_1_intercept, gamma_2_intercept = gamma_2_intercept, 
             gamma_1_linear = gamma_1_linear, gamma_2_linear = gamma_2_linear, 
             gamma_1_semi = gamma_1_semi, gamma_2_semi = gamma_2_semi, 
             naive_var = var_intercept$var0,
             var_intercept1 = var_intercept$var1,
             var_intercept2 = var_intercept$var21,
             var_linear1 = var_linear$var1, 
             var_linear2 = var_linear$var21, 
             var_semi1 = var_semi$var1,
             var_semi2 = var_semi$var21)
  return(res)
}
