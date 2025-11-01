g_x = function(x, beta){
  xbeta = x %*% beta
  return(xbeta)
}
g_dot <- function(x,beta){
  return(1)
}
H <- function(x, beta){
  h = matrix(0, nrow = ncol(x), ncol = ncol(x))
  deriv = g_dot(x,beta)
  for (i in 1:nrow(x)) {
    h = h + deriv*x[i,] %*% t(x[i,])
  }
  return(h/nrow(x))
}
Phi_1 = function(x_obs, y, gamma_obs, H_inv, gamma_miss, mu_miss, cov_miss, obs_index, miss_index) {
  tem_x_obs <- as.matrix(x_obs)
  tem_gamma_obs <- matrix(gamma_obs, ncol = 1)
  tem_mu_miss <- as.matrix(mu_miss)
  tem_gamma_miss <- matrix(gamma_miss, ncol = 1)
  tem1 = as.numeric(tem_x_obs %*% tem_gamma_obs)
  tem2 = as.numeric(tem_mu_miss %*% tem_gamma_miss)
  tem3 = cov_miss %*% tem_gamma_miss
  temy = as.numeric(y)
  ncolx = length(obs_index) + length(miss_index)
  ny = length(temy)
  combine = matrix(0, ncol = ncolx, nrow = ny)
  int1 = tem_x_obs * (temy - tem1 - tem2)
  int2 = tem_mu_miss * (temy - tem1 - tem2) - t(matrix(as.numeric(tem3), nrow = ncol(tem_mu_miss) , ncol = nrow(tem_mu_miss)))
  combine[, obs_index] = int1
  combine[, miss_index] = int2
  res <- H_inv %*% t(combine)
  return(res)
}
## calibration
allocation_real = function(method, phi1, phi2, S, X1_used, X2_used) {
  set.seed(123)
  n_C = nrow(phi1)
  n_IC = nrow(phi2)
  sam_id = sample(1:n_C, size = round(n_C / 2))
  rho = n_IC / n_C
  phi1prime = rho / (1 + rho) * phi1
  S_resi_C = S - phi1prime
  S_resi2 = S / phi1prime - matrix(1, ncol = ncol(S), nrow = nrow(S))
  if (method == "intercept") {
    XX_C = matrix(1, nrow = n_C, ncol = 1)
    XX_IC = matrix(1, nrow = n_IC, ncol = 1)
  } else {
    XX_C = as.matrix(cbind(X1_used, 1))
    XX_IC = as.matrix(cbind(X2_used, 1))
  }
  model_list = vector("list", 2)
  model_list[[1]] = vector("list", ncol(S))
  model_list[[2]] = vector("list", ncol(S))
  for (k in 1:ncol(S)) {
    if (method == "intercept") { # example 1
      temXX_C = phi1prime[, k] * XX_C
      model_list[[1]][[k]] <- solve(var(temXX_C[sam_id, ]) + 1 / rho * var(temXX_C[sam_id, ])) %*% 
        (cov(temXX_C[sam_id, ], S_resi_C[sam_id, k]) - 1 / rho * cov(temXX_C[sam_id, ], phi1prime[sam_id, k]))
      model_list[[2]][[k]] <- solve(var(temXX_C[-sam_id, ]) + 1 / rho * var(temXX_C[-sam_id, ])) %*% 
        (cov(temXX_C[-sam_id, ], S_resi_C[-sam_id, k]) - 1 / rho * cov(temXX_C[-sam_id, ], phi1prime[-sam_id, k]))
    } else if (method == "linear") { # example 2
      temXX_C = phi1[, k] * XX_C
      model_list[[1]][[k]] <- solve(cov(temXX_C[sam_id, ]) + 1 / rho * cov(temXX_C[sam_id, ])) %*% 
        (cov(temXX_C[sam_id, ], S_resi_C[sam_id, k]) - 1 / rho * cov(temXX_C[sam_id, ], phi1prime[sam_id, k]))
      model_list[[2]][[k]] <- solve(cov(temXX_C[-sam_id, ]) + 1 / rho * cov(temXX_C[-sam_id, ])) %*% 
        (cov(temXX_C[-sam_id, ], S_resi_C[-sam_id, k]) - 1 / rho * cov(temXX_C[-sam_id, ], phi1prime[-sam_id, k]))
    } else if (method == "nonlinear") { 
      temXX_C = phi1[, k] * cbind(XX_C, sqrt(abs(XX_C[, 3])), XX_C[, 4] / XX_C[, 5] ^ 2)
      model_list[[1]][[k]] <- solve(cov(temXX_C[sam_id, ]) + 1 / rho * cov(temXX_C[sam_id, ])) %*% 
        (cov(temXX_C[sam_id, ], S_resi_C[sam_id, k]) - 1 / rho * cov(temXX_C[sam_id, ], phi1prime[sam_id, k]))
      model_list[[2]][[k]] <- solve(cov(temXX_C[-sam_id, ]) + 1 / rho * cov(temXX_C[-sam_id, ])) %*% 
        (cov(temXX_C[-sam_id, ], S_resi_C[-sam_id, k]) - 1 / rho * cov(temXX_C[-sam_id, ], phi1prime[-sam_id, k]))
    } else if (method == "KRR") { # example 3
      model_list[[1]][[k]] <- suppressWarnings(ksvm(S_resi2[sam_id, k] ~ XX_C[sam_id, ], kernel = "rbfdot", C = 1, kpar = list(sigma = 0.05)))
      model_list[[2]][[k]] <- suppressWarnings(ksvm(S_resi2[-sam_id, k] ~ XX_C[-sam_id, ], kernel = "rbfdot", C = 1, kpar = list(sigma = 0.05)))
    }
  }
  return(model_list)
}
## asymptotic variance
var_cal2 = function(S, phi_1_C, phi_1_IC, phi_2_C, phi_2_U, correct) {
  n1 = nrow(S)
  n2 = nrow(phi_1_IC)
  n3 = nrow(phi_2_U)
  if (correct == T) {
    var0 = colVars(S)
    var1 = colVars(S) - n2 / (n1 + n2) * colVars(phi_1_C)
    var20 = var1
    var21 = colVars(S) - n2 / (n1 + n2) * colVars(phi_1_C) - n3 / (n1 + n3) * colVars(phi_2_C)
  } else if (correct == "real") {
    var0 = colVars(S)
    var1 = colVars(S - n2 / (n1 + n2) * phi_1_C) + (n1 * n2) / (n1 + n2) ^ 2 * colVars(phi_1_C)
    var20 = var1
    var21 = colVars(S - n2 / (n1 + n2) * phi_1_C - n3 / (n1 + n3) * phi_2_C) + (n1 * n2) / (n1 + n2) ^ 2 * colVars(phi_1_C) + 
      (n1 * n3) / (n1 + n3) ^ 2 * colVars(phi_2_C)
  } else if (correct == "real2") {
    var0 = colVars(S)
    var1 = colVars(S - n2 / (n1 + n2) * phi_1_C) + (n1 * n2) / (n1 + n2) ^ 2 * colVars(phi_1_IC)
    var20 = var1
    var21 = colVars(S - n2 / (n1 + n2) * phi_1_C - n3 / (n1 + n3) * phi_2_C) + (n1 * n2) / (n1 + n2) ^ 2 * colVars(phi_1_IC) + 
      (n1 * n3) / (n1 + n3) ^ 2 * colVars(phi_2_U)
  }
  return(list(var0 = var0, var1 = var1, var20 = var20, var21 = var21))
}
