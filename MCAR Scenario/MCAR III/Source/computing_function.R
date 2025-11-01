# g function for binomial 
g_x = function(x, beta){
  xbeta = x %*% beta
  return(xbeta)
}
# g_dot function
g_dot <- function(x,beta){
  return(1)
}
# hessian function for logistic
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
Phi_2 = function(x_obs, gamma_obs, H_inv, x_miss, gamma_miss, cov_miss, eta, sigma_y, imput_method) {
  x_obs = matrix(x_obs, ncol = 4)
  gamma_obs = matrix(gamma_obs, ncol = 1)
  x_miss = matrix(x_miss, ncol = 2)
  gamma_miss = matrix(gamma_miss, ncol = 1)
  eta = matrix(eta, ncol = 2)
  mu_y = x_obs[, 2:4] %*% gamma_obs + x_miss %*% gamma_miss
  Ey1 = mu_y
  Ey2 = mu_y ^ 2 + sigma_y ^ 2
  Ey3 = mu_y ^ 3 + 3 * mu_y * sigma_y ^ 2
  Ey4 = mu_y ^ 4 + 6 * mu_y ^ 2 * sigma_y ^ 2 + 3 * sigma_y ^ 4 
  higher_complete_data = matrix(cbind(x_obs, mu_y, 
                                      x_obs[, 1] * x_obs[, 2], x_obs[, 1] * x_obs[, 3], x_obs[, 1] * mu_y,
                                      x_obs[, 2] * x_obs[, 3], x_obs[, 2] * mu_y, x_obs[, 3] * mu_y,
                                      x_obs[, 1:3]^2, mu_y^2), ncol = 15)
  complete_data = matrix(cbind(x_obs, mu_y), ncol = 5)
  if (imput_method == "lm") {
    mu_miss = complete_data %*% eta
  } else if (imput_method == "higher_lm") {
    cons = higher_complete_data[,c(1, 2, 3, 4, 6, 7, 9, 12, 13, 14)] %*% eta[c(1, 2, 3, 4, 6, 7, 9, 12, 13, 14),]
    coef_y = higher_complete_data[,c(1, 2, 3, 4)] %*% eta[c(5, 8, 10, 11),]
    coef_yy = higher_complete_data[,1] %*% matrix(eta[15,], nrow = 1)
    mu_miss = higher_complete_data %*% eta
  } else if (imput_method == "ml") {
    yyy = 0
  } else if (imput_method == "enet") {
    yyy = 0
  }
  tem1 = c()
  tem2 = c()
  for (i in 1:nrow(x_obs)) {
    tem_x_obs = matrix(x_obs[i, 2:4], ncol = 1)
    tem_gamma_obs = matrix(gamma_obs, ncol = 1)
    tem_mu_miss = matrix(mu_miss[i, ], ncol = 1)
    tem_gamma_miss = matrix(gamma_miss, ncol = 1)
    if (imput_method == "lm") {
      tem1 = cbind(tem1, tem_x_obs * mu_y[i] - tem_x_obs %*% t(tem_x_obs) %*% tem_gamma_obs - tem_x_obs %*% t(tem_mu_miss) %*% tem_gamma_miss)
      q1 = tem_mu_miss * mu_y[i] + matrix(eta[5,], ncol = 1) * sigma_y ^ 2 
      q2 = (cov_miss + tem_mu_miss %*% t(tem_mu_miss) +  
              matrix(eta[5,], ncol = 1) %*% matrix(eta[5,], nrow = 1) * sigma_y ^ 2 ) %*% tem_gamma_miss
      tem2 = cbind(tem2, q1 - tem_mu_miss %*% t(tem_x_obs) %*% tem_gamma_obs - q2)
    } else if (imput_method == "higher_lm") {
      tem_cons = matrix(cons[i, ], ncol = 1)
      tem_coef_y = matrix(coef_y[i, ], ncol = 1)
      tem_coef_yy = matrix(coef_yy[i, ], ncol = 1)
      tem_higher_mu_miss = (tem_coef_yy * Ey2[i] + tem_coef_y * Ey1[i])
      tem1 = cbind(tem1, tem_x_obs * mu_y[i] - tem_x_obs %*% t(tem_x_obs) %*% tem_gamma_obs - tem_x_obs %*% t(tem_higher_mu_miss) %*% tem_gamma_miss)
      q1 = tem_cons * Ey1[i] + tem_coef_y * Ey2[i] + tem_coef_yy * Ey3[i]
      partial_q2 = tem_coef_yy %*% t(tem_coef_yy) * Ey4[i] + 
        2 * tem_coef_yy %*% t(tem_coef_y) * Ey3[i] + 
        (tem_coef_y %*% t(tem_coef_y) + 2 * tem_coef_yy %*% t(tem_cons)) * Ey2[i] +
        2 * tem_coef_y %*% t(tem_cons) * Ey1[i] + 
        tem_cons %*% t(tem_cons)
      q2 = (cov_miss + partial_q2) %*% tem_gamma_miss 
      mid = tem_higher_mu_miss %*% t(tem_x_obs) %*% tem_gamma_obs 
      tem2 = cbind(tem2, q1 - mid - q2)
    } else if (imput_method == "ml") {
      yyy = 0
    } else if (imput_method == "enet") {
      yyy = 0
    }
  }
  res = H_inv %*% matrix(rbind(tem1, tem2), ncol = nrow(x_obs))
  return(res)
}
