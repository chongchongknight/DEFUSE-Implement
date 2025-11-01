main_func <- function(X, Y, n, N, M, folds, imput_method, gamma.tilt, H_inv, S, sigma_y, gamma.true, nlm, seed){
  set.seed(seed)
  mod4 <- vector("list",5)
  mod5 <- vector("list",5)
  mu_hat_C4 = matrix(0, nrow = n[1], ncol = 1)
  mu_hat_C5 = matrix(0, nrow = n[1], ncol = 1)
  mu_hat_IC4 = matrix(0, nrow = nlm[1], ncol = 5)
  mu_hat_IC5 = matrix(0, nrow = nlm[2], ncol = 5)
  monte_size = 100
  mu_hat_CU4 = matrix(0, nrow = N[1] * monte_size, ncol = 5)
  mu_hat_CU5 = matrix(0, nrow = N[1] * monte_size, ncol = 5)
  ctrl <- trainControl(method = "cv", number = 5)
  complete_data = data.frame(cbind(X[[1]][1:n[1], ], Y[[1]][1:n[1]])) 
  colnames(complete_data) = c("X1", "X2", "X3", "X4", "X5", "Y")
  incomplete_data = data.frame(cbind(X[[2]][1:n[2], ], Y[[2]][1:n[2]])) 
  colnames(incomplete_data) = c("X1", "X2", "X3", "X4", "X5", "Y")
  monte_sample = rnorm(monte_size)
  mu_Y = X[[1]] %*% gamma.tilt
  resi_Y = Y[[1]] - mu_Y
  sigma_Y = sqrt(mean(resi_Y^2))
  monte_Y = rep(mu_Y, each = monte_size) + rep(monte_sample, N[1]) * sigma_Y
  monte_Xobs =  X[[1]][rep(1:nrow(X[[1]]), each = monte_size), ]
  monte_Y_Xobs = data.frame(cbind(monte_Xobs, monte_Y)) 
  colnames(monte_Y_Xobs) = c("X1", "X2", "X3", "X4", "X5", "Y")
  for (i in 1:5) {
    mod4[[i]] = train(X4 ~ X1 + X2 + X3 + X5 + Y, data = complete_data[which(folds != i), ], method = imput_method, trControl = ctrl)
    mod5[[i]] = train(X5 ~ X1 + X2 + X3 + X4 + Y, data = complete_data[which(folds != i), ], method = imput_method, trControl = ctrl)
    mu_hat_C4[which(folds == i), 1] = predict(mod4[[i]], newdata = complete_data[which(folds == i), ])
    mu_hat_C5[which(folds == i), 1] = predict(mod5[[i]], newdata = complete_data[which(folds == i), ])
    mu_hat_IC4[, i] = predict(mod4[[i]], newdata = incomplete_data)[1:nlm[1]]
    mu_hat_IC5[, i] = predict(mod5[[i]], newdata = incomplete_data)[(nlm[1] + 1):(nlm[1] + nlm[2])]
    mu_hat_CU4[, i] = predict(mod4[[i]], newdata = monte_Y_Xobs)
    mu_hat_CU5[, i] = predict(mod5[[i]], newdata = monte_Y_Xobs)
  }
  mu_hat_IC4 = rowMeans(mu_hat_IC4)
  mu_hat_IC5 = rowMeans(mu_hat_IC5)
  mu_hat_CU4 = rowMeans(mu_hat_CU4) 
  mu_hat_CU5 = rowMeans(mu_hat_CU5)
  resi4 = X[[1]][1:n[1],M[[2]][1]] - mu_hat_C4
  resi5 = X[[1]][1:n[1],M[[2]][2]] - mu_hat_C5
  Cov4 = mean(resi4 ^ 2)
  Cov5 = mean(resi5 ^ 2)
  phi_1_C4 = t(Phi_1(X[[1]][1:n[1],c(1:3, 5)], Y[[1]], matrix(gamma.tilt[c(1:3, 5)], ncol = 1), H_inv, gamma.tilt[4], 
                     matrix(c(mu_hat_C4), ncol = 1), Cov4, obs_index = c(1:3, 5), miss_index = 4))
  phi_1_C5 = t(Phi_1(X[[1]][1:n[1],c(1:3, 4)], Y[[1]], matrix(gamma.tilt[c(1:3, 4)], ncol = 1), H_inv, gamma.tilt[5], 
                     matrix(c(mu_hat_C5), ncol = 1), Cov5, obs_index = c(1:3, 4), miss_index = 5))
  phi_1_IC4 = t(Phi_1(X[[2]][1:nlm[1],c(1:3, 5)], Y[[2]][1:nlm[1]], matrix(gamma.tilt[c(1:3, 5)], ncol = 1), H_inv, gamma.tilt[4], 
                      matrix(c(mu_hat_IC4), ncol = 1), Cov4, obs_index = c(1:3, 5), miss_index = 4))
  phi_1_IC5 = t(Phi_1(X[[2]][(nlm[1] + 1):(nlm[1] + nlm[2]),c(1:3, 4)], 
                      Y[[2]][(nlm[1] + 1):(nlm[1] + nlm[2])], 
                      matrix(gamma.tilt[c(1:3, 4)], ncol = 1), H_inv, gamma.tilt[5],
                      matrix(c(mu_hat_IC5), ncol = 1), Cov5, obs_index = c(1:4), miss_index = 5))
  phi_1_C4_0 = phi_1_C4 -  t(matrix(rep(colMeans(phi_1_C4), n[1]), nrow = 5))
  phi_1_C5_0 = phi_1_C5 - t(matrix(rep(colMeans(phi_1_C5), n[1]), nrow = 5))
  phi_1_IC4_0 = phi_1_IC4 - t(matrix(rep(colMeans(phi_1_IC4), nlm[1]), nrow = 5))
  phi_1_IC5_0 = phi_1_IC5 - t(matrix(rep(colMeans(phi_1_IC5), nlm[2]), nrow = 5))
  alpha_1 <- matrix(0, ncol = ncol(S), nrow = 2)
  rho = 1
  for (j in 1:ncol(S)) {
    denominator = matrix(
      c(mean(phi_1_C4_0[, j] ^ 2) + rho * mean(phi_1_IC4_0[, j] ^ 2), 
        mean(phi_1_C4_0[, j] * phi_1_C5_0[, j]),
        mean(phi_1_C5_0[, j] * phi_1_C4_0[, j]),
        mean(phi_1_C5_0[, j] ^ 2) + rho * mean(phi_1_IC5_0[, j] ^ 2)
      ), 
      ncol = 2, nrow = 2)
    numerator = matrix(c(mean(phi_1_C4_0[, j] * S[, j]), 
                         mean(phi_1_C5_0[, j] * S[, j])),
                       ncol = 1)
    alpha_1[, j] = solve(denominator) %*% numerator
  }
  gamma_tilt_1 = gamma.tilt - 
    colMeans(t(diag(alpha_1[1, ]) %*% t(phi_1_C4))) - 
    colMeans(t(diag(alpha_1[2, ]) %*% t(phi_1_C5))) +
    colMeans(t(diag(alpha_1[1, ]) %*% t(phi_1_IC4))) + 
    colMeans(t(diag(alpha_1[2, ]) %*% t(phi_1_IC5))) 
  phi_1_CU4 = t(Phi_1(monte_Xobs[ ,c(1:3, 5)], monte_Y, matrix(gamma.tilt[c(1:3, 5)], ncol = 1), H_inv, gamma.tilt[4], 
                      matrix(c(mu_hat_CU4), ncol = 1), Cov4, obs_index = c(1:3, 5), miss_index = 4))
  phi_1_CU5 = t(Phi_1(monte_Xobs[,c(1:3, 4)], monte_Y, matrix(gamma.tilt[c(1:3, 4)], ncol = 1), H_inv, gamma.tilt[5], 
                      matrix(c(mu_hat_CU5), ncol = 1), Cov5, obs_index = c(1:3, 4), miss_index = 5))
  phi_2_CU4 = apply(phi_1_CU4, 2, function(col) tapply(col, rep(1:N[1], each = nrow(phi_1_CU4)/N[1]), mean))
  phi_2_CU5 = apply(phi_1_CU5, 2, function(col) tapply(col, rep(1:N[1], each = nrow(phi_1_CU5)/N[1]), mean))
  phi_2 = t(diag(alpha_1[1, ]) %*% t(phi_2_CU4)) + t(diag(alpha_1[2, ]) %*% t(phi_2_CU5))
  phi_1 = t(diag(alpha_1[1, ]) %*% t(phi_1_C4)) + t(diag(alpha_1[2, ]) %*% t(phi_1_C4))
  R = S - phi_1
  phi_2_0 = phi_2 -  t(matrix(rep(colMeans(phi_2), N[1]), nrow = 5))
  alpha_2 = matrix(0, ncol = ncol(S), nrow = 1)
  rho2 = 0.5
  for (j in 1:ncol(S)) {
    denominator = mean((phi_2_0[1:n[1],j]) ^ 2) + rho2 * mean((phi_2_0[(n[1] + 1):(N[1]),j]) ^ 2)
    alpha_2[j] = mean(phi_2_0[1:n[1],j] * R[, j]) 
  }
  alpha_2 <- diag(as.numeric(alpha_2))
  gamma_tilt_2 = gamma_tilt_1 - colMeans(t(alpha_2 %*% t(phi_2[1:n[1], ]))) + 
    colMeans(t(alpha_2 %*% t(phi_2[(n[1] + 1):N[1], ])))
  R2 = R - t(alpha_2 %*% t(phi_2[1:n[1], ]))
  var = colMeans(S ^ 2) / n[1]
  var1 = colMeans(R ^ 2) / n[1] + 
    colMeans((t(diag(alpha_1[1, ]) %*% t(phi_1_IC4))) ^ 2) / nlm[1] + 
    colMeans((t(diag(alpha_1[2, ]) %*% t(phi_1_IC5))) ^ 2) / nlm[2]
  var2 = colMeans(R2 ^ 2) / n[1] + 
    colMeans((t(diag(alpha_1[1, ]) %*% t(phi_1_IC4))) ^ 2) / nlm[1] + 
    colMeans((t(diag(alpha_1[2, ]) %*% t(phi_1_IC5))) ^ 2) / nlm[2] +
    colMeans((t(alpha_2 %*% t(phi_2[(n[1] + 1):N[1], ]))) ^ 2) / (N[1] - n[1])
  return(list(gamma.tilt2 = gamma_tilt_2, gamma.tilt1 = gamma_tilt_1, tilt = gamma.tilt, 
              se2 = sqrt(var2), se1 = sqrt(var1), se = sqrt(var)))
}
