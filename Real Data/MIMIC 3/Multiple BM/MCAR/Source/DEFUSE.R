# Get the proposed estimators
sip_dop_get <- function(seed){
  
  set.seed(seed)
  naive = glm(Y_lc_train ~., cbind(Y_lc_train, X_lc_train),family = "binomial")
  gamma.tilt = naive$coefficients
  
  # Step 2: Estimation of Score Function
  # get the score function
  H_inv = solve(H(as.matrix(cbind(1,X_lc_train)), gamma.tilt))
  S = sweep(t(H_inv %*% t(as.matrix(cbind(1,X_lc_train)))), MARGIN=1, (Y_lc_train - g_x(as.matrix(cbind(1,X_lc_train)), gamma.tilt)), '*')
  
  mu_hat_C5 = c()
  mu_hat_C4 = c()
  mu_hat_C45_4 = c()
  mu_hat_C45_5 = c()
  
  mu_hat_IC5 = matrix(0, nrow = nlm[1], ncol = 5)
  mu_hat_IC4 = matrix(0, nrow = nlm[2], ncol = 5)
  mu_hat_IC45_4 = matrix(0, nrow = nlm[3], ncol = 5)
  mu_hat_IC45_5 = matrix(0, nrow = nlm[3], ncol = 5)
  
  mu_hat_IC5 = matrix(0, nrow = nlm[1], ncol = 5)
  mu_hat_IC4 = matrix(0, nrow = nlm[2], ncol = 5)
  mu_hat_IC45_4 = matrix(0, nrow = nlm[3], ncol = 5)
  mu_hat_IC45_5 = matrix(0, nrow = nlm[3], ncol = 5)
  
  mu_hat_CU5_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU4_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_4_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_5_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  
  mu_hat_CU5_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU4_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_4_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_5_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  
  for (i in 1:5) {
    # mean estimation
    mu_hat_C5[which((Y_lc_train == 0) & (folds == i))] = predict(mod5[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[1]]])
    mu_hat_C5[which((Y_lc_train == 1) & (folds == i))] = predict(mod5[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[1]]])
    
    mu_hat_C4[which((Y_lc_train == 0) & (folds == i))] = predict(mod4[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[2]]])
    mu_hat_C4[which((Y_lc_train == 1) & (folds == i))] = predict(mod4[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[2]]])

    mu_hat_C45_4[which((Y_lc_train == 0) & (folds == i))] = predict(mod45_4[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[3]]])
    mu_hat_C45_4[which((Y_lc_train == 1) & (folds == i))] = predict(mod45_4[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[3]]])
    
    mu_hat_C45_5[which((Y_lc_train == 0) & (folds == i))] = predict(mod45_5[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[3]]])
    mu_hat_C45_5[which((Y_lc_train == 1) & (folds == i))] = predict(mod45_5[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[3]]])
    
    mu_hat_IC5[which(Y_lm1_train == 0),i] = predict(mod5[[i]][[1]], newdata = X_lm1_train_withauxi[which(Y_lm1_train == 0),-M[[1]]]) 
    mu_hat_IC5[which(Y_lm1_train == 1),i] = predict(mod5[[i]][[2]], newdata = X_lm1_train_withauxi[which(Y_lm1_train == 1),-M[[1]]])
    
    mu_hat_IC4[which(Y_lm2_train == 0),i] = predict(mod4[[i]][[1]], newdata = X_lm2_train_withauxi[which(Y_lm2_train == 0),-M[[2]]]) 
    mu_hat_IC4[which(Y_lm2_train == 1),i] = predict(mod4[[i]][[2]], newdata = X_lm2_train_withauxi[which(Y_lm2_train == 1),-M[[2]]])
    
    mu_hat_IC45_4[which(Y_lm3_train == 0),i] = predict(mod45_4[[i]][[1]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 0),-M[[3]]]) 
    mu_hat_IC45_4[which(Y_lm3_train == 1),i] = predict(mod45_4[[i]][[2]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 1),-M[[3]]])
    
    mu_hat_IC45_5[which(Y_lm3_train == 0),i] = predict(mod45_5[[i]][[1]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 0),-M[[3]]]) 
    mu_hat_IC45_5[which(Y_lm3_train == 1),i] = predict(mod45_5[[i]][[2]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 1),-M[[3]]])
    
    ## unlabel
    mu_hat_CU5_Y0[ ,i] = predict(mod5[[i]][[1]], newdata = X_CU_train_withauxi[,-M[[1]]]) 
    mu_hat_CU5_Y1[ ,i] = predict(mod5[[i]][[2]], newdata = X_CU_train_withauxi[,-M[[1]]])
    
    mu_hat_CU4_Y0[ ,i] = predict(mod4[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[2]]]) 
    mu_hat_CU4_Y1[ ,i] = predict(mod4[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[2]]])
    
    mu_hat_CU45_4_Y0[ ,i] = predict(mod45_4[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[3]]]) 
    mu_hat_CU45_4_Y1[ ,i] = predict(mod45_4[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[3]]])
    
    mu_hat_CU45_5_Y0[ ,i] = predict(mod45_5[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[3]]]) 
    mu_hat_CU45_5_Y1[ ,i] = predict(mod45_5[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[3]]])
  }
  mu_hat_IC5 = rowMeans(mu_hat_IC5)
  mu_hat_IC4 = rowMeans(mu_hat_IC4) # take the mean values of fitted values from 5 models
  
  mu_hat_IC45_4 = rowMeans(mu_hat_IC45_4)
  mu_hat_IC45_5 = rowMeans(mu_hat_IC45_5)
  
  mu_hat_CU5_Y0 = rowMeans(mu_hat_CU5_Y0)
  mu_hat_CU4_Y0 = rowMeans(mu_hat_CU4_Y0)
  mu_hat_CU45_4_Y0 = rowMeans(mu_hat_CU45_4_Y0)
  mu_hat_CU45_5_Y0 = rowMeans(mu_hat_CU45_5_Y0)
  
  mu_hat_CU5_Y1 = rowMeans(mu_hat_CU5_Y1)
  mu_hat_CU4_Y1 = rowMeans(mu_hat_CU4_Y1)
  mu_hat_CU45_4_Y1 = rowMeans(mu_hat_CU45_4_Y1)
  mu_hat_CU45_5_Y1 = rowMeans(mu_hat_CU45_5_Y1)
  
  # Monte Carlo Integration
  # variance estimation
  resi5 = pull(X_lc_train, M[[1]]) - mu_hat_C5
  resi4 = pull(X_lc_train, M[[2]]) - mu_hat_C4
  
  resi45_4 = pull(X_lc_train, M[[3]][1]) - mu_hat_C45_4
  resi45_5 = pull(X_lc_train, M[[3]][2]) - mu_hat_C45_5
  
  sigma_hat5 = c()
  sigma_hat4 = c()
  
  sigma_hat45_4 = c()
  sigma_hat45_5 = c()
  cor45 = c()
  
  sigma_hat5[which(Y_lc_train==0)] = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5[which(Y_lc_train==1)] = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4[which(Y_lc_train==0)] = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4[which(Y_lc_train==1)] = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4[which(Y_lc_train==0)] = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4[which(Y_lc_train==1)] = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_5[which(Y_lc_train==0)] = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5[which(Y_lc_train==1)] = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45[which(Y_lc_train==0)] = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45[which(Y_lc_train==1)] = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
  sigma_hat5_IC = c() 
  sigma_hat4_IC = c() 
  
  sigma_hat45_4_IC = c()  
  sigma_hat45_5_IC = c()
  cor45_IC = c()
  
  sigma_hat5_IC[which(Y_lm1_train==0)] = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5_IC[which(Y_lm1_train==1)] = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4_IC[which(Y_lm2_train==0)] = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4_IC[which(Y_lm2_train==1)] = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4_IC[which(Y_lm3_train==0)] = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4_IC[which(Y_lm3_train==1)] = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  sigma_hat45_5_IC[which(Y_lm3_train==0)] = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5_IC[which(Y_lm3_train==1)] = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45_IC[which(Y_lm3_train==0)] = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45_IC[which(Y_lm3_train==1)] = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
  # Monte Carlo Integration Process
  monte_sample = 500
  base_sample5 = rnorm(monte_sample, sd = 1)
  base_sample4 = rnorm(monte_sample, sd = 1)
  
  base_sample45_4 = rnorm(monte_sample, sd = 1)
  base_sample45_5 = rnorm(monte_sample, sd = 1)
  
  rows = rep(1:n_lc_train,each = monte_sample)
  obs_C45 = X_lc_train[rows,]
  obs_C4 = X_lc_train[rows,]
  obs_C5 = X_lc_train[rows,]
  
  obs_C45[, M[[3]][1]] = rep(base_sample45_4, n_lc_train) * rep(sigma_hat45_4, each = monte_sample) + 
    rep(mu_hat_C45_4, each=monte_sample)
  obs_C45[, M[[3]][2]] = (rep(base_sample45_4,n_lc_train) * rep(cor45,each = monte_sample) + 
                          rep(base_sample45_5,n_lc_train) * rep(sqrt(1-cor45^2),each=monte_sample)) * 
    rep(sigma_hat45_5,each = monte_sample) + rep(mu_hat_C45_5, each=monte_sample)
  
  obs_C5[, M[[1]]] = rep(base_sample5, n_lc_train) * rep(sigma_hat5, each = monte_sample) + rep(mu_hat_C5, each=monte_sample)
  obs_C4[, M[[2]]] = rep(base_sample4, n_lc_train) * rep(sigma_hat4, each = monte_sample) + rep(mu_hat_C4, each=monte_sample)
  
  obs_C45 <- cbind(1,obs_C45) ## new hessian
  obs_C4 <- cbind(1,obs_C4)
  obs_C5 <- cbind(1,obs_C5)

  S_monte_C45 = t(H_inv %*% t(obs_C45)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C45), gamma.tilt))
  S_monte_C4 = t(H_inv %*% t(obs_C4)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C4), gamma.tilt))
  S_monte_C5 = t(H_inv %*% t(obs_C5)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C5), gamma.tilt))
  
  phi_1_C45 = as.matrix(aggregate(S_monte_C45, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  phi_1_C4 = as.matrix(aggregate(S_monte_C4, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  phi_1_C5 = as.matrix(aggregate(S_monte_C5, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  
  rowsI5 = rep(1:nlm[1],each = monte_sample)
  rowsI4 = rep(1:nlm[2],each = monte_sample)
  rowsI45 = rep(1:nlm[3],each = monte_sample)
  
  obs_IC5 = X_lm1_train[rowsI5,]
  obs_IC4 = X_lm2_train[rowsI4,]
  obs_IC45 = X_lm3_train[rowsI45,]
  
  obs_IC5[, M[[1]]] = rep(base_sample5,nlm[1]) * rep(sigma_hat5_IC, each = monte_sample) + 
    rep(mu_hat_IC5, each=monte_sample)
  obs_IC4[, M[[2]]] = rep(base_sample4,nlm[2]) * rep(sigma_hat4_IC, each = monte_sample) + 
    rep(mu_hat_IC4, each=monte_sample)
  
  obs_IC45[, M[[3]][1]] = rep(base_sample45_4,nlm[3]) * rep(sigma_hat45_4_IC, each = monte_sample) + 
    rep(mu_hat_IC45_4, each=monte_sample) ## same mu but not same variance, so not from true distribution
  obs_IC45[, M[[3]][2]] = (rep(base_sample45_4,nlm[3]) * rep(cor45_IC,each = monte_sample) + 
                             rep(base_sample45_5,nlm[3])*rep(sqrt(1-cor45_IC^2),each=monte_sample)) * rep(sigma_hat45_5_IC, each = monte_sample) + 
    rep(mu_hat_IC45_5, each=monte_sample)

  obs_IC5 <- cbind(1,obs_IC5)
  obs_IC4 <- cbind(1,obs_IC4)
  obs_IC45 <- cbind(1,obs_IC45)
  
  S_monte_IC5 = t(H_inv %*% t(obs_IC5)) * as.numeric(rep(Y_lm1_train,each = monte_sample)-g_x(as.matrix(obs_IC5), gamma.tilt))
  S_monte_IC4 = t(H_inv %*% t(obs_IC4)) * as.numeric(rep(Y_lm2_train,each = monte_sample)-g_x(as.matrix(obs_IC4), gamma.tilt))
  S_monte_IC45 = t(H_inv %*% t(obs_IC45)) * as.numeric(rep(Y_lm3_train,each = monte_sample)-g_x(as.matrix(obs_IC45), gamma.tilt))
  
  phi_1_IC5 = as.matrix(aggregate(S_monte_IC5, list(rep(1:nlm[1], each = monte_sample)), mean)[-1])
  phi_1_IC4 = as.matrix(aggregate(S_monte_IC4, list(rep(1:nlm[2], each = monte_sample)), mean)[-1])
  phi_1_IC45 = as.matrix(aggregate(S_monte_IC45, list(rep(1:nlm[3], each = monte_sample)), mean)[-1])
  
  ### now we allocation the phi_4, phi_4, phi_45
  ## method = "intercept", "linear", "interaction"
  allocation5 = allocation("intercept", phi_1_C5, phi_1_IC5, S, M[[1]], X_lc_train, X_lm1_train, Y_lc_train, Y_lm1_train)
  allocation4 = allocation("intercept", phi_1_C4, phi_1_IC4, S, M[[2]], X_lc_train, X_lm2_train, Y_lc_train, Y_lm2_train)
  allocation45 = allocation("intercept", phi_1_C45, phi_1_IC45, S, M[[3]], X_lc_train, X_lm3_train, Y_lc_train, Y_lm3_train)
  adjusted_phi_1_C5 = allocation5$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_C4 = allocation4$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_C45 = allocation45$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_IC5 = allocation5$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm1_train), ]
  adjusted_phi_1_IC4 = allocation4$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm2_train), ]
  adjusted_phi_1_IC45 = allocation45$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm3_train), ]
  alpha_1 = block_weight(list(phi_1_C5, phi_1_C4, phi_1_C45), list(phi_1_IC5, phi_1_IC4, phi_1_IC45), S)$alpha_1
  
  adjusted_alpha_1 = block_weight(list(adjusted_phi_1_C5, adjusted_phi_1_C4, adjusted_phi_1_C45), 
                                  list(adjusted_phi_1_IC5, adjusted_phi_1_IC4, adjusted_phi_1_IC45), S)$alpha_1
  
  SS_C1 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5))
  SS_C2 = t(diag(alpha_1[2, ]) %*% t(phi_1_C4))
  SS_C3 = t(diag(alpha_1[3, ]) %*% t(phi_1_C45))
  
  SS_IC1 = t(diag(alpha_1[1, ]) %*% t(phi_1_IC5))
  SS_IC2 = t(diag(alpha_1[2, ]) %*% t(phi_1_IC4))
  SS_IC3 = t(diag(alpha_1[3, ]) %*% t(phi_1_IC45))
  
  adjusted_SS_C1 = t(diag(adjusted_alpha_1[1, ]) %*% t(adjusted_phi_1_C5))
  adjusted_SS_C2 = t(diag(adjusted_alpha_1[2, ]) %*% t(adjusted_phi_1_C4))
  adjusted_SS_C3 = t(diag(adjusted_alpha_1[3, ]) %*% t(adjusted_phi_1_C45))
  
  adjusted_SS_IC1 = t(diag(adjusted_alpha_1[1, ]) %*% t(adjusted_phi_1_IC5))
  adjusted_SS_IC2 = t(diag(adjusted_alpha_1[2, ]) %*% t(adjusted_phi_1_IC4))
  adjusted_SS_IC3 = t(diag(adjusted_alpha_1[3, ]) %*% t(adjusted_phi_1_IC45))
  
  gamma_tilt_1 = gamma.tilt - colMeans(SS_C1) - colMeans(SS_C2) - colMeans(SS_C3) + 
    colMeans(SS_IC1) + colMeans(SS_IC2) + colMeans(SS_IC3) 
  
  adjusted_gamma_tilt_1 = gamma.tilt - colMeans(adjusted_SS_C1) - colMeans(adjusted_SS_C2) - colMeans(adjusted_SS_C3) + 
    colMeans(adjusted_SS_IC1) + colMeans(adjusted_SS_IC2) + colMeans(adjusted_SS_IC3) 
  
  var1 = colVars(S - SS_C1 - SS_C2 - SS_C3) / n_lc_train + 
    colVars(SS_IC1) / n_lm1_train + 
    colVars(SS_IC2) / n_lm2_train +
    colVars(SS_IC3) / n_lm3_train
  
  adjusted_var1 = colVars(S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3) / n_lc_train + 
    colVars(adjusted_SS_IC1) / n_lm1_train + 
    colVars(adjusted_SS_IC2) / n_lm2_train +
    colVars(adjusted_SS_IC3) / n_lm3_train
  
  #var1 / adjusted_var1
  ## unlabel
  
  rows_CU = rep(1:(n_uc + n_lc_train),each = monte_sample)
  obs_CU45_Y0 = X_CU_train[rows_CU,]
  obs_CU4_Y0 = X_CU_train[rows_CU,]
  obs_CU5_Y0 = X_CU_train[rows_CU,]
  
  obs_CU45_Y1 = X_CU_train[rows_CU,]
  obs_CU4_Y1 = X_CU_train[rows_CU,]
  obs_CU5_Y1 = X_CU_train[rows_CU,]
  
  sigma_hat5_Y0 = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5_Y1 = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4_Y0 = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4_Y1 = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4_Y0 = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4_Y1 = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_5_Y0 = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5_Y1 = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45_Y0 = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45_Y1 = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
## for y = 0
  obs_CU45_Y0[, M[[3]][1]] = rep(base_sample45_4, (n_uc + n_lc_train)) * rep(sigma_hat5_Y0, monte_sample *  (n_uc + n_lc_train)) + 
    rep(mu_hat_CU45_4_Y0, each=monte_sample)
  obs_CU45_Y0[, M[[3]][2]] = (rep(base_sample45_4,(n_uc + n_lc_train)) * rep(cor45_Y0, monte_sample *  (n_uc + n_lc_train)) + 
                            rep(base_sample45_5,(n_uc + n_lc_train)) * rep(sqrt(1-cor45_Y0^2), monte_sample *  (n_uc + n_lc_train))) * 
    rep(sigma_hat45_5_Y0,each = monte_sample) + rep(mu_hat_CU45_5_Y0, each=monte_sample)
  
  obs_CU5_Y0[, M[[1]]] = rep(base_sample5, (n_uc + n_lc_train)) * rep(sigma_hat5_Y0, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU5_Y0, each=monte_sample)
  obs_CU4_Y0[, M[[2]]] = rep(base_sample4, (n_uc + n_lc_train)) * rep(sigma_hat4_Y0, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU4_Y0, each=monte_sample)
## for y = 1
  
  obs_CU45_Y1[, M[[3]][1]] = rep(base_sample45_4, (n_uc + n_lc_train)) * rep(sigma_hat5_Y1, monte_sample *  (n_uc + n_lc_train)) + 
    rep(mu_hat_CU45_4_Y1, each=monte_sample)
  obs_CU45_Y1[, M[[3]][2]] = (rep(base_sample45_4,(n_uc + n_lc_train)) * rep(cor45_Y1, monte_sample *  (n_uc + n_lc_train)) + 
                                rep(base_sample45_5,(n_uc + n_lc_train)) * rep(sqrt(1-cor45_Y1^2), monte_sample *  (n_uc + n_lc_train))) * 
    rep(sigma_hat45_5_Y1,each = monte_sample) + rep(mu_hat_CU45_5_Y1, each=monte_sample)
  
  obs_CU5_Y1[, M[[1]]] = rep(base_sample5, (n_uc + n_lc_train)) * rep(sigma_hat5_Y1, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU5_Y1, each=monte_sample)
  obs_CU4_Y1[, M[[2]]] = rep(base_sample4, (n_uc + n_lc_train)) * rep(sigma_hat4_Y1, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU4_Y1, each=monte_sample)
  
  obs_CU45_Y0 <- cbind(1,obs_CU45_Y0) ## new hessian
  obs_CU4_Y0 <- cbind(1,obs_CU4_Y0)
  obs_CU5_Y0 <- cbind(1,obs_CU5_Y0)
  
  obs_CU45_Y1 <- cbind(1,obs_CU45_Y1) ## new hessian
  obs_CU4_Y1 <- cbind(1,obs_CU4_Y1)
  obs_CU5_Y1 <- cbind(1,obs_CU5_Y1)
  
  S_monte_CU45_Y0 = t(H_inv %*% t(obs_CU45_Y0)) * as.numeric(-g_x(as.matrix(obs_CU45_Y0), gamma.tilt))
  S_monte_CU4_Y0 = t(H_inv %*% t(obs_CU4_Y0)) * as.numeric(-g_x(as.matrix(obs_CU4_Y0), gamma.tilt))
  S_monte_CU5_Y0 = t(H_inv %*% t(obs_CU5_Y0)) * as.numeric(-g_x(as.matrix(obs_CU5_Y0), gamma.tilt))
  
  S_monte_CU45_Y1 = t(H_inv %*% t(obs_CU45_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU45_Y1), gamma.tilt))
  S_monte_CU4_Y1 = t(H_inv %*% t(obs_CU4_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU4_Y1), gamma.tilt))
  S_monte_CU5_Y1 = t(H_inv %*% t(obs_CU5_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU5_Y1), gamma.tilt))
  
  phi_1_C45_Y0 = as.matrix(aggregate(S_monte_CU45_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C4_Y0 = as.matrix(aggregate(S_monte_CU4_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C5_Y0 = as.matrix(aggregate(S_monte_CU5_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  
  phi_1_C45_Y1 = as.matrix(aggregate(S_monte_CU45_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C4_Y1 = as.matrix(aggregate(S_monte_CU4_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C5_Y1 = as.matrix(aggregate(S_monte_CU5_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  
  phi_1_Y0 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5_Y0)) + 
    t(diag(alpha_1[2, ]) %*% t(phi_1_C4_Y0)) +
    t(diag(alpha_1[3, ]) %*% t(phi_1_C45_Y0))
  
  phi_1_Y1 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5_Y1)) + 
    t(diag(alpha_1[2, ]) %*% t(phi_1_C4_Y1)) +
    t(diag(alpha_1[3, ]) %*% t(phi_1_C45_Y1))
  
  R = (S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3)
  
  prob = predict(naive, newdata = X_CU_train, type = "response")
  
  phi_2 = prob * phi_1_Y1 + (1 - prob) * phi_1_Y0
  
  phi_2_0 = phi_2 -  t(matrix(rep(colMeans(phi_2), (n_uc + n_lc_train)), nrow = ncol(S)))
  
  allocation_uc = allocation2("intercept", phi_2_0[1:n_lc_train, ], phi_2_0[(n_lc_train+1):(n_uc + n_lc_train), ],
                              R, X_lc_train, X_uc_train)
  
  adjusted_gamma_tilt_2 = adjusted_gamma_tilt_1 - colMeans(allocation_uc$adjusted_phi_2[1:n_lc_train, ]) + 
    colMeans(allocation_uc$adjusted_phi_2)
  
  var = colMeans(S^2) / n_lc_train 
  var1 = colMeans((S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3)^2) / n_lc_train + 
    colMeans(adjusted_SS_IC1^2) / n_lm1_train + 
    colMeans(adjusted_SS_IC2^2) / n_lm2_train +
    colMeans(adjusted_SS_IC3^2) / n_lm3_train
  
  var2 = colMeans((S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3 - 
                     allocation_uc$adjusted_phi_2[1:n_lc_train, ])^2) / n_lc_train + 
    colMeans(adjusted_SS_IC1^2) / n_lm1_train + 
    colMeans(adjusted_SS_IC2^2) / n_lm2_train +
    colMeans(adjusted_SS_IC3^2) / n_lm3_train + 
    colMeans(allocation_uc$adjusted_phi_2[(n_lc_train+1):(n_uc + n_lc_train), ] ^ 2) / n_uc 
  return(list(gamma = gamma.tilt, 
              gamma1 = adjusted_gamma_tilt_1, 
              gamma2 = adjusted_gamma_tilt_2, 
              se = sqrt(var),
              se1 = sqrt(var1),
              se2 = sqrt(var2)))
}


# Get the proposed estimators
# Get the proposed estimators
sip_dop_get <- function(seed){
  
  set.seed(seed)
  naive = glm(Y_lc_train ~., cbind(Y_lc_train, X_lc_train),family = "binomial")
  gamma.tilt = naive$coefficients
  
  # Step 2: Estimation of Score Function
  # get the score function
  H_inv = solve(H(as.matrix(cbind(1,X_lc_train)), gamma.tilt))
  S = sweep(t(H_inv %*% t(as.matrix(cbind(1,X_lc_train)))), MARGIN=1, (Y_lc_train - g_x(as.matrix(cbind(1,X_lc_train)), gamma.tilt)), '*')
  
  mu_hat_C5 = c()
  mu_hat_C4 = c()
  mu_hat_C45_4 = c()
  mu_hat_C45_5 = c()
  
  mu_hat_IC5 = matrix(0, nrow = nlm[1], ncol = 5)
  mu_hat_IC4 = matrix(0, nrow = nlm[2], ncol = 5)
  mu_hat_IC45_4 = matrix(0, nrow = nlm[3], ncol = 5)
  mu_hat_IC45_5 = matrix(0, nrow = nlm[3], ncol = 5)
  
  mu_hat_IC5 = matrix(0, nrow = nlm[1], ncol = 5)
  mu_hat_IC4 = matrix(0, nrow = nlm[2], ncol = 5)
  mu_hat_IC45_4 = matrix(0, nrow = nlm[3], ncol = 5)
  mu_hat_IC45_5 = matrix(0, nrow = nlm[3], ncol = 5)
  
  mu_hat_CU5_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU4_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_4_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_5_Y0 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  
  mu_hat_CU5_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU4_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_4_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  mu_hat_CU45_5_Y1 = matrix(0, nrow = (n_uc + n_lc_train), ncol = 5)
  
  for (i in 1:5) {
    # mean estimation
    mu_hat_C5[which((Y_lc_train == 0) & (folds == i))] = predict(mod5[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[1]]])
    mu_hat_C5[which((Y_lc_train == 1) & (folds == i))] = predict(mod5[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[1]]])
    
    mu_hat_C4[which((Y_lc_train == 0) & (folds == i))] = predict(mod4[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[2]]])
    mu_hat_C4[which((Y_lc_train == 1) & (folds == i))] = predict(mod4[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[2]]])
    
    mu_hat_C45_4[which((Y_lc_train == 0) & (folds == i))] = predict(mod45_4[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[3]]])
    mu_hat_C45_4[which((Y_lc_train == 1) & (folds == i))] = predict(mod45_4[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[3]]])
    
    mu_hat_C45_5[which((Y_lc_train == 0) & (folds == i))] = predict(mod45_5[[i]][[1]], newdata = X_lc_train_withauxi[which((Y_lc_train == 0) & (folds == i)),-M[[3]]])
    mu_hat_C45_5[which((Y_lc_train == 1) & (folds == i))] = predict(mod45_5[[i]][[2]], newdata = X_lc_train_withauxi[which((Y_lc_train == 1) & (folds == i)),-M[[3]]])
    
    mu_hat_IC5[which(Y_lm1_train == 0),i] = predict(mod5[[i]][[1]], newdata = X_lm1_train_withauxi[which(Y_lm1_train == 0),-M[[1]]]) 
    mu_hat_IC5[which(Y_lm1_train == 1),i] = predict(mod5[[i]][[2]], newdata = X_lm1_train_withauxi[which(Y_lm1_train == 1),-M[[1]]])
    
    mu_hat_IC4[which(Y_lm2_train == 0),i] = predict(mod4[[i]][[1]], newdata = X_lm2_train_withauxi[which(Y_lm2_train == 0),-M[[2]]]) 
    mu_hat_IC4[which(Y_lm2_train == 1),i] = predict(mod4[[i]][[2]], newdata = X_lm2_train_withauxi[which(Y_lm2_train == 1),-M[[2]]])
    
    mu_hat_IC45_4[which(Y_lm3_train == 0),i] = predict(mod45_4[[i]][[1]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 0),-M[[3]]]) 
    mu_hat_IC45_4[which(Y_lm3_train == 1),i] = predict(mod45_4[[i]][[2]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 1),-M[[3]]])
    
    mu_hat_IC45_5[which(Y_lm3_train == 0),i] = predict(mod45_5[[i]][[1]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 0),-M[[3]]]) 
    mu_hat_IC45_5[which(Y_lm3_train == 1),i] = predict(mod45_5[[i]][[2]], newdata = X_lm3_train_withauxi[which(Y_lm3_train == 1),-M[[3]]])
    
    ## unlabel
    mu_hat_CU5_Y0[ ,i] = predict(mod5[[i]][[1]], newdata = X_CU_train_withauxi[,-M[[1]]]) 
    mu_hat_CU5_Y1[ ,i] = predict(mod5[[i]][[2]], newdata = X_CU_train_withauxi[,-M[[1]]])
    
    mu_hat_CU4_Y0[ ,i] = predict(mod4[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[2]]]) 
    mu_hat_CU4_Y1[ ,i] = predict(mod4[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[2]]])
    
    mu_hat_CU45_4_Y0[ ,i] = predict(mod45_4[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[3]]]) 
    mu_hat_CU45_4_Y1[ ,i] = predict(mod45_4[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[3]]])
    
    mu_hat_CU45_5_Y0[ ,i] = predict(mod45_5[[i]][[1]], newdata = X_CU_train_withauxi[ ,-M[[3]]]) 
    mu_hat_CU45_5_Y1[ ,i] = predict(mod45_5[[i]][[2]], newdata = X_CU_train_withauxi[ ,-M[[3]]])
  }
  mu_hat_IC5 = rowMeans(mu_hat_IC5)
  mu_hat_IC4 = rowMeans(mu_hat_IC4) # take the mean values of fitted values from 5 models
  
  mu_hat_IC45_4 = rowMeans(mu_hat_IC45_4)
  mu_hat_IC45_5 = rowMeans(mu_hat_IC45_5)
  
  mu_hat_CU5_Y0 = rowMeans(mu_hat_CU5_Y0)
  mu_hat_CU4_Y0 = rowMeans(mu_hat_CU4_Y0)
  mu_hat_CU45_4_Y0 = rowMeans(mu_hat_CU45_4_Y0)
  mu_hat_CU45_5_Y0 = rowMeans(mu_hat_CU45_5_Y0)
  
  mu_hat_CU5_Y1 = rowMeans(mu_hat_CU5_Y1)
  mu_hat_CU4_Y1 = rowMeans(mu_hat_CU4_Y1)
  mu_hat_CU45_4_Y1 = rowMeans(mu_hat_CU45_4_Y1)
  mu_hat_CU45_5_Y1 = rowMeans(mu_hat_CU45_5_Y1)
  
  # Monte Carlo Integration
  # variance estimation
  resi5 = pull(X_lc_train, M[[1]]) - mu_hat_C5
  resi4 = pull(X_lc_train, M[[2]]) - mu_hat_C4
  
  resi45_4 = pull(X_lc_train, M[[3]][1]) - mu_hat_C45_4
  resi45_5 = pull(X_lc_train, M[[3]][2]) - mu_hat_C45_5
  
  sigma_hat5 = c()
  sigma_hat4 = c()
  
  sigma_hat45_4 = c()
  sigma_hat45_5 = c()
  cor45 = c()
  
  sigma_hat5[which(Y_lc_train==0)] = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5[which(Y_lc_train==1)] = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4[which(Y_lc_train==0)] = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4[which(Y_lc_train==1)] = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4[which(Y_lc_train==0)] = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4[which(Y_lc_train==1)] = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_5[which(Y_lc_train==0)] = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5[which(Y_lc_train==1)] = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45[which(Y_lc_train==0)] = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45[which(Y_lc_train==1)] = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
  sigma_hat5_IC = c() 
  sigma_hat4_IC = c() 
  
  sigma_hat45_4_IC = c()  
  sigma_hat45_5_IC = c()
  cor45_IC = c()
  
  sigma_hat5_IC[which(Y_lm1_train==0)] = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5_IC[which(Y_lm1_train==1)] = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4_IC[which(Y_lm2_train==0)] = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4_IC[which(Y_lm2_train==1)] = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4_IC[which(Y_lm3_train==0)] = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4_IC[which(Y_lm3_train==1)] = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  sigma_hat45_5_IC[which(Y_lm3_train==0)] = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5_IC[which(Y_lm3_train==1)] = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45_IC[which(Y_lm3_train==0)] = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45_IC[which(Y_lm3_train==1)] = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
  # Monte Carlo Integration Process
  monte_sample = 500
  base_sample5 = rnorm(monte_sample, sd = 1)
  base_sample4 = rnorm(monte_sample, sd = 1)
  
  base_sample45_4 = rnorm(monte_sample, sd = 1)
  base_sample45_5 = rnorm(monte_sample, sd = 1)
  
  rows = rep(1:n_lc_train,each = monte_sample)
  obs_C45 = X_lc_train[rows,]
  obs_C4 = X_lc_train[rows,]
  obs_C5 = X_lc_train[rows,]
  
  obs_C45[, M[[3]][1]] = rep(base_sample45_4, n_lc_train) * rep(sigma_hat45_4, each = monte_sample) + 
    rep(mu_hat_C45_4, each=monte_sample)
  obs_C45[, M[[3]][2]] = (rep(base_sample45_4,n_lc_train) * rep(cor45,each = monte_sample) + 
                            rep(base_sample45_5,n_lc_train) * rep(sqrt(1-cor45^2),each=monte_sample)) * 
    rep(sigma_hat45_5,each = monte_sample) + rep(mu_hat_C45_5, each=monte_sample)
  
  obs_C5[, M[[1]]] = rep(base_sample5, n_lc_train) * rep(sigma_hat5, each = monte_sample) + rep(mu_hat_C5, each=monte_sample)
  obs_C4[, M[[2]]] = rep(base_sample4, n_lc_train) * rep(sigma_hat4, each = monte_sample) + rep(mu_hat_C4, each=monte_sample)
  
  obs_C45 <- cbind(1,obs_C45) ## new hessian
  obs_C4 <- cbind(1,obs_C4)
  obs_C5 <- cbind(1,obs_C5)
  
  S_monte_C45 = t(H_inv %*% t(obs_C45)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C45), gamma.tilt))
  S_monte_C4 = t(H_inv %*% t(obs_C4)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C4), gamma.tilt))
  S_monte_C5 = t(H_inv %*% t(obs_C5)) * as.numeric(rep(Y_lc_train,each = monte_sample)-g_x(as.matrix(obs_C5), gamma.tilt))
  
  phi_1_C45 = as.matrix(aggregate(S_monte_C45, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  phi_1_C4 = as.matrix(aggregate(S_monte_C4, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  phi_1_C5 = as.matrix(aggregate(S_monte_C5, list(rep(1:n_lc_train, each = monte_sample)), mean)[-1])
  
  rowsI5 = rep(1:nlm[1],each = monte_sample)
  rowsI4 = rep(1:nlm[2],each = monte_sample)
  rowsI45 = rep(1:nlm[3],each = monte_sample)
  
  obs_IC5 = X_lm1_train[rowsI5,]
  obs_IC4 = X_lm2_train[rowsI4,]
  obs_IC45 = X_lm3_train[rowsI45,]
  
  obs_IC5[, M[[1]]] = rep(base_sample5,nlm[1]) * rep(sigma_hat5_IC, each = monte_sample) + 
    rep(mu_hat_IC5, each=monte_sample)
  obs_IC4[, M[[2]]] = rep(base_sample4,nlm[2]) * rep(sigma_hat4_IC, each = monte_sample) + 
    rep(mu_hat_IC4, each=monte_sample)
  
  obs_IC45[, M[[3]][1]] = rep(base_sample45_4,nlm[3]) * rep(sigma_hat45_4_IC, each = monte_sample) + 
    rep(mu_hat_IC45_4, each=monte_sample) ## same mu but not same variance, so not from true distribution
  obs_IC45[, M[[3]][2]] = (rep(base_sample45_4,nlm[3]) * rep(cor45_IC,each = monte_sample) + 
                             rep(base_sample45_5,nlm[3])*rep(sqrt(1-cor45_IC^2),each=monte_sample)) * rep(sigma_hat45_5_IC, each = monte_sample) + 
    rep(mu_hat_IC45_5, each=monte_sample)
  
  obs_IC5 <- cbind(1,obs_IC5)
  obs_IC4 <- cbind(1,obs_IC4)
  obs_IC45 <- cbind(1,obs_IC45)
  
  S_monte_IC5 = t(H_inv %*% t(obs_IC5)) * as.numeric(rep(Y_lm1_train,each = monte_sample)-g_x(as.matrix(obs_IC5), gamma.tilt))
  S_monte_IC4 = t(H_inv %*% t(obs_IC4)) * as.numeric(rep(Y_lm2_train,each = monte_sample)-g_x(as.matrix(obs_IC4), gamma.tilt))
  S_monte_IC45 = t(H_inv %*% t(obs_IC45)) * as.numeric(rep(Y_lm3_train,each = monte_sample)-g_x(as.matrix(obs_IC45), gamma.tilt))
  
  phi_1_IC5 = as.matrix(aggregate(S_monte_IC5, list(rep(1:nlm[1], each = monte_sample)), mean)[-1])
  phi_1_IC4 = as.matrix(aggregate(S_monte_IC4, list(rep(1:nlm[2], each = monte_sample)), mean)[-1])
  phi_1_IC45 = as.matrix(aggregate(S_monte_IC45, list(rep(1:nlm[3], each = monte_sample)), mean)[-1])
  
  ### now we allocation the phi_4, phi_4, phi_45
  ## method = "intercept", "linear", "interaction"
  allocation5 = allocation("linear", phi_1_C5, phi_1_IC5, S, M[[1]], X_lc_train, X_lm1_train, Y_lc_train, Y_lm1_train)
  allocation4 = allocation("linear", phi_1_C4, phi_1_IC4, S, M[[2]], X_lc_train, X_lm2_train, Y_lc_train, Y_lm2_train)
  allocation45 = allocation("linear", phi_1_C45, phi_1_IC45, S, M[[3]], X_lc_train, X_lm3_train, Y_lc_train, Y_lm3_train)
  adjusted_phi_1_C5 = allocation5$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_C4 = allocation4$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_C45 = allocation45$adjusted_phi_1[1:n_lc_train, ]
  adjusted_phi_1_IC5 = allocation5$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm1_train), ]
  adjusted_phi_1_IC4 = allocation4$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm2_train), ]
  adjusted_phi_1_IC45 = allocation45$adjusted_phi_1[(1 + n_lc_train):(n_lc_train + n_lm3_train), ]
  alpha_1 = block_weight(list(phi_1_C5, phi_1_C4, phi_1_C45), list(phi_1_IC5, phi_1_IC4, phi_1_IC45), S)$alpha_1
  
  adjusted_alpha_1 = block_weight(list(adjusted_phi_1_C5, adjusted_phi_1_C4, adjusted_phi_1_C45), 
                                  list(adjusted_phi_1_IC5, adjusted_phi_1_IC4, adjusted_phi_1_IC45), S)$alpha_1
  
  SS_C1 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5))
  SS_C2 = t(diag(alpha_1[2, ]) %*% t(phi_1_C4))
  SS_C3 = t(diag(alpha_1[3, ]) %*% t(phi_1_C45))
  
  SS_IC1 = t(diag(alpha_1[1, ]) %*% t(phi_1_IC5))
  SS_IC2 = t(diag(alpha_1[2, ]) %*% t(phi_1_IC4))
  SS_IC3 = t(diag(alpha_1[3, ]) %*% t(phi_1_IC45))
  
  adjusted_SS_C1 = t(diag(adjusted_alpha_1[1, ]) %*% t(adjusted_phi_1_C5))
  adjusted_SS_C2 = t(diag(adjusted_alpha_1[2, ]) %*% t(adjusted_phi_1_C4))
  adjusted_SS_C3 = t(diag(adjusted_alpha_1[3, ]) %*% t(adjusted_phi_1_C45))
  
  adjusted_SS_IC1 = t(diag(adjusted_alpha_1[1, ]) %*% t(adjusted_phi_1_IC5))
  adjusted_SS_IC2 = t(diag(adjusted_alpha_1[2, ]) %*% t(adjusted_phi_1_IC4))
  adjusted_SS_IC3 = t(diag(adjusted_alpha_1[3, ]) %*% t(adjusted_phi_1_IC45))
  
  gamma_tilt_1 = gamma.tilt - colMeans(SS_C1) - colMeans(SS_C2) - colMeans(SS_C3) + 
    colMeans(SS_IC1) + colMeans(SS_IC2) + colMeans(SS_IC3) 
  
  adjusted_gamma_tilt_1 = gamma.tilt - colMeans(adjusted_SS_C1) - colMeans(adjusted_SS_C2) - colMeans(adjusted_SS_C3) + 
    colMeans(adjusted_SS_IC1) + colMeans(adjusted_SS_IC2) + colMeans(adjusted_SS_IC3) 
  
  var1 = colVars(S - SS_C1 - SS_C2 - SS_C3) / n_lc_train + 
    colVars(SS_IC1) / n_lm1_train + 
    colVars(SS_IC2) / n_lm2_train +
    colVars(SS_IC3) / n_lm3_train
  
  adjusted_var1 = colVars(S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3) / n_lc_train + 
    colVars(adjusted_SS_IC1) / n_lm1_train + 
    colVars(adjusted_SS_IC2) / n_lm2_train +
    colVars(adjusted_SS_IC3) / n_lm3_train
  
  #var1 / adjusted_var1
  ## unlabel
  
  rows_CU = rep(1:(n_uc + n_lc_train),each = monte_sample)
  obs_CU45_Y0 = X_CU_train[rows_CU,]
  obs_CU4_Y0 = X_CU_train[rows_CU,]
  obs_CU5_Y0 = X_CU_train[rows_CU,]
  
  obs_CU45_Y1 = X_CU_train[rows_CU,]
  obs_CU4_Y1 = X_CU_train[rows_CU,]
  obs_CU5_Y1 = X_CU_train[rows_CU,]
  
  sigma_hat5_Y0 = sqrt(mean(resi5[which(Y_lc_train==0)]^2))
  sigma_hat5_Y1 = sqrt(mean(resi5[which(Y_lc_train==1)]^2))
  
  sigma_hat4_Y0 = sqrt(mean(resi4[which(Y_lc_train==0)]^2))
  sigma_hat4_Y1 = sqrt(mean(resi4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_4_Y0 = sqrt(mean(resi45_4[which(Y_lc_train==0)]^2))
  sigma_hat45_4_Y1 = sqrt(mean(resi45_4[which(Y_lc_train==1)]^2))
  
  sigma_hat45_5_Y0 = sqrt(mean(resi45_5[which(Y_lc_train==0)]^2))
  sigma_hat45_5_Y1 = sqrt(mean(resi45_5[which(Y_lc_train==1)]^2))
  
  cor45_Y0 = mean(resi45_4[which(Y_lc_train==0)] * resi45_5[which(Y_lc_train==0)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==0)]) * unique(sigma_hat45_4[which(Y_lc_train==0)]))
  cor45_Y1 = mean(resi45_4[which(Y_lc_train==1)] * resi45_5[which(Y_lc_train==1)]) / 
    (unique(sigma_hat45_5[which(Y_lc_train==1)]) * unique(sigma_hat45_4[which(Y_lc_train==1)]))
  
  ## for y = 0
  obs_CU45_Y0[, M[[3]][1]] = rep(base_sample45_4, (n_uc + n_lc_train)) * rep(sigma_hat5_Y0, monte_sample *  (n_uc + n_lc_train)) + 
    rep(mu_hat_CU45_4_Y0, each=monte_sample)
  obs_CU45_Y0[, M[[3]][2]] = (rep(base_sample45_4,(n_uc + n_lc_train)) * rep(cor45_Y0, monte_sample *  (n_uc + n_lc_train)) + 
                                rep(base_sample45_5,(n_uc + n_lc_train)) * rep(sqrt(1-cor45_Y0^2), monte_sample *  (n_uc + n_lc_train))) * 
    rep(sigma_hat45_5_Y0,each = monte_sample) + rep(mu_hat_CU45_5_Y0, each=monte_sample)
  
  obs_CU5_Y0[, M[[1]]] = rep(base_sample5, (n_uc + n_lc_train)) * rep(sigma_hat5_Y0, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU5_Y0, each=monte_sample)
  obs_CU4_Y0[, M[[2]]] = rep(base_sample4, (n_uc + n_lc_train)) * rep(sigma_hat4_Y0, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU4_Y0, each=monte_sample)
  ## for y = 1
  
  obs_CU45_Y1[, M[[3]][1]] = rep(base_sample45_4, (n_uc + n_lc_train)) * rep(sigma_hat5_Y1, monte_sample *  (n_uc + n_lc_train)) + 
    rep(mu_hat_CU45_4_Y1, each=monte_sample)
  obs_CU45_Y1[, M[[3]][2]] = (rep(base_sample45_4,(n_uc + n_lc_train)) * rep(cor45_Y1, monte_sample *  (n_uc + n_lc_train)) + 
                                rep(base_sample45_5,(n_uc + n_lc_train)) * rep(sqrt(1-cor45_Y1^2), monte_sample *  (n_uc + n_lc_train))) * 
    rep(sigma_hat45_5_Y1,each = monte_sample) + rep(mu_hat_CU45_5_Y1, each=monte_sample)
  
  obs_CU5_Y1[, M[[1]]] = rep(base_sample5, (n_uc + n_lc_train)) * rep(sigma_hat5_Y1, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU5_Y1, each=monte_sample)
  obs_CU4_Y1[, M[[2]]] = rep(base_sample4, (n_uc + n_lc_train)) * rep(sigma_hat4_Y1, monte_sample *  (n_uc + n_lc_train)) + rep(mu_hat_CU4_Y1, each=monte_sample)
  
  obs_CU45_Y0 <- cbind(1,obs_CU45_Y0) ## new hessian
  obs_CU4_Y0 <- cbind(1,obs_CU4_Y0)
  obs_CU5_Y0 <- cbind(1,obs_CU5_Y0)
  
  obs_CU45_Y1 <- cbind(1,obs_CU45_Y1) ## new hessian
  obs_CU4_Y1 <- cbind(1,obs_CU4_Y1)
  obs_CU5_Y1 <- cbind(1,obs_CU5_Y1)
  
  S_monte_CU45_Y0 = t(H_inv %*% t(obs_CU45_Y0)) * as.numeric(-g_x(as.matrix(obs_CU45_Y0), gamma.tilt))
  S_monte_CU4_Y0 = t(H_inv %*% t(obs_CU4_Y0)) * as.numeric(-g_x(as.matrix(obs_CU4_Y0), gamma.tilt))
  S_monte_CU5_Y0 = t(H_inv %*% t(obs_CU5_Y0)) * as.numeric(-g_x(as.matrix(obs_CU5_Y0), gamma.tilt))
  
  S_monte_CU45_Y1 = t(H_inv %*% t(obs_CU45_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU45_Y1), gamma.tilt))
  S_monte_CU4_Y1 = t(H_inv %*% t(obs_CU4_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU4_Y1), gamma.tilt))
  S_monte_CU5_Y1 = t(H_inv %*% t(obs_CU5_Y1)) * as.numeric(1-g_x(as.matrix(obs_CU5_Y1), gamma.tilt))
  
  phi_1_C45_Y0 = as.matrix(aggregate(S_monte_CU45_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C4_Y0 = as.matrix(aggregate(S_monte_CU4_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C5_Y0 = as.matrix(aggregate(S_monte_CU5_Y0, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  
  phi_1_C45_Y1 = as.matrix(aggregate(S_monte_CU45_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C4_Y1 = as.matrix(aggregate(S_monte_CU4_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  phi_1_C5_Y1 = as.matrix(aggregate(S_monte_CU5_Y1, list(rep(1:(n_uc + n_lc_train), each = monte_sample)), mean)[-1])
  
  phi_1_Y0 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5_Y0)) + 
    t(diag(alpha_1[2, ]) %*% t(phi_1_C4_Y0)) +
    t(diag(alpha_1[3, ]) %*% t(phi_1_C45_Y0))
  
  phi_1_Y1 = t(diag(alpha_1[1, ]) %*% t(phi_1_C5_Y1)) + 
    t(diag(alpha_1[2, ]) %*% t(phi_1_C4_Y1)) +
    t(diag(alpha_1[3, ]) %*% t(phi_1_C45_Y1))
  
  R = (S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3)
  
  prob = predict(naive, newdata = X_CU_train, type = "response")
  
  phi_2 = prob * phi_1_Y1 + (1 - prob) * phi_1_Y0
  
  phi_2_0 = phi_2 -  t(matrix(rep(colMeans(phi_2), (n_uc + n_lc_train)), nrow = ncol(S)))
  
  allocation_uc = allocation2("linear", phi_2_0[1:n_lc_train, ], phi_2_0[(n_lc_train+1):(n_uc + n_lc_train), ],
                              R, X_lc_train, X_uc_train)
  
  adjusted_gamma_tilt_2 = adjusted_gamma_tilt_1 - colMeans(allocation_uc$adjusted_phi_2[1:n_lc_train, ]) + 
    colMeans(allocation_uc$adjusted_phi_2)
  
  var = colMeans(S^2) / n_lc_train 
  var1 = colMeans((S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3)^2) / n_lc_train + 
    colMeans(adjusted_SS_IC1^2) / n_lm1_train + 
    colMeans(adjusted_SS_IC2^2) / n_lm2_train +
    colMeans(adjusted_SS_IC3^2) / n_lm3_train
  
  var2 = colMeans((S - adjusted_SS_C1 - adjusted_SS_C2 - adjusted_SS_C3 - 
                     allocation_uc$adjusted_phi_2[1:n_lc_train, ])^2) / n_lc_train + 
    colMeans(adjusted_SS_IC1^2) / n_lm1_train + 
    colMeans(adjusted_SS_IC2^2) / n_lm2_train +
    colMeans(adjusted_SS_IC3^2) / n_lm3_train + 
    colMeans(allocation_uc$adjusted_phi_2[(n_lc_train+1):(n_uc + n_lc_train), ] ^ 2) / n_uc 
  return(list(gamma = gamma.tilt, 
              gamma1 = adjusted_gamma_tilt_1, 
              gamma2 = adjusted_gamma_tilt_2, 
              se = sqrt(var),
              se1 = sqrt(var1),
              se2 = sqrt(var2)))
}
