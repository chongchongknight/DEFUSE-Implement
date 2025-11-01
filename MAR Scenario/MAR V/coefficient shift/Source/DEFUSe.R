## data generation process for changed coefficient shift and nonlinearity
mix_gaussian <- function(n_site, n_total, dim, Miss, gamma.bar, seednum, shift, nonlinear){
  #c = 0.2
  #b = 1
  set.seed(seednum)
  error_table = c(-0.3, 0.5, -0.5, -0.2, 0.4, 0)
  x = vector("list", n_site)
  y = vector("list", n_site)
  p = vector("list", n_site)
  for (i in 1:n_site) {
    error = error_table[i]
    x[[i]] = matrix(rnorm(n_total[i] * dim) + error, n_total[i], dim)
    x[[i]][, Miss[1]] = x[[i]][, -Miss] %*% c(0.3, 1.0,  0.5) + rnorm(n_total[i])
    x[[i]][, Miss[2]] = x[[i]][, -Miss] %*% c(0.2, 0.5, -0.3) + rnorm(n_total[i])
    if (i %in% c(3, 5)) {
      p[[i]] = g_x(cbind(1, x[[i]], x[[i]][, 1] ^ 2, x[[i]][, 2] ^ 2),  
                   c(0, -1+shift, -1-shift, 1, 1, 1, nonlinear, nonlinear))
    } else {
      p[[i]] = g_x(cbind(1, x[[i]], x[[i]][, 1] ^ 2, x[[i]][, 2] ^ 2),  
                   c(2 * nonlinear, -1, -1, 1, 1, 1, 0, 0))
    }
    y[[i]] = rnorm(n_total[i]) + p[[i]]
  }
  return(list(X = x, Y = y))
}
## source-wisely calculate the max L(h)
violation = function(Combined_X, Combined_Y, source_index, seednum) {
  ## seednum = 1
  keep_id1 = rep(0, 4)
  keep_id2 = rep(0, 4)
  keep_id3 = rep(0, 4)
  keep_id4 = rep(0, 4)
  X_LC = Combined_X[source_index[[1]], ]
  Y_LC = Combined_Y[source_index[[1]]]
  X_UC = Combined_X[source_index[[6]], ]
  for (lmid in 1:4) {
    error = discrepancy(Y_LC, X_LC, Combined_X[source_index[[lmid + 1]],], Combined_Y[source_index[[lmid + 1]]], X_UC, seednum)
    keep_id1[lmid] = max(error$res1)
    keep_id2[lmid] = max(error$res2)
    keep_id3[lmid] = max(error$res3)
    keep_id4[lmid] = max(error$res4)
  }
  return(list(keep_id1 = keep_id1, keep_id2 = keep_id2, keep_id3 = keep_id3, keep_id4 = keep_id4))
}
## implement code to calculate the max L(h)
discrepancy = function(Y_LC, X_LC, X_LM, Y_LM, X_UC, seednum) {
  set.seed(seednum + 1)
  covariate_z = c(1:3)
  n_fold = 2
  fold_id_LC = sample(rep(1:n_fold, length.out = nrow(X_LC))) 
  fold_id_LM = sample(rep(1:n_fold, length.out = nrow(X_LM))) 
  dml_m_LC_LC = matrix(0, nrow = nrow(X_LC), ncol = 1)
  dml_m_LM_LM = matrix(0, nrow = nrow(X_LM), ncol = 1)
  dml_m_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  dml_m_LM_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  weight_LC = rep(0, nrow(X_LC))
  weight_LM = rep(0, nrow(X_LM))
  ## cross-fitting for DML models
  for (iLC in 1:n_fold) {
    dml_id <- (iLC + 0) %% n_fold + 1
    eva_id <- (iLC + 1) %% n_fold + 1
    weight_lc = density_ratio_xgboost_linear(X_UC[ ,covariate_z], 
                                             X_LC[fold_id_LC == dml_id, covariate_z], 
                                             X_LC[fold_id_LC == eva_id, covariate_z])
    weight_lm = density_ratio_xgboost_linear(X_UC[ ,covariate_z], 
                                             X_LM[fold_id_LM == dml_id, covariate_z], 
                                             X_LM[fold_id_LM == eva_id, covariate_z])
    weight_LC[fold_id_LC == eva_id] = weight_lc$w_apply
    weight_LM[fold_id_LM == eva_id] = weight_lm$w_apply
    best_hp_dml_LC <- tune_gausspr_cv(X_LC[fold_id_LC == dml_id, covariate_z],
                                   Y_LC[fold_id_LC == dml_id])
    invisible(capture.output({
      mod_dml_m_LC <- gausspr(
        x = X_LC[fold_id_LC == dml_id, covariate_z],
        y = Y_LC[fold_id_LC == dml_id],
        kernel = "rbfdot",
        kpar   = list(sigma = best_hp_dml_LC$sigma),
        var    = best_hp_dml_LC$var
      )
    }))
    best_hp_dml_LM <- tune_gausspr_cv(X_LM[fold_id_LM == dml_id, covariate_z], 
                                   Y_LM[fold_id_LM == dml_id])
    invisible(capture.output({
      mod_dml_m_LM <- gausspr(
        x = X_LM[fold_id_LM == dml_id, covariate_z],
        y = Y_LM[fold_id_LM == dml_id],
        kernel = "rbfdot",
        kpar   = list(sigma = best_hp_dml_LM$sigma),
        var    = best_hp_dml_LM$var
      )
    }))
    dml_m_LC_LC[fold_id_LC == eva_id] = predict(mod_dml_m_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z])
    dml_m_LC_UC = dml_m_LC_UC + predict(mod_dml_m_LC, newdata = X_UC[, covariate_z])
    dml_m_LM_LM[fold_id_LM == eva_id] = predict(mod_dml_m_LM, newdata = X_LM[fold_id_LM == eva_id, covariate_z])
    dml_m_LM_UC = dml_m_LM_UC + predict(mod_dml_m_LM, newdata = X_UC[, covariate_z])
  }
  dml_m_LC_UC = dml_m_LC_UC / n_fold
  dml_m_LM_UC = dml_m_LM_UC / n_fold
  ## generate the basis
  V_LC = generate_poly_basis(X_LC[, covariate_z], 2)
  V_LM = generate_poly_basis(X_LM[, covariate_z], 2)
  V_UC = generate_poly_basis(X_UC[, covariate_z], 2)
  res1 = selection_max(dml_m_LM_UC, dml_m_LC_UC,
                       weight_LM, Y_LM, dml_m_LM_LM,
                       weight_LC, Y_LC, dml_m_LC_LC,
                       V_UC, V_LM, V_LC)
  res2 = selection_combination(dml_m_LM_UC, dml_m_LC_UC,
                               weight_LM, Y_LM, dml_m_LM_LM,
                               weight_LC, Y_LC, dml_m_LC_LC,
                               V_UC, V_LM, V_LC)
  final_res1 = rowMaxs(res1) # which basis works best
  final_res2 = res2$value  # which combination of basis works best
  final_res3 = abs(colMeans(dml_m_LM_UC - dml_m_LC_UC)) # without correction
  final_res4 = as.numeric(res1[, 1]) # simple debias with only constant 
  return(list(res1 = final_res1, res2 = final_res2, res3 = final_res3, res4 = final_res4))
}
