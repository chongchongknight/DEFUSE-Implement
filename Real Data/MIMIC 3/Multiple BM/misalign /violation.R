
violation = function(Combined_X, Combined_Y, source_index, seednum) {
  value = rep(0, 3)
  var = rep(0, 3)
  X_LC = Combined_X[source_index[[1]], ]
  Y_LC = Combined_Y[source_index[[1]]]
  X_UC = Combined_X[source_index[[5]], ]
  
  for (lmid in 1:3) {
    error = discrepancy(
      Y_LC, X_LC,
      Combined_X[source_index[[lmid + 1]], ],
      Combined_Y[source_index[[lmid + 1]]],
      X_UC, seednum
    )
    value[lmid] = max(error$res1)
    var[lmid] = max(error$res2)
  }
  return(list(value = value, var = var))
}

discrepancy = function(Y_LC, X_LC, X_LM, Y_LM, X_UC, seednum) {
  set.seed(seednum + 1)
  covariate_z = 1:4
  n_fold = 2
  
  fold_id_LC = sample(rep(1:n_fold, length.out = nrow(X_LC)))
  fold_id_LM = sample(rep(1:n_fold, length.out = nrow(X_LM)))
  
  dml_m_LC_LC = matrix(0, nrow = nrow(X_LC), ncol = 1)
  dml_m_LM_LM = matrix(0, nrow = nrow(X_LM), ncol = 1)
  dml_m_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  dml_m_LM_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  weight_LC = rep(1, nrow(X_LC))
  weight_LM = rep(0, nrow(X_LM))
  
  for (iLC in 1:n_fold) {
    dml_id = (iLC + 0) %% n_fold + 1
    eva_id = (iLC + 1) %% n_fold + 1
    
    weight_lm = density_ratio_xgboost_linear(
      X_UC[, covariate_z],
      X_LM[fold_id_LM == dml_id, covariate_z],
      X_LM[fold_id_LM == eva_id, covariate_z]
    )
    weight_LM[fold_id_LM == eva_id] = weight_lm$w_apply
    
    best_hp_dml_LC = tune_gausspr_cv_auto(
      X_LC[fold_id_LC == dml_id, covariate_z],
      as.factor(Y_LC[fold_id_LC == dml_id])
    )
    invisible(capture.output({
      mod_dml_m_LC = gausspr(
        x = X_LC[fold_id_LC == dml_id, covariate_z],
        y = as.factor(Y_LC[fold_id_LC == dml_id]),
        kernel = "rbfdot",
        kpar = list(sigma = best_hp_dml_LC$sigma),
        var = best_hp_dml_LC$var
      )
    }))
    
    best_hp_dml_LM = tune_gausspr_cv_auto(
      X_LM[fold_id_LM == dml_id, covariate_z],
      as.factor(Y_LM[fold_id_LM == dml_id])
    )
    invisible(capture.output({
      mod_dml_m_LM = gausspr(
        x = X_LM[fold_id_LM == dml_id, covariate_z],
        y = as.factor(Y_LM[fold_id_LM == dml_id]),
        kernel = "rbfdot",
        kpar = list(sigma = best_hp_dml_LM$sigma),
        var = best_hp_dml_LM$var
      )
    }))
    
    dml_m_LC_LC[fold_id_LC == eva_id] =
      predict(mod_dml_m_LC, newdata = X_LC[fold_id_LC == eva_id, covariate_z],
              type = "probabilities")[, 2]
    dml_m_LC_UC = dml_m_LC_UC +
      matrix(predict(mod_dml_m_LC, newdata = X_UC[, covariate_z],
                     type = "probabilities")[, 2])
    dml_m_LM_LM[fold_id_LM == eva_id] =
      predict(mod_dml_m_LM, newdata = X_LM[fold_id_LM == eva_id, covariate_z],
              type = "probabilities")[, 2]
    dml_m_LM_UC = dml_m_LM_UC +
      matrix(predict(mod_dml_m_LM, newdata = X_UC[, covariate_z],
                     type = "probabilities")[, 2])
  }
  
  dml_m_LC_UC = dml_m_LC_UC / n_fold
  dml_m_LM_UC = dml_m_LM_UC / n_fold
  
  V_LC = generate_poly_basis(X_LC[, covariate_z], 2)
  V_LM = generate_poly_basis(X_LM[, covariate_z], 2)
  V_UC = generate_poly_basis(X_UC[, covariate_z], 2)
  
  train_id_LC = sample(rep(0:1, length.out = nrow(V_LC)))
  train_id_LM = sample(rep(0:1, length.out = nrow(V_LM)))
  train_id_UC = sample(rep(0:1, length.out = nrow(V_UC)))
  
  res1 = selection_combination(
    as.matrix(dml_m_LM_UC[train_id_UC == 0, ]),
    as.matrix(dml_m_LC_UC[train_id_UC == 0, ]),
    weight_LM[train_id_LM == 0], Y_LM[train_id_LM == 0],
    as.matrix(dml_m_LM_LM[train_id_LM == 0, ]),
    weight_LC[train_id_LC == 0], Y_LC[train_id_LC == 0],
    as.matrix(dml_m_LC_LC[train_id_LC == 0, ]),
    V_UC[train_id_UC == 0, ], V_LM[train_id_LM == 0, ], V_LC[train_id_LC == 0, ]
  )
  
  res2 = selection_combination(
    as.matrix(dml_m_LM_UC), as.matrix(dml_m_LC_UC),
    weight_LM, Y_LM, as.matrix(dml_m_LM_LM),
    weight_LC, Y_LC, as.matrix(dml_m_LC_LC),
    V_UC, V_LM, V_LC
  )
  
  res3 = selection_eval(
    as.matrix(dml_m_LM_UC[train_id_UC == 1, ]),
    as.matrix(dml_m_LC_UC[train_id_UC == 1, ]),
    weight_LM[train_id_LM == 1], Y_LM[train_id_LM == 1],
    as.matrix(dml_m_LM_LM[train_id_LM == 1, ]),
    weight_LC[train_id_LC == 1], Y_LC[train_id_LC == 1],
    as.matrix(dml_m_LC_LC[train_id_LC == 1, ]),
    V_UC[train_id_UC == 1, ], V_LM[train_id_LM == 1, ],
    V_LC[train_id_LC == 1, ], res1$alpha
  )
  
  return(list(res1 = res3$value, res2 = res3$var))
}

violation2 = function(Combined_X, Combined_Y, source_index, seednum) {
  value = rep(0, 3)
  var = rep(0, 3)
  X_LC = Combined_X[source_index[[1]], ]
  Y_LC = Combined_Y[source_index[[1]]]
  X_UC = Combined_X[source_index[[5]], ]
  
  for (lmid in 1:3) {
    error = discrepancy2(
      Y_LC, X_LC,
      Combined_X[source_index[[lmid + 1]], ],
      Combined_Y[source_index[[lmid + 1]]],
      X_UC, seednum
    )
    value[lmid] = max(error$res1)
    var[lmid] = max(error$res2)
  }
  return(list(value = value, var = var))
}

discrepancy2 = function(Y_LC, X_LC, X_LM, Y_LM, X_UC, seednum) {
  set.seed(seednum + 1)
  covariate_z = 1:4
  n_fold = 2
  
  fold_id_LC = sample(rep(1:n_fold, length.out = nrow(X_LC)))
  fold_id_LM = sample(rep(1:n_fold, length.out = nrow(X_LM)))
  
  dml_m_LC_LC = matrix(0, nrow = nrow(X_LC), ncol = 1)
  dml_m_LM_LM = matrix(0, nrow = nrow(X_LM), ncol = 1)
  dml_m_LC_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  dml_m_LM_UC = matrix(0, nrow = nrow(X_UC), ncol = 1)
  weight_LC = rep(1, nrow(X_LC))
  weight_LM = rep(1, nrow(X_LM))
  
  dml_m_LC_UC = dml_m_LC_UC / n_fold
  dml_m_LM_UC = dml_m_LM_UC / n_fold
  
  V_LC = generate_poly_basis(X_LC[, covariate_z], 2)
  V_LM = generate_poly_basis(X_LM[, covariate_z], 2)
  V_UC = generate_poly_basis(X_UC[, covariate_z], 2)
  
  train_id_LC = sample(rep(0:1, length.out = nrow(V_LC)))
  train_id_LM = sample(rep(0:1, length.out = nrow(V_LM)))
  train_id_UC = sample(rep(0:1, length.out = nrow(V_UC)))
  
  res1 = selection_combination(
    as.matrix(dml_m_LM_UC[train_id_UC == 0, ]),
    as.matrix(dml_m_LC_UC[train_id_UC == 0, ]),
    weight_LM[train_id_LM == 0], Y_LM[train_id_LM == 0],
    as.matrix(dml_m_LM_LM[train_id_LM == 0, ]),
    weight_LC[train_id_LC == 0], Y_LC[train_id_LC == 0],
    as.matrix(dml_m_LC_LC[train_id_LC == 0, ]),
    V_UC[train_id_UC == 0, ], V_LM[train_id_LM == 0, ], V_LC[train_id_LC == 0, ]
  )
  
  res2 = selection_combination(
    as.matrix(dml_m_LM_UC), as.matrix(dml_m_LC_UC),
    weight_LM, Y_LM, as.matrix(dml_m_LM_LM),
    weight_LC, Y_LC, as.matrix(dml_m_LC_LC),
    V_UC, V_LM, V_LC
  )
  
  res3 = selection_eval(
    as.matrix(dml_m_LM_UC[train_id_UC == 1, ]),
    as.matrix(dml_m_LC_UC[train_id_UC == 1, ]),
    weight_LM[train_id_LM == 1], Y_LM[train_id_LM == 1],
    as.matrix(dml_m_LM_LM[train_id_LM == 1, ]),
    weight_LC[train_id_LC == 1], Y_LC[train_id_LC == 1],
    as.matrix(dml_m_LC_LC[train_id_LC == 1, ]),
    V_UC[train_id_UC == 1, ], V_LM[train_id_LM == 1, ],
    V_LC[train_id_LC == 1, ], res1$alpha
  )
  
  return(list(res1 = res3$value, res2 = res3$var))
}
