### Data generation procedure for continuous response
mix_gaussian <- function(n_site, n_total, A, XOmgea_mimusA, XOmgea_nonsense, XGamma_use,
                         XGamma, XGamma_nonsense, Xunobs, seednum, mean_shift, dim,
                         XOmgea_mimusA_coef, XGamma_coef) {
  # A = c(1:3)
  # XOmgea_mimusA = c(4:5)
  # XOmgea_nonsense = c(6:15)
  # XGamma = c(16:17)
  # XGamma_nonsense = c(18:25)
  # Xunobs = c(26:30)
  # seednum = 1
  # mean_shift = 0.2
  # dim = 30
  set.seed(seednum)
  error_table = c(-mean_shift, 0, mean_shift)
  # K is 2, 
  x = vector("list", n_site)
  y = vector("list", n_site)
  n_XOmgea_mimusA = length(XOmgea_mimusA)
  n_XGamma = length(XGamma_use)
  for (i in 1:n_site) {
    error = error_table[i]
    x[[i]] = matrix(0, n_total[i], dim)
    x[[i]][, A] = matrix(rnorm(n_total[i] * length(A)) + error, n_total[i], length(A))
    for (j in 1:n_XOmgea_mimusA) {
      x[[i]][, XOmgea_mimusA[j]] = x[[i]][, A] %*% XOmgea_mimusA_coef[, j]+ rnorm(n_total[i],sd = 1)
    }
    x[[i]][, XOmgea_nonsense] = matrix(rnorm(n_total[i] * length(XOmgea_nonsense)), n_total[i], length(XOmgea_nonsense))
    x[[i]][, XGamma_nonsense] = matrix(rnorm(n_total[i] * length(XGamma_nonsense)), n_total[i], length(XGamma_nonsense))
    x[[i]][, Xunobs] = matrix(rnorm(n_total[i] * length(Xunobs)), n_total[i], length(Xunobs))
    y[[i]] = rnorm(n_total[i], g_x(x[[i]][, c(A, XOmgea_mimusA)], gamma.bar[i, 1:(length(c(A, XOmgea_mimusA)))]), 1)
    for (k in 1:n_XGamma) {
      x[[i]][, XGamma_use[k]] =  XGamma_coef[1, k] + XGamma_coef[2, k] * y[[i]] + rnorm(n_total[i],sd = 1)
    }
  }
  return(list(X = x, Y = y))
}