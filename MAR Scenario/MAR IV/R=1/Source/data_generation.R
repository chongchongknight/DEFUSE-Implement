### Data generation procedure for continuous response
mix_gaussian <- function(n_site, n_total, n_label, mu, sigma, imput_linear, y_correct, M, gamma.bar, seednum,
                         mean_shift, nonlinear){
  set.seed(seednum)
  error_table = c(-mean_shift, 0, mean_shift)
  x = vector("list", n_site)
  y = vector("list", n_site)
  p = vector("list", n_site)
  for (i in 1:n_site) {
    error = error_table[i]
    x[[i]] = matrix(rnorm(n_total[i] * nrow(sigma)) + error, n_total[i], nrow(sigma))
    if (imput_linear == 1) {
      x[[i]][,M[[2]][1]] = x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5) + rnorm(n_total[i],sd = 1)
      x[[i]][,M[[2]][2]] = x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5) + rnorm(n_total[i],sd = 1)
    } else {
      x[[i]][,M[[2]][1]] = rnorm(n_total[i]) + sin(x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5))
      x[[i]][,M[[2]][2]] = rnorm(n_total[i]) + sin(x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5))
    }
    if (y_correct == 1) {
      p[[i]] = g_x(x[[i]], gamma.bar[i, ])
    } else {
      p[[i]] = g_x(cbind(x[[i]],x[[i]][,1]^2, x[[i]][,2]^2), c(gamma.bar[i, ], nonlinear, nonlinear))
    }
    for (j in 1:n_label[i]) {
      y[[i]][j] = rnorm(1, p[[i]][j], 1)
    }
  }
  return(list(X = x, Y = y))
}