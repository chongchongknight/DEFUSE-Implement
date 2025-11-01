mix_gaussian <- function(n_site, n_total, n_label, mu, sigma, imput_linear, y_correct, M, gamma.bar, seednum){
  set.seed(seednum)
  # K is 2, 
  x = vector("list", n_site)
  y = vector("list", n_site)
  p = vector("list", n_site)
  for (i in 1:n_site) {
    # x[[i]] = mvrnorm(n_total[i], Mu, Sigma)
    x[[i]] = matrix(runif(n_total[i] * length(sigma[,1]), -1, 1), n_total[i], length(sigma[,1]))
    # linear or non_linear imputation of missing data
    if (imput_linear == 1) {
      x[[i]][,M[[2]][1]] = x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5) + rnorm(n_total[i],sd = 1)
      x[[i]][,M[[2]][2]] = x[[i]][,-M[[2]]] %*% c(0.5, 0.5, 0.5) + rnorm(n_total[i],sd = 1)
    } else {
      x[[i]][,M[[2]][1]] = rnorm(n_total[i]) + sin(x[[i]][,1])
      x[[i]][,M[[2]][2]] = rnorm(n_total[i]) + sin(x[[i]][,2])
    }
    # misspecified or correct-specified of y model
    if (y_correct == 1) {
      p[[i]] = g_x(x[[i]], gamma.bar)
    } else {
      p[[i]] = g_x(cbind(x[[i]],x[[i]][,1]^2,x[[i]][,2]^2), c(gamma.bar,1,1))
    }
    for (j in 1:n_label[i]) {
      y[[i]][j] = rnorm(1, p[[i]][j], 1)
    }
  }
  return(list(X = x, Y = y))
}





