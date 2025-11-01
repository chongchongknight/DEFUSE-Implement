# generate the distribution for binary outcomes
mix_gaussian <- function(n_site, n_total, n_label, mu, sigma, imput_linear, y_correct, M, gamma.bar, seednum,  generate_mu, prop){
  prop = 0.3
  set.seed(seednum)
  # K is 2, 
  x = vector("list", n_site)
  y = vector("list", n_site)
  p = vector("list", n_site)
  for (i in 1:n_site) {
    x[[i]] = matrix(0, n_total[i], length(sigma[, 1]))
    y[[i]] = rbinom(n_total[i], 1, 0.5)
    Mu_obs1 = generate_mu[1:3]
    Sigma_obs1 = diag(0.5, 3)
    Mu_obs0 = generate_mu[4:6]
    Sigma_obs0 = diag(0.5, 3)
    x1 = rmvnorm(n_total[i], mean = Mu_obs1, sigma = Sigma_obs1)
    x0 = rmvnorm(n_total[i], mean = Mu_obs0, sigma = Sigma_obs0)
    x[[i]][which(y[[i]] == 1), 1:3] = x1[which(y[[i]] == 1), ]
    x[[i]][which(y[[i]] == 0), 1:3] = x1[which(y[[i]] == 0), ]
    # linear or non_linear imputation of missing data
    if (imput_linear == 1) {
      x[[i]][which(y[[i]] == 1), M[[2]][1]] = x[[i]][which(y[[i]] == 1), -M[[2]]] %*% c(1, -2, 3) * prop + rnorm(length(which(y[[i]] == 1)), sd = 1)
      x[[i]][which(y[[i]] == 0), M[[2]][1]] = x[[i]][which(y[[i]] == 0), -M[[2]]] %*% c(-1, 3, -2) * prop + rnorm(length(which(y[[i]] == 0)), sd = 1)
      x[[i]][which(y[[i]] == 1), M[[2]][2]] = x[[i]][which(y[[i]] == 1), -M[[2]]] %*% c(2, -2, -1) * prop + rnorm(length(which(y[[i]] == 1)), sd = 1)
      x[[i]][which(y[[i]] == 0), M[[2]][2]] = x[[i]][which(y[[i]] == 0), -M[[2]]] %*% c(-3, 2, 2) * prop + rnorm(length(which(y[[i]] == 0)), sd = 1)
    } else {
      x[[i]][which(y[[i]] == 1), M[[2]][1]] = x[[i]][which(y[[i]] == 1), -M[[2]]] ^ 2 %*% c(1, -2, 3) * prop + rnorm(length(which(y[[i]] == 1)), sd = 1)
      x[[i]][which(y[[i]] == 0), M[[2]][1]] = x[[i]][which(y[[i]] == 0), -M[[2]]] ^ 2 %*% c(-1, 3, -2) * prop + rnorm(length(which(y[[i]] == 0)), sd = 1)
      x[[i]][which(y[[i]] == 1), M[[2]][2]] = x[[i]][which(y[[i]] == 1), -M[[2]]] ^ 2 %*% c(2, -2, -1) * prop + rnorm(length(which(y[[i]] == 1)), sd = 1)
      x[[i]][which(y[[i]] == 0), M[[2]][2]] = x[[i]][which(y[[i]] == 0), -M[[2]]] ^ 2 %*% c(-3, 2, 2) * prop + rnorm(length(which(y[[i]] == 0)), sd = 1)
    }
    y[[i]] = y[[i]][1:n_label[i]]
  }
  return(list(X = x, Y = y))
}
