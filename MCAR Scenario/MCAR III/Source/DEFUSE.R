double_projection = function(y_form, imput, sam, gamma.bar, seed) {
  #y_form = 1
  #imput = 0
  #sam = c(sample[[1]])
  #gamma.bar = c(coefficient[[j]])
  #seed = 6
  set.seed(seed)
  K = 2
  N = c(rep(sam[1],K))
  n = c(sam[2], sam[3])
  dim = 5
  M = vector("list", length = K)
  M[[1]] = NA
  M[[2]] = c(4,5)
  nlm = c(500, 500)
  Mu = c(0,0,0,0,0)
  rho = 0
  cor_1 = matrix(rho, nrow = dim, ncol = dim)
  var = 1
  cor_1[dim,] = 0
  cor_1[,dim] = 0
  for (i in 1:dim) {
    cor_1[i,i] = var
  }
  Sigma = cor_1
  gamma.bar = as.matrix(gamma.bar, ncol = 1)
  y_form = y_form
  imput = imput
  gamma.true = gamma.bar
  dat = mix_gaussian(n_site = K, n_total = N, n_label = n, mu = Mu, 
                     sigma = Sigma, imput_linear = imput, y_correct = y_form, 
                     M = M, gamma.bar = gamma.bar, seed)
  X = dat$X
  Y = dat$Y
  mod = lm(Y[[1]] ~ X[[1]][1:n[1],] - 1)
  coef = matrix(0, 3, dim)
  
  
  gamma.tilt = mod$coefficients
  sigma_y = sqrt(mean(mod$residuals ^ 2))
  H_inv = solve(H(X[[1]][1:n[1],], gamma.tilt))
  S = t(H_inv %*% t(X[[1]][1:n[1],])) * (Y[[1]] - g_x(X[[1]], gamma.tilt)[1:n[1]])
  folds <- sample(cut(seq(1,n[[1]]),breaks = 5,labels = FALSE))
  linear_result <- main_func(X, Y, n, N, M, folds, imput_method = "lm", gamma.tilt, H_inv, S, sigma_y, gamma.true, nlm, seed)
  coef[2, ] = linear_result$gamma.tilt1
  coef[3, ] = linear_result$gamma.tilt2
  coef[1, ] = gamma.tilt
  res = list(coef = coef)
  return(res)
}
