set.seed(1)

library(mvtnorm)
beta_int <- cbind(seq(3.7, 3.8, length.out = 20), seq(2.9, 3.1, length.out = 20))
vim_mult <- rep(NA, 2)
corrs <- c(0, 0.8)
n <- 5 * 10^6
int_vims <- matrix(NA, dim(beta_int)[1], 2)
for (l in 1:dim(beta_int)[1]) {
  cbeta <- beta_int[l, ]
  for (i in 1:length(corrs)) {
    sigma <- matrix(corrs[i], 2, 2)
    diag(sigma) <- 1
    mu <- c(0, 0)
    bisamp <- rmvnorm(n, mu, sigma)
    fluc <- rnorm(n)
    samp <- cbind(bisamp, fluc)
    int_vim <- apply(samp, 1, function(x) {
      (x[1] - x[3])^2 * (1 + cbeta[i] * x[2])^2
    })
    int_vims[l, i] <- mean(int_vim)
  }
}
save(int_vims, file = "int_vims.RData")
beta_int_f <- c(beta_int[, 1][which.min(abs(30 - int_vims[, 1]))], beta_int[, 2][which.min(abs(30 - 
  int_vims[, 2]))])

n <- 10^8
cube_coef <- seq(0.995, 1.05, length.out = 20)
cube_vims <- rep(NA, length(cube_coef))
for (i in 1:length(cube_vims)) {
  vim_dat <- matrix(rnorm(2 * n), n, 2)
  ccoef <- cube_coef[i]
  vim_dat <- vim_dat^3
  vim_dat <- (vim_dat[, 1] - vim_dat[, 2])^2
  mvim_dat <- mean(vim_dat)
  ctrue_vim <- ccoef^2 * mvim_dat
  cube_vims[i] <- ctrue_vim
  print(i)
}
save(cube_vims, file = "cube_vims.RData")
cube_coef_f <- cube_coef[which.min(abs(30 - cube_vims))]
