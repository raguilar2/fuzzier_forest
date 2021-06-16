# In contrast to the linear simulations, the nonlinear simulations generate
# results from conditional inference forests, random forests, and fuzzy forests
# all at once.

# Note: requires results from VIM simulations which are expected to reside in
# ../../VimSimulations/.

# id goes from 1 to 400 for simulations run on a cluster of workstations.
# This code should work for any cluster that uses SGE.

id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(id)) id <- 1 # if run manually only do for id = 1

set.seed(id)

library(WGCNA)
library(mvtnorm)
library(Matrix)
library(fuzzyforest)
library(party)
library(randomForest)
sim_number <- 5
sim_results <- list()
# Each row corresponds to one sample.
gen_mod <- function(n, p, corr) {
  sigma <- matrix(corr, nrow = p, ncol = p)
  diag(sigma) <- 1
  X <- rmvnorm(n, sigma = sigma)
  return(X)
}

beta_int <- cbind(seq(3.7, 3.8, length.out = 20), seq(2.9, 3.1, length.out = 20))
load("../../VimSimulations/int_vims.RData")
beta_int_f <- c(beta_int[, 1][which.min(abs(30 - int_vims[, 1]))], beta_int[, 2][which.min(abs(30 - 
  int_vims[, 2]))])

cube_coef <- seq(0.995, 1.05, length.out = 20)
load("../../VimSimulations/cube_vims.RData")
cube_coef_f <- cube_coef[which.min(abs(30 - cube_vims))]

beta <- c(1, beta_int_f, sqrt(15), cube_coef_f)
# note that the interaction term requires different coefficients for the
# correlated group and the uncorrelated group.
beta_list <- list(beta1 = beta[-2], beta2 = beta[-3])

# f_calc calculates f(x) for a vector, x var_block_ind is where block of
# important variables begins
f_calc <- function(x, beta_list, var_block_ind) {
  lp <- 0
  for (i in 1:length(beta_list)) {
    cbeta <- beta_list[[i]]
    cx <- x[var_block_ind[i]:(var_block_ind[i] + length(cbeta) - 1)]
    lp <- lp + cx[1] * cbeta[1] + cx[2] * cbeta[1] + cx[1] * cx[2] * cbeta[2] + 
      cx[3] * beta[3] + cx[4]^3 * cbeta[4]
  }
  lp
}

# this function generates the data mod_sizes is vector giving the size of each of
# m modules
gen_reg <- function(n, beta_list, var_block_ind, mod_sizes, corr, sd_err) {
  m <- length(mod_sizes)
  X_list <- vector("list", length = m)
  for (i in 1:m) {
    X_list[[i]] <- gen_mod(n, mod_sizes[i], corr[i])
  }
  X <- do.call("cbind", X_list)
  f <- apply(X, 1, f_calc, beta_list = beta_list, var_block_ind = var_block_ind)
  y <- f + rnorm(n, sd_err)
  dat <- cbind(y, X)
  return(dat)
}

## Test Fuzzy Forests
var_block_ind <- c(1, 76)
mod_sizes <- rep(25, 4)
corr <- rep(0.8, 4)
corr[4] <- 0
sim_num <- 5
cif_mat <- matrix(NA, sim_num, 10)
ff_mat <- matrix(NA, sim_num, 10)
rf_mat <- matrix(NA, sim_num, 10)
n <- 250
if (id > 200) {
  n <- 500
}

for (k in 1:sim_num) {
  sd_err <- 0.5
  dat <- gen_reg(n, beta_list, var_block_ind, mod_sizes, corr, sd_err)
  dat <- as.data.frame(dat)
  names(dat) <- c("y", paste("V", 1:100, sep = ""))
  X <- dat[, 2:101]
  y <- dat[, 1]
  
  rf <- randomForest(X, y, importance = T, mtry = 100)
  rf_list <- importance(rf, type = 1, scale = F)
  rf_list <- rf_list[order(rf_list[, 1], decreasing = T), ]
  rf_list <- rf_list[1:10]
  rf_mat[k, ] <- names(rf_list)
  
  screen_params <- screen_control(keep_fraction = 0.25, ntree_factor = 1, mtry_factor = 15, 
    min_ntree = 500)
  select_params <- select_control(number_selected = 10, drop_fraction = 0.1, ntree_factor = 1, 
    mtry_factor = 15, min_ntree = 500)
  wff_fit <- wff(X, y, select_params = select_params, screen_params = screen_params)
  feat_list <- paste("V", c(1:4, 76:79), sep = "")
  ff_list <- wff_fit$feature_list[, 1]
  ff_mat[k, ] <- ff_list
  
  if (n < 500) {
    cif <- cforest(y ~ ., data = dat, controls = cforest_unbiased(ntree = 100, 
      mtry = 100))
    cif_list <- varimp(cif, conditional = T)
    
    cif_list <- sort(cif_list, decreasing = T)
    cif_list <- cif_list[1:10]
    cif_mat[k, ] <- names(cif_list)
  }
  
  # summarize the results
  true_feat <- paste("V", c(1:5, 76:80), sep = "")
  # get pct of times each feature is selected
  feat_pct <- function(v, x) {
    sim_num <- dim(x)[1]
    feat_count <- apply(x, 1, function(l) {
      v %in% l
    })
  }
  rf_feat_pct <- sapply(true_feat, feat_pct, x = rf_mat)
  ff_feat_pct <- sapply(true_feat, feat_pct, x = ff_mat)
  if (n < 500) {
    cif_feat_pct <- sapply(true_feat, feat_pct, x = cif_mat)
  }
}

if (n >= 500) {
  cif_feat_pct <- NULL
}

dir.create("out", showWarnings = FALSE)
out_filename <- paste("out/out", id, sep = "")
out <- list(rf = rf_feat_pct, ff = ff_feat_pct, cif = cif_feat_pct)
save(out, file = out_filename)
