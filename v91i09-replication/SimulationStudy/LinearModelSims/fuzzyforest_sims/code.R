# id goes from 1 to 12000 for simulations run on a cluster of workstations.
# This code should work for any cluster that uses SGE.

id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(id)) id <- 1 # if run manually only do for id = 1

set.seed(id)

library(WGCNA)
library(randomForest)
library(mvtnorm)
library(fuzzyforest)

rep_num <- 100
keep_frac <- c(0.01, 0.05, 0.1, 0.15, 0.25)
drop_frac <- c(0.05, 0.1, 0.25, 0.5)
mtry_factor <- c(0.5, 1, 2)
p <- c(100, 1000)
n <- c(100)

param_list <- list(keep_frac, drop_frac, mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 5:1]
names(param_settings) <- c("n", "p", "mtry_factor", "drop_fraction", "keep_fraction")

param_settings
current_sim_params <- param_settings[ceiling(id/rep_num), ]

sim_number <- 5
sim_results <- list()
sim_mod <- function(n, p, corr) {
  sigma <- matrix(corr, nrow = p, ncol = p)
  diag(sigma) <- 1
  X <- rmvnorm(n, sigma = sigma)
  return(X)
}

n <- as.numeric(current_sim_params[1])
p <- as.numeric(current_sim_params[2])
mtry_factor <- as.numeric(current_sim_params[3])
keep_fraction <- as.numeric(current_sim_params[4])
drop_fraction <- as.numeric(current_sim_params[5])
corr <- 0.8
if (p == 100) {
  number_of_groups <- 4
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:3, 76:78)
  vim_interest <- c(1:4, 76:79)
  beta_list <- rep(c(5, 5, 2), 2)
}
if (p == 1000) {
  number_of_groups <- 10
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:3, 901:903)
  vim_interest <- c(1:4, 901:904)
  beta_list <- rep(c(5, 5, 2), 2)
}

for (l in 1:sim_number) {
  all_modules <- lapply(1:number_of_mods, function(j) sim_mod(n, p_per_group, corr))
  all_modules[[number_of_groups]] <- matrix(rnorm(p_per_group * n), nrow = n, ncol = p_per_group)
  X <- do.call(cbind, all_modules)
  beta <- rep(0, p_per_group * (number_of_mods + 1))
  beta[vim_list] <- beta_list
  y <- X %*% beta + rnorm(n, sd = 0.1)
  X <- as.data.frame(X)
  names(X) <- paste("V", 1:p, sep = "")
  mtry_factor <- 1
  screen_params <- screen_control(drop_fraction = drop_fraction, keep_fraction = keep_fraction, 
    mtry_factor = mtry_factor)
  select_params <- select_control(number_selected = 10, drop_fraction = drop_fraction, 
    mtry_factor = mtry_factor)
  y <- as.numeric(y)
  ff <- wff(X, y, screen_params = screen_params, select_params = select_params, 
    num_processors = 1, nodesize = 1)
  vim <- ff$feature_list
  sim_results[[l]] <- vim
}
summarize_sim <- function(sim_results, vim_interest) {
  num_interesting_feat <- length(vim_interest)
  vims <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  selected <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  for (r in 1:num_interesting_feat) {
    current_feature <- paste("V", vim_interest[r], sep = "")
    for (t in 1:length(sim_results)) {
      current_vims <- sim_results[[t]]
      if (current_feature %in% current_vims[, 1]) {
        current_rank <- which(current_vims[, 1] == current_feature)
        vims[r, t] <- current_vims[current_rank, 2]
        selected[r, t] <- 1
      } else {
        vims[r, t] <- 0
      }
    }
  }
  mean_vims <- apply(vims, 1, mean)
  selected_props <- apply(selected, 1, mean)
  sim_out <- cbind(mean_vims, selected_props)
  return(sim_out)
}

out <- list(summarize_sim(sim_results, vim_interest), vim_interest, current_sim_params)

dir.create("out", showWarnings = FALSE)
out_filename <- paste("out/out", id, sep = "")
save(out, file = out_filename)
