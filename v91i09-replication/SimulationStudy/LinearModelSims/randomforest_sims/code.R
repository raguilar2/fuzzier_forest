# id goes from 1 to 4000 for simulations run on a cluster of workstations.
# This code should work for any cluster that uses SGE.

id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(id)) id <- 3900 # if run manually only do for id = 1

set.seed(id)

library(randomForest)
library(mvtnorm)

rep_num <- 500
n <- c(100)
p <- c(100, 1000)
mtry_factor <- c(1, 3/4, 1/2, 1/3)

param_list <- list(mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 3:1]
names(param_settings) <- c("n", "p", "mtry_factor")

param_settings
current_sim_params <- param_settings[ceiling(id/rep_num), ]

sim_number <- 1
sim_results <- list()
sim_mod <- function(n, p, corr) {
  sigma <- matrix(corr, nrow = p, ncol = p)
  diag(sigma) <- 1
  X <- rmvnorm(n, sigma = sigma)
  return(X)
}

corr <- 0.8
n <- as.numeric(current_sim_params[1])
p <- as.numeric(current_sim_params[2])
mtry_factor <- as.numeric(current_sim_params[3])

if (p == 10) {
  number_of_groups <- 5
  number_of_mods <- number_of_groups - 1
  p_per_group <- p/number_of_groups
  vim_list <- c(1:2, 9:10)
  vim_interest <- c(1:3, 8:10)
  beta_list <- rep(c(5, 5), 2)
  keep_fraction <- 0.25
  keep_fraction <- 0.5
  keep_fraction <- 0.75
  keep_fraction <- 1
}
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
  y <- as.numeric(y)
  df <- cbind(y, X)
  fit <- randomForest(X, y, mtry = floor(p * mtry_factor), ntree = 5000, importance = T)
  vim <- importance(fit, scale = F, type = 1)
  vim <- as.data.frame(cbind(rownames(vim), vim))
  vim[, 2] <- as.numeric(as.character(vim[, 2]))
  names(vim) <- c("feature_name", "variable_importance")
  row.names(vim) <- NULL
  vim <- vim[order(vim[, 2], decreasing = T), ]
  sim_results[[l]] <- vim
}

summarize_sim <- function(sim_results, vim_interest, num_selected = 10) {
  num_interesting_feat <- length(vim_interest)
  vims <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  selected <- matrix(0, nrow = num_interesting_feat, ncol = length(sim_results))
  for (r in 1:num_interesting_feat) {
    current_feature <- paste("V", vim_interest[r], sep = "")
    for (t in 1:length(sim_results)) {
      current_vims <- sim_results[[t]]
      if (current_feature %in% current_vims[1:num_selected, 1]) {
        selected[r, t] <- 1
      }
      current_rank <- which(current_vims[, 1] == current_feature)
      vims[r, t] <- current_vims[current_rank, 2]
    }
  }
  mean_vims <- apply(vims, 1, mean)
  selected_props <- apply(selected, 1, mean)
  sim_out <- cbind(mean_vims, selected_props)
  return(sim_out)
}

out <- summarize_sim(sim_results, vim_interest, num_selected = 10)

dir.create("out", showWarnings = FALSE)
out_filename <- paste("out/out", id, sep = "")
save(out, file = out_filename)
