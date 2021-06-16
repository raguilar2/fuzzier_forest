#new version
files <- list.files("out", full.names = TRUE)
get_indices <- function(files) {
  split <- strsplit(files, "out/out")
  indices <- unlist(lapply(split, function(x) {
    x[2]
  }))
  indices <- as.numeric(indices)
  indices <- indices[is.na(indices) == FALSE]
  indicies <- sort(indices)
}

reorder_files <- function(files, n){
  bad_ordering <- sort(as.character(1:n))
  files <- files[order(as.numeric(bad_ordering))]
  files
}
files <- reorder_files(files, length(files))

rep_num <- 100
n <- c(100)
p <- c(100, 1000)
mtry_factor <- c(1, 3/4, 1/2, 1/3)

param_list <- list(mtry_factor, p, n)
param_settings <- expand.grid(param_list)
param_settings <- param_settings[, 3:1]
names(param_settings) <- c("n", "p", "mtry_factor")

total_sim_settings <- dim(param_settings)[1]
total_sims <- rep_num * total_sim_settings

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


model_index <- cbind(seq(1, total_sims, by = rep_num), seq(rep_num, total_sims, by = rep_num))
find_model_index <- function(id) {
  which(apply(model_index, 1, function(x) {
    (id >= x[1]) && (id <= x[2])
  } == TRUE))
}

results <- rep(list(matrix(0, 8, 2)), dim(param_settings)[1])
indices <- get_indices(files)
out_list <- vector("list", length(files))
for (i in 1:length(files)) {
  load(files[i])
  if(!all(unique(out[, "selected_props"]) %in% c(1, 0))){
    browser()
  }
  id <- indices[i]
  current_model <- find_model_index(id)
  results[[current_model]] <- results[[current_model]] + out
}

results <- lapply(results, function(x) {
  x/100
})

cif_results <- list(param_settings, results)
save(cif_results, file = "cif_results.RData")
