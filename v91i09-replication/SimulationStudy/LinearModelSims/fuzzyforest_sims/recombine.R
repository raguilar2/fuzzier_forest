files <- list.files("out", full.names = TRUE)
get_indices <- function(files) {
  split <- strsplit(files, "out/out")
  indices <- unlist(lapply(split, function(x) {
    x[2]
  }))
  indices <- as.numeric(indices)
  indices <- indices[is.na(indices) == FALSE]
  indices <- sort(indices)
}

reorder_files <- function(files, n){
  bad_ordering <- sort(as.character(1:n))
  files <- files[order(as.numeric(bad_ordering))]
  files
}
files <- reorder_files(files, length(files))

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

total_sim_settings <- dim(param_settings)[1]
total_sims <- rep_num * total_sim_settings

indices <- get_indices(files)
model_index <- cbind(seq(1, total_sims, by = rep_num), seq(rep_num, total_sims, by = rep_num))
find_model_index <- function(id) {
  which(apply(model_index, 1, function(x) {
    (id >= x[1]) && (id <= x[2])
  } == TRUE))
}

ff_results <- vector("list", total_sim_settings)
for (i in 1:length(ff_results)) {
  ff_results[[i]] <- matrix(0, 8, 2)
}
model_counts <- rep(0, total_sim_settings)

for (i in 1:length(files)) {
  load(files[i])
  id <- indices[i]
  current_model <- find_model_index(id)
  model_counts[current_model] <- model_counts[current_model] + 1
  ff_results[[current_model]] <- out[[1]] + ff_results[[current_model]]
}
results <- vector("list", 120)
for (i in 1:120) {
  results[[i]] <- ff_results[[i]]/model_counts[i]
}
fuzzyforest_results <- list(param_settings, results)
save(fuzzyforest_results, file = "fuzzyforest_results.RData")
