library(ggplot2)

load("LinearModelSims/fuzzyforest_sims/fuzzyforest_results.RData")
ff_results <- fuzzyforest_results[[2]]
ff_param_settings <- fuzzyforest_results[[1]]

load("LinearModelSims/conditional_vim_sims/cif_results.RData")
cif_param_settings <- cif_results[[1]]
cif_param_settings <- data.frame(cif_param_settings, mtry = floor(cif_param_settings$p * 
  cif_param_settings$mtry_factor))
cif_results <- cif_results[[2]]

load("LinearModelSims/randomforest_sims/rf_results.RData")
rf_param_settings <- rf_results[[1]]
rf_param_settings <- data.frame(rf_param_settings, mtry = floor(rf_param_settings$p * 
  rf_param_settings$mtry_factor))
rf_results <- rf_results[[2]]

# choose the best result from each result
parameter_settings <- vector("list", 3)

ranking_error <- function(est_prop) {
  true_membership <- c(rep(1, 3), 0, rep(1, 3), 0)
  out <- mean((true_membership - est_prop)^2)
}
get_best_params <- function(param_settings, results, p) {
  p_ind <- which(param_settings$p == p)
  rankings <- rep(NA, length(p_ind))
  for (i in 1:length(p_ind)) {
    rankings[i] <- ranking_error(results[[p_ind[i]]][, "selected_props"])
  }
  best_rank_ind <- p_ind[which.min(rankings)]
  best_param_setting <- param_settings[best_rank_ind, ]
  best_result <- results[[best_rank_ind]]
  out <- list(param_setting = best_param_setting, result = best_result)
}

best_ff_results <- get_best_params(ff_param_settings, ff_results, p = 100)
ff_best <- best_ff_results[[2]]
parameter_settings[[1]] <- best_ff_results[[1]]

best_cif_results <- get_best_params(cif_param_settings, cif_results, p = 100)
cif_best <- best_cif_results[[2]]
parameter_settings[[2]] <- best_cif_results[[1]]

best_rf_results <- get_best_params(rf_param_settings, rf_results, p = 100)
rf_best <- best_rf_results[[2]]
parameter_settings[[3]] <- best_rf_results[[1]]

names(parameter_settings) <- c("Fuzzy Forests", "Conditional Inference Forests", 
  "Random Forests")
results <- as.data.frame(rbind(ff_best, cif_best, rf_best))
results <- cbind(rep(paste("X", c(1:4, 76:79), sep = ""), 3), results)
results <- cbind(rep(names(parameter_settings), each = 8), results)
names(results)[1:2] <- c("Method", "feature")
results[, 3:4] <- round(results[, 3:4], 2)
results[, 2] <- as.factor(results[, 2])
results[, 1] <- factor(results[, 1], levels(results[, 1])[c(2, 3, 1)])
p <- ggplot(results, aes(x = feature, y = selected_props, colour = Method, group = Method)) + 
  geom_line(size = 1.1)
p <- p + geom_point(size = 5, alpha = 0.8) + labs(x = "Feature", y = "Proportion of Times Feature was Selected", 
  title = "Feature Selection Performance n=100, p=100") + geom_line(size = 2, alpha = 0.8)
p <- p + theme(plot.title = element_text(hjust = 0.5))
pdf("lm_p100plot.pdf", width = 8.2)
plot(p)
dev.off()
parameter_settings <- vector("list", 3)


best_ff_results <- get_best_params(ff_param_settings, ff_results, p = 1000)
ff_best <- best_ff_results[[2]]
parameter_settings[[1]] <- best_ff_results[[1]]

best_cif_results <- get_best_params(cif_param_settings, cif_results, p = 1000)
cif_best <- best_cif_results[[2]]
parameter_settings[[2]] <- best_cif_results[[1]]

best_rf_results <- get_best_params(rf_param_settings, rf_results, p = 1000)
rf_best <- best_rf_results[[2]]
parameter_settings[[3]] <- best_rf_results[[1]]

names(parameter_settings) <- c("Fuzzy Forests", "Conditional Inference Forests", 
  "Random Forests")
results <- as.data.frame(rbind(ff_best, cif_best, rf_best))
results <- cbind(rep(paste("X", c(1:4, 901:904), sep = ""), 3), results)
results <- cbind(rep(names(parameter_settings), each = 8), results)
names(results)[1:2] <- c("Method", "feature")
results[, 3:4] <- round(results[, 3:4], 2)
results[, 2] <- as.factor(results[, 2])
results[, 1] <- factor(results[, 1], levels(results[, 1])[c(2, 3, 1)])
p <- ggplot(results, aes(x = feature, y = selected_props, colour = Method, group = Method)) + 
  geom_line(size = 2, alpha = 0.7) + geom_point(size = 5, alpha = 0.7)
p <- p + labs(x = "Feature", y = "Proportion of Times Feature was Selected", title = "Feature Selection Performance n=100, p=1,000")
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + ylim(0, 1)
pdf("lm_p1000plot.pdf", width = 8.2)
plot(p)
dev.off()

# plot nonlinear simulations
load("NonLinearModelSims/nl_simspaper/nl_sims.RData")
n250dat <- rep(paste("X", c(1:5, 76:80), sep = ""), 3)
n250dat <- cbind(rep(c("Fuzzy Forests", "Random Forests", "Conditional Inference Forests"), 
  each = 10), n250dat)
n250sim <- nl_sims[[1]]
n250props <- do.call("c", n250sim)
n250dat <- cbind(n250dat, n250props)
row.names(n250dat) <- NULL
n250dat <- as.data.frame(n250dat, stringsAsFactors = TRUE)
n250dat[, 1] <- factor(n250dat[, 1], levels(n250dat[, 1])[c(2, 3, 1)])
names(n250dat) <- c("Method", "feature", "selected_props")
n250dat[, 3] <- as.numeric(as.character(n250dat[, 3]))
p <- ggplot(n250dat, aes(x = feature, y = selected_props, colour = Method, group = Method)) + 
  geom_line(size = 2, alpha = 0.7) + geom_point(size = 5, alpha = 0.7)
p <- p + labs(x = "Feature", y = "Proportion of Times Feature was Selected", title = "Feature Selection Performance n=250, p=100")
p <- p + ylim(0, 1)
p <- p + theme(plot.title = element_text(hjust = 0.5))
pdf("nlm_n250plot.pdf", width = 8.2)
plot(p)
dev.off()

library(scales)
ggpal <- c("#F8766D", "#00BA38")
n500dat <- rep(paste("X", c(1:5, 76:80), sep = ""), 2)
n500dat <- cbind(rep(c("Fuzzy Forests", "Random Forests"), each = 10), n500dat)
n500sim <- nl_sims[[2]]
n500props <- do.call("c", n500sim)
n500dat <- cbind(n500dat, n500props)
row.names(n500dat) <- NULL
n500dat <- as.data.frame(n500dat, stringsAsFactors = TRUE)
names(n500dat) <- c("Method", "feature", "selected_props")
n500dat[, 3] <- as.numeric(as.character(n500dat[, 3]))
p <- ggplot(n500dat, aes(x = feature, y = selected_props, colour = Method, group = Method)) + 
  geom_line(size = 2, alpha = 0.7) + geom_point(size = 5, alpha = 0.7)
p <- p + labs(x = "Feature", y = "Proportion of Times Feature was Selected", title = "Feature Selection Performance n=500, p=100")
p <- p + ylim(0, 1)
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + scale_color_manual(values = ggpal)
pdf("nlm_n500plot.pdf", width = 8.2)
plot(p)
dev.off()
