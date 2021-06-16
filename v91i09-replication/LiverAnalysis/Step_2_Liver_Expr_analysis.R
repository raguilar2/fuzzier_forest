set.seed(1)
library(randomForest)
library(WGCNA)
library(fuzzyforest)
library(doParallel)
load("final_Liver_Expr.RData")
Liver_Expr <- final_Liver_Expr
weight <- Liver_Expr[, 1]
expression_levels <- Liver_Expr[, -1]

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function Insurance over pickSoftThreshold
registerDoParallel(cores = 1)
sft <- pickSoftThreshold(expression_levels, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)

pdf("scale_free_plot.pdf")
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers,
  cex = cex1, col = "red")

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.9, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
  type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

net <- blockwiseModules(expression_levels, power = 9, TOMType = "unsigned", minModuleSize = 15,
  reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = FALSE, pamRespectsDendro = FALSE,
  verbose = 3)
module_membership <- net$colors
mtry_factor <- 1
drop_fraction <- 0.25
number_selected <- 10
keep_fraction <- 0.05
min_ntree <- 5000
ntree_factor <- 5
final_ntree <- 5000
screen_params <- screen_control(drop_fraction = drop_fraction, keep_fraction = keep_fraction,
  min_ntree = min_ntree, mtry_factor = mtry_factor, ntree_factor = ntree_factor)
select_params <- select_control(drop_fraction = drop_fraction, number_selected = number_selected,
  min_ntree = min_ntree, mtry_factor = mtry_factor, ntree_factor = ntree_factor)

# WGCNA fuzzy forest
WGCNA_params <- WGCNA_control(power = 6, minModuleSize = 30, TOMType = "unsigned", reassignThreshold = 0,
  mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE)

wff_fit <- wff(expression_levels, weight, WGCNA_params = WGCNA_params, screen_params = screen_params,
  select_params = select_params, final_ntree = final_ntree, num_processors = 1)

modplot(wff_fit)

pdf("FemaleLiverVimPlot.pdf")
modplot(wff_fit)
dev.off()

cex1 <- 0.9
wff_fit$feature_list
