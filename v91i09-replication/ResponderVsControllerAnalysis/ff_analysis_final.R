rm(list=ls())
set.seed(1)
library(doParallel)
library(randomForest)
library(WGCNA)
library(fuzzyforest)

# First we extract the data for responders, non-responders, and controllers, Then
# we impute the missing data and write the imputed data set to a .csv file.  This
# section can be commented out after running once.
dat <- read.csv("ResponderVsController.csv")

# rename outcome to 'y'
names(dat)[which(names(dat) == "Group_Discrete")] <- "y"

# Responder:'ART (High CD4)' Non-Responder:'ART (Low CD4)' Controller: 'Elite
# Controller'
RNRC <- dat[dat$y %in% c("ART (High CD4)", "ART (Low CD4)", "Elite Controller"), ]
X <- RNRC
X <- X[, -1]
y <- droplevels(RNRC$y)

# remove features unrelated to flow cytometry and features that were used to
# construct the outcome
dat_imputed <- rfImpute(X, y, iter = 25, ntree = 10000, data = X)

write.csv(dat_imputed, file = "Imputed_ResponderVsController.csv")
dat_imputed <- read.csv("Imputed_ResponderVsController.csv")

y <- dat_imputed[, 2]
X <- dat_imputed[, -c(1:2)]

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
# this is insurance to get the code working in case pickSoftThreshold fails
registerDoParallel(cores = 1)
sft <- pickSoftThreshold(X, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "scale_free_plot.pdf", width = 6.69)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers,
  cex = cex1, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
  type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# calculate modules
net <- blockwiseModules(X, power = 8, TOMtype = "unsigned", deepSplit = 4, numbericLabels = TRUE,
  saveTOMs = FALSE, minModuleSize = 10)

RespVsController <- subset(dat_imputed, y %in% c("ART (High CD4)", "Elite Controller"))
y_resp_controller <- droplevels(RespVsController[, 2])
X_resp_controller <- RespVsController[, -c(1, 2)]
table(y_resp_controller)

# fit fuzzy forests
module_membership <- net$colors
fit <- ff(X = X_resp_controller, y = y_resp_controller, module_membership = module_membership,
  screen_params = screen_control(min_ntree = 6000, ntree_factor = 1, keep_fraction = 0.25),
  select_params = select_control(min_ntree = 6000, ntree_factor = 1, number_selected = 10),
  num_processors = 4)

pdf(file = "modplot.pdf", width = 6.69)
modplot(fit, main = "Module Membership for Controllers versus Responders")
dev.off()
fit_rf <- fit$final_rf

library(lattice)

# re-calculate dot plot with the feature names re-labeled so that they are more
# understandable.
X_select <- subset(X_resp_controller, select = fit$feature_list[, 1])

get_direction <- function(x) {
  grp_means <- aggregate(x ~ y_resp_controller, FUN = mean)
  return(as.numeric(grp_means[2, 2] > grp_means[1, 2]))
}

feature_names <- c("RA+R7+PD1+ (MFI on CD4+ naive cells)",
                   "RA-CCR7+CCR5-CD38+HLADR+PD1- (% of CD4+ T cells)",
                   "IL2+ (% HIV gag-specific CD8+ T cells)",
                   "PD1+ (MFI on CD4+ T cells)",
                   "TNF+ (% HIV gag-specific CD8+ cells)",
                   "CD107+IFN+IL2+TNF+ (% of HIV gag specific CD4+ cells)",
                   "CD107-IFN+IL2+TNF+ (% of HIV gag specific CD4+ cells)",
                    "RA-CCR7+PD1+ (MFI on CD8+ T cells)",
                    "IFN+ (% of CD4+ cells)",
                   "RA+R7-PD1+ (MFI on CD4+ T cells)")

directions <- apply(X_select, 2, get_direction)
direction_cols <- rep(NA, length(directions))
direction_cols[directions == 1] <- "red"
direction_cols[directions == 0] <- "black"
vims <- fit$feature_list$variable_importance
dotchart(rev(vims), labels = rev(feature_names), main = "Variable Importance Chart",
  xlab = "Variable Importance", color = rev(direction_cols), pch = 19, cex = 0.75)
pdf(file = "vimplot.pdf", width=12, height=8)
dotchart(rev(vims), labels = rev(feature_names), main = "Variable Importance Chart",
  xlab = "Variable Importance", color = rev(direction_cols), pch = 19, cex = 0.75)
dev.off()

