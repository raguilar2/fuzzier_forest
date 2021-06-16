# This file is designed to demonstrate how the data set used in the vignette was
# created from the raw data set. For information about the data set go to the the
# following site:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# Go to the following part of the page: I. Network analysis of liver expression
# data from female mice: finding modules related to body weight

set.seed(323516)
library(randomForest)
library(WGCNA)
library(fuzzyforest)
library(flashClust)
disableWGCNAThreads()

# Read in the female liver data set
femData <- read.csv("LiverFemale3600.csv", stringsAsFactors = FALSE)

# Take a quick look at what is in the data set:
dim(femData)
names(femData)

# remove columns that do not gene expression.
datExpr0 <- as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) <- femData$substanceBXH
rownames(datExpr0) <- names(femData)[-c(1:8)]

# remove genes and samples of low quality
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
      collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples],
      collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- flashClust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches The
# user should change the dimensions if the window is too large or too small
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
  cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red")

# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# read in clinical trait data set
traitData <- read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits <- traitData[, -c(31, 16)]
allTraits <- allTraits[, c(2, 11:36)]
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical
# traits
femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Mice)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# n is 134
# datExpr holds the expression level data, there are 3600 features in it
# datTraits contains data on traits

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
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

net <- blockwiseModules(datExpr, power = 6, TOMType = "unsigned", minModuleSize = 30)
  #reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = FALSE, pamRespectsDendro = FALSE,
  #saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM", verbose = 3)

y <- datTraits$weight_g
X <- datExpr
missing_y <- which(is.na(y))
X <- X[-missing_y, ]
y <- y[-missing_y]
X_impute <- rfImpute(X, y, iter = 5, ntree = 5000)
final_Liver_Expr <- X_impute
save(final_Liver_Expr, file = "final_Liver_Expr.RData")
