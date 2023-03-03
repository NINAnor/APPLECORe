# Import the pollen deposition data

library(here)       # Library for relative links
library(readxl)     # Library to read in Excel files
library(tidyverse)  # Import tidyverse functions
library(brms)       # Import the brms library for diverse model specification
library(DHARMa)     # Import the model checking library
library(ggplot2)    # Import plotting functions
library(parallel)   # Import the parallel libraries
# Import the pollen deposition data from the file
pollenDepData <- read_excel("Work/APPLECORe/PollenDeposition/data_deposition.xlsx", "Ark1", na = c("", "NA", "missing"))
# Create a Genus column so that Genus can be analysed
pollenDepData$Genus <- sapply(X = strsplit(pollenDepData$Species, " ", fixed = TRUE), FUN = function(curVec) { curVec[1] })

# Set the MCMC parameters
numChains <- 4
numIterations <- 10000
numBurnIn <- 5000
numThin <- 1
numCores <- min(detectCores(), numChains)

# Set the colour palette for the bee types/genera
beeTypeColour <- setNames(c(
  rgb(255, 170, 125, maxColorValue = 255),
  rgb(255, 236, 139, maxColorValue = 255),
  rgb(189, 183, 107, maxColorValue = 255)),
c("Bumblebee", "Honeybee", "Solitary Bee"))
beeGenusColour<- setNames(
  c(beeTypeColour[c("Honeybee", "Bumblebee")],
    rgb(238, 233, 233, maxColorValue = 255),
    rgb(154, 205, 50, maxColorValue = 255)
  ),
c("Apis", "Bombus", "Lassioglossum", "Andrena"))

####################
#  DATA WRANGLING  #
####################

# Number of stigma per flower
numStigma <- 5
# Create a long-form version of the deposition data (one row per stigma)
longDepData <- do.call(rbind, lapply(X = 1:nrow(pollenDepData), FUN = function(curRowInd, inData, numStigma) {
  tempFrame <- inData[rep(curRowInd, numStigma),
    !grepl("^Stigma", colnames(inData), perl = TRUE) &
    !grepl("^Tubes", colnames(inData), perl = TRUE) &
    !grepl("^Grains", colnames(inData), perl = TRUE)]
  colSelector <- function(pistolIndex, inData, prefix, curRowInd) {
    inData[curRowInd, paste(prefix, pistolIndex, sep = " ")]
  }
  tempFrame <- cbind(tempFrame, data.frame(
    Grains = sapply(X = 1:numStigma, FUN = colSelector, inData = inData, prefix = "Stigma", curRowInd = curRowInd),
    Tubes = sapply(X = 1:numStigma, FUN = colSelector, inData = inData, prefix = "Tubes", curRowInd = curRowInd)
  ))
  rownames(tempFrame) <- paste("Flower", tempFrame[, "Flower ID"], "_", "Stigma", 1:numStigma, "_Orchard", tempFrame[, "Orchard"], sep = "")
  tempFrame
}, inData = as.data.frame(pollenDepData), numStigma = numStigma))
# Reformat the column names so that they use camel case
colnames(longDepData) <- sapply(X = strsplit(colnames(longDepData), "\\s+", perl = TRUE), FUN = function(curCol) {
  tempCol <- sapply(X = curCol, FUN = function(colDel) {
    paste(toupper(substring(colDel, 1, 1)), substring(colDel, 2), sep = "")
  })
  tempCol[1] <- tolower(tempCol[1])
  paste(tempCol, collapse = "")
})
longDepData$pollinator <- factor(
  as.character(longDepData$pollinator),
  levels = c("Bumblebee", "Honeybee", "Solitary"),
  labels = c("Bumblebee", "Honeybee", "Solitary Bee"))
longDepData$genus <- factor(
  as.character(longDepData$genus),
  levels = unique(longDepData$genus[!is.na(longDepData$genus)])
)
longDepData <- cbind(longDepData, data.frame(
  timeSinceFirstVisit = as.double(difftime(longDepData$date, sort(unique(longDepData$date))[1], units = "days"))
))

# Function to make an error bar table from model outputs
makeErrorBarTable <- function(modMean, modSE) {
  newModMean <- modMean + c(0, rep(modMean[1], length(modMean) - 1))
  data.frame(
    pred = names(modMean),
    mean = exp(newModMean),
    uppConf = exp(newModMean + 1.96 * modSE),
    lowConf = exp(newModMean - 1.96 * modSE)
  )
}

############################
#  VISIT LENGTH SUB-MODEL  #
############################

# Look at visit length sub-model
#   1. Do a model using the pollinator group
visitLengthModel_pollinator <- glm(visitLength ~ pollinator, data = longDepData, family = gaussian(link = "log"))
summary(visitLengthModel_pollinator)
#   Create a table of confidence intervals for the coefficients
levelLabel <- levels(longDepData$pollinator)
visitLengthModel_pollinator_confTab <- makeErrorBarTable(
  setNames(summary(visitLengthModel_pollinator)$coefficients[, "Estimate"], levelLabel),
  setNames(summary(visitLengthModel_pollinator)$coefficients[, "Std. Error"], levelLabel)
)
# Plot violin plot of the data and the model coefficient predictions
ggplot(data = longDepData[!is.na(longDepData$pollinator), ], aes(x = pollinator, y = visitLength)) +
  geom_violin(aes(fill = pollinator), col = NA, alpha = 0.6) +
  geom_segment(aes(x = pred, xend = pred, y = lowConf, yend = uppConf), data = visitLengthModel_pollinator_confTab) +
  geom_point(aes(x = pred, y = mean), data = visitLengthModel_pollinator_confTab) +
  scale_y_continuous(trans = "log2") + scale_fill_manual(values = beeTypeColour) +
  ylab("Visit length (seconds)") +
  xlab("") +
  theme_classic() + theme(legend.position = "none")
#   2. Do a model using the pollinator genus
visitLengthModel_genus <- glm(visitLength ~ genus, data = longDepData, family = gaussian(link = "log"))
summary(visitLengthModel_genus)
#   Create a table of confidence intervals for the coefficients
levelLabel <- levels(longDepData$genus)
visitLengthModel_genus_confTab <- makeErrorBarTable(
  setNames(summary(visitLengthModel_genus)$coefficients[, "Estimate"], levelLabel),
  setNames(summary(visitLengthModel_genus)$coefficients[, "Std. Error"], levelLabel)
)
# Plot violin plot of the data and the model coefficient predictions
ggplot(data = longDepData[!is.na(longDepData$genus), ], aes(x = genus, y = visitLength)) +
  geom_violin(aes(fill = genus), col = NA, alpha = 0.6) +
  geom_segment(aes(x = pred, xend = pred, y = lowConf, yend = uppConf), data = visitLengthModel_genus_confTab) +
  geom_point(aes(x = pred, y = mean), data = visitLengthModel_genus_confTab) +
  scale_y_continuous(trans = "log2") + scale_fill_manual(values = beeGenusColour) +
  ylab("Visit length (seconds)") +
  xlab("") +
  theme_classic() + theme(legend.position = "none")

################################
#  GRAIN DEPOSITION SUB-MODEL  #
################################

# Retrieve the average visit length
avVisitLength <- mean(longDepData$visitLength, na.rm = TRUE)
# 1. Version of the model using pollinator group
# Full model for pollen grain deposition
grainDep_pollinator <- brm(grains ~ visitLength + pollinator + (1 | date:orchard),
  data = longDepData, family = zero_inflated_negbinomial(),
  chains = numChains,
  iter = numBurnIn + numIterations,
  warmup = numBurnIn,
  thin = numThin,
  backend = "cmdstanr",
  cores = numCores)
# Produce a summary of the grain deposition parameter
summary(grainDep_pollinator)
# Plot to assess convergence of the parameters
# >1.1 too high - no convergence
# 1.05-1.1 OK
# <1.05 good
mcmc_plot(grainDep_pollinator, type = "rhat", variable = "^b_", regex = TRUE)
# Plot to look at effect of bee type
grainDep_pollinator_postSamples <- do.call(rbind, lapply(X = as_draws(grainDep_pollinator), FUN = as.data.frame))
bumbleBeeAv <- grainDep_pollinator_postSamples$b_Intercept + avVisitLength * grainDep_pollinator_postSamples$b_visitLength
grainDep_pollinator_beeTypeFrame <- data.frame(
  grainDep = exp(c(
    bumbleBeeAv,
    bumbleBeeAv + grainDep_pollinator_postSamples$b_pollinatorHoneybee,
    bumbleBeeAv + grainDep_pollinator_postSamples$b_pollinatorSolitaryBee)),
  relGrainDep = c(
    rep(NA, nrow(grainDep_pollinator_postSamples)),
    exp(bumbleBeeAv + grainDep_pollinator_postSamples$b_pollinatorHoneybee) - exp(bumbleBeeAv),
    exp(bumbleBeeAv + grainDep_pollinator_postSamples$b_pollinatorSolitaryBee) - exp(bumbleBeeAv)),
  beeType = c(
    rep("Bumblebee", nrow(grainDep_pollinator_postSamples)),
    rep("Honeybee", nrow(grainDep_pollinator_postSamples)),
    rep("Solitary Bee", nrow(grainDep_pollinator_postSamples))
  )
)
credIntFunc <- function(inValues, probs) {
  inQuants <- setNames(c(quantile(inValues, probs), mean(inValues)), c("ymin", "lower", "upper", "ymax", "middle"))
}
# Plot of grain deposition per bee type
ggplot(grainDep_pollinator_beeTypeFrame, aes(y = grainDep, x = beeType)) + geom_violin(aes(fill = beeType), alpha = 0.6, col = NA) +
  stat_summary(geom = "boxplot", fun.data = credIntFunc, fun.args = list(probs = c(0.025, 0.25, 0.75, 0.975)), width = 0.2) +
  scale_fill_manual(values = beeTypeColour) + scale_y_continuous(trans = "log2") +
  ylab("Predicted Grain Deposition at Mean Visitiation Time") + xlab("") +
  theme_classic() + theme(legend.position = "none") + coord_flip()
# Plot of relative grain deposition per bee type
ggplot(grainDep_pollinator_beeTypeFrame[grainDep_pollinator_beeTypeFrame$beeType != "Bumblebee", ], aes(y = relGrainDep, x = beeType)) + geom_violin(aes(fill = beeType), alpha = 0.6, col = NA) +
  stat_summary(geom = "boxplot", fun.data = credIntFunc, fun.args = list(probs = c(0.025, 0.25, 0.75, 0.975)), width = 0.2) +
  geom_hline(yintercept = 0.0, col = "grey", linetype = "dashed") +
  scale_fill_manual(values = beeTypeColour) +
  ylab("Predicted Grain Deposition Relative to Bumblebees at Mean Visitiation Time") + xlab("") +
  scale_y_continuous(limits = c(quantile(grainDep_pollinator_beeTypeFrame$relGrainDep, probs = 0.01, na.rm = TRUE), max(grainDep_pollinator_beeTypeFrame$relGrainDep))) +
  theme_classic() + theme(legend.position = "none") + coord_flip()
# Plot of effect of visitation time

# Plot of bee effect
coeffsToPlot <- c("Intercept", "pollinatorHoneybee", "pollinatorSolitaryBee")
labelsToPlot <- c("Bumblebee", "Honeybee", "Solitary Bee")
mcmc_plot(grainDep_pollinator, type = "areas", variable = paste("b", coeffsToPlot, sep = "_"), prob = 0.50, prob_outer = 0.95, point_est = "mean", area_method = "scaled height") + scale_y_discrete(labels = labelsToPlot)
# 2. Version of the model using genus
grainDep_genus <- brm(grains ~ visitLength + genus + (1 | data:orchard),
  data = longDepData, family = zero_inflated_negbinomial(),
  chains = numChains,
  iter = numBurnIn + numIterations,
  warmup = numBurnIn,
  thin = numThin,
  backend = "cmdstanr",
  cores = numCores)

##############################
#  TUBE FORMATION SUB-MODEL  #
##############################

# Model for the tube formation
tubeForm_pollinator <- glm(cbind(tubes, grains - tubes) ~ pollinator, data = longDepData, family = binomial)
summary(tubeForm_pollinator)
