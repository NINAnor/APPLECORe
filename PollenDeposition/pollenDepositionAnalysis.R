# Import the pollen deposition data

library(here)       # Library for relative links
library(readxl)     # Library to read in Excel files
library(tidyverse)  # Import tidyverse functions
library(brms)       # Import the brms library for diverse model specification
library(DHARMa)     # Import the model checking library
library(ggplot2)    # Import plotting functions
library(parallel)   # Import the parallel libraries
# Import the pollen deposition data from the file
pollenDepData <- read_excel("PollenDeposition/data_deposition.xlsx", "Ark1", na = c("", "NA", "missing"))
# Create a Genus column so that Genus can be analysed
pollenDepData$Genus <- sapply(X = strsplit(pollenDepData$Species, " ", fixed = TRUE), FUN = function(curVec) { curVec[1] })

# Set the MCMC 

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
  rownames(tempFrame) <- paste("Flower", tempFrame[, "Flower ID"], "_", "Stigma", 1:numStigma, sep = "")
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
  geom_violin(fill = "grey", col = NA) +
  geom_segment(aes(x = pred, xend = pred, y = lowConf, yend = uppConf), data = visitLengthModel_pollinator_confTab) +
  geom_point(aes(x = pred, y = mean), data = visitLengthModel_pollinator_confTab) +
  scale_y_continuous(trans = "log2") +
  ylab("Visit length (seconds)") +
  xlab("") +
  theme_classic()
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
  geom_violin(fill = "grey", col = NA) +
  geom_segment(aes(x = pred, xend = pred, y = lowConf, yend = uppConf), data = visitLengthModel_genus_confTab) +
  geom_point(aes(x = pred, y = mean), data = visitLengthModel_genus_confTab) +
  scale_y_continuous(trans = "log2") +
  ylab("Visit length (seconds)") +
  xlab("") +
  theme_classic()

# MCMC parameters
numChains <- 4
numIterations <- 100000
numBurnIn <- 20000
numThin <- 1
numCores <- min(detectCores(), numChains)
# options(cmdstanr_write_stan_file_dir = tempdir())
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
coeffsToPlot <- c("visitLength", "pollinatorHoneybee", "pollinatorSolitaryBee")
labelsToPlot <- c("Visit Length", "Honeybee", "Solitary Bee")
# Plot to assess convergence of the parameters
# >1.1 too high - no convergence
# 1.05-1.1 OK
# <1.05 good
mcmc_plot(grainDep_pollinator, type = "rhat", variable = "^b_", regex = TRUE)
# Plot to look at model effects
mcmc_plot(grainDep_pollinator, type = "areas", variable = paste("b", coeffsToPlot, sep = "_"), prob = 0.50, prob_outer = 0.95, point_est = "mean", area_method = "scaled height") + scale_y_discrete(labels = labelsToPlot)

# Model for the tube formation
tubeForm_pollinator <- glm(cbind(grains, tubes + grains) ~ pollinator, data = longDepData, family = binomial)
summary(tubeForm_pollinator)
