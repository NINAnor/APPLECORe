##############################
# BEFORE YOU RUN THIS SCRIPT #
##############################
#   You must set an environment variable to where the APPLECORe Analysis file is
#   located:
#     1. Run usethis::edit_r_environ() - this will load up an environmental variables file
#       that will set the environmental variables when R starts up
#     2. Add a line to the file with the text: APPLECORE_SHAREPOINTANALYSIS_LOC=[file location here] where
#       you replace the text [file location here] to the location of directory that contains
#       the "Analysis" folder
#     3. Save the file after modifying it
#     4. Restart R
#   You must also install CmdStanR which involves doing the following:
#     1. If you are on a Windows machine then you must install RTools.  Easiest way to do
#       this is to install the 'installr' package (run the command 'install.packages("installr"))
#       and then run the command 'library(installr)' and then 'install.Rtools()'
#     2. Run the command 'install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))'
#     3. Run the command 'library(cmdstanr)' then 'check_cmdstan_toolchain(fix = TRUE)'
#       finally followed by 'install_cmdstan()'
######
# Import the neccessary libraries to run the analysis
library(readxl)     # Library to read in Excel files
library(tidyverse)  # Import tidyverse functions
library(cmdstanr)   # Import the STAN libraries
library(brms)       # Import the brms library for diverse model specification
# library(DHARMa)     # Import the model checking library
library(ggplot2)    # Import plotting functions
library(parallel)   # Import the parallel libraries
# Import the pollen deposition data from the file
pollenDepData <- read_excel(file.path(Sys.getenv("APPLECORE_SHAREPOINTANALYSIS_LOC"), "PollenDeposition", "data_deposition.xlsx"), "Ark1", na = c("", "NA", "missing"))
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

# Set the density of prediction when plotting marginal effects of continuous variables
predDens <- 100

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
makeErrorBarTable <- function(modMean, modSE, link = "log") {
  newModMean <- modMean + c(0, rep(modMean[1], length(modMean) - 1))
  outFrame <- data.frame(
    pred = names(modMean),
    mean = exp(newModMean),
    uppConf = exp(newModMean + 1.96 * modSE),
    lowConf = exp(newModMean - 1.96 * modSE)
  )
  if(link != "log") {
    outFrame$mean <- outFrame$mean / (1.0 + outFrame$mean)
    outFrame$uppConf <- outFrame$uppConf / (1.0 + outFrame$uppConf)
    outFrame$lowConf <- outFrame$lowConf / (1.0 + outFrame$lowConf)
  }
  outFrame
}
# Retrieve a table of dates where depositions happened at each of the sites
orchardDateFrame <- do.call(rbind, lapply(X = unique(longDepData$orchard), FUN = function(curOrchard, longDepData) {
  # Retrieve the data taken from the current orchard
  tempData <- longDepData[longDepData$orchard == curOrchard, ]
  # Get a set of the unique dates with samples taken at that orchard
  uniqueDates <- unique(tempData$date[!is.na(tempData$date)])
  outFrame <- data.frame()
  if(length(uniqueDates) > 0) {
    # Count the number of samples at each of the unique dates
    samplesAtDate <- sapply(X = uniqueDates, FUN = function(curDate, tempData) {
      sum(tempData$date == curDate, na.rm = TRUE)
    }, tempData = tempData)
    outFrame <- data.frame(
      orchard = rep(curOrchard, length(uniqueDates)),
      date = strptime(uniqueDates, format = "%Y-%m-%d"),
      count = samplesAtDate
    )
  }
  outFrame
}, longDepData = longDepData))
rownames(orchardDateFrame) <- paste(orchardDateFrame$date, orchardDateFrame$orchard, sep = "_")

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
#  backend = "cmdstanr",
  cores = numCores)
# Produce a summary of the grain deposition model
summary(grainDep_pollinator)
# Plot to assess convergence of the parameters
# >1.1 too high - no convergence
# 1.05-1.1 OK
# <1.05 good
mcmc_plot(grainDep_pollinator, type = "rhat", variable = "^b_", regex = TRUE)
# Create a new data frame to plot the effect of bee type on pollen deposition
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
# Make a data frame to assess the effect of visitation length on pollen deposition
grainDep_pollinator_visitLenFrame <- data.frame(
  visitLength = rep(seq(min(longDepData$visitLength, na.rm = TRUE), max(longDepData$visitLength, na.rm = TRUE), length.out = predDens), 3),
  pollinator = rep(c("Bumblebee", "Honeybee", "Solitary Bee"), rep(predDens, 3))
)
grainDep_pollinator_visitLenFrame <- cbind(grainDep_pollinator_visitLenFrame, predict(grainDep_pollinator, newdata = grainDep_pollinator_visitLenFrame, re_formula = ~ visitLength + pollinator, summary = TRUE, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975)))
# Plot of the marginal effect of visit length
ggplot(grainDep_pollinator_visitLenFrame, aes(x = visitLength, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = gray(0.97)) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = gray(0.95)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = gray(0.90)) +
  geom_line(linewidth = 0.7, color = gray(0.7)) +
  geom_point(aes(x = visitLength, y = grains, colour = pollinator), data = longDepData[!is.na(longDepData$pollinator), ]) +
  facet_wrap(~ pollinator) +
  scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2") +
  xlab("Visit Length (s)") + ylab("Pollen Grain Deposition") +
  theme_classic() + theme(legend.position = "none") + scale_colour_manual(values = beeTypeColour)
waic(grainDep_pollinator)
# Create a data frame of the MCMC samples
grainDep_pollinator_randEffFrame <- posterior_summary(grainDep_pollinator, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
grainDep_pollinator_randEffFrame <- grainDep_pollinator_randEffFrame[grepl("^r_date\\:orchard\\[", row.names(grainDep_pollinator_randEffFrame), perl = TRUE), ]
row.names(grainDep_pollinator_randEffFrame) <- gsub(",Intercept\\]$", "",
  gsub("^r_date\\:orchard\\[", "", row.names(grainDep_pollinator_randEffFrame), perl = TRUE), perl = TRUE)
grainDep_pollinator_randEffFrame <- cbind(orchardDateFrame, grainDep_pollinator_randEffFrame)
# Plot of the random effects
ggplot(grainDep_pollinator_randEffFrame, aes(x = date, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = gray(0.97)) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = gray(0.95)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = gray(0.90)) +
  geom_line(linewidth = 0.7, color = gray(0.7)) +
  geom_point() +
  geom_hline(yintercept = 0.0, linetype = "dotted") +
  facet_wrap(~ orchard) +
  xlab("Date") + ylab("Random Effect") +
  theme_classic() + theme(legend.position = "none") + scale_colour_manual(values = beeTypeColour)

# 2. Version of the model using genus
# Full model for pollen grain deposition
longDepData$genus <- factor(as.character(longDepData$genus),levels = c("Apis", "Andrena", "Bombus", "Lasioglossum"))
grainDep_genus <- brm(grains ~ visitLength + genus + (1 | date:orchard),
  data = longDepData, family = zero_inflated_negbinomial(),
  chains = numChains,
  iter = numBurnIn + numIterations,
  warmup = numBurnIn,
  thin = numThin,
  #backend = "cmdstanr",
  cores = numCores)
# Produce a summary of the grain deposition parameter
summary(grainDep_genus)
# Plot to assess convergence of the parameters
# >1.1 too high - no convergence
# 1.05-1.1 OK
# <1.05 good
mcmc_plot(grainDep_genus, type = "rhat", variable = "^b_", regex = TRUE)
# Create a new data frame to plot the effect of bee type on pollen deposition
grainDep_genus_postSamples <- do.call(rbind, lapply(X = as_draws(grainDep_genus), FUN = as.data.frame))
apisAv <- grainDep_genus_postSamples$b_Intercept + avVisitLength * grainDep_pollinator_postSamples$b_visitLength
grainDep_genus_beeTypeFrame <- data.frame(
  grainDep = exp(c(
    apisAv,
    apisAv + grainDep_genus_postSamples$b_genusAndrena,
    apisAv + grainDep_genus_postSamples$b_genusBombus,
    apisAv + grainDep_genus_postSamples$b_genusLasioglossum)),
  relGrainDep = c(
    rep(NA, nrow(grainDep_genus_postSamples)),
    exp(apisAv + grainDep_genus_postSamples$b_genusAndrena) - exp(apisAv),
    exp(apisAv + grainDep_genus_postSamples$b_genusBombus) - exp(apisAv),
    exp(apisAv + grainDep_genus_postSamples$b_genusLasioglossum) - exp(apisAv)),
  beeType = c(
    rep("Apis", nrow(grainDep_genus_postSamples)),
    rep("Andrena", nrow(grainDep_genus_postSamples)),
    rep("Bombus", nrow(grainDep_genus_postSamples)),
    rep("Lasioglossum", nrow(grainDep_genus_postSamples)))
)
# Make a data frame of difference summaries
diffSummary <- do.call(rbind, lapply(X = 2:length(levels(longDepData$genus)), FUN = function(fromIndex, outSamples) {
  # Calculate the base deposition
  baseDep <- exp(outSamples[outSamples$beeType == unique(outSamples$beeType)[fromIndex], "grainDep"])
  do.call(rbind, lapply(X = 1:(fromIndex - 1), FUN = function(toIndex, fromIndex, baseDep, outSamples) {
    toDep <- exp(outSamples[outSamples$beeType == unique(outSamples$beeType)[toIndex], "grainDep"])
    diffSamples <- toDep - baseDep
    probCalc <- sum(diffSamples > 0) / length(diffSamples)
    data.frame(
      genusOne = unique(outSamples$beeType)[fromIndex],
      genusTwo = unique(outSamples$beeType)[toIndex],
      probTwoGreaterThanOne = probCalc,
      probOneGreaterThanTwo = 1.0 - probCalc
    )
  }, fromIndex = fromIndex, baseDep = baseDep, outSamples = outSamples))
}, outSamples = grainDep_genus_beeTypeFrame))
# Plot of grain deposition per bee type
ggplot(grainDep_genus_beeTypeFrame, aes(y = grainDep, x = beeType)) + geom_violin(aes(fill = beeType), alpha = 0.6, col = NA) +
  stat_summary(geom = "boxplot", fun.data = credIntFunc, fun.args = list(probs = c(0.025, 0.25, 0.75, 0.975)), width = 0.2) +
  scale_fill_manual(values = beeGenusColour) + scale_y_continuous(trans = "log2") +
  ylab("Predicted Grain Deposition at Mean Visitiation Time") + xlab("") +
  theme_classic() + theme(legend.position = "none") + coord_flip()
# Plot of relative grain deposition per bee type
ggplot(grainDep_genus_beeTypeFrame[grainDep_genus_beeTypeFrame$beeType != "Apis", ], aes(y = relGrainDep, x = beeType)) + geom_violin(aes(fill = beeType), alpha = 0.6, col = NA) +
  stat_summary(geom = "boxplot", fun.data = credIntFunc, fun.args = list(probs = c(0.025, 0.25, 0.75, 0.975)), width = 0.2) +
  geom_hline(yintercept = 0.0, col = "grey", linetype = "dashed") +
  scale_fill_manual(values = beeGenusColour) +
  ylab("Predicted Grain Deposition Relative to Apis at Mean Visitiation Time") + xlab("") +
  scale_y_continuous(limits = c(quantile(grainDep_pollinator_beeTypeFrame$relGrainDep, probs = 0.01, na.rm = TRUE), max(grainDep_pollinator_beeTypeFrame$relGrainDep))) +
  theme_classic() + theme(legend.position = "none") + coord_flip()
# Make a data frame to assess the effect of visitation length on pollen deposition
grainDep_genus_visitLenFrame <- data.frame(
  visitLength = rep(seq(min(longDepData$visitLength, na.rm = TRUE), max(longDepData$visitLength, na.rm = TRUE), length.out = predDens), 4),
  genus = rep(c("Apis", "Andrena", "Bombus", "Lasioglossum"), rep(predDens, 4))
)
grainDep_genus_visitLenFrame <- cbind(grainDep_genus_visitLenFrame, predict(grainDep_genus,newdata = grainDep_genus_visitLenFrame, re_formula = ~ visitLength + genus, summary = TRUE, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975)))
# Plot of the marginal effect of visit length
ggplot(grainDep_genus_visitLenFrame, aes(x = visitLength, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = gray(0.97)) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = gray(0.95)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = gray(0.90)) +
  geom_line(linewidth = 0.7, color = gray(0.7)) +
  geom_point(aes(x = visitLength, y = grains, colour = genus), data = longDepData[!is.na(longDepData$genus), ]) +
  facet_wrap(~ genus) +
  scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2") +
  xlab("Visit Length (s)") + ylab("Pollen Grain Deposition") +
  theme_classic() + theme(legend.position = "none") + scale_colour_manual(values = beeGenusColour)
waic(grainDep_genus)
# Create a data frame of the MCMC samples
grainDep_genus_randEffFrame <- posterior_summary(grainDep_genus, probs = c(0.025, 0.05, 0.25, 0.75, 0.95, 0.975))
grainDep_genus_randEffFrame <- grainDep_genus_randEffFrame[grepl("^r_date\\:orchard\\[", row.names(grainDep_genus_randEffFrame), perl = TRUE), ]
row.names(grainDep_genus_randEffFrame) <- gsub(",Intercept\\]$", "",
  gsub("^r_date\\:orchard\\[", "", row.names(grainDep_genus_randEffFrame), perl = TRUE), perl = TRUE)
grainDep_genus_randEffFrame <- cbind(orchardDateFrame, grainDep_genus_randEffFrame)
# Plot of the random effects
ggplot(grainDep_genus_randEffFrame, aes(x = date, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = gray(0.97)) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = gray(0.95)) +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), fill = gray(0.90)) +
  geom_line(linewidth = 0.7, color = gray(0.7)) +
  geom_point() +
  geom_hline(yintercept = 0.0, linetype = "dotted") +
  facet_wrap(~ orchard) +
  xlab("Date") + ylab("Random Effect") +
  theme_classic() + theme(legend.position = "none") + scale_colour_manual(values = beeGenusColour)

##############################
#  TUBE FORMATION SUB-MODEL  #
##############################

# Model for the tube formation
tubeForm_genus <- glm(cbind(tubes, pmax(grains - tubes, 0)) ~ genus, data = longDepData, family = binomial)
summary(tubeForm_genus)
# Create a table of confidence intervals for the coefficients
levelLabel <- levels(longDepData$genus)
tubeFormModel_genus_confTab <- makeErrorBarTable(
  setNames(summary(tubeForm_genus)$coefficients[, "Estimate"], levelLabel),
  setNames(summary(tubeForm_genus)$coefficients[, "Std. Error"], levelLabel),
  link = "logit")
# Plot violin plot of the data and the model coefficient predictions
longDepData$propTubes <- pmin(longDepData$tubes / longDepData$grains, 1)
ggplot(data = longDepData[!is.na(longDepData$genus), ], aes(x = genus, y = propTubes)) +
  geom_violin(aes(fill = genus), col = NA, alpha = 0.6) +
  geom_segment(aes(x = pred, xend = pred, y = lowConf, yend = uppConf), data = tubeFormModel_genus_confTab) +
  geom_point(aes(x = pred, y = mean), data = tubeFormModel_genus_confTab) +
  scale_fill_manual(values = beeGenusColour) +
  ylab("Proportion of Tubes Formed From Grains") +
  xlab("") +
  theme_classic() + theme(legend.position = "none")
