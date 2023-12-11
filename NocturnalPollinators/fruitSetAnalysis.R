library(openxlsx2)
library(nimble)

# Replace the following line with the location to where you have stored the "FruitSet.xlsx" data file
fruitSetDataLoc <- file.path(Sys.getenv("APPLECORE_SHAREPOINTDATA_LOC"), "kristiMasters", "FruitSet.xlsx")

# Import the fruit set data and tidy it up
fruitSetData <- read_xlsx(file = fruitSetDataLoc, sheet = "FruitSet", start_row = 3)
fruitSetData <- fruitSetData[, !is.na(names(fruitSetData))]
# Import the flowers per cluster data and tidy it up
flowersData <- read_xlsx(file = fruitSetDataLoc, sheet = "FlowerFormation")

# Combine the datasets into one combined data frame
combinedData <- data.frame(
  id = c(fruitSetData$ID, paste("F", substring(flowersData$Site, 1, 1), gsub("^S$", "SR", substring(flowersData$Cultivar, 1, 1), perl = TRUE), 1:nrow(flowersData), sep = "-")),
  experiment = factor(rep(c("fruitSet", "flowersPerCluster"), c(nrow(fruitSetData), nrow(flowersData))), c("flowersPerCluster", "fruitSet")),
  site = tolower(as.character(c(fruitSetData$Site, flowersData$Site))),
  cultivar = tolower(as.character(c(fruitSetData$Cultivar, flowersData$Cultivar))),
  treatment = factor(c(tolower(as.character(fruitSetData$Treatment)), rep("open", nrow(flowersData))), c("open", "closed", "day", "night")),
  napples = c(fruitSetData$n_apples, rep(NA, nrow(flowersData))),
  nflowers = c(rep(NA, nrow(fruitSetData)), flowersData$n_flowers),
  nclusters = c(fruitSetData$n_clusters, rep(1, nrow(flowersData))),
  nexpcount = c(fruitSetData$n_apples, flowersData$n_flowers),
  stringsAsFactors = TRUE
)
# Add extra treatment columns
#combinedData <- cbind(combinedData, data.frame(
#  treatment.closed = ifelse(as.character(combinedData$treatment) == "closed", 1, 0),
#  treatment.day = ifelse(as.character(combinedData$treatment) == "day", 1, 0),
#  treatment.night = ifelse(as.character(combinedData$treatment) == "night", 1, 0)
#))

modelOut <- glm(
  nexpcount ~ offset(log(nclusters)) + site + experiment + cultivar * treatment,
  data = combinedData[combinedData$nclusters > 0, ],
  family = poisson
)
