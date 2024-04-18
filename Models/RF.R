library(Rglimclim)
library(RandomForestDist)
library(dplyr)
library(tidyr)
future::plan(future::multisession, workers = parallel::detectCores())

SimData <- read.GLCdata("Simulation_from_GLM.dat")
station.data <- read.table("stations.dat",header=TRUE,
                           stringsAsFactors = FALSE)
load("DailyPredictors.rda")

#################################################################
# Temperature
#################################################################
# Extract only temperature data first
Rhine.temp <- SimData[-5]

sites <- unique(Rhine.temp$Site)
#################################################################
# Fold 1
sitesimulations1_T <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_T <- Rhine.temp[Rhine.temp$Site == site,]
  site_T <- merge(site_T, DailyPreds, by = c("Year", "Month", "Day"))
  site_T <- site_T[order(as.Date(paste(site_T$Year, site_T$Month, site_T$Day, sep = "-"))),]
  colnames(site_T) <- NULL
  sitetrain1.x_T <- as.matrix(site_T[site_T[,1] > 1974,][,c(-1:-5)])
  sitetrain1.y_T <- as.matrix(site_T[site_T[,1] > 1974,][,5])
  sitetest1.x_T <- as.matrix(site_T[site_T[,1] <= 1974,][,c(-1:-5)])
  siterf1_T <- randomForestTrain(x = sitetrain1.x_T, y = sitetrain1.y_T, minbucket=10, parallel.plan = "auto")
  sitepr.distr1_T <- randomForestPredict(siterf1_T, newdata = sitetest1.x_T, method = "mle")
  sitesimulations1_T <- data.frame(sitesimulations1_T,randomForestSimulate(sitepr.distr1_T))
}
colnames(sitesimulations1_T) <- sites
save(sitesimulations1_T,file = "sitesimulations1_T.rda")

# Fold 2
sitesimulations2_T <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_T <- Rhine.temp[Rhine.temp$Site == site,]
  site_T <- merge(site_T, DailyPreds, by = c("Year", "Month", "Day"))
  site_T <- site_T[order(as.Date(paste(site_T$Year, site_T$Month, site_T$Day, sep = "-"))),]
  colnames(site_T) <- NULL
  sitetrain2.x_T <- as.matrix(site_T[site_T[,1] < 1975 | site_T[,1] > 1990,][,c(-1:-5)])
  sitetrain2.y_T <- as.matrix(site_T[site_T[,1] < 1975 | site_T[,1] > 1990,][,5])
  sitetest2.x_T <- as.matrix(site_T[site_T[,1] >= 1975 & site_T[,1] <= 1990,][,c(-1:-5)])
  siterf2_T <- randomForestTrain(x = sitetrain2.x_T, y = sitetrain2.y_T, minbucket=10, parallel.plan = "auto")
  sitepr.distr2_T <- randomForestPredict(siterf2_T, newdata = sitetest2.x_T, method = "mle")
  sitesimulations2_T <- data.frame(sitesimulations2_T,randomForestSimulate(sitepr.distr2_T))
}
colnames(sitesimulations2_T) <- sites

# Fold 3
sitesimulations3_T <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_T <- Rhine.temp[Rhine.temp$Site == site,]
  site_T <- merge(site_T, DailyPreds, by = c("Year", "Month", "Day"))
  site_T <- site_T[order(as.Date(paste(site_T$Year, site_T$Month, site_T$Day, sep = "-"))),]
  colnames(site_T) <- NULL
  sitetrain3.x_T <- as.matrix(site_T[site_T[,1] < 1991 | site_T[,1] > 2006,][,c(-1:-5)])
  sitetrain3.y_T <- as.matrix(site_T[site_T[,1] < 1991 | site_T[,1] > 2006,][,5])
  sitetest3.x_T <- as.matrix(site_T[site_T[,1] >= 1991 & site_T[,1] <= 2006,][,c(-1:-5)])
  siterf3_T <- randomForestTrain(x = sitetrain3.x_T, y = sitetrain3.y_T, minbucket=10, parallel.plan = "auto")
  sitepr.distr3_T <- randomForestPredict(siterf3_T, newdata = sitetest3.x_T, method = "mle")
  sitesimulations3_T <- data.frame(sitesimulations3_T,randomForestSimulate(sitepr.distr3_T))
}
colnames(sitesimulations3_T) <- sites

# Fold 4
sitesimulations4_T <- data.frame(matrix(nrow = 5479, ncol = 0))
for (site in sites){
  site_T <- Rhine.temp[Rhine.temp$Site == site,]
  site_T <- merge(site_T, DailyPreds, by = c("Year", "Month", "Day"))
  site_T <- site_T[order(as.Date(paste(site_T$Year, site_T$Month, site_T$Day, sep = "-"))),]
  colnames(site_T) <- NULL
  sitetrain4.x_T <- as.matrix(site_T[site_T[,1] < 2007,][,c(-1:-5)])
  sitetrain4.y_T <- as.matrix(site_T[site_T[,1] < 2007,][,5])
  sitetest4.x_T <- as.matrix(site_T[site_T[,1] >= 2007,][,c(-1:-5)])
  siterf4_T <- randomForestTrain(x = sitetrain4.x_T, y = sitetrain4.y_T, minbucket=10, parallel.plan = "auto")
  sitepr.distr4_T <- randomForestPredict(siterf4_T, newdata = sitetest4.x_T, method = "mle")
  sitesimulations4_T <- data.frame(sitesimulations4_T,randomForestSimulate(sitepr.distr4_T))
}
colnames(sitesimulations4_T) <- sites

###################################################################
#Precipitation
###################################################################
#this is precipitation for both wet and dry days
Rhine.full <- SimData[1:5]

#extract only wet days (Var1 >= 0.95)
Rhine.wet <- subset(SimData, SimData$Var1 > 0)[1:5]

#####################################################################
# Fold 1
sitesimulations1_P <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_P <- Rhine.full[Rhine.full$Site == site,]
  site_P <- merge(site_P, DailyPreds, by = c("Year", "Month", "Day"))
  site_P <- site_P[order(as.Date(paste(site_P$Year, site_P$Month, site_P$Day, sep = "-"))),]
  site_P.wet <- subset(site_P, site_P$Var1 > 0)
  colnames(site_P) <- NULL
  colnames(site_P.wet) <- NULL
  sitetrain1.x_P <- as.matrix(site_P.wet[site_P.wet[,1] > 1974,][,c(-1:-5)])
  sitetrain1.y_P <- as.matrix(site_P.wet[site_P.wet[,1] > 1974,][,5])
  sitetest1.x_P <- as.matrix(site_P[site_P[,1] <= 1974,][,c(-1:-5)])
  siterf1_P <- randomForestTrain(x = sitetrain1.x_P, y = sitetrain1.y_P - 0.949, method = "gammaDeviation", minbucket=10, parallel.plan = "auto")
  sitepr.distr1_P <- randomForestPredict(siterf1_P, newdata = sitetest1.x_P, method = "bc3")
  sitesimulations1_P <- data.frame(sitesimulations1_P, (randomForestSimulate(sitepr.distr1_P) + 0.949))
}
colnames(sitesimulations1_P) <- sites


# Fold 2
sitesimulations2_P <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_P <- Rhine.full[Rhine.full$Site == site,]
  site_P <- merge(site_P, DailyPreds, by = c("Year", "Month", "Day"))
  site_P <- site_P[order(as.Date(paste(site_P$Year, site_P$Month, site_P$Day, sep = "-"))),]
  site_P.wet <- subset(site_P, site_P$Var1 > 0)
  colnames(site_P) <- NULL
  colnames(site_P.wet) <- NULL
  sitetrain2.x_P <- as.matrix(site_P.wet[site_P.wet[,1] < 1975 | site_P.wet[,1] > 1990,][,c(-1:-5)])
  sitetrain2.y_P <- as.matrix(site_P.wet[site_P.wet[,1] < 1975 | site_P.wet[,1] > 1990,][,5])
  sitetest2.x_P <- as.matrix(site_P[site_P[,1] >= 1975 & site_P[,1] <= 1990,][,c(-1:-5)])
  siterf2_P <- randomForestTrain(x = sitetrain2.x_P, y = sitetrain2.y_P - 0.949, method = "gammaDeviation", minbucket=10, parallel.plan = "auto")
  sitepr.distr2_P <- randomForestPredict(siterf2_P, newdata = sitetest2.x_P, method = "bc3")
  sitesimulations2_P <- data.frame(sitesimulations2_P, (randomForestSimulate(sitepr.distr2_P) + 0.949))
}
colnames(sitesimulations2_P) <- sites


# Fold 3
sitesimulations3_P <- data.frame(matrix(nrow = 5844, ncol = 0))
for (site in sites){
  site_P <- Rhine.full[Rhine.full$Site == site,]
  site_P <- merge(site_P, DailyPreds, by = c("Year", "Month", "Day"))
  site_P <- site_P[order(as.Date(paste(site_P$Year, site_P$Month, site_P$Day, sep = "-"))),]
  site_P.wet <- subset(site_P, site_P$Var1 > 0)
  colnames(site_P) <- NULL
  colnames(site_P.wet) <- NULL
  sitetrain3.x_P <- as.matrix(site_P.wet[site_P.wet[,1] < 1991 | site_P.wet[,1] > 2006,][,c(-1:-5)])
  sitetrain3.y_P <- as.matrix(site_P.wet[site_P.wet[,1] < 1991 | site_P.wet[,1] > 2006,][,5])
  sitetest3.x_P <- as.matrix(site_P[site_P[,1] >= 1991 & site_P[,1] <= 2006,][,c(-1:-5)])
  siterf3_P <- randomForestTrain(x = sitetrain3.x_P, y = sitetrain3.y_P - 0.949, method = "gammaDeviation", minbucket=10, parallel.plan = "auto")
  sitepr.distr3_P <- randomForestPredict(siterf3_P, newdata = sitetest3.x_P, method = "bc3")
  sitesimulations3_P <- data.frame(sitesimulations3_P, (randomForestSimulate(sitepr.distr3_P) + 0.949))
}
colnames(sitesimulations3_P) <- sites


# Fold 4
sitesimulations4_P <- data.frame(matrix(nrow = 5479, ncol = 0))
for (site in sites){
  site_P <- Rhine.full[Rhine.full$Site == site,]
  site_P <- merge(site_P, DailyPreds, by = c("Year", "Month", "Day"))
  site_P <- site_P[order(as.Date(paste(site_P$Year, site_P$Month, site_P$Day, sep = "-"))),]
  site_P.wet <- subset(site_P, site_P$Var1 > 0)
  colnames(site_P) <- NULL
  colnames(site_P.wet) <- NULL
  sitetrain4.x_P <- as.matrix(site_P.wet[site_P.wet[,1] < 2007,][,c(-1:-5)])
  sitetrain4.y_P <- as.matrix(site_P.wet[site_P.wet[,1] < 2007,][,5])
  sitetest4.x_P <- as.matrix(site_P[site_P[,1] >= 2007,][,c(-1:-5)])
  siterf4_P <- randomForestTrain(x = sitetrain4.x_P, y = sitetrain4.y_P - 0.949, method = "gammaDeviation", minbucket=10, parallel.plan = "auto")
  sitepr.distr4_P <- randomForestPredict(siterf4_P, newdata = sitetest4.x_P, method = "bc3")
  sitesimulations4_P <- data.frame(sitesimulations4_P, (randomForestSimulate(sitepr.distr4_P) + 0.949))
}
colnames(sitesimulations4_P) <- sites