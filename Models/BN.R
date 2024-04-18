library(Rglimclim)
library(loadeR)
library(loadeR.java)
library(climate4R.UDG)
library(dplyr)
library(BNWeatherGen,  quietly = T)
library(devtools)

SimData <- read.GLCdata("Simulation_from_GLM.dat")
station.data <- read.table("stations.dat",header=TRUE,
                           stringsAsFactors = FALSE)
sites <- unique(SimData$Site)
#Modify precipitation to occurence first
Rhine.threshold <- SimData
Rhine.threshold$Var1 <- ifelse(test = Rhine.threshold$Var1 > 0,yes = 1, no = 0)

#VALUE (text) format for station data
names(station.data)[names(station.data) == 'Scode'] <- 'station_id'
names(station.data)[names(station.data) == 'Longitude'] <- 'longitude'
names(station.data)[names(station.data) == 'Latitude'] <- 'latitude'

#VALUE (text) format for variables
variables.VALUE <- data.frame(variable_id = "precip",
                              name = "Total_precipitation_accumulated_in_24 hours",
                              unit = "mm",
                              missing_code = "NaN")

##########################################################
#fold1
fold1_Rhine.threshold <- Rhine.threshold[Rhine.threshold[,1] <= 1974,]

#VALUE (text) format for station data
#stations.txt
write.table(station.data, file = "stations.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#variables.txt
write.table(variables.VALUE, file = "variables.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#precip.txt
#VALUE (text) format for precip (OccOnly)
fold1_Rhine.OccOnly <- fold1_Rhine.threshold
fold1_Rhine.OccOnly$YYYYMMDD <- (1e4 * fold1_Rhine.OccOnly$Year) + (100 * fold1_Rhine.OccOnly$Month) + fold1_Rhine.OccOnly$Day
fold1_Rhine.OccOnly <- fold1_Rhine.OccOnly[,c(7,4,5)]
fold1_Rhine.OccOnly <- reshape(fold1_Rhine.OccOnly, direction = "wide", idvar = "YYYYMMDD", timevar = "Site")
colnames(fold1_Rhine.OccOnly) <- c(gsub("Var1.", "", colnames(fold1_Rhine.OccOnly))) #delete "Var1." in column that generates from reshape()
colnames(fold1_Rhine.OccOnly) #sites of stations.txt and variables.txt should have the same order
sorted_colnames <- sort(colnames(fold1_Rhine.OccOnly))
sorted_colnames <- c(colnames(fold1_Rhine.OccOnly)[1], sort(colnames(fold1_Rhine.OccOnly)[2:88]))
fold1_Rhine.OccOnly <- fold1_Rhine.OccOnly %>% select(sorted_colnames)

write.table(fold1_Rhine.OccOnly, file = "precip.txt", sep = ", ", row.names = FALSE, quote = FALSE, na = "NaN")

# Check
#a <- read.csv("precip.txt")
#b <- read.csv("stations.txt")

# Then, put three files into one folder
fold1_di <-dataInventory('fold1', return.stats= TRUE)
fold1_input <- loadStationData(dataset = 'fold1', var = 'precip')
fold1_input[[1]] <- fold1_input[[1]][-2]

fold1_BNmodel <- buildWeatherGenerator(y = fold1_input)

##################################################################
#fold2
fold2_Rhine.threshold <- Rhine.threshold[Rhine.threshold$Year >= 1975 & Rhine.threshold$Year <= 1990,]

#VALUE (text) format for station data
#stations.txt
write.table(station.data, file = "stations.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#variables.txt
write.table(variables.VALUE, file = "variables.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#precip.txt
#VALUE (text) format for precip (OccOnly)
fold2_Rhine.OccOnly <- fold2_Rhine.threshold
fold2_Rhine.OccOnly$YYYYMMDD <- (1e4 * fold2_Rhine.OccOnly$Year) + (100 * fold2_Rhine.OccOnly$Month) + fold2_Rhine.OccOnly$Day
fold2_Rhine.OccOnly <- fold2_Rhine.OccOnly[,c(7,4,5)]
fold2_Rhine.OccOnly <- reshape(fold2_Rhine.OccOnly, direction = "wide", idvar = "YYYYMMDD", timevar = "Site")
colnames(fold2_Rhine.OccOnly) <- c(gsub("Var1.", "", colnames(fold2_Rhine.OccOnly))) #delete "Var1." in column that generates from reshape()
colnames(fold2_Rhine.OccOnly) #sites of stations.txt and variables.txt should have the same order
sorted_colnames <- sort(colnames(fold2_Rhine.OccOnly))
sorted_colnames <- c(colnames(fold2_Rhine.OccOnly)[1], sort(colnames(fold2_Rhine.OccOnly)[2:88]))
fold2_Rhine.OccOnly <- fold2_Rhine.OccOnly %>% select(sorted_colnames)

write.table(fold2_Rhine.OccOnly, file = "precip.txt", sep = ", ", row.names = FALSE, quote = FALSE, na = "NaN")

# Then, put three files into one folder
fold2_di <-dataInventory('fold2', return.stats= TRUE)
fold2_input <- loadStationData(dataset = 'fold2', var = 'precip')
fold2_input[[1]] <- fold2_input[[1]][-2]

fold2_BNmodel <- buildWeatherGenerator(y = fold2_input)

############################################################
#fold3
fold3_Rhine.threshold <- Rhine.threshold[Rhine.threshold$Year >= 1991 & Rhine.threshold$Year <= 2006,]

#VALUE (text) format for station data
#stations.txt
write.table(station.data, file = "stations.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#variables.txt
write.table(variables.VALUE, file = "variables.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#precip.txt
#VALUE (text) format for precip (OccOnly)
fold3_Rhine.OccOnly <- fold3_Rhine.threshold
fold3_Rhine.OccOnly$YYYYMMDD <- (1e4 * fold3_Rhine.OccOnly$Year) + (100 * fold3_Rhine.OccOnly$Month) + fold3_Rhine.OccOnly$Day
fold3_Rhine.OccOnly <- fold3_Rhine.OccOnly[,c(7,4,5)]
fold3_Rhine.OccOnly <- reshape(fold3_Rhine.OccOnly, direction = "wide", idvar = "YYYYMMDD", timevar = "Site")
colnames(fold3_Rhine.OccOnly) <- c(gsub("Var1.", "", colnames(fold3_Rhine.OccOnly))) #delete "Var1." in column that generates from reshape()
colnames(fold3_Rhine.OccOnly) #sites of stations.txt and variables.txt should have the same order
sorted_colnames <- sort(colnames(fold3_Rhine.OccOnly))
sorted_colnames <- c(colnames(fold3_Rhine.OccOnly)[1], sort(colnames(fold3_Rhine.OccOnly)[2:88]))
fold3_Rhine.OccOnly <- fold3_Rhine.OccOnly %>% select(sorted_colnames)

write.table(fold3_Rhine.OccOnly, file = "precip.txt", sep = ", ", row.names = FALSE, quote = FALSE, na = "NaN")

# Then, put three files into one folder
fold3_di <-dataInventory('fold3', return.stats= TRUE)
fold3_input <- loadStationData(dataset = 'fold3', var = 'precip')
fold3_input[[1]] <- fold3_input[[1]][-2]

fold3_BNmodel <- buildWeatherGenerator(y = fold3_input)

###################################################################
#fold4
fold4_Rhine.threshold <- Rhine.threshold[Rhine.threshold[,1] >= 2007,]

#VALUE (text) format for station data
#stations.txt
write.table(station.data, file = "stations.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#variables.txt
write.table(variables.VALUE, file = "variables.txt", sep = ", ", row.names = FALSE, quote = FALSE)

#precip.txt
#VALUE (text) format for precip (OccOnly)
fold4_Rhine.OccOnly <- fold4_Rhine.threshold
fold4_Rhine.OccOnly$YYYYMMDD <- (1e4 * fold4_Rhine.OccOnly$Year) + (100 * fold4_Rhine.OccOnly$Month) + fold4_Rhine.OccOnly$Day
fold4_Rhine.OccOnly <- fold4_Rhine.OccOnly[,c(7,4,5)]
fold4_Rhine.OccOnly <- reshape(fold4_Rhine.OccOnly, direction = "wide", idvar = "YYYYMMDD", timevar = "Site")
colnames(fold4_Rhine.OccOnly) <- c(gsub("Var1.", "", colnames(fold4_Rhine.OccOnly))) #delete "Var1." in column that generates from reshape()
colnames(fold4_Rhine.OccOnly) #sites of stations.txt and variables.txt should have the same order
sorted_colnames <- sort(colnames(fold4_Rhine.OccOnly))
sorted_colnames <- c(colnames(fold4_Rhine.OccOnly)[1], sort(colnames(fold4_Rhine.OccOnly)[2:88]))
fold4_Rhine.OccOnly <- fold4_Rhine.OccOnly %>% select(sorted_colnames)

write.table(fold4_Rhine.OccOnly, file = "precip.txt", sep = ", ", row.names = FALSE, quote = FALSE, na = "NaN")

# Then, put three files into one folder
fold4_di <-dataInventory('fold4', return.stats= TRUE)
fold4_input <- loadStationData(dataset = 'fold4', var = 'precip')
fold4_input[[1]] <- fold4_input[[1]][-2]

fold4_BNmodel <- buildWeatherGenerator(y = fold4_input)

#################################################################
# Simulate the results
fold1_weather.series <- generateWeather(wg = fold1_BNmodel, n = 5844)
fold1_weather.series <- fold1_weather.series[-1,]
colnames(fold1_weather.series) <- sites

fold2_weather.series <- generateWeather(wg = fold2_BNmodel, n = 5844)
fold2_weather.series <- fold2_weather.series[-1,]
colnames(fold2_weather.series) <- sites

fold3_weather.series <- generateWeather(wg = fold3_BNmodel, n = 5844)
fold3_weather.series <- fold3_weather.series[-1,]
colnames(fold3_weather.series) <- sites

fold4_weather.series <- generateWeather(wg = fold4_BNmodel, n = 5479)
fold4_weather.series <- fold4_weather.series[-1,]
colnames(fold4_weather.series) <- sites