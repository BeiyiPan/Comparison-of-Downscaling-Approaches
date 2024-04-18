library(Rglimclim)
library(maps)
library(lattice)

# 4-fold cross validation
#PrecipTemp
PrecipTemp <- read.GLCdata("PrecipTemp_Processed.dat")

fold1.PrecipTemp <- subset(PrecipTemp, PrecipTemp$Year <= 1974)
fold2.PrecipTemp <- subset(PrecipTemp, PrecipTemp$Year >= 1975 & PrecipTemp$Year <= 1990)
fold3.PrecipTemp <- subset(PrecipTemp, PrecipTemp$Year >= 1991 & PrecipTemp$Year <= 2006)
fold4.PrecipTemp <- subset(PrecipTemp, PrecipTemp$Year >= 2007)

write.GLCdata(fold1.PrecipTemp, file = "fold1.PrecipTemp.dat")
write.GLCdata(fold2.PrecipTemp, file = "fold2.PrecipTemp.dat")
write.GLCdata(fold3.PrecipTemp, file = "fold3.PrecipTemp.dat")
write.GLCdata(fold4.PrecipTemp, file = "fold4.PrecipTemp.dat")

train1.PrecipTemp <- rbind(fold2.PrecipTemp, fold3.PrecipTemp, fold4.PrecipTemp)
train2.PrecipTemp <- rbind(fold1.PrecipTemp, fold3.PrecipTemp, fold4.PrecipTemp)
train3.PrecipTemp <- rbind(fold1.PrecipTemp, fold2.PrecipTemp, fold4.PrecipTemp)
train4.PrecipTemp <- rbind(fold1.PrecipTemp, fold2.PrecipTemp, fold3.PrecipTemp)

write.GLCdata(train1.PrecipTemp, file = "train1.PrecipTemp.dat")
write.GLCdata(train2.PrecipTemp, file = "train2.PrecipTemp.dat")
write.GLCdata(train3.PrecipTemp, file = "train3.PrecipTemp.dat")
write.GLCdata(train4.PrecipTemp, file = "train4.PrecipTemp.dat")

###

fold4.DailyPreds <- subset(DailyPreds, DailyPreds$Year >= 2007)
write.GLCexternal(fold4.DailyPreds, file = "fold4.DailyPreds.dat")
###
Rver <- paste(R.version$major, R.version$minor, sep=".")
if (Rver >= "3.6.1") {
  warn.settings <- options()$warn
  options(warn=-1)
}

#####
set.seed(2000)
#####

##########################################################################
#train1
#TempModel
train1.TempModel.init <-
  read.modeldef("DEFINITIONS/TempModel12Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train1.TempModelDisp.init <-
  read.modeldef("DEFINITIONS/TempModel12Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train1.TempModel <- GLCfit(train1.TempModel.init, dispersion.def=train1.TempModelDisp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="train1.PrecipTemp.dat",
                      external.files=c(NA,NA,"DailyPreds.dat"),
                      diagnostics=2,nprev.required=-1,
                      cor.file="train1.TempModel_Corr.dat",
                      resid.file="train1.TempModel_Resids.dat")

print(train1.TempModel)

summary(train1.TempModel, tables=NULL)

save(siteinfo,train1.TempModel,file="train1.Temperature.rda")


#OccModel
train1.OccModel.init <-
  read.modeldef("DEFINITIONS/OccModel13Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train1.OccModel <- GLCfit(train1.OccModel.init,which.response=1,siteinfo=siteinfo,
                     data.file="train1.PrecipTemp.dat",
                     diagnostics=1, nprev.required=-1, verbosity=1,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "train1.OccModel_Corr.dat")

summary(train1.OccModel, tables=NULL)

save(siteinfo,train1.OccModel,file="train1.Occurrence.rda")


#IntModel
train1.IntModel.init <-
  read.modeldef("DEFINITIONS/IntModel12Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train1.IntModel <- GLCfit(train1.IntModel.init,which.response=1,siteinfo=siteinfo,
                     data.file="train1.PrecipTemp.dat",
                     diagnostics=2, nprev.required=-1, verbosity=2,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "train1.IntModel_Corr.dat",
                     resid.file = "train1.IntModel_Resids.dat")

summary(train1.IntModel, tables=NULL)

save(siteinfo,train1.IntModel,file="train1.Intensity.rda")



sim.fold1 <-
  GLCsim(list(list(Occurrence=train1.OccModel,Intensity=train1.IntModel),train1.TempModel),
         nsims=100, start=195901, end=197412, impute.until=195812,
         data.file="fold1.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold1", file.prefix="Simfold1",
         external.files = c(NA, NA, "fold1.DailyPreds.dat"))

obs.fold1 <-
  GLCsim(list(list(Occurrence=train1.OccModel,Intensity=train1.IntModel),train1.TempModel),
         nsims=39,start=195901, end=197412,
         data.file="fold1.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold1", file.prefix="Obsfold1",
         external.files = c(NA, NA, "fold1.DailyPreds.dat"))

seasons <- list(3:5,6:8,9:11,c(12,1,2)) 

sim.fold1.summary <- summary(sim.fold1, season.defs=seasons, thresholds=c(0,NA))

obs.fold1.summary <- summary(obs.fold1, season.defs=seasons, thresholds=c(0,NA))


########################################################################
#train2
#TempModel
train2.TempModel.init <-
  read.modeldef("DEFINITIONS/TempModel12Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train2.TempModelDisp.init <-
  read.modeldef("DEFINITIONS/TempModel12Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train2.TempModel <- GLCfit(train2.TempModel.init, dispersion.def=train2.TempModelDisp.init,
                           model.type="normal-heteroscedastic", which.response=2,
                           siteinfo=siteinfo,
                           data.file="train2.PrecipTemp.dat",
                           external.files=c(NA,NA,"DailyPreds.dat"),
                           diagnostics=2,nprev.required=-1,
                           cor.file="train2.TempModel_Corr.dat",
                           resid.file="train2.TempModel_Resids.dat")

print(train2.TempModel)

summary(train2.TempModel, tables=NULL)

save(siteinfo,train2.TempModel,file="train2.Temperature.rda")


#OccModel
train2.OccModel.init <-
  read.modeldef("DEFINITIONS/OccModel13Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train2.OccModel <- GLCfit(train2.OccModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train2.PrecipTemp.dat",
                          diagnostics=1, nprev.required=-1, verbosity=1,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train2.OccModel_Corr.dat")

summary(train2.OccModel, tables=NULL)

save(siteinfo,train2.OccModel,file="train2.Occurrence.rda")


#IntModel
train2.IntModel.init <-
  read.modeldef("DEFINITIONS/IntModel12Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train2.IntModel <- GLCfit(train2.IntModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train2.PrecipTemp.dat",
                          diagnostics=2, nprev.required=-1, verbosity=2,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train2.IntModel_Corr.dat",
                          resid.file = "train2.IntModel_Resids.dat")

summary(train2.IntModel, tables=NULL)

sim.fold2 <-
  GLCsim(list(list(Occurrence=train2.OccModel,Intensity=train2.IntModel),train2.TempModel),
         nsims=100, start=197501, end=199012, impute.until=197412,
         data.file="fold2.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold2", file.prefix="Simfold2",
         external.files = c(NA, NA, "fold2.DailyPreds.dat"))

obs.fold2 <-
  GLCsim(list(list(Occurrence=train2.OccModel,Intensity=train2.IntModel),train2.TempModel),
         nsims=39,start=197501, end=199012,
         data.file="fold2.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold2", file.prefix="Obsfold2",
         external.files = c(NA, NA, "fold2.DailyPreds.dat"))

save(siteinfo,train2.IntModel,file="train2.Intensity.rda")

sim.fold2.summary <- summary(sim.fold2, season.defs=seasons, thresholds=c(0,NA))

obs.fold2.summary <- summary(obs.fold2, season.defs=seasons, thresholds=c(0,NA))

##########################################################################
#train3
#TempModel
train3.TempModel.init <-
  read.modeldef("DEFINITIONS/TempModel12Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train3.TempModelDisp.init <-
  read.modeldef("DEFINITIONS/TempModel12Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train3.TempModel <- GLCfit(train3.TempModel.init, dispersion.def=train3.TempModelDisp.init,
                           model.type="normal-heteroscedastic", which.response=2,
                           siteinfo=siteinfo,
                           data.file="train3.PrecipTemp.dat",
                           external.files=c(NA,NA,"DailyPreds.dat"),
                           diagnostics=2,nprev.required=-1,
                           cor.file="train3.TempModel_Corr.dat",
                           resid.file="train3.TempModel_Resids.dat")

print(train3.TempModel)

summary(train3.TempModel, tables=NULL)

save(siteinfo,train3.TempModel,file="train3.Temperature.rda")


#OccModel
train3.OccModel.init <-
  read.modeldef("DEFINITIONS/OccModel13Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train3.OccModel <- GLCfit(train3.OccModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train3.PrecipTemp.dat",
                          diagnostics=1, nprev.required=-1, verbosity=1,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train3.OccModel_Corr.dat")

summary(train3.OccModel, tables=NULL)

save(siteinfo,train3.OccModel,file="train3.Occurrence.rda")


#IntModel
train3.IntModel.init <-
  read.modeldef("DEFINITIONS/IntModel12Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train3.IntModel <- GLCfit(train3.IntModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train3.PrecipTemp.dat",
                          diagnostics=2, nprev.required=-1, verbosity=2,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train3.IntModel_Corr.dat",
                          resid.file = "train3.IntModel_Resids.dat")

summary(train3.IntModel, tables=NULL)

save(siteinfo,train3.IntModel,file="train3.Intensity.rda")

sim.fold3 <-
  GLCsim(list(list(Occurrence=train3.OccModel,Intensity=train3.IntModel),train3.TempModel),
         nsims=100, start=199101, end=200612, impute.until=199012,
         data.file="fold3.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold3", file.prefix="Simfold3",
         external.files = c(NA, NA, "fold3.DailyPreds.dat"))

obs.fold3 <-
  GLCsim(list(list(Occurrence=train3.OccModel,Intensity=train3.IntModel),train3.TempModel),
         nsims=39,start=199101, end=200612,
         data.file="fold3.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold3", file.prefix="Obsfold3",
         external.files = c(NA, NA, "fold3.DailyPreds.dat"))


sim.fold3.summary <- summary(sim.fold3, season.defs=seasons, thresholds=c(0,NA))

obs.fold3.summary <- summary(obs.fold3, season.defs=seasons, thresholds=c(0,NA))

########################################################################

#train4
#TempModel
train4.TempModel.init <-
  read.modeldef("DEFINITIONS/TempModel12Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train4.TempModelDisp.init <-
  read.modeldef("DEFINITIONS/TempModel12Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)

train4.TempModel <- GLCfit(train4.TempModel.init, dispersion.def=train4.TempModelDisp.init,
                           model.type="normal-heteroscedastic", which.response=2,
                           siteinfo=siteinfo,
                           data.file="train4.PrecipTemp.dat",
                           external.files=c(NA,NA,"DailyPreds.dat"),
                           diagnostics=2,nprev.required=-1,
                           cor.file="train4.TempModel_Corr.dat",
                           resid.file="train4.TempModel_Resids.dat")

print(train4.TempModel)

summary(train4.TempModel, tables=NULL)

save(siteinfo,train4.TempModel,file="train4.Temperature.rda")


#OccModel
train4.OccModel.init <-
  read.modeldef("DEFINITIONS/OccModel13Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train4.OccModel <- GLCfit(train4.OccModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train4.PrecipTemp.dat",
                          diagnostics=1, nprev.required=-1, verbosity=1,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train4.OccModel_Corr.dat")

summary(train4.OccModel, tables=NULL)

save(siteinfo,train4.OccModel,file="train4.Occurrence.rda")


#IntModel
train4.IntModel.init <-
  read.modeldef("DEFINITIONS/IntModel12Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))

train4.IntModel <- GLCfit(train4.IntModel.init,which.response=1,siteinfo=siteinfo,
                          data.file="train4.PrecipTemp.dat",
                          diagnostics=2, nprev.required=-1, verbosity=2,
                          external.files = c(NA, NA, "DailyPreds.dat"),
                          cor.file = "train4.IntModel_Corr.dat",
                          resid.file = "train4.IntModel_Resids.dat")

summary(train4.IntModel, tables=NULL)

save(siteinfo,train4.IntModel,file="train4.Intensity.rda")
#######
sim.fold4 <-
  GLCsim(list(list(Occurrence=train4.OccModel,Intensity=train4.IntModel),train4.TempModel),
         nsims=100, start=200701, end=202112, impute.until=200612,
         data.file="fold4.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold4", file.prefix="Simfold4",
         external.files = c(NA, NA, "fold4.DailyPreds.try.dat"))
##########
obs.fold4 <-
  GLCsim(list(list(Occurrence=train4.OccModel,Intensity=train4.IntModel),train4.TempModel),
         nsims=39,start=200701, end=202112,
         data.file="fold4.PrecipTemp.dat", which.regions=0:1,
         simdir="./SimFiles.fold4", file.prefix="Obsfold4",
         external.files = c(NA, NA, "fold4.DailyPreds.dat"))

sim.fold4.summary <- summary(sim.fold4, season.defs=seasons, thresholds=c(0,NA))

obs.fold4.summary <- summary(obs.fold4, season.defs=seasons, thresholds=c(0,NA))





########################################################################
seasons <- list(3:5,6:8,9:11,c(12,1,2))

#Precipitation
if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=12,height=8)
par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold1.summary, imputation=obs.fold1.summary, which.sites=NULL,
     which.regions="Rhine", which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold2.summary, imputation=obs.fold2.summary, which.sites=NULL,
     which.regions="Rhine", which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold3.summary, imputation=obs.fold1.summary, which.sites=NULL,
     which.regions="Rhine", which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold4.summary, imputation=obs.fold1.summary, which.sites=NULL,
     which.regions="Rhine", which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

#Temperature
par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold1.summary, imputation=obs.fold1.summary, which.sites=NULL,
     which.timescales="daily", which.variables="Temperature",
     which.regions="Rhine",
     which.stats=stat.names[c(1:7,11)],
     quantiles=quantiles, colours.sim=plotcols,
     plot.titles=stat.names[c(1:7,11)],
     ylabs=c(rep(expression(degree*C), 4), rep("Correlation", 4)))

par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold2.summary, imputation=obs.fold2.summary, which.sites=NULL,
     which.timescales="daily", which.variables="Temperature",
     which.regions="Rhine",
     which.stats=stat.names[c(1:7,11)],
     quantiles=quantiles, colours.sim=plotcols,
     plot.titles=stat.names[c(1:7,11)],
     ylabs=c(rep(expression(degree*C), 4), rep("Correlation", 4)))

par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold3.summary, imputation=obs.fold3.summary, which.sites=NULL,
     which.timescales="daily", which.variables="Temperature",
     which.regions="Rhine",
     which.stats=stat.names[c(1:7,11)],
     quantiles=quantiles, colours.sim=plotcols,
     plot.titles=stat.names[c(1:7,11)],
     ylabs=c(rep(expression(degree*C), 4), rep("Correlation", 4)))

par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold4.summary, imputation=obs.fold4.summary, which.sites=NULL,
     which.timescales="daily", which.variables="Temperature",
     which.regions="Rhine",
     which.stats=stat.names[c(1:7,11)],
     quantiles=quantiles, colours.sim=plotcols,
     plot.titles=stat.names[c(1:7,11)],
     ylabs=c(rep(expression(degree*C), 4), rep("Correlation", 4)))



#########
# Annual graphs

# Precipitation
# Fold 1
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2,0.75,0))
season.names <-
  c("March, April, May", "June, July, August",
    "September, October, November", "December, January, February")
plot(sim.fold1.summary, imputation=obs.fold1.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Precipitation",
     quantiles=quantiles, colours.sim="colour",
     ylabs=rep("mm",4), plot.titles=season.names)

# Fold 2
plot(sim.fold2.summary, imputation=obs.fold2.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Precipitation",
     quantiles=quantiles, colours.sim="colour",
     ylabs=rep("mm",4), plot.titles=season.names)

# Fold 3
plot(sim.fold3.summary, imputation=obs.fold3.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Precipitation",
     quantiles=quantiles, colours.sim="colour",
     ylabs=rep("mm",4), plot.titles=season.names)

# Fold 4
plot(sim.fold4.summary, imputation=obs.fold4.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Precipitation",
     quantiles=quantiles, colours.sim="colour",
     ylabs=rep("mm",4), plot.titles=season.names)



##################

# Temperature

# Fold 1
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2,0.75,0))
plot(sim.fold1.summary, imputation=obs.fold1.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Temperature",
     quantiles=quantiles, colours.sim=plotcols,
     ylabs=rep(expression(degree*C),4), plot.titles=season.names)

# Fold 2
plot(sim.fold2.summary, imputation=obs.fold2.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Temperature",
     quantiles=quantiles, colours.sim=plotcols,
     ylabs=rep(expression(degree*C),4), plot.titles=season.names)

# Fold 3
plot(sim.fold3.summary, imputation=obs.fold3.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Temperature",
     quantiles=quantiles, colours.sim=plotcols,
     ylabs=rep(expression(degree*C),4), plot.titles=season.names)

# Fold 4
plot(sim.fold4.summary, imputation=obs.fold4.summary, which.regions="Rhine",
     which.timescales="monthly", which.variable="Temperature",
     quantiles=quantiles, colours.sim=plotcols,
     ylabs=rep(expression(degree*C),4), plot.titles=season.names)



###########################################################################
# Precipitation summary statistics for a single site
######
# Low altitude: S096
if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(2,0.75,0))


plot(sim.fold1.summary,imputation=obs.fold1.summary,which.sites="S096",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold2.summary,imputation=obs.fold2.summary,which.sites="S096",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold3.summary,imputation=obs.fold3.summary,which.sites="S096",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold4.summary,imputation=obs.fold4.summary,which.sites="S096",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))


# High altitude: S028

plot(sim.fold1.summary,imputation=obs.fold1.summary,which.sites="S028",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold2.summary,imputation=obs.fold2.summary,which.sites="S028",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold3.summary,imputation=obs.fold3.summary,which.sites="S028",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))

plot(sim.fold4.summary,imputation=obs.fold4.summary,which.sites="S028",
     which.regions=NULL, which.timescales="daily",
     which.variables="Precipitation", which.stats=stat.names[c(1:3,5:10)],
     colours.sim="colour", plot.titles=stat.names[c(1:3,5:10)],
     ylabs=c(rep("mm",3), rep("Correlation", 3), "Probability", rep("mm", 2)))