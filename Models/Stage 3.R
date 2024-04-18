library(Rglimclim)
library(maps)
library(lattice)

# save it as DailyPreds dataframe
load("DailyPredictors.rda")
write.GLCexternal(DailyPreds, file="DailyPreds.dat")

#Temperature model 9: including large-scale monthly temperature effects
#choose the predictors with the maximum likelihood
write.modeldef(TempModel6, file=c("TempModel9Init.def",
                                  "TempModel9Init_Dispersion.def")) # Copied to DEFINITIONS and edited
TempModel9.init <-
  read.modeldef("DEFINITIONS/TempModel9Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel9.init)

TempModel9Disp.init <-
  read.modeldef("DEFINITIONS/TempModel9Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel9Disp.init)

TempModel9 <- GLCfit(TempModel9.init, dispersion.def=TempModel9Disp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      external.files=c(NA,NA,"DailyPreds.dat"),
                      diagnostics=2,nprev.required=-1,
                      cor.file="TempModel9_Corr.dat",
                      resid.file="TempModel9_Resids.dat")

print(TempModel9)
summary(TempModel9, tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel9,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel9,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=4))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel9,which.plots=4:5,plot.cols=c("blue","purple"))
#QQ plot becomes worse


#main: delete the last interaction
#dispersion: delete 2nd interaction
#Temperature model 10: delete insignificant interactions
write.modeldef(TempModel9, file=c("TempModel10Init.def",
                                  "TempModel10Init_Dispersion.def")) # Copied to DEFINITIONS and edited
TempModel10.init <-
  read.modeldef("DEFINITIONS/TempModel10Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel10.init)

TempModel10Disp.init <-
  read.modeldef("DEFINITIONS/TempModel10Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel10Disp.init)

TempModel10 <- GLCfit(TempModel10.init, dispersion.def=TempModel10Disp.init,
                     model.type="normal-heteroscedastic", which.response=2,
                     siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     external.files=c(NA,NA,"DailyPreds.dat"),
                     diagnostics=2,nprev.required=-1,
                     cor.file="TempModel10_Corr.dat",
                     resid.file="TempModel10_Resids.dat")

anova(TempModel10, TempModel9)
print(TempModel10)
summary(TempModel10, tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel10,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel10,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=4))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel10,which.plots=4:5,plot.cols=c("blue","purple"))

#Temperature model 11: include other covariates
write.modeldef(TempModel10, file=c("TempModel11Init.def",
                                  "TempModel11Init_Dispersion.def")) # Copied to DEFINITIONS and edited
TempModel11.init <-
  read.modeldef("DEFINITIONS/TempModel11Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel11.init)

TempModel11Disp.init <-
  read.modeldef("DEFINITIONS/TempModel11Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel11Disp.init)

TempModel11 <- GLCfit(TempModel11.init, dispersion.def=TempModel11Disp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      external.files=c(NA,NA,"DailyPreds.dat"),
                      diagnostics=2,nprev.required=-1,
                      cor.file="TempModel11_Corr.dat",
                      resid.file="TempModel11_Resids.dat")

print(TempModel11)
summary(TempModel11, tables=NULL)


par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel11,which.plots=4:5,plot.cols=c("blue","purple"))

#Temperature model 12:delete atmospheric predictors in dispersion model
write.modeldef(TempModel11, file=c("TempModel12Init.def",
                                   "TempModel12Init_Dispersion.def")) # Copied to DEFINITIONS and edited
TempModel12.init <-
  read.modeldef("DEFINITIONS/TempModel12Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo,which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel12.init)

TempModel12Disp.init <-
  read.modeldef("DEFINITIONS/TempModel12Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                external.files=c(NA,NA,"DailyPreds.dat"),
                oldGlimClim.warning = FALSE)
print(TempModel12Disp.init)

TempModel12 <- GLCfit(TempModel12.init, dispersion.def=TempModel12Disp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      external.files=c(NA,NA,"DailyPreds.dat"),
                      diagnostics=2,nprev.required=-1,
                      cor.file="TempModel12_Corr.dat",
                      resid.file="TempModel12_Resids.dat")

print(TempModel12)
summary(TempModel12, tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel12,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel12,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=4))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel12,which.plots=4:5,plot.cols=c("blue","purple"))

save(siteinfo,TempModel12,file="Temperature.rda")




#Precipitation models including large-scale covariates
#Precipitation occurrence
#Occurrence model 11: incorporating large-scale covariates, plus interactions
write.modeldef(OccModel10,file="OccModel11Init.def") # Copied to DEFINITIONS/ and edited
OccModel11.init <-
  read.modeldef("DEFINITIONS/OccModel11Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(OccModel11.init)

OccModel11 <- GLCfit(OccModel11.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1, nprev.required=-1, verbosity=2,
                    external.files = c(NA, NA, "DailyPreds.dat"))

print(OccModel11)

summary(OccModel11, tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel11,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(OccModel11,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=2))


#Occurrence model 12: delete tag 3
write.modeldef(OccModel11,file="OccModel12Init.def")

OccModel12.init <-
  read.modeldef("DEFINITIONS/OccModel12Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(OccModel12.init)

OccModel12 <- GLCfit(OccModel12.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1, verbosity=2,
                     external.files = c(NA, NA, "DailyPreds.dat"))

print(OccModel12)

summary(OccModel12, tables=NULL)

anova(OccModel12, OccModel11)

#OccModel13: remove interaction LP:Temp700
write.modeldef(OccModel12,file="OccModel13Init.def")

OccModel13.init <-
  read.modeldef("DEFINITIONS/OccModel13Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(OccModel13.init)

OccModel13 <- GLCfit(OccModel13.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1, verbosity=1,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "OccModel13_Corr.dat")

print(OccModel13)

summary(OccModel13, tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel13,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(OccModel13,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=2))

plot(OccModel13,which.plots=5,plot.cols=c("blue","purple"), distance.units="degrees")

save(siteinfo,OccModel13,file="Occurrence.rda")



###
#Precipitation intensity
#Intensity model 10: incorporating large-scale covariates, plus interactions
write.modeldef(IntModel9,file="IntModel10Init.def") # Copied to DEFINITIONS/ and edited
IntModel10.init <-
  read.modeldef("DEFINITIONS/IntModel10Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(IntModel10.init)

IntModel10 <- GLCfit(IntModel10.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1, verbosity=2,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "IntModel10_Corr.dat")

print(IntModel10)

summary(IntModel10, tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel10,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel10,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=2))

#Intensity model 11: delete 850hPa temperature
write.modeldef(IntModel10,file="IntModel11Init.def") # Copied to DEFINITIONS/ and edited
IntModel11.init <-
  read.modeldef("DEFINITIONS/IntModel11Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(IntModel11.init)

IntModel11 <- GLCfit(IntModel11.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1, verbosity=2,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "IntModel11_Corr.dat")

print(IntModel11)

summary(IntModel11, tables=NULL)


# 4&5, 6&7, LP:GeoPot_1000
# IntModel11: Final(delete insignificant interactions)
write.modeldef(IntModel11,file="IntModel12Init.def") # Copied to DEFINITIONS/ and edited
IntModel12.init <-
  read.modeldef("DEFINITIONS/IntModel12Init.def",model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE,
                external.files = c(NA, NA, "DailyPreds.dat"))
print(IntModel12.init)

IntModel12 <- GLCfit(IntModel12.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=2, nprev.required=-1, verbosity=2,
                     external.files = c(NA, NA, "DailyPreds.dat"),
                     cor.file = "IntModel12_Corr.dat",
                     resid.file = "IntModel12_Resids.dat")

print(IntModel12)

summary(IntModel12, tables=NULL)
summary(IntModel12, tables="site")

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel12,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel12,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=2))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(IntModel12,which.plots=4:5,plot.cols=c("blue","purple"))

save(siteinfo,IntModel12,file="Intensity.rda")

file.remove(list.files(pattern="def$"))
file.remove(list.files(pattern="Corr.dat$"))
file.remove(list.files(pattern="Resids.dat$"))
if (Rver >= "3.6.1") options(warn=warn.settings) # And reset warnings if necessary