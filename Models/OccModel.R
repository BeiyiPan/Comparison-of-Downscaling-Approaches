###Package loading and station definition

library(Rglimclim)
library(maps)
library(lattice)

load("Topography.rda")
RhineData <- read.GLCdata("PrecipTemp_Processed.dat")

summary(topo.lats)
summary(topo.longs)
#Limits of region
LatLims <- c(50.1,52)
LongLims <- c(5.6,8.9)

station.data <- read.table("stations.dat",header=TRUE,
                           stringsAsFactors = FALSE)

# Identify the stations that don't appear in the "fitting" dataset
station.data$TestStation <-
  as.numeric(!(station.data$Scode %in%
                 read.GLCdata("PrecipTemp_Processed.dat")$Site))

site.attributes <- c("Longitude","Latitude","Altitude (100m)",
                     "Mapped Altitude (100m)",
                     "10km^2 mean altitude",
                     "100km^2 mean altitude",
                     "1000km^2 mean altitude",
                     "10km^2 altitude std dev",
                     "100km^2 altitude std dev",
                     "1000km^2 altitude std dev",
                     "10km^2 E-W slope (100m / km)",
                     "100km^2 E-W slope (100m / km)",
                     "1000km^2 E-W slope (100m / km)",
                     "10km^2 N-S slope (100m / km)",
                     "100km^2 N-S slope (100m / km)",
                     "1000km^2 N-S slope (100m / km)")


siteinfo <-
  make.siteinfo(station.data,coord.cols=c(3,4),
                site.codes=1,site.names=2,
                attr.names=site.attributes,
                region.col=which(names(station.data)=="TestStation"),
                regions=define.regions(c("Rhine",
                                         "Validation stations"),
                                       codes=0:1))

print(siteinfo)


#Model-building: stage 1(marginal baseline models)
data(GLCdemo)
print(ConstantModel)                       

write.modeldef(ConstantModel,file = "OccModel0Init.def")

#Prevent warnings
Rver <- paste(R.version$major, R.version$minor, sep=".")
if (Rver >= "3.6.1") {
  warn.settings <- options()$warn
  options(warn=-1)
}

OccModel0.init <-
  read.modeldef("DEFINITIONS/OccModel0Init.def", model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel0.init)

OccModel0 <- GLCfit(OccModel0.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel0)

summary(OccModel0,tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel0,which.plots=1:2,lwd=2)


#Occurrence model 1: model with simple seasonal cycle
write.modeldef(OccModel0,file="OccModel1Init.def")
OccModel1.init <-read.modeldef("DEFINITIONS/OccModel1Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel1.init)

OccModel1 <- GLCfit(OccModel1.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel1)

summary(OccModel1,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel1,which.plots=1:2,lwd=2)


#Occurrence model 2: improved seasonal cycle (first harmonic)
write.modeldef(OccModel1,file="OccModel2Init.def")
OccModel2.init <-read.modeldef("DEFINITIONS/OccModel2Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel2.init)

OccModel2 <- GLCfit(OccModel2.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel2)

summary(OccModel2,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel2,which.plots=1:2,lwd=2)
#Spike in March

#Occurrence model 3: regional variation
if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))

image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel2,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))



write.modeldef(OccModel2,file="OccModel3Init.def")
OccModel3.init <-read.modeldef("DEFINITIONS/OccModel3Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel3.init)

OccModel3 <- GLCfit(OccModel3.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel3)

summary(OccModel3,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel3,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))

anova(OccModel3,OccModel2)



if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
site.resids <- OccModel3$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.1,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs"),2)))
}
#altitude highest

siteinfo$Attribute.names

########################################################
#Occurrence model 4: add 100 km^2 mean altitude
write.modeldef(OccModel3,file="OccModel4Init.def")
OccModel4.init <-read.modeldef("DEFINITIONS/OccModel4Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel4.init)

OccModel4 <- GLCfit(OccModel4.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel4)

summary(OccModel4,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel4,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))

anova(OccModel4,OccModel3)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
site.resids <- OccModel4$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.05,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs"),2)))
}



#Occurrence model 5: add 1000 km^2 E-W slope and 1000 km^2 N-S slope
write.modeldef(OccModel4,file="OccModel5Init.def")
OccModel5.init <-read.modeldef("DEFINITIONS/OccModel5Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel5.init)

OccModel5 <- GLCfit(OccModel5.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(OccModel5)

summary(OccModel5,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel5,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))

anova(OccModel5,OccModel4)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
site.resids <- OccModel5$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.05,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs"),2)))
}


#Occurrence model 6: Adding autocorrelation
write.modeldef(OccModel5,file="OccModel6Init.def")
OccModel6.init <-read.modeldef("DEFINITIONS/OccModel6Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(OccModel6.init)

#Modify nprev.required to 4
OccModel6 <- GLCfit(OccModel6.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=4)

print(OccModel6)
summary(OccModel6,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel6,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


#OccModel6b:  proportion of wet sites yesterday #13
write.modeldef(OccModel6,file="OccModel6bInit.def")
OccModel6b.init <-read.modeldef("DEFINITIONS/OccModel6bInit.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                               oldGlimClim.warning=FALSE)
print(OccModel6b.init)

OccModel6b <- GLCfit(OccModel6b.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=4)

summary(OccModel6b,tables=NULL)
print(OccModel6b)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel6b,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


#OccModel6c:  proportion of wet sites yesterday #23
write.modeldef(OccModel6,file="OccModel6cInit.def")
OccModel6c.init <-read.modeldef("DEFINITIONS/OccModel6cInit.def", model.type="logistic",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                oldGlimClim.warning=FALSE)
print(OccModel6c.init)

OccModel6c <- GLCfit(OccModel6c.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

summary(OccModel6c,tables=NULL)
print(OccModel6c)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel6c,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


#OccModel6d:#33
write.modeldef(OccModel6c,file="OccModel6dInit.def")
OccModel6d.init <-read.modeldef("DEFINITIONS/OccModel6dInit.def", model.type="logistic",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                oldGlimClim.warning=FALSE)
print(OccModel6d.init)

OccModel6d <- GLCfit(OccModel6d.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

summary(OccModel6d,tables=NULL)

print(OccModel6d)
#standard error too high of two parameters are too high
if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel6c,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


#Model6.2days: 2 previous days
write.modeldef(OccModel6c,file="OccModel6.2days.Init.def")

OccModel6.2days.init <-read.modeldef("DEFINITIONS/OccModel6.2days.Init.def", model.type="logistic",
                                      siteinfo=siteinfo, var.names=c("Precipitation","Temperature")
                                      ,oldGlimClim.warning=FALSE)
print(OccModel6.2days.init)

OccModel6.2days <- GLCfit(OccModel6.2days.init,which.response=1,siteinfo=siteinfo,
                           data.file="PrecipTemp_Processed.dat",
                           diagnostics=1,nprev.required=4)

summary(OccModel6.2days,tables=NULL)



#Model11.3days: 3 previous days
write.modeldef(OccModel6.2days,file="OccModel6.3days.Init.def")

OccModel6.3days.init <-read.modeldef("DEFINITIONS/OccModel6.3days.Init.def", model.type="logistic",
                                      siteinfo=siteinfo, var.names=c("Precipitation","Temperature")
                                      ,oldGlimClim.warning=FALSE)
print(OccModel6.3days.init)

OccModel6.3days <- GLCfit(OccModel6.3days.init,which.response=1,siteinfo=siteinfo,
                           data.file="PrecipTemp_Processed.dat",
                           diagnostics=1,nprev.required=4)

summary(OccModel6.3days,tables=NULL)



#Model6.4days: 4 previous days
write.modeldef(OccModel6.3days,file="OccModel6.4days.Init.def")

OccModel6.4days.init <-read.modeldef("DEFINITIONS/OccModel6.4days.Init.def", model.type="logistic",
                                      siteinfo=siteinfo, var.names=c("Precipitation","Temperature")
                                      ,oldGlimClim.warning=FALSE)
print(OccModel6.4days.init)

OccModel6.4days <- GLCfit(OccModel6.4days.init,which.response=1,siteinfo=siteinfo,
                           data.file="PrecipTemp_Processed.dat",
                           diagnostics=1,nprev.required=4)

summary(OccModel6.4days,tables=NULL)

anova(OccModel6c, OccModel6.2days, OccModel6.3days, OccModel6.4days)
#3days is the best
print(OccModel6.3days)

summary(OccModel6.3days,tables="site")


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel6.3days,which.plots=1:2,lwd=2)


par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel6.3days,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))

#seasonality in residual standard deviation


#Occurence Model7: Interaction
write.modeldef(OccModel6.3days,file="OccModel7Init.def")
OccModel7.init <-read.modeldef("DEFINITIONS/OccModel7Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                               oldGlimClim.warning=FALSE)
print(OccModel7.init)

OccModel7 <- GLCfit(OccModel7.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=4)

summary(OccModel7,tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)

plot(OccModel7,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))

# 4 and 5,6 and 7,????1 and 2
#Occurrence Model 8: Final
write.modeldef(OccModel7,file="OccModel8Init.def")
OccModel8.init <-read.modeldef("DEFINITIONS/OccModel8Init.def", model.type="logistic",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                               oldGlimClim.warning=FALSE)
print(OccModel8.init)

OccModel8 <- GLCfit(OccModel8.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1, verbosity=1,
                    cor.file="OccModel8_Corr.dat")

summary(OccModel8,tables=NULL)
summary(OccModel8,tables="site")
print(OccModel8)



if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel8,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="Longitude",ylab="Latitude",xlim=LongLims,ylim=LatLims)
plot(OccModel8,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=6,height=4)


pdf("Occ_Inter-site.pdf",width = 8, height = 6)
par(mfrow=c(1,1),mar=c(4,3,3,3))
plot(OccModel8,which.plots=5,plot.cols=c("blue","purple"), distance.units="degrees")
dev.off()

save(siteinfo, OccModel8, file="Occurrence_marginal.rda")
