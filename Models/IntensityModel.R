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

####################################################################
IntModel0.init <- read.modeldef("DEFINITIONS/OccModel0Init.def",
                                model.type="gamma", siteinfo=siteinfo,
                                var.names=c("Precipitation","Temperature"))
IntModel0 <- GLCfit(IntModel0.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

summary(IntModel0,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel0,which.plots=1:2,lwd=2)


#Intensity model 1: model with simple seasonal cycle
write.modeldef(IntModel0,file="IntModel1Init.def")
IntModel1.init <-read.modeldef("DEFINITIONS/IntModel1Init.def", model.type="gamma",
                               siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel1.init)

IntModel1 <- GLCfit(IntModel1.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(IntModel1)

summary(IntModel1,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel1,which.plots=1:2,lwd=2)

#Intensity model 2: model with enhanced seasonal cycle
write.modeldef(IntModel1,file="IntModel2Init.def")
IntModel2.init <- read.modeldef("DEFINITIONS/IntModel2Init.def", model.type="gamma",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel2.init)

IntModel2 <- GLCfit(IntModel2.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(IntModel2)

summary(IntModel2,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel2,which.plots=1:2,lwd=2)



#Intensity model 3: initial representation of regional variation
write.modeldef(IntModel2,file="IntModel3Init.def")
IntModel3.init <- read.modeldef("DEFINITIONS/IntModel3Init.def", model.type="gamma",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel3.init)

IntModel3 <- GLCfit(IntModel3.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(IntModel3)

summary(IntModel3,tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel3,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
site.resids <- IntModel3$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.2,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs", method="spearman"),2)))
}


#Intensity model 4: add 100 km^2 mean altitude
siteinfo$Attribute.names

write.modeldef(IntModel3,file="IntModel4Init.def")
IntModel4.init <- read.modeldef("DEFINITIONS/IntModel4Init.def", model.type="gamma",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel4.init)

IntModel4 <- GLCfit(IntModel4.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=-1)

print(IntModel4)

summary(IntModel4,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
site.resids <- IntModel4$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.15,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs", method="spearman"),2)))
}


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel4,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=3))


#Intensity model 5a: representing autocorrelation
write.modeldef(IntModel4,file="IntModel5aInit.def")
IntModel5a.init <- read.modeldef("DEFINITIONS/IntModel5aInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel5a.init)

IntModel5a <- GLCfit(IntModel5a.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5a)

summary(IntModel5a,tables=NULL)


#Intensity model 5b: adding interaction between transformation and seasonality
write.modeldef(IntModel5a,file="IntModel5bInit.def")
IntModel5b.init <- read.modeldef("DEFINITIONS/IntModel5bInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"))
print(IntModel5b.init)

IntModel5b <- GLCfit(IntModel5b.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5b)

summary(IntModel5b,tables=NULL)


#Intensity model 5c: transformations #22
write.modeldef(IntModel5b,file="IntModel5cInit.def")
IntModel5c.init <- read.modeldef("DEFINITIONS/IntModel5cInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5c.init)

IntModel5c <- GLCfit(IntModel5c.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5c)

summary(IntModel5c,tables=NULL)
#The standard error of exponential decay rate is too large

#Intensity model 5d: transformations #112
write.modeldef(IntModel5b,file="IntModel5dInit.def")
IntModel5d.init <- read.modeldef("DEFINITIONS/IntModel5dInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5d.init)

IntModel5d <- GLCfit(IntModel5d.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5d)

summary(IntModel5d,tables=NULL)


#Intensity model 5e: transformations #122
write.modeldef(IntModel5d,file="IntModel5eInit.def")
IntModel5e.init <- read.modeldef("DEFINITIONS/IntModel5eInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5e.init)

IntModel5e <- GLCfit(IntModel5e.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5e)

summary(IntModel5e,tables=NULL)
#Standard error of exponential decay rate is too large


#Models 5f to 5h: lags 2 through 4
#Model 5f: lags 2
write.modeldef(IntModel5d,file="IntModel5fInit.def")
IntModel5f.init <- read.modeldef("DEFINITIONS/IntModel5fInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5f.init)

IntModel5f <- GLCfit(IntModel5f.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5f)

summary(IntModel5f,tables=NULL)


#Model 5g: lags 3
write.modeldef(IntModel5f,file="IntModel5gInit.def")
IntModel5g.init <- read.modeldef("DEFINITIONS/IntModel5gInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5g.init)

IntModel5g <- GLCfit(IntModel5g.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5g)

summary(IntModel5g,tables=NULL)


#Model 5h: lags 4
write.modeldef(IntModel5g,file="IntModel5hInit.def")
IntModel5h.init <- read.modeldef("DEFINITIONS/IntModel5hInit.def", model.type="gamma",
                                 siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(IntModel5h.init)

IntModel5h <- GLCfit(IntModel5h.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=4)

print(IntModel5h)

summary(IntModel5h,tables=NULL)

anova(IntModel5d, IntModel5f, IntModel5g, IntModel5h)
#Lag 1 is the best

print(IntModel5d)

summary(IntModel5d, tables="site")

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel5d,which.plots=1:2,lwd=2)


#Intensity Model 6: Adding interaction
write.modeldef(IntModel5d,file="IntModel6Init.def")
IntModel6.init <- read.modeldef("DEFINITIONS/IntModel6Init.def", model.type="gamma",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                oldGlimClim.warning=FALSE)
print(IntModel6.init)

IntModel6 <- GLCfit(IntModel6.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1,nprev.required=4)

print(IntModel6)

summary(IntModel5d,tables = NULL)
summary(IntModel6,tables = NULL)

anova(IntModel5d, IntModel6)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel6,which.plots=1:2,lwd=2)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel6,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=2))


#Intensity Model 7: Final
write.modeldef(IntModel6,file="IntModel7Init.def")
IntModel7.init <- read.modeldef("DEFINITIONS/IntModel7Init.def", model.type="gamma",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                oldGlimClim.warning=FALSE)
print(IntModel7.init)

IntModel7 <- GLCfit(IntModel7.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=2,nprev.required=-1,
                    cor.file="IntModel7_Corr.dat",
                    resid.file="IntModel7_Resids.dat")


print(IntModel7)

summary(IntModel7, tables=NULL)
summary(IntModel7, tables="site")

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel7,which.plots=1:2,lwd=2)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="Longitude",ylab="Latitude",xlim=LongLims,ylim=LatLims)
plot(IntModel7,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(IntModel7,which.plots=4:5,plot.cols=c("blue","purple"))

save(siteinfo, IntModel7, file="Intensity_marginal.rda")