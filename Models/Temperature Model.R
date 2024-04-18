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


###############################################################################
#Temperature model 0: homoscedastic null model with constant only
write.modeldef(OccModel0,file="TempModel0Init.def") # Copied to DEFINITIONS/ and edited
TempModel0.init <- read.modeldef("DEFINITIONS/TempModel0Init.def",
                                 model.type="normal", siteinfo=siteinfo,
                                 which.response=2,
                                 var.names=c("Precipitation","Temperature"))
print(TempModel0.init)

TempModel0 <- GLCfit(TempModel0.init, which.response=2,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1)

summary(TempModel0,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel0,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel0,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))


#Temperature model 1: simple seasonality
write.modeldef(TempModel0,file="TempModel1Init.def") # Copied to DEFINITIONS/ and edited
TempModel1.init <- read.modeldef("DEFINITIONS/TempModel1Init.def",
                                 model.type="normal", siteinfo=siteinfo,
                                 which.response=2,
                                 var.names=c("Precipitation","Temperature"))
print(TempModel1.init)

TempModel1 <- GLCfit(TempModel1.init, which.response=2,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1)

summary(TempModel1,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel1,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel1,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))

#Temperature model 2: seasonal harmonic and altitude effect
write.modeldef(TempModel1,file="TempModel2Init.def") # Copied to DEFINITIONS/ and edited
TempModel2.init <- read.modeldef("DEFINITIONS/TempModel2Init.def",
                                 model.type="normal",
                                 siteinfo=siteinfo,
                                 which.response=2,
                                 var.names=c("Precipitation","Temperature"))
print(TempModel2.init)

TempModel2 <- GLCfit(TempModel2.init, which.response=2,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1)

print(TempModel2)

summary(TempModel2,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel2,which.plots=1:2,lwd=2)
#Feb, Aug, Sept significant

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel2,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))


#Temperature model 3: seasonal harmonic and altitude effect, interaction of Legendre polynomial, limits for polynomial
write.modeldef(TempModel2,file="TempModel3Init.def") # Copied to DEFINITIONS/ and edited
TempModel3.init <- read.modeldef("DEFINITIONS/TempModel3Init.def",
                                 model.type="normal",
                                 siteinfo=siteinfo,
                                 which.response=2,
                                 var.names=c("Precipitation","Temperature"))
print(TempModel3.init)

TempModel3 <- GLCfit(TempModel3.init, which.response=2,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1)

print(TempModel3)

summary(TempModel3,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel3,which.plots=1:2,lwd=2)
#Sept, Oct significant

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel3,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))

anova(TempModel3, TempModel2)


#Temperature model 4a: adding previous day??s temperature, interaction between previous days and daily annual cycle
write.modeldef(TempModel3,file="TempModel4aInit.def") # Copied to DEFINITIONS/ and edited
TempModel4a.init <- read.modeldef("DEFINITIONS/TempModel4aInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"))
print(TempModel4a.init)

TempModel4a <- GLCfit(TempModel4a.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)

print(TempModel4a)

summary(TempModel4a,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel4a,which.plots=1:2,lwd=2)
#No significant month, but regular spikes

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel4a,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))

#Different transformations

#Temperature model 4b:#10
write.modeldef(TempModel4a,file="TempModel4bInit.def") # Copied to DEFINITIONS/ and edited
TempModel4b.init <- read.modeldef("DEFINITIONS/TempModel4bInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4b.init)

TempModel4b <- GLCfit(TempModel4b.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)

print(TempModel4b)

summary(TempModel4b,tables=NULL)

#Temperature model 4c:#20
write.modeldef(TempModel4b,file="TempModel4cInit.def") # Copied to DEFINITIONS/ and edited
TempModel4c.init <- read.modeldef("DEFINITIONS/TempModel4cInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4c.init)

TempModel4c <- GLCfit(TempModel4c.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)

print(TempModel4c)

summary(TempModel4c,tables=NULL)
#log-likelihood improves, exponential decay rate and its standard error is not very large,
#use this transformation

#Temperature model 4d:#30
write.modeldef(TempModel4c,file="TempModel4dInit.def") # Copied to DEFINITIONS/ and edited
TempModel4d.init <- read.modeldef("DEFINITIONS/TempModel4dInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4d.init)

TempModel4d <- GLCfit(TempModel4d.init, which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)

print(TempModel4d)

summary(TempModel4d,tables=NULL)
#the estimates of the offsets are very close to zero, especially considering their standard errors. 
#This suggests that don??t need to use the ??shifted?? version of the model.

#Temperature model 4e-g: lag 2-4
#Temperature model 4e: lag 2
write.modeldef(TempModel4c,file="TempModel4eInit.def") # Copied to DEFINITIONS/ and edited
TempModel4e.init <- read.modeldef("DEFINITIONS/TempModel4eInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4e.init)

TempModel4e <- GLCfit(TempModel4e.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)


#Temperature model 4f: lag 3
write.modeldef(TempModel4e,file="TempModel4fInit.def") # Copied to DEFINITIONS/ and edited
TempModel4f.init <- read.modeldef("DEFINITIONS/TempModel4fInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4f.init)

TempModel4f <- GLCfit(TempModel4f.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)


#Temperature model 4g: lag 4
write.modeldef(TempModel4f,file="TempModel4gInit.def") # Copied to DEFINITIONS/ and edited
TempModel4g.init <- read.modeldef("DEFINITIONS/TempModel4gInit.def",
                                  model.type="normal", siteinfo=siteinfo,
                                  which.response=2,
                                  var.names=c("Precipitation","Temperature"),
                                  oldGlimClim.warning=FALSE)
print(TempModel4g.init)

TempModel4g <- GLCfit(TempModel4g.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=4)

anova(TempModel4c, TempModel4e, TempModel4f, TempModel4g)
#Lag 3 best

TempModel4f <- GLCfit(TempModel4f.init, which.response=2,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=-1)

print(TempModel4f)

summary(TempModel4f,tables=NULL)

if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel4f,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel4f,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))

site.resids <- TempModel4f$Residuals$Pearson$Site.table$Mean
sites.wanted <- !is.na(site.resids)
par(mfrow=c(4,4),mar=c(2,2,3,0.5),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
for (i in 3:16) {
  tmp <- siteinfo$Attribute.values[,i]
  plot(tmp[sites.wanted],site.resids[sites.wanted],xlab="",ylab="",
       pch=15,col="lightblue",main=siteinfo$Attribute.names[i])
  text(quantile(tmp[sites.wanted],0.98,na.rm=TRUE),0.11,adj=1,
       paste("Corr:",round(cor(tmp,site.resids,use="complete.obs"),2)))
}



#Temperature model 5: add additional interaction
write.modeldef(TempModel4f,file="TempModel5Init.def") # Copied to DEFINITIONS/ and edited
TempModel5.init <- read.modeldef("DEFINITIONS/TempModel5Init.def",
                                 model.type="normal", siteinfo=siteinfo,
                                 which.response=2,
                                 var.names=c("Precipitation","Temperature"),
                                 oldGlimClim.warning=FALSE)
print(TempModel5.init)

TempModel5 <- GLCfit(TempModel5.init, which.response=2,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1)

summary(TempModel5,tables=NULL)

print(TempModel5)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel5,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel5,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="none",scale=4))

# Temperature model 6: combine mean model 6 with seasonal and regional variation in residual variance
write.modeldef(TempModel5,file="TempModel6Init.def")

TempModel6.init <-
  read.modeldef("DEFINITIONS/TempModel6Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel6.init)

write.modeldef(TempModel3,file="TempModel6Init_Dispersion.def")

TempModel6Disp.init <-
  read.modeldef("DEFINITIONS/TempModel6Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"))
print(TempModel6Disp.init)

TempModel6 <- GLCfit(TempModel6.init, dispersion.def=TempModel6Disp.init,
                     model.type="normal-heteroscedastic", which.response=2,
                     siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=2,nprev.required=-1,
                     cor.file="TempModel6_Corr.dat",
                     resid.file="TempModel6_Resids.dat")

print(TempModel6)

summary(TempModel6,tables="site")
summary(TempModel6,tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel6,which.plots=1:2,lwd=2)


par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="Longitude",ylab="Latitude",xlim=LongLims,ylim=LatLims)
plot(TempModel6,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel6,which.plots=4:5,plot.cols=c("blue","purple"))

par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
par(mfrow=c(2,2),mar=c(3,4,3,1),mgp=c(2,0.75,0),lwd=2,ask=FALSE)
plot(1:12, TempModel7$Residuals$Pearson$Month.table$Mean, type="l",
     xlab="Month", ylab="Mean\n", main="Monthly residual means",
     ylim=c(-0.06, 0.06), las=1)
lines(1:12, 1.96*TempModel7$Residuals$Pearson$Month.table$`S.E. mean`, lty=2)
lines(1:12, -1.96*TempModel7$Residuals$Pearson$Month.table$`S.E. mean`, lty=2)
abline(h=0,lwd=3)
plot(1:12, TempModel7$Residuals$Pearson$Month.table$`Std Dev`, type="l",
     xlab="Month", ylab="Std Dev\n", main="Monthly residual standard deviations",
     ylim=c(0, 1.2), las=1)
abline(h=1,lwd=3)
Years <- rownames(TempModel6$Residuals$Pearson$Year.table)
plot(Years, TempModel6$Residuals$Pearson$Year.table$Mean, type="l",
     xlab="Year", ylab="Mean\n", main="Annual residual means",
     ylim=c(-0.15, 0.15), las=1)
lines(Years, 1.96*TempModel6$Residuals$Pearson$Year.table$`S.E. mean`, lty=2)
lines(Years, -1.96*TempModel6$Residuals$Pearson$Year.table$`S.E. mean`, lty=2)
abline(h=0,lwd=3)
plot(Years, TempModel6$Residuals$Pearson$Year.table$`Std Dev`, type="l",
     xlab="Year", ylab="Std Dev\n", main="Annual residual standard deviations",
     ylim=c(0, 1.2), las=1)
abline(h=1,lwd=3)