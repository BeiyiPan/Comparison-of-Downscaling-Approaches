library(Rglimclim)
library(maps)
library(lattice)

#OccModel9a: including temperature
write.modeldef(OccModel8,file="OccModel9aInit.def") # Copied to DEFINITIONS/ and edited
OccModel9a.init <-read.modeldef("DEFINITIONS/OccModel9aInit.def",model.type="logistic",
                                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                                oldGlimClim.warning=FALSE)

print(OccModel9a.init)

OccModel9a <- GLCfit(OccModel9a.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1,verbosity=1,
                     cor.file="OccModel9a_Corr.dat")

print(OccModel9a)

summary(OccModel9a, tables=NULL)


#OccModel9b: transformation of temperature #10
write.modeldef(OccModel9a,file="OccModel9bInit.def") # Copied to DEFINITIONS/ and edited
OccModel9b.init <-
  read.modeldef("DEFINITIONS/OccModel9bInit.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)

print(OccModel9b.init)

OccModel9b <- GLCfit(OccModel9b.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1,verbosity=1,
                     cor.file="OccModel9b_Corr.dat")

print(OccModel9b)

summary(OccModel9b, tables=NULL)


#OccModel9c: transformation of temperature #20
write.modeldef(OccModel9b,file="OccModel9cInit.def") # Copied to DEFINITIONS/ and edited
OccModel9c.init <-
  read.modeldef("DEFINITIONS/OccModel9cInit.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)

print(OccModel9c.init)

OccModel9c <- GLCfit(OccModel9c.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1,nprev.required=-1,verbosity=1,
                     cor.file="OccModel9c_Corr.dat")

print(OccModel9c)

summary(OccModel9c, tables=NULL)
# Decay rate and std err too large

print(OccModel9b)

#Final baseline occurrence model conditioned on temperature
#Delete interaction 2 and 3,daoshu 5 and 6, 6 and 7
write.modeldef(OccModel9b,file="OccModel10Init.def") # Copied to DEFINITIONS/ and edited
OccModel10.init <-
  read.modeldef("DEFINITIONS/OccModel10Init.def",model.type="logistic",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)
print(OccModel10.init)
OccModel10 <- GLCfit(OccModel10.init,which.response=1,siteinfo=siteinfo,
                    data.file="PrecipTemp_Processed.dat",
                    diagnostics=1, nprev.required=-1, verbosity=1,
                    cor.file="OccModel10_Corr.dat")

anova(OccModel9b, OccModel10)

print(OccModel10)

summary(OccModel10, tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(OccModel10,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(OccModel10,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

plot(OccModel10,which.plots=5,plot.cols=c("blue","purple"), distance.units="degrees")


###########################################################################
#Precipitation intensity
#IntModel8a:incorporating temperature as a covariate
write.modeldef(IntModel7,file="IntModel8aInit.def") # Copied to DEFINITIONS/ and edited
IntModel8a.init <-
  read.modeldef("DEFINITIONS/IntModel8aInit.def", model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)
print(IntModel8a.init)

IntModel8a <- GLCfit(IntModel8a.init,which.response=1,siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1, nprev.required=-1,
                      cor.file="IntModel8a_Corr.dat")

print(IntModel8a)

summary(IntModel8a, tables=NULL)


#IntModel8b:#10
write.modeldef(IntModel7,file="IntModel8bInit.def") # Copied to DEFINITIONS/ and edited
IntModel8b.init <-
  read.modeldef("DEFINITIONS/IntModel8bInit.def", model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)
print(IntModel8b.init)

IntModel8b <- GLCfit(IntModel8b.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1,
                     cor.file="IntModel8b_Corr.dat")

print(IntModel8b)

summary(IntModel8b, tables=NULL)


#IntModel8c:#20
write.modeldef(IntModel8b,file="IntModel8cInit.def") # Copied to DEFINITIONS/ and edited
IntModel8c.init <-
  read.modeldef("DEFINITIONS/IntModel8cInit.def", model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)
print(IntModel8c.init)

IntModel8c <- GLCfit(IntModel8c.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=1, nprev.required=-1,
                     cor.file="IntModel8c_Corr.dat")

print(IntModel8c)

summary(IntModel8c, tables=NULL)
#standard error too large


print(IntModel8b)
#IntModel9:Final
#TRY Delete interaction2 and 3 daoshu 1and 2
write.modeldef(IntModel8b,file="IntModel9Init.def") # Copied to DEFINITIONS/ and edited
IntModel9.init <-
  read.modeldef("DEFINITIONS/IntModel9Init.def", model.type="gamma",
                siteinfo=siteinfo, var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning=FALSE)
print(IntModel9.init)

IntModel9 <- GLCfit(IntModel9.init,which.response=1,siteinfo=siteinfo,
                     data.file="PrecipTemp_Processed.dat",
                     diagnostics=2, nprev.required=-1,
                     cor.file="IntModel9_Corr.dat",
                     resid.file="IntModel9_Resids.dat")

print(IntModel9)

summary(IntModel9, tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=10,height=6)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(IntModel9,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(IntModel9,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(IntModel9,which.plots=4:5,plot.cols=c("blue","purple"), distance.units="degrees")
#The diagnostics improves slightly and the predictive performance improves