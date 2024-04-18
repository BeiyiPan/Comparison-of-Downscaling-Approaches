library(Rglimclim)
library(maps)
library(lattice)

#Temperature model 7a: incorporating wet / dry precipitation indicator and interactions
#3
write.modeldef(TempModel6, file=c("TempModel7aInit.def",
                                  "TempModel7aInit_Dispersion.def")) # Copied to DEFINITIONS/ and edite

TempModel7a.init <-
  read.modeldef("DEFINITIONS/TempModel7aInit.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7a.init)

TempModel7aDisp.init <-
  read.modeldef("DEFINITIONS/TempModel7aInit_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7aDisp.init)

TempModel7a <- GLCfit(TempModel7a.init, dispersion.def=TempModel7aDisp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=-1,
                      cor.file="TempModel7a_Corr.dat")

print(TempModel7a)

summary(TempModel7a,tables=NULL)

####################################################################
#Different transformations
#Temperature model 7b: #13 (#2 worse
write.modeldef(TempModel7a, file=c("TempModel7bInit.def",
                                  "TempModel7bInit_Dispersion.def")) # Copied to DEFINITIONS/ and edite

TempModel7b.init <-
  read.modeldef("DEFINITIONS/TempModel7bInit.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7b.init)

TempModel7bDisp.init <-
  read.modeldef("DEFINITIONS/TempModel7bInit_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7bDisp.init)

TempModel7b <- GLCfit(TempModel7b.init, dispersion.def=TempModel7bDisp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=-1,
                      cor.file="TempModel7b_Corr.dat")

print(TempModel7b)

summary(TempModel7b,tables=NULL)


#Temperature model 7c: #23
write.modeldef(TempModel7b, file=c("TempModel7cInit.def",
                                   "TempModel7cInit_Dispersion.def")) # Copied to DEFINITIONS/ and edite

TempModel7c.init <-
  read.modeldef("DEFINITIONS/TempModel7cInit.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7c.init)

TempModel7cDisp.init <-
  read.modeldef("DEFINITIONS/TempModel7cInit_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7cDisp.init)

TempModel7c <- GLCfit(TempModel7c.init, dispersion.def=TempModel7cDisp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=-1,
                      cor.file="TempModel7c_Corr.dat")

print(TempModel7c)

summary(TempModel7c,tables=NULL)
# Improves very slightly, and Std Err:1.1077

#Temperature model 7d:#113(worse)
write.modeldef(TempModel7b, file=c("TempModel7dInit.def",
                                   "TempModel7dInit_Dispersion.def")) # Copied to DEFINITIONS/ and edite

TempModel7d.init <-
  read.modeldef("DEFINITIONS/TempModel7dInit.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7d.init)

TempModel7dDisp.init <-
  read.modeldef("DEFINITIONS/TempModel7dInit_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel7dDisp.init)

TempModel7d <- GLCfit(TempModel7d.init, dispersion.def=TempModel7dDisp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=1,nprev.required=-1,
                      cor.file="TempModel7d_Corr.dat")

print(TempModel7d)

summary(TempModel7d,tables=NULL)


print(TempModel7b)
#Temperature model 8:
write.modeldef(TempModel7b, file=c("TempModel8Init.def",
                                   "TempModel8Init_Dispersion.def")) # Copied to DEFINITIONS/ and edite

TempModel8.init <-
  read.modeldef("DEFINITIONS/TempModel8Init.def",
                model.type="normal-heteroscedastic", which.part="mean",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel8.init)

TempModel8Disp.init <-
  read.modeldef("DEFINITIONS/TempModel8Init_Dispersion.def",
                model.type="normal-heteroscedastic", which.part="dispersion",
                siteinfo=siteinfo, which.response=2,
                var.names=c("Precipitation","Temperature"),
                oldGlimClim.warning = FALSE)
print(TempModel8Disp.init)

TempModel8 <- GLCfit(TempModel8.init, dispersion.def=TempModel8Disp.init,
                      model.type="normal-heteroscedastic", which.response=2,
                      siteinfo=siteinfo,
                      data.file="PrecipTemp_Processed.dat",
                      diagnostics=2,nprev.required=-1,
                      cor.file="TempModel8_Corr.dat",
                      resid.file="TempModel8_Resids.dat")

print(TempModel8)

summary(TempModel8,tables=NULL)


if (!(names(dev.cur()) %in% c("windows","X11cairo"))) x11(width=7,height=12)
par(mfrow=c(2,2),mar=c(3,3,3,2),mgp=c(2,0.75,0),lwd=2,ask=TRUE)
plot(TempModel8,which.plots=1:2,lwd=2)

par(mfrow=c(1,1),mar=c(5,3,5,3))
image(topo.longs,topo.lats,topo.data,col=terrain.colors(100),
      xlab="",ylab="",xlim=LongLims,ylim=LatLims)
plot(TempModel8,which.plots=3,
     site.options=list(add.to.map=TRUE,site.labels="significant",scale=1))

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(TempModel8,which.plots=4:5,plot.cols=c("blue","purple"))