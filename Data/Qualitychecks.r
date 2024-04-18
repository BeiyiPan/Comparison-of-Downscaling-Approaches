#######################################################################
#######################################################################
#######################################################################
######                                                           ######
######  Checking the quality of the precipitation data used in   ######
######  the study. The reason for doing this is that modelling   ######
######  results can be particularly sensitive to inhomogeneities ######
######  in precipitation data, which are themselves quite        ######
######  common.                                                  ######
######                                                           ######
######  Run makedata.r first, to generate the file PrecipTemp.dat######
######                                                           ######
#######################################################################
#######################################################################
#######################################################################
library(lattice); library(RColorBrewer); library(Rglimclim)
#######################################################################
#
#	Read data from file
#
#######################################################################
all.data <- read.GLCdata("PrecipTemp.dat")
names(all.data)[5:6] <- c("Rain", "Temp")
#######################################################################
#
#	First check: plot annual proportions of wet days at each site, for
#	different wet day thresholds. It's also worth looking at outlying
# residuals from a two-way ANOVA, try and locate any anomalies in
# these plots. 
#
#######################################################################
pwet <- function(x,tau=0) { mean(x > tau, na.rm=TRUE) }
# colscl <- grey((rev(0:100)/100)^(0.5))
colscl <- rainbow(100,start=0,end=5/6)
if (dev.cur() == 1) x11(width=11,height=8)
tau <- c(0,0.5,1,2,5)
for (i in 1:length(tau)) {
 annual.pwet <- tapply(all.data$Rain,INDEX=list(all.data$Year,all.data$Site),
                       FUN=pwet,tau=tau[i])
 plot(levelplot(annual.pwet, 
           xlab="Year",ylab="Station",aspect="fill",
           main=paste("Proportions of days with rain > ",tau[i],"mm",sep=""),
           scales=list(cex=0.5,x=list(rot=90)),col.regions=colscl))
 flnmroot <- paste("annual_pwet_",tau[i],"mm",sep="")
 dev.copy(pdf,paste(flnmroot,".pdf",sep=""),width=11,height=8)
 dev.off()
 annual.pwet <- data.frame(Pwet=as.numeric(annual.pwet),
                           Year=rownames(annual.pwet)[row(annual.pwet)],
                           Site=colnames(annual.pwet)[col(annual.pwet)])
 cat(paste("\nAnomalous proportions of wet sites: threshold ", 
           tau[i],"mm:\n", sep=""))
 e <- resid(lm(Pwet ~ Year + Site, data=annual.pwet, na.action=na.exclude))
 annual.pwet$Resid <- e
 print(annual.pwet[as.numeric(names(boxplot(e, plot=FALSE)$out)),])
}
#######################################################################
#
#	At zero threshold, visual inspection identifies the following stations
# that consistently have unrealistically high proportions of wet days
# (very close to 1 in some instances): S046, S048, S058, S073, S084,
# S086, S088, S091, S101. These are completely unrealistic and will 
# be removed. 
#
#######################################################################
wanted.rows <- !(all.data$Site %in% c("S046", "S048", "S058", "S073", 
                                      "S084", "S086", "S088", "S091", "S101"))
all.data <- all.data[wanted.rows,]
#
#   Now looking at the ANOVA residuals together with the plots, the
#   following additional queries are flagged at zero threshold (the
#   list below excludes cases that are probably associated with 
#   incomplete years at the start or end of a station's record, 
#   as well as instances where there's just one slightly anomalous
#   year within a station's record):
#
# - S028, 1962 and several years during the 1970s: much wetter than 
#   expected. This needs further investigation. 
# - S037: some very dry years in the 1960s, then unusually wet in 1994.
#   Further investigation needed. 
# - S049: seems unusually wet from 1966 to 1985. 
# - S051: a handful of unusual values, no obvious pattern.
# - S055: some erratic values in the 1990s before a break in the record,
#   which possibly suggests problems with the operation of the station
#   during this period. 
# - S060: slightly wetter than expected in 1967 and 2015. 
# - S061: looks a bit low until about 1964. 
# - S065: seems a bit dry until mid-1960s, and possibly a bit wet in
#   1989 and 1994. 
# - S074: unusually wet in 2003 when other sites indicate a dry year. 
#   It's a short record as well: can probably remove without losing
#   much information. 
# - S075: relatively dry in 1971 and 1974. 
# - S079: looks too dry during the early 1980s and too wet from 
#   2016-2019. 
# - S087: looks a bit too wet from about 2000-2005. 
# - S092: a bit too wet in 1959 and 1960. 
# - S097: unusually dry in 1971 and 1976. 
# - S108: suggestion of being wetter than expected until mid-1970s. 
# - S113: looks rather dry until mid-1970s. 
#
#	This overall pattern persists at higher thresholds, with 
# a few additional sites showing anomalies; and S092 
# shows some very odd behaviour after 2000 with a 5mm threshold. 
# The additional sites are: 
#
# - S041: seems particularly dry in 2000-2005 (and it's very dry overall)
# - S072: seems unusually wet in 1968 and 1984 (but it's a very dry site
#   overall)
# - S085: has relatively few values above 5mm in early 1970s. 
# - S102 seems a bit wet in 1967-1968, and unusually dry in 1994. 
#
# Let's now look at the frequency of observations per year for the
# "potentially suspect" sites identified above. 
#
#######################################################################
wanted.rows <- (all.data$Site %in% c("S028","S037","S041","S049","S051",
                                     "S055","S060","S061","S065","S072",
                                     "S074","S075","S079","S085","S087",
                                     "S092","S097","S102","S108","S113"))
dev.off(); x11(width=11,height=6)
Nobs.table <- table(all.data[wanted.rows,c("Year","Site")])
Nobs.table[Nobs.table==0] <- NA
plot(levelplot(Nobs.table,xlab="Year",ylab="Station",aspect="fill",
               main="Annual numbers of non-missing days",
               scales=list(cex=0.7,x=list(rot=90)),col.regions=colscl))
dev.copy(pdf,"AnnualNObs_bySite.pdf",width=11,height=6)
dev.off()
#
#   Not much there in the way of incomplete records. The only 
#   noteworthy points are:
#
#   - S041: several missing observations in 2003 and 2004, suggesting
#     that there may have ben some operational difficulties during
#     this period (although this wasn't the period that was flagged
#     from the previous analysis). DECISION: remove data for these
#     two years from this station's record. 
#   - S055: missing values in both years (1993 and 1994) prior to
#     a break in the record starting in 1995. Again, perhaps indicates
#     operational difficulties. DECISION: remove data for 1993, 1994
#     and 1995 for this station. 
#   - S060: some missing values in 1991 and 1992, again not the
#     period flagged in the previous analysis but treat in the
#     same way as for S041. 
#   - S072: looks like a couple of missing months in 1990, probably
#     nothing to worry about. 
#   - S074: short record is a bit patchy, particularly towards the
#     end. DECISION: discard this site.
#   - S079: missing observations in 2017 and in 2021, suggesting
#     potential operational difficulties. DECISION: discard data
#     from this station from 2016 onwards (NB this period was
#     flagged previously for this station).
# 
#   There's nothing particularly unusual at any of the other
#   stations: occasional years with missing data, nothing more. 
#
all.data <- 
  all.data[!(all.data$Site=="S041" & all.data$Year %in% 2003:2004),]
all.data <- 
  all.data[!(all.data$Site=="S055" & all.data$Year %in% 1993:1995),]
all.data <- 
  all.data[!(all.data$Site=="S060" & all.data$Year %in% 1991:1992),]
all.data <- all.data[!(all.data$Site=="S074"),]
all.data <- all.data[!(all.data$Site=="S079" & all.data$Year >= 2016),]
#######################################################################
#
# Next check is for differences in recording resolution between
# sites. 
#
#######################################################################
decimal <- round(10*(all.data$Rain %% 1),1)
decimal.tab <- t(table(decimal,all.data$Site))
decimal.tab <- decimal.tab / rowSums(decimal.tab)
plot(levelplot(decimal.tab[,-1], 
           xlab="Station",ylab="Decimal",aspect="fill",
           main="Frequencies of first decimals at each site (zeroes excluded)",
           scales=list(cex=0.75,x=list(rot=90)),col.regions=colscl))
dev.copy(pdf,"Decimals_bySite.pdf",width=11,height=6)
dev.off()
#######################################################################
#
#	Most seem fairly homogeneous EXCEPT: 
# 
# - Stations S037 and S054 which sem to have wildly different
#   recording reesolutions to anything else. DECISION: get
#   rid of these sites since they may not be comparable with 
#   the others (S037 was picked up before).
# - Stations S028, S034, S045, S049, S063, S068, S079, S081, 
#   S085, S087, S099, S105 and S108: these tend to have a
#   concentration of values ending x.1 to x.5 (except S079
#   where the distribution is almost uniform). Some of these
#   may relate to the recording of small values - particularly
#   for some of the stations that have already featured in 
#   earlier checks. 
#
#######################################################################
all.data <- all.data[!(all.data$Site %in% c("S037", "S054")),]
#######################################################################
# 
# It's probably also worth checking whether the recording 
# resolution has changed over time ... 
#
#######################################################################
decimal <- round(10*(all.data$Rain %% 1),1)
decimal.tab <- t(table(decimal,all.data$Year))
decimal.tab <- decimal.tab / rowSums(decimal.tab)
n.years <- nrow(decimal.tab)
year.numbers <- as.numeric(row.names(decimal.tab))
which.years <- (year.numbers / 5 == floor(year.numbers/5))
xlabs <- rep("",n.years)
xlabs[which.years] <- as.character(year.numbers[which.years])
plot(levelplot(decimal.tab[,-1], xlab="Year",ylab="Decimal",aspect="fill",
               main=list(paste("Relative frequencies of first decimal digits in",
                               "each year\n(zeroes excluded)"),cex=1.3),
               scales=list(x=list(cex=1,rot=90,labels=xlabs),
                           y=list(cex=1)),col.regions=colscl))
dev.copy(pdf,"Decimals_byYear.pdf",width=11,height=6)
dev.off()
#######################################################################
#
#	That's weird! There seem to be some clear but unsystematic differences
# between years, in particular with an increasing frequency of x.1
# records in the last decade or so. To try and remove the effects
# of any such inhomogeneities, I'm going to round all the precipitation
# values to the nearest 0.5mm (this is done later).
#
#	Next check: tabulate proportion of wet days for each day of the 
#	week at each site.
#
#######################################################################
dayofweek <- factor(weekdays(ISOdate(all.data$Year,all.data$Month,
                                     all.data$Day)),
                    levels=c("Monday","Tuesday","Wednesday",
                             "Thursday","Friday","Saturday","Sunday"))
for (i in 1:length(tau)) {
 weekday.pwet <- tapply(all.data$Rain,INDEX=list(all.data$Site,dayofweek),
                       FUN=pwet,tau=tau[i])
 plot(levelplot(weekday.pwet, 
           xlab="Station",ylab="Day of week",aspect="fill",
           main=paste("Proportions of days with rain > ",tau[i],"mm",sep=""),
           scales=list(x=list(cex=0.75,rot=90)),col.regions=colscl))
 flnmroot <- paste("weekday_pwet_",tau[i],"mm",sep="")
 dev.copy(pdf,paste(flnmroot,".pdf",sep=""),width=11,height=8)
 dev.off()
}
#######################################################################
#
# The main thing that leaps out here is that some sites are substantially 
# drier than the others - no obvious weekend effects anywhere. 
#
# It is nonetheless possible that dry days are recorded 
# preferentially. To look at this, plot pwet against proportion of 
# observations 
#
#######################################################################
prop.obs <- aggregate(all.data$Rain,by=list(all.data$Year,all.data$Site),
                   FUN=length)
names(prop.obs) <- c("Year","Site","Proportion")
prop.obs[,3] <- prop.obs[,3] / 366
annual.pwet <- aggregate(all.data$Rain,by=list(all.data$Year,all.data$Site),
                       FUN=pwet,tau=0)
names(annual.pwet) <- c("Year","Site","Pwet")
plot(xyplot(annual.pwet$Pwet ~ prop.obs$Proportion | annual.pwet$Site,pch=16,
       xlab="Proportion of observations",ylab="Proportion of wet days",
       scales=list(cex=0.5,tck=0.5),as.table=TRUE,
       main="Annual proportions of wet days vs proportions of available observations",
       par.strip.text=list(cex=0.7,font=2)))
dev.copy(jpeg,"PwetvsPobs_byYear.jpg",width=11*72,height=8*72,quality=99)
dev.off()
prop.obs <- aggregate(prop.obs$Proportion,by=list(prop.obs$Site), FUN=mean)
names(prop.obs) <- c("Site","Proportion")
site.pwet <- aggregate(annual.pwet$Pwet,by=list(annual.pwet$Site),FUN=mean)
names(site.pwet) <- c("Site","Pwet")
plot(prop.obs$Proportion,site.pwet$Pwet,type="n",
     xlab="Proportion of observations",
     ylab="Proportion of wet days",
     main="Proportion of wet days vs proportion of available observations, by site")
text(prop.obs$Proportion,site.pwet$Pwet,labels=site.pwet$Site,cex=0.8,font=2)
cor.res <- cor.test(prop.obs$Proportion,site.pwet$Pwet,method="kendall")
legend("topleft",bty="n",
       legend=c(paste("Kendall tau:",round(cor.res$estimate,2)),
                paste("p:",round(cor.res$p.value,3))))
box(lwd=2)
dev.copy(jpeg,"PwetvsPobs_bySite.jpg",width=11*72,height=8*72,quality=99)
dev.off()
#######################################################################
#
#	All looks fine: there's no relationship (unsurprising perhaps, 
# given that there are few missing values in any station's record)
#
#	Putting that lot together: there's nothing obviously wrong with
# the precipitation data, except for apparent differences in 
# recording resolution between different sites. To reduce the, 
# effect of this in the modelling, round all precipitation
# values to the nearest 0.5mm and take 1mm as the threshold for 
# defining a "wet" day. 
#
#	Here is the revised gaugvals.dat, with all values rounded to
#	0.5mm
#
#######################################################################
all.data$Rain <- round(2*all.data$Rain)/2
cat("Writing data files ...\n")
write.GLCdata(all.data, file="PrecipTemp_Processed.dat")
cat("Done.\n")
