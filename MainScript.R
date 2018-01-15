##############################################################################################
# Script for "Extreme mid-winter drought weakens tree hydraulic and carbohydrate systems"
# Earles et al.
# last updated: Oct. 27th, 2017
###############################################################################################

######################
# Load data and libs #
######################

#Load libraries
library('reshape');library('lme4');library('ggplot2');library('pbkrtest');library('zoo');library('MuMIn');library('merTools');library('gridExtra'); library('reshape'); library('grid'); library('dplR')

#Read VPD and Preip table generation script
source('subScripts/VPD and Precip Anomaly Table Generation.R')

#read relative water content (RWC) and weather data tables
d.RWC <- melt(read.table('Script Input Data/RWC.txt',header=T))
d.weather <- read.table('Script Input Data/weather.txt',header=T)
d.VPDPrecipAnom <- read.table('Script Input Data/VPD and Precip Anomaly.txt')
d.weather <- cbind(d.weather,d.VPDPrecipAnom[,-1])
d.weather <- d.weather[d.weather$day!=24,]

#Calculate precipitation deviation from 22 year normal
P30.norm <- mean(d.precip[d.precip$Year %in% 1993:2014,]$precip.ma30, na.rm=TRUE)
P45.norm <- mean(d.precip[d.precip$Year %in% 1993:2014,]$precip.ma45, na.rm=TRUE)
P60.norm <- mean(d.precip[d.precip$Year %in% 1993:2014,]$precip.ma60, na.rm=TRUE)
d.weather$P30.sml <- d.weather$P30*0.01
d.weather$P45.sml <- d.weather$P45*0.01
d.weather$P60.sml <- d.weather$P60*0.01

#Add day variable to RWC table
d.RWC$day <- as.numeric(substr(d.RWC$variable,2,4))
d.RWC <- d.RWC[d.RWC$day!=24,]

#name columns
names(d.RWC) <- c("ID","","RWC","day")

#merge RWC and weather data
d.all <- merge(d.RWC,d.weather,by=c("day"),all.x=T,all.y=T)

#change day to factor
d.all$day <- as.factor(d.all$day)

#add species variable
d.all$species <- substr(d.all$ID,3,4)

#add plot variable
d.all$plot <- as.factor(substr(d.all$ID,1,1))

#add variable for July 2014 RWC
d.all$RWC.July14 <- d.all[d.all$day == 205,]$RWC

#Add scaled precip variable
d.all$P60.scale <- scale(d.all$P60)

#For dredging: Remove columns that have VPD45 and VPD60 to avoid incomplete cases
d.all <- d.all[,-c(11,16,17,22,23,28,29,34,35,40,41)]

#remove incomplete cases
d.all.complete <- d.all[complete.cases(d.all),] #remove incomplete cases

#################################
#Fig 1: Plot precip and VPD data#
#################################

png(file="Script Output/Fig1_VPDandPrecipMonthly.png", units='px', width=1000, height=700, pointsize=12, res=150, bg='white')
par(mfrow=c(2,1), mar=c(1.5,1,1,1), oma=c(1,2,0,0))
boxplot(vpd.norecent$VPD~vpd.norecent$Month,ylim=c(0,3),xlim=c(0,12),frame.plot=FALSE,axes=FALSE,xlab="")
axis(1,at=1:12,labels=FALSE) 
axis(2,at=seq(0,3,by=1.5),pos=0.25, ylab="VPD (kPa)")
par(new=T)
plot(vpd.recent1$Month,vpd.recent1$VPD,col="turquoise",bg="turquoise",ylim=c(0,2.5),xlim=c(0,12),pch=21,cex=1.5,ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
par(new=T)
plot(vpd.recent2$Month,vpd.recent2$VPD,col="red",bg="red",ylim=c(0,2.5),xlim=c(0,12),pch=21,cex=1.5,ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
legend(1,2.5,c('1993-2013','2014','2015'), pch=c(0,21,21), col=c('black','turquoise','red'), pt.bg=c('black','turquoise','red'), bty='n', y.intersp=1.4)
mtext("VPD (kPa)", side=2, line=1)

boxplot(precip.norecent$Precip~precip.norecent$Month,ylim=c(0,1000),xlim=c(0,12),frame.plot=FALSE,axes=FALSE,ylab='',xlab="")
axis(1,at=1:12,labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) 
axis(2,at=seq(0,1000,by=500), pos=0.25)
par(new=T)
plot(precip.recent1$Month,precip.recent1$Precip,col="turquoise",bg="turquoise",ylim=c(0,1000),xlim=c(0,12),pch=21,cex=1.5,ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
par(new=T)
plot(precip.recent2$Month,precip.recent2$Precip,col="red",bg="red",ylim=c(0,1000),xlim=c(0,12),pch=21,cex=1.5,ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
mtext('precip (mm)', side=2, line=1)
dev.off()

##########################################
## Read RWC and NSCs time series script ##
## Generate Figs. 2 and S1-S3           ##
##########################################
source('subScripts/RWCandNSC_TimeSeries_Analysis.R')

################################################################################
##Statistical Model to explore relationship between RWC and climatic variables##
################################################################################

# Option 1: Run exhaustive model search to identify best climatic variables to predict RWC
# Run script separately to systematically examine search process and output
# options(na.action = "na.fail")
# source('subScripts/ModelSearchFunction.R') #This takes ~5 mins to run
# options(na.action = "na.omit")
 
# Option 2: Assume best model is described below
best.mod <- lmer(RWC ~  P60.sml + VPD.zscore5 + (P60.sml | ID) + (1 | species/plot), data=d.all.complete) #This model removes correlation of the random effects, minimizes AIC and results in most significant fixed effect precipitation response

#Summarize best model
summary(best.mod)

#Get p-values of best model
GetME_PVals=function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs #Everything is significantly different from UB
  return(coefs)
}
GetME_PVals(best.mod)

#Get R^2 equivalent for best model; see reference in SI text
1-var(residuals(best.mod))/(var(model.response(model.frame(best.mod)))) #from Ronghui Xu et al.


##############################################
#Fig 3: Observed vs Predicted RWC by species#
##############################################

#Plot observed and predicted RWC for all species combined over time
d.all.pred <- data.frame(d.all.complete$day,d.all.complete$ID,d.all.complete$RWC,predict(best.mod))
names(d.all.pred) <- c("day","ID","RWC","RWCpred")
d.all.pred$species <- substr(d.all.pred$ID,3,4)
RWCobs.plot <- ddply(d.all.pred, .(day,species), RWCmed = median(RWC,na.rm=TRUE), RWC5th = quantile(RWC,probs=0.025,na.rm=TRUE), RWC25th = quantile(RWC,probs=0.25,na.rm=TRUE),RWC75th = quantile(RWC,probs=0.75,na.rm=TRUE),RWC95th = quantile(RWC,probs=0.975,na.rm=TRUE),summarize)
RWCobs.plot$OP <- "observed"
RWCpred.plot <- ddply(d.all.pred, .(day,species), RWCmed = median(RWCpred,na.rm=TRUE), RWC5th = quantile(RWCpred,probs=0.025,na.rm=TRUE), RWC25th = quantile(RWCpred,probs=0.25,na.rm=TRUE),RWC75th = quantile(RWCpred,probs=0.75,na.rm=TRUE),RWC95th = quantile(RWCpred,probs=0.975,na.rm=TRUE),summarize)
RWCpred.plot$OP <- "predicted"
RWC.plot <- rbind(RWCobs.plot,RWCpred.plot)
RWCplot <- melt(RWC.plot)
names(RWCplot) <- c("day","species","OP","RWC.quant","RWC")
RWCplot$day <- as.numeric(as.character(RWCplot$day))
RWCplot$species <- as.factor(RWCplot$species)
RWCplot$group.unq <- as.factor(paste(RWCplot$RWC.quant,RWCplot$species,RWCplot$OP))
species_names <- c('DF' = "P. menziesii",'PP' = "P. ponderosa",'SP' = "P. lambertiana") #change labels for plot
p1 <- ggplot(RWCplot) + geom_line(aes(x=day,y=RWC,colour=species,group=RWC.quant,linetype=RWC.quant)) + facet_grid(species~OP,labeller = labeller(species = species_names)) + ylab(expression(paste('RWC [',m^3,' ', m^-3,']'))) + xlab('days since Jan. 1, 2014') + theme_bw() + scale_colour_manual(values=c('cadetblue','grey50','coral2')) + theme(strip.text.y = element_text(face = "italic", size=11), strip.text.x = element_text(size = 11), axis.text = element_text(size=10), strip.background = element_blank(), legend.title = element_text(size=11),legend.text = element_text(size=11), legend.position = 'none') + annotate("rect",xmin=380, xmax=460, ymin=0.2, ymax=1, alpha=0.3) + annotate("segment",x=365, xend=365, y=0.2, yend=1.0, ymin=0.2, ymax=1, linetype='dashed')

# Save plot output
png(file="Script Output/Fig5_RWC_PredObs.png", units='in', width=4, height=5, pointsize=12, res=300, bg='white')
p1
dev.off()

##########################
#Fig 4: RWC vs. Kh model #
##########################

# Read RWC vs hydraulic conductivity (Kh) calibration script
source('subScripts/RWC_Kh_calibration.R')

############################################################################
#Fig 5a: Historic/recent climate data over hydraulic conductivity heat map #
############################################################################

#Plot historic climate data colored by PLC
VPD5.zscore.sim <- seq(-3,3,by=0.1)
P60.sml.sim <- seq(-5,14,by=0.1)
Sim.grid <- expand.grid(P60.sml.sim,VPD5.zscore.sim) 
names(Sim.grid) <- c("P60.sml.sim","VPD5.zscore.sim")
Sim.grid$RWC.sim.DF <- fixef(best.mod)[1]+Sim.grid[,1]*fixef(best.mod)[2]+Sim.grid[,2]*fixef(best.mod)[3]+ ranef(best.mod)$species[1,1]
Sim.grid$Kh.est.DF <- predict(fit.DF, newdata=list(RWCXYL = Sim.grid$RWC.sim.DF))
Sim.grid$RWC.sim.PP <- fixef(best.mod)[1]+Sim.grid[,1]*fixef(best.mod)[2]+Sim.grid[,2]*fixef(best.mod)[3] + ranef(best.mod)$species[2,1]
Sim.grid$Kh.est.PP <- predict(fit.PP, newdata=list(RWCXYL = Sim.grid$RWC.sim.PP))
Sim.grid$RWC.sim.SP <- fixef(best.mod)[1]+Sim.grid[,1]*fixef(best.mod)[2]+Sim.grid[,2]*fixef(best.mod)[3] + ranef(best.mod)$species[3,1] 
Sim.grid$Kh.est.SP <- predict(fit.SP, newdata=list(RWCXYL = Sim.grid$RWC.sim.SP))
Sim.grid$RWC.sim <- rowMeans(cbind(Sim.grid$RWC.sim.DF,Sim.grid$RWC.sim.PP,Sim.grid$RWC.sim.SP))
Sim.grid$Kh.est <- rowMeans(cbind(Sim.grid$Kh.est.DF,Sim.grid$Kh.est.PP,Sim.grid$Kh.est.SP))


#Calculate historic z-scores by month for VPD5.zscore and P60.sml
d.clim.hist <- data.frame(d.precip,d.vpd)[,-c(9,10,11)]
d.clim.hist$MoDay <- paste(d.clim.hist$Month,"_",d.clim.hist$Day)
d.clim.hist.wn <- d.clim.hist[d.clim.hist$Month %in% c(1,2,3,11,12) & d.clim.hist$Year %in% 1992:2014,]
d.clim.hist.su <- d.clim.hist[d.clim.hist$Month %in% 6:9 & d.clim.hist$Year %in% 1992:2014,]
d.clim.recent.wn <- d.clim.hist[d.clim.hist$Month %in% c(1,2,3) & d.clim.hist$Year %in% 2015,]
d.clim.recent.su <- d.clim.hist[d.clim.hist$Month %in% 6:9 & d.clim.hist$Year %in% 2014:2015,]
clim.mo.musd <- ddply(d.clim.hist,.(Year,Month),P60.mean = mean(precip.ma60,na.rm=TRUE),VPD5.mean = mean(VPD.ma5,na.rm=TRUE),P60.sd = sd(precip.ma60,na.rm=TRUE),VPD5.sd = sd(VPD.ma5,na.rm=TRUE),summarize)
clim.yr.musd <- ddply(d.clim.hist,.(Month,Day),P60.mean = mean(precip.ma60,na.rm=TRUE),VPD5.mean = mean(VPD.ma5,na.rm=TRUE),P60.sd = sd(precip.ma60,na.rm=TRUE),VPD5.sd = sd(VPD.ma5,na.rm=TRUE),summarize)
clim.yr.musd$MoDay <- paste(clim.yr.musd$Month,"_",clim.yr.musd$Day)
clim.musd <- merge(d.clim.hist,clim.yr.musd,by="MoDay")
clim.musd <- clim.musd[order(clim.musd$Year,clim.musd$Month.x,clim.musd$Day.x),]
clim.musd$P60.sml <- clim.musd$precip.ma60*0.01
clim.musd$VPD5.hist.z <- (clim.musd$VPD.ma5-clim.musd$VPD5.mean)/clim.musd$VPD5.sd
clim.musd <- clim.musd[,c(2,3,4,22,23)]
names(clim.musd) <- c("Month","Day","Year","P60.sml","VPD5.hist.z")
clim.musd.wn.recent <- clim.musd[clim.musd$Month %in% c(1,2,3) & clim.musd$Year %in% 2015,]
clim.musd.su.recent <- clim.musd[clim.musd$Month %in% 6:9 & clim.musd$Year %in% 2014,]
clim.musd.wn.historic <- clim.musd[clim.musd$Month %in% c(1,2,3) & clim.musd$Year %in% 1992:2014,]
clim.musd.su.historic <- clim.musd[clim.musd$Month %in% 6:9 & clim.musd$Year %in% 1992:2014,]

#Calculate historic monthly averages of VPD5.zscore and P60.sml
clim.musd.moavg <- ddply(clim.musd,.(Year,Month), VPD5.hist.z = mean(VPD5.hist.z, na.rm=TRUE), P60.sml = mean(P60.sml, na.rm=TRUE),summarize)
clim.musd.moavg.wn.recent <- clim.musd.moavg[clim.musd.moavg$Month %in% c(1,2,3) & clim.musd.moavg$Year %in% 2015,]
clim.musd.moavg.su.recent <- clim.musd.moavg[clim.musd.moavg$Month %in% 6:9 & clim.musd.moavg$Year %in% 2014,]
clim.musd.moavg.wn.historic <- clim.musd.moavg[clim.musd.moavg$Month %in% c(1,2,3) & clim.musd.moavg$Year %in% 1993:2014,]
clim.musd.moavg.su.historic <- clim.musd.moavg[clim.musd.moavg$Month %in% 6:9 & clim.musd.moavg$Year %in% 1993:2014,]

#Plot monthly average P60 vs. VPD5.zscore with historic Kh values overlaid as points
p0 <- ggplot() + geom_tile(data=Sim.grid,aes(x=VPD5.zscore.sim,y=P60.sml.sim*100,z=Sim.grid$Kh.est,fill=Sim.grid$Kh.est)) + 
  scale_fill_gradientn(colours=c("blue","yellow","red"), limit=c(0.4,2.8), name = expression(K[h])) + theme_bw() + 
  geom_point(data=clim.musd.moavg.wn.historic,aes(x=VPD5.hist.z,y=P60.sml*100),color="white", size=2) + 
  geom_point(data=clim.musd.moavg.su.historic,aes(x=VPD5.hist.z,y=P60.sml*100),color="black", size=2) + 
  geom_point(data=clim.musd.moavg.wn.recent[1,],aes(x=VPD5.hist.z,y=P60.sml*100),color="turquoise", size=2.1) + 
  geom_point(data=clim.musd.moavg.wn.recent[2,],aes(x=VPD5.hist.z,y=P60.sml*100),color="turquoise", size=2.1) + 
  geom_point(data=clim.musd.moavg.wn.recent[3,],aes(x=VPD5.hist.z,y=P60.sml*100),color="turquoise", size=2.1) + 
  annotate('text', x = clim.musd.moavg.wn.recent[1:3,3]*c(1.4,1.35,1.4), y = clim.musd.moavg.wn.recent[1:3,4]*100*c(1.15,1.1,0.9),
           label=c('J15','F15','M15'), size=3.5) +
  ylim(-50,1300) + xlim(-2.5,2) + xlab(expression(paste(VPD['5d'], ' (z-score)'))) + 
  ylab(expression(paste(P['60d'],' (mm)'))) + 
  theme(legend.position = c(0.15,0.62), legend.text = element_text(size=10),
        legend.title = element_text(size=10), axis.text = element_text(size=10),
        axis.title = element_text(size=10), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot as png file
png(file="Script Output/Fig8a_KhHeatMap.png", units='in', width=3, height=2.5, pointsize=10, res=300, bg='white')
p0
dev.off()

########################################################
#Fig 5b: Density plots and quantiles for monthly PLCest#
########################################################

#Density plot of PLCest values
#Historic winter PLCest
day.unq <- paste(clim.musd.wn.historic[,1],clim.musd.wn.historic[,2],clim.musd.wn.historic[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.wn.historic$day.unq <- paste(clim.musd.wn.historic[,1],clim.musd.wn.historic[,2],clim.musd.wn.historic[,3],sep="_")
RWC.wn.hist.est <- merge(Tree.Sim.grid, clim.musd.wn.historic, by="day.unq")[,-1]
RWC.wn.hist.est$plot <- substr(RWC.wn.hist.est$ID,1,1)
RWC.wn.hist.est$species <- substr(RWC.wn.hist.est$ID,3,4)
names(RWC.wn.hist.est) <- c("ID","Month","Day","Year","P60.sml","VPD.zscore5","plot","species")
RWC.wn.hist.est$RWC.est <- predict(best.mod, RWC.wn.hist.est)
Kh.est <- vector(length=nrow(RWC.wn.hist.est))
Kh.est[which(RWC.wn.hist.est$species=="DF")] <- predict(fit.DF, newdata = list(RWCXYL = RWC.wn.hist.est[which(RWC.wn.hist.est$species=="DF"),9]))
Kh.est[which(RWC.wn.hist.est$species=="PP")] <- predict(fit.PP, newdata = list(RWCXYL = RWC.wn.hist.est[which(RWC.wn.hist.est$species=="PP"),9]))
Kh.est[which(RWC.wn.hist.est$species=="SP")] <- predict(fit.SP, newdata = list(RWCXYL = RWC.wn.hist.est[which(RWC.wn.hist.est$species=="SP"),9]))
RWC.wn.hist.est$Kh.est <- Kh.est
Khest.wn.hist.mo <- ddply(RWC.wn.hist.est,.(Year,Month,ID),Kh.est.mo = median(Kh.est,na.rm=TRUE),summarize)

#Historic summer Kh.est
day.unq <- paste(clim.musd.su.historic[,1],clim.musd.su.historic[,2],clim.musd.su.historic[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.su.historic$day.unq <- paste(clim.musd.su.historic[,1],clim.musd.su.historic[,2],clim.musd.su.historic[,3],sep="_")
RWC.su.hist.est <- merge(Tree.Sim.grid, clim.musd.su.historic, by="day.unq")[,-1]
RWC.su.hist.est$plot <- substr(RWC.su.hist.est$ID,1,1)
RWC.su.hist.est$species <- substr(RWC.su.hist.est$ID,3,4)
names(RWC.su.hist.est) <- c("ID","Month","Day","Year","P60.sml","VPD.zscore5","plot","species")
RWC.su.hist.est$RWC.est <- predict(best.mod, RWC.su.hist.est)
Kh.est <- vector(length=nrow(RWC.su.hist.est))
Kh.est[which(RWC.su.hist.est$species=="DF")] <- predict(fit.DF, newdata = list(RWCXYL = RWC.su.hist.est[which(RWC.su.hist.est$species=="DF"),9]))
Kh.est[which(RWC.su.hist.est$species=="PP")] <- predict(fit.PP, newdata = list(RWCXYL = RWC.su.hist.est[which(RWC.su.hist.est$species=="PP"),9]))
Kh.est[which(RWC.su.hist.est$species=="SP")] <- predict(fit.SP, newdata = list(RWCXYL = RWC.su.hist.est[which(RWC.su.hist.est$species=="SP"),9]))
RWC.su.hist.est$Kh.est <- Kh.est
Khest.su.hist.mo <- ddply(RWC.su.hist.est,.(Year,Month,ID),Kh.est.mo = median(Kh.est,na.rm=TRUE),summarize)

#Recent winter Jan PLCest
clim.musd.jan.recent <- clim.musd.wn.recent[clim.musd.wn.recent$Month==1,]
day.unq <- paste(clim.musd.jan.recent[,1],clim.musd.jan.recent[,2],clim.musd.jan.recent[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.jan.recent$day.unq <- paste(clim.musd.jan.recent[,1],clim.musd.jan.recent[,2],clim.musd.jan.recent[,3],sep="_")
RWC.jan.recent.est <- merge(Tree.Sim.grid, clim.musd.jan.recent, by="day.unq")[,c(-1,-3,-4,-5)]
RWC.jan.recent.est$plot <- substr(RWC.jan.recent.est$ID,1,1)
RWC.jan.recent.est$species <- substr(RWC.jan.recent.est$ID,3,4)
names(RWC.jan.recent.est) <- c("ID","P60.sml","VPD.zscore5","plot","species")
RWC.jan.recent.est$RWC.est <- predict(best.mod, RWC.jan.recent.est)
Kh.est <- vector(length=nrow(RWC.jan.recent.est))
Kh.est[which(RWC.jan.recent.est$species=="DF")] <- predict(fit.DF, newdata = list(RWCXYL = RWC.jan.recent.est[which(RWC.jan.recent.est$species=="DF"),6]))
Kh.est[which(RWC.jan.recent.est$species=="PP")] <- predict(fit.PP, newdata = list(RWCXYL = RWC.jan.recent.est[which(RWC.jan.recent.est$species=="PP"),6]))
Kh.est[which(RWC.jan.recent.est$species=="SP")] <- predict(fit.SP, newdata = list(RWCXYL = RWC.jan.recent.est[which(RWC.jan.recent.est$species=="SP"),6]))
RWC.jan.recent.est$Kh.est <- Kh.est
Kh.jan.recent.est <- Kh.est

#Recent winter feb PLCest
clim.musd.feb.recent <- clim.musd.wn.recent[clim.musd.wn.recent$Month==2,]
day.unq <- paste(clim.musd.feb.recent[,1],clim.musd.feb.recent[,2],clim.musd.feb.recent[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.feb.recent$day.unq <- paste(clim.musd.feb.recent[,1],clim.musd.feb.recent[,2],clim.musd.feb.recent[,3],sep="_")
RWC.feb.recent.est <- merge(Tree.Sim.grid, clim.musd.feb.recent, by="day.unq")[,c(-1,-3,-4,-5)]
RWC.feb.recent.est$plot <- substr(RWC.feb.recent.est$ID,1,1)
RWC.feb.recent.est$species <- substr(RWC.feb.recent.est$ID,3,4)
names(RWC.feb.recent.est) <- c("ID","P60.sml","VPD.zscore5","plot","species")
RWC.feb.recent.est$RWC.est <- predict(best.mod, RWC.feb.recent.est)
Kh.est <- vector(length=nrow(RWC.feb.recent.est))
Kh.est[which(RWC.feb.recent.est$species=="DF")] <- predict(fit.DF, newdata = list(RWCXYL = RWC.feb.recent.est[which(RWC.feb.recent.est$species=="DF"),6]))
Kh.est[which(RWC.feb.recent.est$species=="PP")] <- predict(fit.PP, newdata = list(RWCXYL = RWC.feb.recent.est[which(RWC.feb.recent.est$species=="PP"),6]))
Kh.est[which(RWC.feb.recent.est$species=="SP")] <- predict(fit.SP, newdata = list(RWCXYL = RWC.feb.recent.est[which(RWC.feb.recent.est$species=="SP"),6]))
RWC.feb.recent.est$Kh.est <- Kh.est
Kh.feb.recent.est <- Kh.est

#Recent winter mar PLCest
clim.musd.mar.recent <- clim.musd.wn.recent[clim.musd.wn.recent$Month==3,]
day.unq <- paste(clim.musd.mar.recent[,1],clim.musd.mar.recent[,2],clim.musd.mar.recent[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.mar.recent$day.unq <- paste(clim.musd.mar.recent[,1],clim.musd.mar.recent[,2],clim.musd.mar.recent[,3],sep="_")
RWC.mar.recent.est <- merge(Tree.Sim.grid, clim.musd.mar.recent, by="day.unq")[,c(-1,-3,-4,-5)]
RWC.mar.recent.est$plot <- substr(RWC.mar.recent.est$ID,1,1)
RWC.mar.recent.est$species <- substr(RWC.mar.recent.est$ID,3,4)
names(RWC.mar.recent.est) <- c("ID","P60.sml","VPD.zscore5","plot","species")
RWC.mar.recent.est$RWC.est <- predict(best.mod, RWC.mar.recent.est)
Kh.est <- vector(length=nrow(RWC.mar.recent.est))
Kh.est[which(RWC.mar.recent.est$species=="DF")] <- predict(fit.DF, newdata = list(RWCXYL = RWC.mar.recent.est[which(RWC.mar.recent.est$species=="DF"),6]))
Kh.est[which(RWC.mar.recent.est$species=="PP")] <- predict(fit.PP, newdata = list(RWCXYL = RWC.mar.recent.est[which(RWC.mar.recent.est$species=="PP"),6]))
Kh.est[which(RWC.mar.recent.est$species=="SP")] <- predict(fit.SP, newdata = list(RWCXYL = RWC.mar.recent.est[which(RWC.mar.recent.est$species=="SP"),6]))
RWC.mar.recent.est$Kh.est <- Kh.est
Kh.mar.recent.est <- Kh.est

#Recent winter mar PLCest
clim.musd.recent <- clim.musd.wn.recent[clim.musd.wn.recent$Month %in% c(1,2,3),]
day.unq <- paste(clim.musd.recent[,1],clim.musd.recent[,2],clim.musd.recent[,3],sep="_")
Tree.Sim.grid <- expand.grid(unique(d.all$ID),day.unq)
names(Tree.Sim.grid) <- c("ID","day.unq")
clim.musd.recent$day.unq <- paste(clim.musd.recent[,1],clim.musd.recent[,2],clim.musd.recent[,3],sep="_")
RWC.recent.est <- merge(Tree.Sim.grid, clim.musd.recent, by="day.unq")[,c(-1,-3,-4,-5)]
RWC.recent.est$plot <- substr(RWC.recent.est$ID,1,1)
RWC.recent.est$species <- substr(RWC.recent.est$ID,3,4)
names(RWC.recent.est) <- c("ID","P60.sml","VPD.zscore5","plot","species")
RWC.recent.est$RWC.est <- predict(best.mod, RWC.recent.est)
Kh.est <- vector(length=nrow(RWC.recent.est))
RWC.recent.est$Kh = NA
RWC.recent.est[which(RWC.recent.est$species=="DF"),7] <- predict(fit.DF, newdata = list(RWCXYL = RWC.recent.est[which(RWC.recent.est$species=="DF"),6]))
RWC.recent.est[which(RWC.recent.est$species=="PP"),7] <- predict(fit.PP, newdata = list(RWCXYL = RWC.recent.est[which(RWC.recent.est$species=="PP"),6]))
RWC.recent.est[which(RWC.recent.est$species=="SP"),7] <- predict(fit.SP, newdata = list(RWCXYL = RWC.recent.est[which(RWC.recent.est$species=="SP"),6]))


# Plot probability density and overlaid box plots for predicted historic and recent Kh
png(file="Script Output/Fig8b_KhProbDensity_ThreeSpecies.png", units='in', width=2.5, height=2.5, pointsize=10, res=300, bg='white')
source('subScripts/alphafn.R') # Load plotting transparency function

# Probability density plots
par(mfrow=c(3,1), oma = c(4,1,1,1) + 0.1, mar = c(1,1,1,1) + 0.1, cex.lab=1.5, cex.axis=1.5)

# Douglas fir
plot(0,0,type='n',xlim=c(0,3),ylim=c(0,3.5), xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
polygon(density(Khest.su.hist.mo[grep("DF", Khest.su.hist.mo$ID),4],na.rm=TRUE), col='grey50', border='grey50')
par(new=TRUE)
polygon(density(Khest.wn.hist.mo[grep("DF", Khest.wn.hist.mo$ID),4],na.rm=TRUE), col = add.alpha('white',alpha=0.2))
par(new=TRUE)
polygon(density(RWC.recent.est[grep("DF", RWC.recent.est$ID),7],na.rm=TRUE), col = add.alpha('turquoise',alpha=0.3), border = 'turquoise')
axis(side=1, at=c(0,1,2,3,4))
text(x = 0.02, y = 3, pos = 4, cex = 1.2, labels = expression(italic("Pseudotsuga menziesii")))
#mtext(expression(paste(K[h],' [kg ',s^-1,' ',m^-2,' ',MPa^-1,']')), side = 1, outer = TRUE, cex = 1.15, line = 2, col = "grey20", at=0.52)

# Ponderosa pine
plot(0,0,type='n',xlim=c(0,3),ylim=c(0,3.5), xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
polygon(density(Khest.su.hist.mo[grep("PP", Khest.su.hist.mo$ID),4],na.rm=TRUE), col='grey50', border='grey50')
par(new=TRUE)
polygon(density(Khest.wn.hist.mo[grep("PP", Khest.wn.hist.mo$ID),4],na.rm=TRUE), col = add.alpha('white',alpha=0.2))
par(new=TRUE)
polygon(density(RWC.recent.est[grep("PP", RWC.recent.est$ID),7],na.rm=TRUE), col = add.alpha('turquoise',alpha=0.3), border = 'turquoise')
axis(side=1, at=c(0,1,2,3,4))
text(x = 0.02, y = 3, pos = 4, cex = 1.2, labels = expression(italic("Pinus ponderosa")))
#mtext(expression(paste(K[h],' [kg ',s^-1,' ',m^-2,' ',MPa^-1,']')), side = 1, outer = TRUE, cex = 1.15, line = 2, col = "grey20", at=0.52)

# Add legend to plot
points(x = rep(1.6, 3), y = c(0.7, 1.6, 2.5), col = point.col[c(1,4,5)], bg = point.bg[c(1,4,5)], pch = 22, cex = 1.75)
text(x = rep(1.7, 3), y = c(0.7, 1.6, 2.5), labels = c('winter drought 2015','summer 1992-2014','winter 1992-2014'), pos=4, cex = 1)

# Sugar pine
plot(0,0,type='n',xlim=c(0,3),ylim=c(0,3.5), xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
polygon(density(Khest.su.hist.mo[grep("SP", Khest.su.hist.mo$ID),4],na.rm=TRUE), col='grey50', border='grey50')
par(new=TRUE)
polygon(density(Khest.wn.hist.mo[grep("SP", Khest.wn.hist.mo$ID),4],na.rm=TRUE), col = add.alpha('white',alpha=0.2))
par(new=TRUE)
polygon(density(RWC.recent.est[grep("SP", RWC.recent.est$ID),7],na.rm=TRUE), col = add.alpha('turquoise',alpha=0.3), border = 'turquoise')
axis(side=1, at=c(0,1,2,3,4))
text(x = 0.02, y = 3, pos = 4, cex = 1.2, labels = expression(italic("Pinus lambertiana")))
mtext(expression(paste(K[h],' [kg ',s^-1,' ',m^-2,' ',MPa^-1,']')), side = 1, outer = TRUE, cex = 0.9, line = 2, col = "grey20", at=0.52)

dev.off()


############################################
#Fig 5c: Historic/recent PLCest time series#
############################################

#Plot VPD, Precip anom, and PLCest since 1993 for winter and summer
clim.musd.moavg <- clim.musd.moavg[clim.musd.moavg$Year %in% 1993:2014,]
clim.musd.moavg <- clim.musd.moavg[clim.musd.moavg$Month %in% c(1,2,3,6,7,8,9),]
clim.musd.moavg$group <- NA
clim.musd.moavg[clim.musd.moavg$Month %in% c(1,2,3),]$group <- "winter"
clim.musd.moavg[clim.musd.moavg$Month %in% c(6,7,8,9),]$group <- "summer"

#Predict RWC and PLCest based on historic weather data
RWC.hist.pred.DF <- fixef(best.mod)[1]+clim.musd.moavg$P60.sml*fixef(best.mod)[2]+clim.musd.moavg$VPD5.hist*fixef(best.mod)[3]+ ranef(best.mod)$species[1,1]
Kh.hist.pred.DF <- predict(fit.DF, newdata = list(RWCXYL = RWC.hist.pred.DF))
RWC.hist.pred.PP <- fixef(best.mod)[1]+clim.musd.moavg$P60.sml*fixef(best.mod)[2]+clim.musd.moavg$VPD5.hist*fixef(best.mod)[3]+ ranef(best.mod)$species[2,1]
Kh.hist.pred.PP <-  predict(fit.PP, newdata = list(RWCXYL = RWC.hist.pred.PP))
RWC.hist.pred.SP <- fixef(best.mod)[1]+clim.musd.moavg$P60.sml*fixef(best.mod)[2]+clim.musd.moavg$VPD5.hist*fixef(best.mod)[3]+ ranef(best.mod)$species[3,1]
Kh.hist.pred.SP <-  predict(fit.SP, newdata = list(RWCXYL = RWC.hist.pred.SP))
Kh.hist.pred.stand <- rowMeans(data.frame(Kh.hist.pred.DF,Kh.hist.pred.PP,Kh.hist.pred.SP)) #average predictions for three species
clim.musd.moavg$Kh.hist.pred.stand <- Kh.hist.pred.stand

#Plot VPD5 anomaly time series
p1 <- ggplot(clim.musd.moavg)  + stat_smooth(aes(x=Year,y=VPD5.hist.z,group=group,colour=group,fill=group),method="lm") + theme_bw() + theme(legend.position="none", plot.margin = unit(c(.5,.5,.5,.5), 'cm')) + scale_color_manual(values=c(summer='black', winter='white')) + scale_fill_manual(values=c(summer='black', winter='grey65')) + xlab('') + ylab(expression(paste(VPD['5d'], ' (z-score)'))) + scale_x_continuous(breaks=c(1995,2005,2015)) + theme(axis.text = element_text(size=10), axis.title= element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#Plot P60.sml time series
p2 <- ggplot(clim.musd.moavg)  + stat_smooth(aes(x=Year,y=P60.sml,group=group,colour=group,fill=group),method="lm") + theme_bw() + theme(legend.position="none", plot.margin = unit(c(.5,.5,.5,.5), 'cm')) + scale_color_manual(values=c(summer='black', winter='white')) + scale_fill_manual(values=c(summer='black', winter='grey65')) + xlab('Year') + ylab(expression(paste(P['60d'], ' (mm)'))) + scale_x_continuous(breaks=c(1995,2005,2015)) + theme(axis.text = element_text(size = 10), axis.title= element_text(size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#Plot PLCest time series
p3 <- ggplot(clim.musd.moavg)  + stat_smooth(aes(x = Year,y = Kh.hist.pred.stand, group=group, colour=group, fill=group), method="lm") + theme_bw() + theme(legend.position="none", plot.margin = unit(c(.5,.5,.5,.5), 'cm')) + scale_color_manual(values=c(summer='black', winter='white')) + scale_fill_manual(values=c(summer='black', winter='grey65')) + xlab('') + ylab(expression(K[h])) + scale_x_continuous(breaks=c(1995,2005,2015)) + theme(axis.text = element_text(size = 10), axis.title= element_text(size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + annotate('text', x = c(1996, 2007), y = c(1,1.25), label=c('JJAS', 'JFM'), size = 3.5)

#Arrange plots on single panel and save as png file
png(file="Script Output/Fig8c_HistoricTrendsTimeseries.png", units='in', width=6, height=2, pointsize=10, res=300, bg='white')
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

#Check for significance of temporal trend
mmm1.sm <- lm(VPD5.hist.z ~ Year, , data=clim.musd.moavg[clim.musd.moavg$group=="summer",])
mmm1.wn <- lm(VPD5.hist.z ~ Year, data=clim.musd.moavg[clim.musd.moavg$group=="winter",])
mmm2.sm <- lm(P60.sml ~ Year, data=clim.musd.moavg[clim.musd.moavg$group=="summer",])
mmm2.wn <- lm(P60.sml ~ Year, data=clim.musd.moavg[clim.musd.moavg$group=="winter",])
mmm3.sm <- lm(Kh.hist.pred.stand ~ Year, data=clim.musd.moavg[clim.musd.moavg$group=="summer",])
mmm3.wn <- lm(Kh.hist.pred.stand ~ Year, data=clim.musd.moavg[clim.musd.moavg$group=="winter",])
summary(mmm1.sm)
summary(mmm1.wn)
summary(mmm2.sm)
summary(mmm2.wn)
summary(mmm3.sm)
summary(mmm3.wn)


########################################
# Figs. 6 and 7: Ring width index data #
########################################

#Load ring width series as a list
rwi_path = "Script Input Data/RWLfiles/"
rw.series <- dir(rwi_path)[1:52]
rw.list <- lapply(paste0(rwi_path,rw.series), read.rwl)

#Function for extracting ring width and binding it with year
nameCols <- function(i) substr(rw.series[[i]],1,5) #Extract unique tree ID
tree.ID <- sapply(1:52, nameCols)
rwi.fn <- function(i) {
  tmp <- cbind(as.numeric(row.names(rw.list[[i]])),
               rw.list[[i]][,1])
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c('year','ring.width')
  tmp$ID <- rep(tree.ID[i],nrow(tmp))
  return(tmp)
}

#Run functions
rwi_list <- lapply(1:52, rwi.fn)
rwi_df <- melt(rwi_list, id.vars = c('ID','year','ring.width'))[,-4]
rwi_df$species <- substr(rwi_df$ID,3,4)
rwi_df$plot <- substr(rwi_df$ID,1,1)
rwi_df$trt <- substr(rwi_df$ID,2,2)

# Detrend ring width record for each tree using negative exponential assumption
rwi_negexp <- lapply(1:length(rw.list), function(i) detrend(rw.list[[i]])[[1]][,2])
rwi_mat <- matrix(nrow=max(ldply(1:52,function(i) length(rwi_negexp[[i]]))),ncol=52)
for(i in 1:52){
  rwi_mat[1:length(rwi_negexp[[i]]),i] <- rwi_negexp[[i]]
}
rwi_df <- as.data.frame(rwi_mat)

# Calculate stand-level median and 25th/75th quantiles of 
pct75 <- function(x) quantile(x, 0.75, na.rm=TRUE)
pct25 <- function(x) quantile(x, 0.25, na.rm=TRUE)
rwi_med <- apply(rwi_df,1,median, na.rm=TRUE)
rwi_75 <- apply(rwi_df,1,pct75)
rwi_25 <- apply(rwi_df,1,pct25)

# Calculate mean VPD_5d and P_60d for winter (JFM) and summer (JJAS) 
clim_wntr = rbind(aggregate(clim.musd.moavg.wn.historic, by = list(Year = clim.musd.moavg.wn.historic$Year), mean)[,c(1,4,5)], apply(clim.musd.moavg.wn.recent,2,mean)[c(1,3,4)])
clim_sumr = rbind(aggregate(clim.musd.moavg.su.historic, by = list(Year = clim.musd.moavg.su.historic$Year), mean)[,c(1,4,5)], apply(clim.musd.moavg.su.recent,2,mean)[c(1,3,4)])


# Examine statistical relationships
summary(lm(rwi_med[-1] ~ clim_wntr[-1,2]*clim_wntr[-1,3]))
summary(lm(rwi_med[-1] ~ clim_sumr[-1,2]*clim_wntr[-1,3]))
summary(lm(rwi_med[-1] ~ clim_wntr[-1,2]))
summary(lm(rwi_med[-1] ~ clim_sumr[-1,2]))
summary(lm(rwi_med[-1] ~ clim_wntr[-1,3]))
summary(lm(rwi_med[-1] ~ clim_sumr[-1,3]))


#######
# Fig. 6 - Ring width versus climate time series
#######

# Ring width time series
q0 <- ggplot() +
  geom_line(aes(x = 1994:2015, y = rwi_med[-1])) +
  geom_ribbon(aes(x = 1994:2015, ymin = rwi_25[-1], ymax = rwi_75[-1]), alpha = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  xlab('year') +
  ylab('ring width index')

# VPD time series
q1 <- ggplot() +
  geom_line(aes(x = 1994:2015, y = clim_sumr[-1,2]), col = 'coral3') +
  geom_line(aes(x = 1994:2015, y = clim_wntr[-1,2]), col = 'cadetblue3') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  ylab(expression(paste(VPD['5d'], ' (z-score)'))) +
  xlab('') + 
  annotate("text", label="JFM", x = 2012.5, y = 1, col = 'cadetblue3', size = 2.1) + 
  annotate("text", label="JJAS", x = 2014.2, y = 0.39, col = 'coral3', size = 2.1)

# VPD time series
q2 <- ggplot() +
  geom_line(aes(x = 1994:2015, y = clim_sumr[-1,3]*100), col = 'coral3') +
  geom_line(aes(x = 1994:2015, y = clim_wntr[-1,3]*100), col = 'cadetblue3') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  xlab('') + 
  ylab(expression(paste(P['60d'], ' (mm)'))) + 
  annotate("text", label="JFM", x = 2014, y = 390, col = 'cadetblue3', size = 2.1) +
  annotate("text", label="JJAS", x = 2014, y = 90, col = 'coral3', size = 2.1)

# Get the grobs
gA <- ggplotGrob(q0)
gB <- ggplotGrob(q1)
gC <- ggplotGrob(q2)

# Combine the plots
g <- arrangeGrob(grobs = list(rbind(gB, gC, gA, size = "last")),
                 ncol = 1)

png(file="/Users/masonearles/Dropbox/Yale/Projects/PRIORITY_1/Mid-winter Drought/NewPhyt/Fig6_RWI_Clim_Time.png", units='in', width=2.2, height=4.5, pointsize=11, res=300, bg='white')
grid.newpage()
grid.draw(g)
dev.off()

#######
# Fig. 7 - Ring width versus climate scatter plots
#######
library(gridExtra)
p0 <- ggplot() + geom_point(aes(y = rwi_med[-1], x = clim_wntr[-1,2]), col = 'black', size = 1.75) + stat_smooth(aes(y = rwi_med[-1], x = clim_wntr[-1,2]), method = "lm", col = 'black', linetype = 'dashed', size = 0.5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(-.85,1.1) + ylab('ring width index') + xlab(expression(paste(VPD['5d,winter'], ' (z-score)'))) + ylim(0.6,1.2)

p1 <- ggplot() + geom_point(aes(y = rwi_med[-1], x = clim_sumr[-1,2]), col = 'black', size = 1.75) + stat_smooth(aes(y = rwi_med[-1], x = clim_sumr[-1,2]), method = "lm", col = 'black', linetype = 'dashed', size = 0.5)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank()) + xlim(-.85,1) + ylab('') + xlab(expression(paste(VPD['5d,summer'], ' (z-score)'))) + ylim(0.6,1.2)

p2 <- ggplot() + geom_point(aes(y = rwi_med[-1], x = clim_wntr[-1,3]*100), col = 'black', size = 1.75) + stat_smooth(aes(y = rwi_med[-1], x = clim_wntr[-1,3]*100), method = "lm", col = 'black', linetype = 'dashed', size = 0.5)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank())  + xlim(0,1000) + ylab('') + xlab(expression(paste(P['60d,winter'], ' (mm)'))) + ylim(0.6,1.2)

p3 <- ggplot() + geom_point(aes(y = rwi_med[-1], x = clim_sumr[-1,3]*100), col = 'black', size = 1.75) + stat_smooth(aes(y = rwi_med[-1], x = clim_sumr[-1,3]*100), method = "lm", col = 'black', linetype = 'dashed', size = 0.5)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank()) + xlim(0,1000) + ylab('') + xlab(expression(paste(P['60d,summer'], ' (mm)'))) + ylim(0.6,1.2)

# Get the grobs
gA <- ggplotGrob(p0)
gB <- ggplotGrob(p1)
gC <- ggplotGrob(p2)
gD <- ggplotGrob(p3)

# Combine the plots
g <- arrangeGrob(grobs = list(cbind(gA, gB, gC, gD, size = "last")),
                 ncol = 1)

# Save plot
png(file="Fig7_RWIvsClim.png", units='in', width=7, height=1.8, pointsize=11, res=300, bg='white')
grid.newpage()
grid.draw(g)
dev.off()



#######################################################
#Fig S4: Plot fixed and random effects from best model#
#######################################################
#Plot simulated fixed and random effects with confidence bands
p1 <- plotFEsim(FEsim(best.mod))
p2 <- plotREsim(REsim(best.mod))

png(file="Script Output/FigS1_ModelFixandRanEffs.png", units='in', width=8, height=3.5, pointsize=10, res=300, bg='white')
grid.arrange(p1,p2,nrow=1)
dev.off()

#########################################################
#Fig. S5 Histogram of days below freezing in JFM of 2015#
#########################################################
png(file="Script Output/FigS2_DaysBelowFreezing.png", units='px', width=500, height=500, pointsize=12, res=150, bg='white')
par(mfrow=c(1,1), mar=c(4.5,4.5,2,2), oma=c(2,2,2,2))
hist((d[d$Year==2015 & d$Month %in% c(1,2,3),]$T_min),
     xlab="Tmin (deg C)",
     main="") #only one day hits 0.0...all other Tmin's above freezing
dev.off()
