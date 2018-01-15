#Read ring series
library('dplR')
library('plyr')
library('reshape2')
library('lme4')
library('ggplot2')
library('gridExtra')

# Set wd
setwd('/Users/masonearles/Dropbox/UC Davis/Explorations/Experimental Drought/Data/Ring Width data/RWLfiles')
#Load RWC data
RWC.max <- apply(read.table('RWC.txt',header=T)[,-1], 1, function(x) max(x,na.rm=TRUE))
RWC.dMax <- 100*(RWC.max-read.table('RWC.txt',header=T)[,-1])
RWC.dMax <- melt(RWC.dMax)


d.RWC <- melt(read.table('RWC.txt',header=T))
d.RWC$dmax.RWC <- RWC.dMax[,2]
d.RWC$day <- as.numeric(substr(d.RWC$variable,2,4))
d.RWC <- d.RWC[d.RWC$day!=24,]
names(d.RWC) <- c("ID","","RWC","dmax.RWC","day")
d.RWC <- d.RWC[,-2]
d.RWC

#Load ring width series as a list
rw.series <- dir()[1:52]
rw.list <- lapply(rw.series, read.rwl)

#Look at detrended ring width values
plot(detrend(rw.list[[1]])[[1]][,2], type = "l", ylim=c(0,2))
lapply(2:length(rw.list), function(i) lines(detrend(rw.list[[i]])[[1]][,2], type = "l", ylim=c(0,2)))

#Function for extracting ring width and binding it with year
nameCols <- function(i) substr(rw.series[[i]],1,5) #Extract unique tree ID
tree.ID <- sapply(1:52, nameCols)

BAI.fn <- function(i) {
  tmp <- cbind(as.numeric(row.names(rw.list[[i]])),
               rw.list[[i]][,1])
  colnames(tmp) <- c('year','ring.width')
  tmp <- as.data.frame(tmp)
  BA <- 3.14*cumsum(tmp$ring.width)^2
  tmp$BA <- BA
  tmp$BAI <- c(NA,diff(BA))
  tmp$sig.BAI <- c(NA,cumsum(diff(BA)))
  tmp$RGR <- c(NA,tmp$BAI[-1])/BA
  #tmp$RGR <- tmp$BAI/max(BA,na.rm=TRUE)
  tmp$maxBA <- max(BA,na.rm=TRUE)
  tmp$dRGR <- c(NA,diff(tmp$RGR))
  tmp$ID <- rep(tree.ID[i],nrow(tmp))
  return(tmp)
}

#Run functions
BAI.trimmed <- lapply(1:52, BAI.fn)

#Trim data to include only 1997 forward
BAI.trimmed <- lapply(1:52, function(i) BAI[[i]][(nrow(BAI[[i]])-18):nrow(BAI[[i]]),])

#Collapse list to dataframe
BAI.trim.melt <- melt(BAI.trimmed, id.vars=c('ID','year','BAI','sig.BAI','RGR','ring.width','maxBA','dRGR','BA'))
BAI.trim.melt$species <- substr(BAI.trim.melt$ID,3,4)
BAI.trim.melt$plot <- substr(BAI.trim.melt$ID,1,1)
BAI.trim.melt$trt <- substr(BAI.trim.melt$ID,2,2)

#Merge BAI and RWC dataframes
BAI.RWC.merge <- merge(BAI.trim.melt, d.RWC, by=c('ID'))

#Add column to dataframe describing RWC during mid-winter drought
RWC.feb15 <- d.RWC[d.RWC$day==420,-4]
BAI.RWC.merge <- merge(BAI.trim.melt, d.RWC, by=c('ID'))
BAI.RWC.merge <- merge(BAI.RWC.merge, RWC.feb15, by=c('ID'))

#Add PLCest to dataframe
BAI.RWC.merge[BAI.RWC.merge$species=='DF',]$dmax.RWC.y <- 100/(1+10^((0.346-(1-BAI.RWC.merge[BAI.RWC.merge$species=='DF',]$RWC.y))*5.15))
BAI.RWC.merge[BAI.RWC.merge$species=='PP',]$dmax.RWC.y <- 100/(1+10^((0.268-(1-BAI.RWC.merge[BAI.RWC.merge$species=='PP',]$RWC.y))*5.98))
BAI.RWC.merge[BAI.RWC.merge$species=='SP',]$dmax.RWC.y <- 100/(1+10^((0.268-(1-BAI.RWC.merge[BAI.RWC.merge$species=='SP',]$RWC.y))*5.98))


ggplot(BAI.RWC.merge[BAI.RWC.merge$species=='DF',]) + geom_smooth(aes(x=year, y=BAI, color=dmax.RWC.y, group=ID)) + theme_bw() + scale_colour_gradient2(low = 'yellow', mid = 'red', high = 'blue', midpoint=25)
ggplot(BAI.RWC.merge[BAI.RWC.merge$species=='PP',]) + geom_smooth(aes(x=year, y=BAI, color=dmax.RWC.y, group=ID)) + theme_bw() + scale_colour_gradient2(low = 'yellow', mid = 'red', high = 'purple', midpoint=50) + xlim(2009,2015)
ggplot(BAI.RWC.merge[BAI.RWC.merge$species=='SP',]) + geom_smooth(aes(x=year, y=BAI, color=dmax.RWC.y, group=ID)) + theme_bw() + scale_colour_gradient2(low = 'yellow', mid = 'red', high = 'blue', midpoint=35)

#Testing if trees that initially had higher relative growth rates ended up with lower relative growth rates
#PP trees with low RGR early on in development (ca. 2002-2006) show lower RWC during mid-winter drought; PP trees with low recent RGR (ca. 2010-2014) show lower RWC during previous 2014 summer drought
#DF trees with low BAI early on in development (ca. 2002-2006) show lower RWC during mid-winter drought
#SP trees with low recent BAI (ca. 2010-2014) show lower RWC during 2015 summer drought
dev.off()
  i=2011
  par(new=TRUE)
  a <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==i,]
  b <- BAI.RWC.merge[BAI.RWC.merge$species=='PP' & BAI.RWC.merge$year==2011,]
  f <- merge(a,b,by=c('ID'))
  plot(f$RWC.y.x~f$RGR.x, ylim=c(0.6,0.9), cex=f$maxBAI.y/15000, col=i)

m0 <- lm(f$RWC.y.x ~ f$RGR.x + f$BA.x)
m1 <- lm(f$RWC.y.x ~ f$RGR.x*f$BA.x)
summary(m1)
AIC(m0,m1)
summary(lm(f$RWC.y.x ~ f$RGR.x*f$BA.x))

RGR.minmax <- seq(round(min(f$RGR.x),2)-0.01,round(max(f$RGR.x),2)+0.01,0.002)
BA.minmax <- seq(round(min(f$BA.x),2),round(max(f$BA.x),2),20)
RGR.BA.grid <- expand.grid(RGR.minmax, BA.minmax)
RGR.BA.grid$RWC.pred <- coef(m1)[1] + coef(m1)[2]*RGR.BA.grid[,1] + coef(m1)[3]*RGR.BA.grid[,2] + coef(m1)[4]*RGR.BA.grid[,1]*RGR.BA.grid[,2]
names(RGR.BA.grid) <- c('RGR','BA','RWC.pred')
ggplot(RGR.BA.grid) + geom_tile(aes(x=RGR, y=BA, fill=RWC.pred)) + scale_fill_gradientn(colours=c("red","yellow","blue"), limits=c(0.7,0.85)) + geom_point(data=f, aes(x=RGR.x, y=BA.x, color=RWC.y.x), size=1.5) + theme_bw()

#Testing out averaging BAI and RGR during several year windows
species <- 'PP'
a.yr <- 2004
b.yr <- 2007
ff <- aggregate(BAI.RWC.merge[BAI.RWC.merge$year>a.yr & BAI.RWC.merge$year<b.yr & BAI.RWC.merge$species==species,], by = list(Year = BAI.RWC.merge[BAI.RWC.merge$year>a.yr & BAI.RWC.merge$year<b.yr & BAI.RWC.merge$species==species,]$ID), function(x) mean(x, na.rm=TRUE))
summary(lm(RWC.y ~ RGR, data=ff))

#Parameter estimates for DF in 2014
.2935+4.253*0.07+.000023*5000-.000235*0.07*5000
.2935+4.253*0.03+.000023*5000-.000235*0.03*5000 #low BA and slow BAI leads to lowest water content in stem during mid-winter drought in PP
.2935+4.253*0.07+.000023*12000-.000235*0.07*12000
.2935+4.253*0.03+.000023*12000-.000235*0.03*12000
.2935+4.253*0.07+.000023*30000-.000235*0.07*30000
.2935+4.253*0.03+.000023*30000-.000235*0.03*30000 #high BA and slow BAI leads to greatest water content in stem during mid-winter drought in PP

#Parameter estimates for DF in 2011
#Signs are different in a logical way
#In 2011, DF trees with high RGR and low BA are most vulnerable
#In 2014 (3 yrs into drought), DF trees with low RGR and low BA are most vulnerable
#In 2011 and 2014, DF trees with low RGR and high BA are least vulnerable

.6483-2.385*0.11+.0000087*5000+.000085*0.11*5000
.6483-2.385*0.045+.0000087*5000+.000085*0.045*5000
.6483-2.385*0.11+.0000087*20000+.000085*0.11*20000
.6483-2.385*0.045+.0000087*25000+.000085*0.045*20000



names(BAI.RWC.merge)
m0 <- lmer(dmax.RWC.y ~ scale(BAI)*scale(maxBA) + (1|species) + (1|plot), data=BAI.RWC.merge)
m1 <- lmer(dmax.RWC.y ~ scale(BAI)*scale(maxBA) + (1|species) + (1|plot), data=BAI.RWC.merge)
m2 <- lmer(dmax.RWC.y ~ scale(BAI)*scale(BA) + (1|species) + (1|plot), data=BAI.RWC.merge)
m3 <- lmer(dmax.RWC.y ~ scale(BAI) + (scale(BAI)|species) + (1|plot), data=BAI.RWC.merge)
summary(m2)

#Get R^2 equivalent for best model; see reference in SI text
mod <- m2
1-var(residuals(mod))/(var(model.response(model.frame(mod))))


summary(m3)
AIC(m0,m1)
fixef(m3)
ranef(m3)


dev.off()
a <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==2000,]
b <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==2014,]
f <- merge(a,b,by=c('ID'))
plot(f$RGR.x~f$RGR.y)
summary(lm(f$RGR.x~f$RGR.y))


#Examining the effect of previous years' RGR on 2015 RGR


#Douglas fir
year.seq <- 1997:2015
DF.output.list <- list()
for(i in 1:length(year.seq)){
  a <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==year.seq[i],]
  b <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==2014,]
  f <- merge(a,b,by=c('ID'))
  plot(f$BAI.x~f$dmax.RWC.y.y)
  DF.output.list[[i]] <- summary(lm(log(f$dmax.RWC.y.y)~f$BAI.x))
}

ff <- f[f$dmax.RWC.y.y<50,]
summary(lm(ff$BAI.x~exp(ff$dmax.RWC.y.y)))
summary(lm(log(ff$BAI.x)~ff$dmax.RWC.y.y))


plot(ff$BAI.x~ff$dmax.RWC.y.y)

DF.BAI.corr <- data.frame(year.seq,
                          unlist(lapply(DF.output.list, function(x) coef(x)[1])),
                          unlist(lapply(DF.output.list, function(x) coef(x)[2])),
                          unlist(lapply(DF.output.list, function(x) x$r.squared)))

colnames(DF.BAI.corr) <- c('year','intercept','slope','r.squared')

#Sugar pine
year.seq <- 1997:2015
SP.output.list <- list()
for(i in 1:length(year.seq)){
  a <- BAI.RWC.merge[BAI.RWC.merge$species=='SP' & BAI.RWC.merge$year==year.seq[i],]
  b <- BAI.RWC.merge[BAI.RWC.merge$species=='SP' & BAI.RWC.merge$year==2014,]
  f <- merge(a,b,by=c('ID'))
  plot(f$BAI.x, f$dmax.RWC.y.y)
  SP.output.list[[i]] <- summary(lm(log(f$dmax.RWC.y.y)~f$RGR.x))
}

SP.BAI.corr <- data.frame(year.seq,
  unlist(lapply(SP.output.list, function(x) coef(x)[1])),
  unlist(lapply(SP.output.list, function(x) coef(x)[2])),
  unlist(lapply(SP.output.list, function(x) x$r.squared)))

colnames(SP.BAI.corr) <- c('year','intercept','slope','r.squared')


#Ponderosa pine
year.seq <- 1997:2015
PP.output.list <- list()
for(i in 1:length(year.seq)){
  a <- BAI.RWC.merge[BAI.RWC.merge$species=='PP' & BAI.RWC.merge$year==year.seq[i],]
  b <- BAI.RWC.merge[BAI.RWC.merge$species=='PP' & BAI.RWC.merge$year==2009,]
  f <- merge(a,b,by=c('ID'))
  plot(f$BAI.x~f$dmax.RWC.y.y)
  PP.output.list[[i]] <- summary(lm(log(f$dmax.RWC.y.y) ~ f$RGR.x))
}

PP.BAI.corr <- data.frame(year.seq,
                          unlist(lapply(PP.output.list, function(x) coef(x)[1])),
                          unlist(lapply(PP.output.list, function(x) coef(x)[2])),
                          unlist(lapply(PP.output.list, function(x) x$r.squared)))

colnames(PP.BAI.corr) <- c('year','intercept','slope','r.squared')

#summarize data for all three species
DF.BAI.corr$species <- 'DF'
PP.BAI.corr$species <- 'PP'
SP.BAI.corr$species <- 'SP'
BAI.corr <- rbind(DF.BAI.corr,PP.BAI.corr,SP.BAI.corr)
p1 <- ggplot(BAI.corr, aes(x=year, y=slope)) + geom_point(aes(color=species, group=species)) + geom_smooth(aes(color=species, group=species), method=lm) + theme_bw()
p2 <- ggplot(BAI.corr, aes(x=year, y=intercept)) + geom_point(aes(color=species, group=species)) + geom_smooth(aes(color=species, group=species)) + theme_bw()
p3 <- ggplot(BAI.corr, aes(x=year, y=r.squared)) + geom_point(aes(color=species, group=species)) + geom_smooth(aes(color=species, group=species)) + theme_bw()
grid.arrange(p1,p2,p3)
p1
p2
p3



d.DF <- BAI.RWC.merge[BAI.RWC.merge$species=='DF' & BAI.RWC.merge$year==2015 & BAI.RWC.merge$day==420,]
plot(d.DF$RWC.x, d.DF$sig.BAI)
ggplot(d.DF, aes(y=RGR.x, x=ye)) + geom_point() + theme_bw() + stat_function(fun = function(x) -1919.3 + 1360.9*exp(x)) + ylim(-200,2500)
  
  
geom_abline(intercept=461.9, slope=-289.5, linetype='dashed') + geom_abline(intercept=-5842, slope=8822, linetype='dashed') + xlim(0.2,0.95)

summary(lm(BAI~RWC.x + I(RWC.x^2), data=d.DF))
summary(lm(BAI~RWC.x, data=d.DF[d.DF$RWC.x<0.7,]))

d.DF <- BAI.RWC.merge[BAI.RWC.merge$species=='SP' & BAI.RWC.merge$year==2014 & BAI.RWC.merge$day==420,]
plot(d.DF$RWC.x, d.DF$RGR)
summary(lm(log(d.DF$RGR)~d.DF$RWC.x))

d.SP <- BAI.RWC.merge[BAI.RWC.merge$species=='SP' & BAI.RWC.merge$year==2015 & BAI.RWC.merge$day==584,]
plot(d.SP$RWC.x, d.SP$sig.BAI)
ggplot() + geom_point(aes(x=d.SP$RWC.x, y=d.SP$sig.BAI))
summary(lm(log(d.SP$sig.BAI)~d.SP$RWC.x))

d.PP <- BAI.RWC.merge[BAI.RWC.merge$species=='PP' & BAI.RWC.merge$year==2015 & BAI.RWC.merge$day==205,]
plot(d.PP$RWC.x, d.PP$sig.BAI)
summary(lm(log(d.PP$RGR)~d.PP$RWC.x))


#Douglas fir
#high vuln: 2CDF1 2CDF2 1TDF3 2TDF2 2TDF3
#low vuln: 1CDF2 1CDF3 3CDF3 2TDF1 3TDF3
#Extract and summarize mean/sd of high vulnerability trees
BAI.highvul <- ldply(lapply(c(11,17:19,27,28), function(i) BAI.trimmed[[i]]))
BAI.highvul.mu <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), mean)
BAI.highvul.se <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), sd)

#Extract and summarize mean/sd of low vulnerability trees
BAI.lowvul <- ldply(lapply(c(2,3,9,37,26,46), function(i) BAI.trimmed[[i]]))
BAI.lowvul.mu <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), mean)
BAI.lowvul.se <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), sd)

#Plot stuff
par(new=FALSE)
for(i in c(11,17:19,27,28)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red', lt=2)
  par(new=TRUE)
}
for(i in c(2,3,9,37,26,46)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue', lt=2)
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red')
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue')


par(new=FALSE)
for(i in c(11,17:19,27,28)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(1995,2015), ylim=c(0,0.12), col='red')
  par(new=TRUE)
}
for(i in c(2,3,9,37,26,46)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(1995,2015), ylim=c(0,0.12), col='blue')
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$RGR, type='l', xlim=c(1995,2015), ylim=c(0,1), col='red')
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$RGR, type='l', xlim=c(1995,2015), ylim=c(0,1), col='blue')


#Ponderosa pine
#low vuln = 1CPP1 2CPP3 2TPP2 2TPP3 3TPP2
#high vuln = 1CPP2 3CPP1 1TPP3 2TPP1 3TPP3
#Extract and summarize mean/sd of high vulnerability trees
BAI.highvul <- ldply(lapply(c(5,38,14,29,49), function(i) BAI.trimmed[[i]]))
BAI.highvul.mu <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), mean, na.rm=TRUE)
BAI.highvul.se <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), sd, na.rm=TRUE)

#Extract and summarize mean/sd of low vulnerability trees
BAI.lowvul <- ldply(lapply(c(4,22,30,31,48), function(i) BAI.trimmed[[i]]))
BAI.lowvul.mu <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), mean, na.rm=TRUE)
BAI.lowvul.se <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), sd, na.rm=TRUE)

#Plot stuff
par(new=FALSE)
for(i in c(5,38,14,29,49)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red', lt=2)
  par(new=TRUE)
}
for(i in c(4,22,30,31,48)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue', lt=2)
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red')
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue')


par(new=FALSE)
for(i in c(5,38,14,29,49)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='red')
  par(new=TRUE)
}
for(i in c(4,22,30,31,48)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='blue')
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='red', lw=2)
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='blue', lw=2)

#Sugar pine
#high vuln: 2CSP1 2CSP2 3CSP2 3CSP3 1TSP2
#low vuln: 3CSP1 1TSP1 2TSP1 2TSP3 3TSP3
#Extract and summarize mean/sd of high vulnerability trees
BAI.highvul <- ldply(lapply(c(23,24,42,43,16), function(i) BAI.trimmed[[i]]))
BAI.highvul.mu <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), mean, na.rm=TRUE)
BAI.highvul.se <- aggregate(BAI.highvul, by=list(Year=BAI.highvul$year), sd, na.rm=TRUE)

#Extract and summarize mean/sd of low vulnerability trees
BAI.lowvul <- ldply(lapply(c(41,15,32,34,52), function(i) BAI.trimmed[[i]]))
BAI.lowvul.mu <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), mean, na.rm=TRUE)
BAI.lowvul.se <- aggregate(BAI.lowvul, by=list(Year=BAI.lowvul$year), sd, na.rm=TRUE)

#Plot stuff
par(new=FALSE)
for(i in c(23,24,42,43,16)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red', lt=2)
  par(new=TRUE)
}
for(i in c(41,15,32,34,52)) {
  plot(BAI[[i]]$year, BAI[[i]]$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue', lt=2)
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='red')
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$BAI, type='l', xlim=c(1995,2015), ylim=c(0,3000), col='blue')





par(new=FALSE)
for(i in c(23,24,42,43,16)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='red')
  par(new=TRUE)
}
for(i in c(41,15,32,34,52)) {
  plot(BAI[[i]]$year, BAI[[i]]$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='blue')
  par(new=TRUE)
}

plot(BAI.highvul.mu$year,BAI.highvul.mu$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='red', lw=2)
par(new=TRUE)
plot(BAI.lowvul.mu$year,BAI.lowvul.mu$RGR, type='l', xlim=c(2005,2015), ylim=c(0,0.2), col='blue', lw=2)



##########
##########

names(BAI.trim.melt)
years = unique(BAI.trim.melt$year)
plot(years, sapply(years, function(i) median(na.omit(BAI.trim.melt[BAI.trim.melt$year == i & BAI.trim.melt$species == 'PP',]$dRGR))), type = 'l', col = 'cadetblue')
lines(years, sapply(years, function(i) median(na.omit(BAI.trim.melt[BAI.trim.melt$year == i & BAI.trim.melt$species == 'DF',]$dRGR))), type = 'l')
lines(years, sapply(years, function(i) median(na.omit(BAI.trim.melt[BAI.trim.melt$year == i & BAI.trim.melt$species == 'SP',]$dRGR))), type = 'l', col = 'coral2')


mean(BAI.trim.melt[BAI.trim.melt$year=i])
aggregate(BAI.trim.melt, by = list(year='year'), mean)



#Look at detrended ring width values
plot(detrend(rw.list[[1]])[[1]][,3], type = "l", ylim=c(0,2))
lapply(2:length(rw.list), function(i) lines(detrend(rw.list[[i]])[[1]][,2], type = "l", ylim=c(0,2)))

rwi_negexp <- lapply(1:length(rw.list), function(i) detrend(rw.list[[i]])[[1]][,2])
rwi_mat <- matrix(nrow=max(ldply(1:52,function(i) length(rwi_negexp[[i]]))),ncol=52)


for(i in 1:52){
  rwi_mat[1:length(rwi_negexp[[i]]),i] <- rwi_negexp[[i]]
}

# Define percentile function
pct90 <- function(x) quantile(x, 0.75, na.rm=TRUE)
pct10 <- function(x) quantile(x, 0.25, na.rm=TRUE)
rwi_df <- as.data.frame(rwi_mat)
rwi_med <- apply(rwi_df,1,median, na.rm=TRUE)
plot(1993:2015, rwi_med, type='l', ylim=c(0.5,1.5))
plot(1993:2015, apply(rwi_df,1,median, na.rm=TRUE), type='l', lwd = 2, ylim=c(0.5,1.5), xlab = "year", ylab = "ring width index")
lines(1993:2015, apply(rwi_df,1,pct90), type='l', ylim=c(0.2,1.7))
lines(1993:2015, apply(rwi_df,1,pct10), type='l', ylim=c(0.2,1.7))
abline(h=1, lty=2)
