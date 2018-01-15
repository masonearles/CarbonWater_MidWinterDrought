#####Generate tables of VPD and Precip anomalies for different windows of time
#####J. Mason Earles
#####last updated: Jan. 2nd, 2016

#Load libraries
library('plyr')
library('reshape')
library('zoo')

#load weather data (separate table)
d.weather <- read.table('Script Input Data/weather.txt',header=T)

#Julian days of interest
days.of.interest <- c(24,54,84,129,167,205,244,306,360,420,466,498,566,624)

#load data
d <- read.table('Script Input Data/ChallengeVPD.txt',head=TRUE)
d[d$Precip<0,4:13] <- 0
d[d$AVP>5000,4:13] <- 0
vpd.month <- aggregate(VPD~Month,data=d,mean)
vpd.month.yr <- aggregate(VPD~Month+Year,data=d,mean)
vpd.month.yr[261,3] <- NA
vpd.month <- aggregate(VPD~Month,data=d[d$Year==2014,],mean)
vpd.norecent <- vpd.month.yr[vpd.month.yr$Year!="2014",]
vpd.norecent <- vpd.norecent[vpd.norecent$Year!="2015",]
vpd.recent1 <- vpd.month.yr[vpd.month.yr$Year=="2014" & vpd.month.yr$Month>1,]
vpd.recent2 <- vpd.month.yr[vpd.month.yr$Year=="2015" & vpd.month.yr$Month<=9,]

precip.month <- aggregate(Precip~Month,data=d,sum)
precip.month.yr <- aggregate(Precip~Month+Year,data=d,sum)
precip.month.mean <- aggregate(Precip~Month,data=precip.month.yr,mean)
precip.norecent <- precip.month.yr[precip.month.yr$Year!="2014",]
precip.norecent <- precip.norecent[precip.norecent$Year!="2015",]
precip.recent1 <- precip.month.yr[precip.month.yr$Year=="2014",]
precip.recent2 <- precip.month.yr[precip.month.yr$Year=="2015" & precip.month.yr$Month<=9,]

########################
##VPD table generation##
########################

#function for calculating moving average of VPD
ma <- function(n){stats::filter(d[,13],rep(1/n,n), sides=1)}

#define sequence of window lengths for VPD moving average
n <- c(5,10,15,30,45,60)

#calculate moving average for window lengths defined above
d.vpd <- sapply(n, FUN=ma)
colnames(d.vpd) <- c("VPD.ma5","VPD.ma10","VPD.ma15","VPD.ma30","VPD.ma45","VPD.ma60") #add column names
d.vpd <- data.frame(d[,1:3],d.vpd) #add d/m/y columns

#calculate VPD normal for all moving windows
#d.vpd.ma.mean <- aggregate(x=d.vpd,by=list(Month=d.vpd$Month,Day=d.vpd$Day), mean, na.rm=TRUE)[,-c(1,2,5)]
d.vpd.ma.mean <- ddply(d.vpd,.(Month,Day),VPD5.mean = mean(VPD.ma5,na.rm=TRUE),VPD10.mean = mean(VPD.ma10,na.rm=TRUE),VPD15.mean = mean(VPD.ma15,na.rm=TRUE),VPD30.mean = mean(VPD.ma30,na.rm=TRUE),VPD45.mean = mean(VPD.ma45,na.rm=TRUE),VPD60.mean = mean(VPD.ma60,na.rm=TRUE),summarize)
d.vpd.ma.mean <- d.vpd.ma.mean[order(d.vpd.ma.mean$Month),]
d.vpd.ma.mean <- rbind(d.vpd.ma.mean,d.vpd.ma.mean) #temporary data frame
d.vpd.ma.mean$DayJulian <- 1:732
d.vpd.mean <- d.vpd.ma.mean[d.vpd.ma.mean$DayJulian %in% days.of.interest,] #select only days of interest

#calculate VPD sd for all moving windows
#d.vpd.ma.sd <- aggregate(x=d.vpd,by=list(Month=d.vpd$Month,Day=d.vpd$Day), sd, na.rm=TRUE)[,-c(1,2,5)]
d.vpd.ma.sd <- ddply(d.vpd,.(Month,Day),VPD5.sd = sd(VPD.ma5,na.rm=TRUE),VPD10.sd = sd(VPD.ma10,na.rm=TRUE),VPD15.sd = sd(VPD.ma15,na.rm=TRUE),VPD30.sd = sd(VPD.ma30,na.rm=TRUE),VPD45.sd = sd(VPD.ma45,na.rm=TRUE),VPD60.sd = sd(VPD.ma60,na.rm=TRUE),summarize)
d.vpd.ma.sd <- d.vpd.ma.sd[order(d.vpd.ma.sd$Month),]
d.vpd.ma.sd <- rbind(d.vpd.ma.sd,d.vpd.ma.sd)
d.vpd.ma.sd$DayJulian <- 1:732
d.vpd.sd <- d.vpd.ma.sd[d.vpd.ma.sd$DayJulian %in% days.of.interest,]

#select only 2014 and 2015 VPD data
d.vpd.20142015 <- d.vpd[d.vpd$Year %in% c(2014,2015),]
d.vpd.20142015$DayJulian <- 1:nrow(d.vpd.20142015) #add column with julian day
d.vpd.mean.20142015 <- d.vpd.20142015[d.vpd.20142015$DayJulian %in% days.of.interest,]

#Combine recent/historic VPD mean and sd into dataframe, and calculate z-scores; re-name columns
d.vpd.mean.recent <- data.frame(d.vpd.mean.20142015$DayJulian,d.vpd.mean[,3:8],d.vpd.sd[,3:8],d.vpd.mean.20142015[,4:9],((d.vpd.mean.20142015[,4:9]-d.vpd.mean[,3:8])/d.vpd.sd[,3:8]))

names(d.vpd.mean.recent) <- c("day","VPD.hist.mean5","VPD.hist.mean10","VPD.hist.mean15","VPD.hist.mean30","VPD.hist.mean45","VPD.hist.mean60","VPD.hist.sd5","VPD.hist.sd10","VPD.hist.sd15","VPD.hist.sd30","VPD.hist.sd45","VPD.hist.sd60","VPD.recent5","VPD.recent10","VPD.recent15","VPD.recent30","VPD.recent45","VPD.recent60","VPD.zscore5","VPD.zscore10","VPD.zscore15","VPD.zscore30","VPD.zscore45","VPD.zscore60")


###########################
##Precip table generation##
###########################

#function for calculating moving average of VPD
roll <- function(x,n) rollapplyr(x, n+1, function(x) sum(x[-n-1]), fill = NA)

#define sequence of window lengths for precip moving average
n <- c(10,15,30,45,60)

#calculate moving average for window lengths defined above
d.precip <- sapply(n,function(n) roll(d[,10],n))
colnames(d.precip) <- c("precip.ma10","precip.ma15","precip.ma30","precip.ma45","precip.ma60") #add column names
d.precip <- data.frame(d[,1:3],d.precip) #add d/m/y columns

#calculate precip normal for all moving windows
#d.precip.ma.mean <- aggregate(x=d.precip,by=list(Month=d.precip$Month,Day=d.precip$Day), mean, na.rm=TRUE)[,-c(1,2,5)] 
d.precip.ma.mean <- ddply(d.precip,.(Month,Day),P10.mean = mean(precip.ma10,na.rm=TRUE),P15.mean = mean(precip.ma15,na.rm=TRUE),P30.mean = mean(precip.ma30,na.rm=TRUE),P45.mean = mean(precip.ma45,na.rm=TRUE),P60.mean = mean(precip.ma60,na.rm=TRUE),summarize)
d.precip.ma.mean <- d.precip.ma.mean[order(d.precip.ma.mean$Month),]
d.precip.ma.mean <- rbind(d.precip.ma.mean,d.precip.ma.mean) #temporary data frame
d.precip.ma.mean$DayJulian <- 1:732
d.precip.mean <- d.precip.ma.mean[d.precip.ma.mean$DayJulian %in% days.of.interest,] #select only days of interest

#calculate precip sd for all moving windows
#d.precip.ma.sd <- aggregate(x=d.precip,by=list(Month=d.precip$Month,Day=d.precip$Day), sd, na.rm=TRUE)[,-c(1,2,5)]
d.precip.ma.sd <- ddply(d.precip,.(Month,Day),P10.sd = sd(precip.ma10,na.rm=TRUE),P15.sd = sd(precip.ma15,na.rm=TRUE),P30.sd = sd(precip.ma30,na.rm=TRUE),P45.sd = sd(precip.ma45,na.rm=TRUE),P60.sd = sd(precip.ma60,na.rm=TRUE),summarize)
d.precip.ma.sd <- d.precip.ma.sd[order(d.precip.ma.sd$Month),]
d.precip.ma.sd <- rbind(d.precip.ma.sd,d.precip.ma.sd)
d.precip.ma.sd$DayJulian <- 1:732
d.precip.sd <- d.precip.ma.sd[d.precip.ma.sd$DayJulian %in% days.of.interest,]

#select only 2014 and 2015 precip data
d.precip.20142015 <- d.precip[d.precip$Year %in% c(2014,2015),]
d.precip.20142015$DayJulian <- 1:nrow(d.precip.20142015) #add column with julian day
d.precip.mean.20142015 <- d.precip.20142015[d.precip.20142015$DayJulian %in% days.of.interest,]

#Combine recent/historic precip mean and sd into dataframe, and calculate z-scores; re-name columns
d.precip.mean.recent <- data.frame(d.precip.mean.20142015$DayJulian,d.precip.mean[,3:7],d.precip.sd[,3:7],d.precip.mean.20142015[,4:8],(d.precip.mean.20142015[,4:8]-d.precip.mean[,3:7])/d.precip.sd[,3:7])

names(d.precip.mean.recent) <- c("day","precip.hist.mean10","precip.hist.mean15","precip.hist.mean30","precip.hist.mean45","precip.hist.mean60","precip.hist.sd10","precip.hist.sd15","precip.hist.sd30","precip.hist.sd45","precip.hist.sd60","precip.recent10","precip.recent15","precip.recent30","precip.recent45","precip.recent60","precip.zscore10","precip.zscore15","precip.zscore30","precip.zscore45","precip.zscore60")

#Combine VPD and precip tables
d.mean.recent <- data.frame(d.vpd.mean.recent,d.precip.mean.recent[,-1])

#Write table for statistical model
write.table(d.mean.recent,"Script Input Data/VPD and Precip Anomaly.txt",col.names=TRUE)
