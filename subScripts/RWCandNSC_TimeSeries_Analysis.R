#####################################################################################################
# Script for analyzing temporal trends in relative water content, soluble soluble sugar, and starch #
# J. Mason Earles                                                                                   #
# Last updated: Jan. 4th, 2017                                                                      #
#####################################################################################################

# Load libraries
library('agricolae')
library('ggplot2')
library('gridExtra')
library('grid')
library('multcomp')
library('reshape')
library('lme4')
library('merTools')
library('plyr')

# Relative water content analysis
d.RWC <- melt(read.table('Script Input Data/RWC.txt',header=T))
d.RWC$day <- as.numeric(substr(d.RWC$variable,2,4))
d.RWC <- d.RWC[,-2]
names(d.RWC) <- c('ID','RWC','DAY')
d.RWC$DAYf <- as.factor(as.character(d.RWC$DAY))
d.RWC$species <- substr(d.RWC$ID,3,4) #add species variable
d.RWC$plot <- as.factor(substr(d.RWC$ID,1,1)) #add plot variable
d.RWC$DAYf.species <- interaction(d.RWC$DAYf, d.RWC$species)
d.RWC.om <- na.omit(d.RWC)
d.RWC.om$IDunq <- paste(d.RWC.om$ID,d.RWC.om$DAY,sep="_")

# Repeated measures mixed effects model with ID as a random effect. Tukey HSD for multiple comparison.
# All species
d.RWC.om.mu <- aggregate(d.RWC.om, by=list(day=d.RWC.om$DAY), mean)[,c(1,3)]
mm0 <- lmer(RWC ~ DAYf + (1|species/ID), data=d.RWC.om)
mm0.glht <- glht(mm0, linfct=mcp(DAYf = "Tukey"))
mm0.tuk <- cld(summary(mm0.glht))
mm0.tuk <- data.frame(as.numeric(names(mm0.tuk$mcletters$Letters)),mm0.tuk$mcletters$Letters)
names(mm0.tuk) <- c('day','M')
mm0.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm0))[2:13,1]), 5, nchar(names(coef(summary(mm0))[2:13,1]))))), coef(summary(mm0))[1:13,2])
names(mm0.se) <- c('day','se')
mm0 <- merge(merge(d.RWC.om.mu, mm0.se, by='day'), mm0.tuk, by='day')

# Douglas fir
d.RWC.om.mu <- aggregate(d.RWC.om[d.RWC.om$species=='DF',], by=list(day=d.RWC.om[d.RWC.om$species=='DF',]$DAY), mean)[,c(1,3)]
mm1 = lmer(RWC ~ DAYf + (1|ID), data=d.RWC.om[d.RWC.om$species=='DF',])
mm1.tuk <- cld(summary(glht(mm1, linfct=mcp(DAYf = "Tukey"))))
mm1.tuk <- data.frame(as.numeric(names(mm1.tuk$mcletters$Letters)),mm1.tuk$mcletters$Letters)
names(mm1.tuk) <- c('day','M')
mm1.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm1))[2:13,1]), 5, nchar(names(coef(summary(mm1))[2:13,1]))))), coef(summary(mm1))[1:13,2])
names(mm1.se) <- c('day','se')
mm1 <- merge(merge(d.RWC.om.mu, mm1.se, by='day'), mm1.tuk, by='day')

# Ponderosa pine
d.RWC.om.mu <- aggregate(d.RWC.om[d.RWC.om$species=='PP',], by=list(day=d.RWC.om[d.RWC.om$species=='PP',]$DAY), mean)[,c(1,3)]
mm2 = lmer(RWC ~ DAYf + (1|ID), data=d.RWC.om[d.RWC.om$species=='PP',])
mm2.tuk <- cld(summary(glht(mm2, linfct=mcp(DAYf = "Tukey"))))
mm2.tuk <- data.frame(as.numeric(names(mm2.tuk$mcletters$Letters)),mm2.tuk$mcletters$Letters)
names(mm2.tuk) <- c('day','M')
mm2.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm2))[2:13,1]), 5, nchar(names(coef(summary(mm2))[2:13,1]))))), coef(summary(mm2))[1:13,2])
names(mm2.se) <- c('day','se')
mm2 <- merge(merge(d.RWC.om.mu, mm2.se, by='day'), mm2.tuk, by='day')

# Sugar pine
d.RWC.om.mu <- aggregate(d.RWC.om[d.RWC.om$species=='SP',], by=list(day=d.RWC.om[d.RWC.om$species=='SP',]$DAY), mean)[,c(1,3)]
mm3 = lmer(RWC ~ DAYf + (1|ID), data=d.RWC.om[d.RWC.om$species=='SP',])
mm3.tuk <- cld(summary(glht(mm3, linfct=mcp(DAYf = "Tukey"))))
mm3.tuk <- data.frame(as.numeric(names(mm3.tuk$mcletters$Letters)),mm3.tuk$mcletters$Letters)
names(mm3.tuk) <- c('day','M')
mm3.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm3))[2:13,1]), 5, nchar(names(coef(summary(mm3))[2:13,1]))))), coef(summary(mm3))[1:13,2])
names(mm3.se) <- c('day','se')
mm3 <- merge(merge(d.RWC.om.mu, mm3.se, by='day'), mm3.tuk, by='day')

# Plot RWC time series with Tukey HSD comparison 
scaleFUN <- function(x) sprintf("%.2f", x)

p0 <- ggplot(mm0, aes(y=RWC, x=day)) + geom_point(colour='grey50', size = 3) + geom_errorbar(aes(ymin=RWC-se, ymax=RWC+se), colour='grey50', width=0) + geom_line(colour='grey50', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 0.88, label=M), size=3.5) + annotate("text", x=100, y=0.72, label='all species', size=4, fontface="italic") + ylim(0.7,0.88) + xlim(0,650)

p1 <- ggplot(mm1, aes(y=RWC, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=RWC-se, ymax=RWC+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 0.88, label=M), size=3.5) + annotate("text", x=170, y=0.65, label='Pseudotsuga menziesii', size=4, fontface="italic") + scale_y_continuous(limits=c(0.6,0.88), labels=scaleFUN) + xlim(0,650)

p2 <- ggplot(mm2, aes(y=RWC, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=RWC-se, ymax=RWC+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 0.86, label=M), size=3.5) + annotate("text", x=140, y=0.72, label='Pinus ponderosa', size=4, fontface="italic") + ylim(0.68,0.86) + xlim(0,650)

p3 <- ggplot(mm3, aes(y=RWC, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=RWC-se, ymax=RWC+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('days since Jan. 1, 2014') + theme_bw() + theme(legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 0.94, label=M), size=3.5) + annotate("text", x=140, y=0.81, label='Pinus lambertiana', size=4, fontface="italic") + ylim(0.79,0.94) + xlim(0,650)

# Save plot to file
png(file="Script Output/Fig2_RWC_timeseries.png", units='px', width=1000, height=1100, pointsize=12, res=150, bg='white')
q0 <- grid.arrange(p0, p1, p2, p3, ncol=1, heights=c(1,1,1,1.3))
grid.arrange(textGrob(expression(paste('RWC [',m^3,' ', m^-3,']')), rot=90), q0, ncol=2, widths=c(0.05,0.95))
dev.off()


#########################
# Soluble sugar analysis#
#########################

# Soluble carbohydrate analysis
d.sol <- read.table('Script Input Data/NSC_soluble.txt', header=TRUE)
d.sol$IDunq <- paste(d.sol$ID,d.sol$DAY,sep="_")
d.sol$DAYf <- as.factor(as.character(d.sol$DAY))
d.sol$species <- substr(d.sol$ID,3,4) #add species variable
d.sol$plot <- as.factor(substr(d.sol$ID,1,1)) #add plot variable
d.sol <- na.omit(d.sol)

# Repeated measures mixed effects model with ID as a random effect. Tukey HSD for multiple comparison.
# All species
d.sol.mu <- aggregate(d.sol, by=list(day=d.sol$DAY), mean)[,c(1,6)]
mm0 = lmer(CONCAVG ~ DAYf*species + (1|ID), data=d.sol)
mm0.tuk <- cld(summary(glht(mm0, linfct=mcp(DAYf = "Tukey"))))
mm0.tuk <- data.frame(as.numeric(names(mm0.tuk$mcletters$Letters)),mm0.tuk$mcletters$Letters)
names(mm0.tuk) <- c('day','M')
mm0.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm0))[2:12,1]), 5, nchar(names(coef(summary(mm0))[2:12,1]))))), coef(summary(mm0))[1:12,2])
names(mm0.se) <- c('day','se')
mm0 <- merge(merge(d.sol.mu, mm0.se, by='day'), mm0.tuk, by='day')

# Douglas fir
d.sol.mu <- aggregate(d.sol[d.sol$species=='DF',], by=list(day=d.sol[d.sol$species=='DF',]$DAY), mean)[,c(1,6)]
mm1 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.sol[d.sol$species=='DF',])
mm1.tuk <- cld(summary(glht(mm1, linfct=mcp(DAYf = "Tukey"))))
mm1.tuk <- data.frame(as.numeric(names(mm1.tuk$mcletters$Letters)),mm1.tuk$mcletters$Letters)
names(mm1.tuk) <- c('day','M')
mm1.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm1))[2:12,1]), 5, nchar(names(coef(summary(mm1))[2:12,1]))))), coef(summary(mm1))[1:12,2])
names(mm1.se) <- c('day','se')
mm1 <- merge(merge(d.sol.mu, mm1.se, by='day'), mm1.tuk, by='day')

# Ponderosa pine
d.sol.mu <- aggregate(d.sol[d.sol$species=='PP',], by=list(day=d.sol[d.sol$species=='PP',]$DAY), mean)[,c(1,6)]
mm2 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.sol[d.sol$species=='PP',])
mm2.tuk <- cld(summary(glht(mm2, linfct=mcp(DAYf = "Tukey"))))
mm2.tuk <- data.frame(as.numeric(names(mm2.tuk$mcletters$Letters)),mm2.tuk$mcletters$Letters)
names(mm2.tuk) <- c('day','M')
mm2.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm2))[2:12,1]), 5, nchar(names(coef(summary(mm2))[2:12,1]))))), coef(summary(mm2))[1:12,2])
names(mm2.se) <- c('day','se')
mm2 <- merge(merge(d.sol.mu, mm2.se, by='day'), mm2.tuk, by='day')

# Sugar pine
d.sol.mu <- aggregate(d.sol[d.sol$species=='SP',], by=list(day=d.sol[d.sol$species=='SP',]$DAY), mean)[,c(1,6)]
mm3 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.sol[d.sol$species=='SP',])
mm3.tuk <- cld(summary(glht(mm3, linfct=mcp(DAYf = "Tukey"))))
mm3.tuk <- data.frame(as.numeric(names(mm3.tuk$mcletters$Letters)),mm3.tuk$mcletters$Letters)
names(mm3.tuk) <- c('day','M')
mm3.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm3))[2:12,1]), 5, nchar(names(coef(summary(mm3))[2:12,1]))))), coef(summary(mm3))[1:12,2])
names(mm3.se) <- c('day','se')
mm3 <- merge(merge(d.sol.mu, mm3.se, by='day'), mm3.tuk, by='day')

# Plot RWC time series with Tukey HSD comparison 
scaleFUN <- function(x) sprintf("%.2f", x)

p0 <- ggplot(mm0, aes(y=CONCAVG, x=day)) + geom_point(colour='grey50', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='grey50', width=0) + geom_line(colour='grey50', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 42, label=M), size=3.5) + annotate("text", x=500, y=36, label='all species', size=4, fontface="italic") + ylim(14,43) + xlim(0,650)

p1 <- ggplot(mm1, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 46, label=M), size=3.5) + annotate("text", x=500, y=38, label='Pseudotsuga menziesii', size=4, fontface="italic") + ylim(13,47) + xlim(0,650)

p2 <- ggplot(mm2, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 38, label=M), size=3.5) + annotate("text", x=500, y=31, label='Pinus ponderosa', size=4, fontface="italic") + ylim(15,39) + xlim(0,650)

p3 <- ggplot(mm3, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('days since Jan. 1, 2014') + theme_bw() + theme(legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 44, label=M), size=3.5) + annotate("text", x=500, y=36, label='Pinus lambertiana', size=4, fontface="italic") + ylim(11,45) + xlim(0,650)

# Save plots to file
png(file="Script Output/Fig3_SolubleSugar_timeseries.png", units='px', width=1000, height=1100, pointsize=12, res=150, bg='white')
q0 <- grid.arrange(p0, p1, p2, p3, ncol=1, heights=c(1,1,1,1.25))
grid.arrange(textGrob(expression(paste('soluble carbohydrate [mg ',g^-1,']')), rot=90), q0, ncol=2, widths=c(0.05,0.95))
dev.off()


##################
# Starch analysis#
##################

# Starch analysis
d.insol <- read.table('Script Input Data/NSC_insoluble.txt', header=TRUE)
d.insol$IDunq <- paste(d.insol$ID,d.insol$DAY,sep="_")
d.insol$DAYf <- as.factor(as.character(d.insol$DAY))
d.insol$species <- substr(d.insol$ID,3,4) #add species variable
d.insol$plot <- as.factor(substr(d.insol$ID,1,1)) #add plot variable
d.insol <- na.omit(d.insol)

# Repeated measures mixed effects model with ID as a random effect. Tukey HSD for multiple comparison.
# All species
d.insol.mu <- aggregate(d.insol, by=list(day=d.insol$DAY), mean)[,c(1,6)]
mm0 = lmer(CONCAVG ~ DAYf*species + (1|ID), data=d.insol)
mm0.tuk <- cld(summary(glht(mm0, linfct=mcp(DAYf = "Tukey")), test=adjusted('BH')))
mm0.tuk <- data.frame(as.numeric(names(mm0.tuk$mcletters$Letters)),mm0.tuk$mcletters$Letters)
names(mm0.tuk) <- c('day','M')
mm0.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm0))[2:12,1]), 5, nchar(names(coef(summary(mm0))[2:12,1]))))), coef(summary(mm0))[1:12,2])
names(mm0.se) <- c('day','se')
mm0 <- merge(merge(d.insol.mu, mm0.se, by='day'), mm0.tuk, by='day')

# Douglas fir
d.insol.mu <- aggregate(d.insol[d.insol$species=='DF',], by=list(day=d.insol[d.insol$species=='DF',]$DAY), mean)[,c(1,6)]
mm1 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.insol[d.insol$species=='DF',])
mm1.tuk <- cld(summary(glht(mm1, linfct=mcp(DAYf = "Tukey")), test=adjusted('BH')))
mm1.tuk <- data.frame(as.numeric(names(mm1.tuk$mcletters$Letters)),mm1.tuk$mcletters$Letters)
names(mm1.tuk) <- c('day','M')
mm1.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm1))[2:12,1]), 5, nchar(names(coef(summary(mm1))[2:12,1]))))), coef(summary(mm1))[1:12,2])
names(mm1.se) <- c('day','se')
mm1 <- merge(merge(d.insol.mu, mm1.se, by='day'), mm1.tuk, by='day')

# Ponderosa pine
d.insol.mu <- aggregate(d.insol[d.insol$species=='PP',], by=list(day=d.insol[d.insol$species=='PP',]$DAY), mean)[,c(1,6)]
mm2 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.insol[d.insol$species=='PP',])
mm2.tuk <- cld(summary(glht(mm2, linfct=mcp(DAYf = "Tukey")), test=adjusted('BH')))
mm2.tuk <- data.frame(as.numeric(names(mm2.tuk$mcletters$Letters)),mm2.tuk$mcletters$Letters)
names(mm2.tuk) <- c('day','M')
mm2.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm2))[2:12,1]), 5, nchar(names(coef(summary(mm2))[2:12,1]))))), coef(summary(mm2))[1:12,2])
names(mm2.se) <- c('day','se')
mm2 <- merge(merge(d.insol.mu, mm2.se, by='day'), mm2.tuk, by='day')

# Sugar pine
d.insol.mu <- aggregate(d.insol[d.insol$species=='SP',], by=list(day=d.insol[d.insol$species=='SP',]$DAY), mean)[,c(1,6)]
mm3 = lmer(CONCAVG ~ DAYf + (1|ID), data=d.insol[d.insol$species=='SP',])
mm3.tuk <- cld(summary(glht(mm3, linfct=mcp(DAYf = "Tukey")), test=adjusted('BH')))
mm3.tuk <- data.frame(as.numeric(names(mm3.tuk$mcletters$Letters)),mm3.tuk$mcletters$Letters)
names(mm3.tuk) <- c('day','M')
mm3.se <- data.frame(c(129,as.numeric(substr(names(coef(summary(mm3))[2:12,1]), 5, nchar(names(coef(summary(mm3))[2:12,1]))))), coef(summary(mm3))[1:12,2])
names(mm3.se) <- c('day','se')
mm3 <- merge(merge(d.insol.mu, mm3.se, by='day'), mm3.tuk, by='day')

# Plot RWC time series with Tukey HSD comparison 
scaleFUN <- function(x) sprintf("%.2f", x)

p0 <- ggplot(mm0, aes(y=CONCAVG, x=day)) + geom_point(colour='grey50', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='grey50', width=0) + geom_line(colour='grey50', size = 0.75) + xlab('days since Jan. 1, 2014') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 25, label=M), size=3.5) + annotate("text", x=500, y=19, label='all species', size=4, fontface="italic") + ylim(-5,26) + xlim(0,650)

p1 <- ggplot(mm1, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 24, label=M), size=3.5) + annotate("text", x=500, y=15, label='Pseudotsuga menziesii', size=4, fontface="italic") + ylim(-3,25) + xlim(0,650)

p2 <- ggplot(mm2, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 23, label=M), size=3.5) + annotate("text", x=500, y=14, label='Pinus ponderosa', size=4, fontface="italic") + ylim(-3,24) + xlim(0,650)

p3 <- ggplot(mm3, aes(y=CONCAVG, x=day)) + geom_point(colour='red', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='red', width=0) + geom_line(colour='red', size = 0.75) + xlab('days since Jan. 1, 2014') + theme_bw() + theme(legend.position='none',text = element_text(size=13)) + ylab('') + geom_text(aes(y = 28, label=M), size=3.5) + annotate("text", x=500, y=21, label='Pinus lambertiana', size=4, fontface="italic") + ylim(-5,29) + xlim(0,650)

# Save plots to file
png(file="Script Output/Fig4_Starch_timeseries.png", units='px', width=1000, height=1100, pointsize=12, res=150, bg='white')
q0 <- grid.arrange(p0, p1, p2, p3, ncol=1, heights=c(1,1,1,1.25))
grid.arrange(textGrob(expression(paste('starch [mg ',g^-1,']')), rot=90), q0, ncol=2, widths=c(0.05,0.95))
dev.off()

######################
# Species-level only #
######################

# Save models
mm0_sp = mm0
mm1_sp = mm0
mm2_sp = mm0

# Plot time series for RWC, soluble sugar, and strach
scaleFUN <- function(x) sprintf("%.2f", x)

p0 <- ggplot(mm0_sp, aes(y=RWC, x=day)) +
  geom_point(colour='grey50', size = 3) +
  geom_errorbar(aes(ymin=RWC-se, ymax=RWC+se), colour='grey50', width=0) +
  geom_line(colour='grey50', size = 0.75) +
  annotate('rect', xmin=mm0_sp[1,1]-5, xmax=mm0_sp[1,1]+5, ymin = 0.7, ymax = mm0_sp[1,2], fill = 'cadetblue3', alpha=0.6) +
  annotate('rect', xmin=mm0_sp[9,1]-5, xmax=mm0_sp[9,1]+5, ymin = 0.7, ymax = mm0_sp[9,2], fill = 'cadetblue3', alpha=0.6) + xlab('') +
  annotate('segment', x=365, xend=365, y = 0.7, yend = 0.86, col = 'black', linetype='dashed') +
  annotate('text', x = mm0_sp[1,1]+45, y = 0.72, label = 'F14', col = 'cadetblue4', size = 4.5, fontface = 'bold') +
  annotate('text', x = mm0_sp[9,1]+45, y = 0.72, label = 'F15', col = 'cadetblue4', size = 4.5, fontface = 'bold') +
  xlab('') +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y = element_text(size = 12), legend.position='none',text = element_text(size=12)) +
  ylab(expression(paste('RWC [', m^3, ' ', m^-3,']'))) +
  geom_text(aes(y = 0.88, label=M), size=4) +
  ylim(0.7,0.88) + xlim(0,650)

p1 <- ggplot(mm1_sp, aes(y=CONCAVG, x=day)) + geom_point(colour='grey50', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='grey50', width=0) + geom_line(colour='grey50', size = 0.75) + xlab('') + annotate('rect', xmin=mm1_sp[3,1]-5, xmax=mm1_sp[3,1]+5, ymin = 14, ymax = mm1_sp[3,2], fill = 'darkorchid3', alpha=0.4) + annotate('rect', xmin=mm1_sp[11,1]-5, xmax=mm1_sp[11,1]+5, ymin = 14, ymax = mm1_sp[11,2], fill = 'darkorchid3', alpha=0.4) + annotate('segment', x=365, xend=365, y = 14, yend = 39, col = 'black', linetype='dashed') +
  annotate('text', x = mm1_sp[3,1]-45, y = 17, label = 'M14', col = 'darkorchid3', size = 4.5, fontface = 'bold') +
  annotate('text', x = mm1_sp[11,1]-45, y = 17, label = 'M15', col = 'darkorchid3', size = 4.5, fontface = 'bold') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y = element_text(size = 12), legend.position='none',text = element_text(size=13)) + ylab(expression(paste('SC [mg', ' ', g^-1,']'))) + geom_text(aes(y = 42, label=M), size=4) + ylim(14,43) + xlim(0,650)

p2 <- ggplot(mm2_sp, aes(y=CONCAVG, x=day)) + geom_point(colour='grey50', size = 3) + geom_errorbar(aes(ymin=CONCAVG-se, ymax=CONCAVG+se), colour='grey50', width=0) + geom_line(colour='grey50', size = 0.75) + annotate('rect', xmin=mm1_sp[3,1]-5, xmax=mm1_sp[3,1]+5, ymin = -5, ymax = mm2_sp[3,2], fill = 'darkorchid3', alpha=0.4) + annotate('rect', xmin=mm1_sp[11,1]-5, xmax=mm1_sp[11,1]+5, ymin = -5, ymax = mm2_sp[11,2], fill = 'darkorchid3', alpha=0.4) + annotate('segment', x=365, xend=365, y = -5, yend = 22, col = 'black', linetype='dashed') +
  annotate('text', x = mm1_sp[3,1]-45, y = -2, label = 'M14', col = 'darkorchid3', size = 4.5, fontface = 'bold') +
  annotate('text', x = mm1_sp[11,1]-45, y = -2, label = 'M15', col = 'darkorchid3', size = 4.5, fontface = 'bold') +
  theme_bw() + theme( legend.position='none',text = element_text(size=12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0))) + ylab(expression(paste('St [mg', ' ', g^-1,']'))) + geom_text(aes(y = 25, label=M), size=4) + ylim(-5,26) + xlim(0,650) + xlab('\ndays since Jan. 1, 2014')

pdf('Script Output/Fig2_Species_Timeseries.pdf', width = 6, height = 5, pointsize = 12)
grid.draw(rbind(ggplotGrob(p0),ggplotGrob(p1),ggplotGrob(p2), size = 'first'))
dev.off()
