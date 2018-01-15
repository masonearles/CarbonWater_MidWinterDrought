###########################################
# Examine relationship between RWC and KH #
# J. Mason Earles                         #
# Updated: Jan. 18th, 2017                #
###########################################

# Load libraries
library('nls2')
library('nlstools')
library('investr')

# Load data
d = read.table('Script Input Data/RWC_KH_calibration.txt', header=TRUE)
d$WP = -d$WP

# Examine relationship between Kh and RWC
# Fit power function for DF and PP
p.fn = function(RWCXYL,b,z) b*RWCXYL^z
fit.DF = nls(Kh ~ p.fn(RWCXYL,b,z), start = list(b=1.25, z=3), data = d[d$SPECIES=='DF',])
summary(fit.DF)
fit.PP = nls(Kh ~ p.fn(RWCXYL,b,z), start = list(b=1.25, z=3), data = d[d$SPECIES=='PP',])
summary(fit.PP)
fit.SP = nls(Kh ~ p.fn(RWCXYL,b,z), start = list(b=1.25, z=3), data = d[d$SPECIES=='SP',])
summary(fit.SP)

# Calculate RMSE
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

RMSE(d[d$SPECIES=='DF',6], predict(fit.DF))
RMSE(na.omit(d[d$SPECIES=='PP',-5])[,5], predict(fit.PP))
RMSE(d[d$SPECIES=='SP',6], predict(fit.SP))


# Plot relationship between Kh and RWC
png(file="Script Output/Fig6_KhvsRWC.png", units='in', width=6, height=2.5, pointsize=11, res=300, bg='white')
par(mfrow=c(1,3), oma = c(4,4,2,2) + 0.1, mar = c(1,1,1,1) + 0.1, cex.lab=1.5, cex.axis=1.25)
# Plot for DF
plotFit(fit.DF, interval = 'confidence', xlim = c(0.2,1), ylim = c(0,2.5), yaxt='n',xlab='', ylab='', cex=1.5, pch=21, col=2, bg=2, shade = TRUE, col.conf = 'grey80')
text(x=0.2, y=2.2, labels=expression(italic('P. menziesii')), cex=1.5, pos=4)
axis(side=2, at = c(0,1,2))

# Plot for PP
plotFit(fit.PP, interval = 'confidence', xlim = c(0.2,1), ylim = c(0,2.5), xlab='', ylab='', yaxt='n', cex=1.5, pch=21, col=2, bg=2, shade = TRUE, col.conf = 'grey80')
text(x=0.2, y=2.2, labels=expression(italic('P. ponderosa')), cex=1.5, pos=4)

# Plot for SP
plotFit(fit.SP, interval = 'confidence', xlim = c(0.2,1), ylim = c(0,2.5), xlab='', ylab='', yaxt='n', cex=1.5, pch=21, col=2, bg=2, shade = TRUE, col.conf = 'grey80')
text(x=0.2, y=2.2, labels=expression(italic('P. lambertiana')), cex=1.5, pos=4)

mtext(expression(paste('RWC [',m^3,' ',m^-3,']')), side = 1, outer = TRUE, cex = 1, line = 2.5, col = "grey20", at=0.52)
mtext(expression(paste(K[h],' [kg ',s^-1,' ',m^-2,' ',MPa^-1,']')), side = 2, outer = TRUE, cex = 1, line = 2, col = "grey20", at=0.52)

dev.off()