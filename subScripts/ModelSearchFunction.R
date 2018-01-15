#Model search sub-script
#Earles et al.
#last updated: Feb. 8th, 2017

#Test precip*VPD models with random intercept only; set up global models for systematic AIC model comparison
m1.pv.RI <- lmer(RWC~  precip.zscore15*VPD.zscore5 + precip.zscore15*VPD.zscore15 + precip.zscore15*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m2.pv.RI <- lmer(RWC~  precip.zscore30*VPD.zscore5 + precip.zscore30*VPD.zscore15 + precip.zscore30*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m3.pv.RI <- lmer(RWC~  precip.zscore45*VPD.zscore5 + precip.zscore45*VPD.zscore15 + precip.zscore45*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m4.pv.RI <- lmer(RWC~  precip.zscore60*VPD.zscore5 + precip.zscore60*VPD.zscore15 + precip.zscore60*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m5.pv.RI <- lmer(RWC~  P30.sml*VPD.zscore5 + P30.sml*VPD.zscore15 + P30.sml*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m6.pv.RI <- lmer(RWC~  P45.sml*VPD.zscore5 + P45.sml*VPD.zscore15 + P45.sml*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m7.pv.RI <- lmer(RWC~  P60.sml*VPD.zscore5 + P60.sml*VPD.zscore15 + P60.sml*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)

#Test precip*VPD models with random intercept only; set up global models for systematic AIC model comparison
m1.vp.RI <- lmer(RWC ~ precip.zscore15*VPD.zscore5 + precip.zscore30*VPD.zscore5 + precip.zscore45*VPD.zscore5 + precip.zscore60*VPD.zscore5 +(1|ID) + (1|species/plot), data=d.all.complete)
m2.vp.RI <- lmer(RWC ~ precip.zscore15*VPD.zscore15 + precip.zscore30*VPD.zscore15 + precip.zscore45*VPD.zscore15 + precip.zscore60*VPD.zscore15 +(1|ID) + (1|species/plot), data=d.all.complete)
m3.vp.RI <- lmer(RWC ~ precip.zscore15*VPD.zscore30 + precip.zscore30*VPD.zscore30 + precip.zscore45*VPD.zscore30 + precip.zscore60*VPD.zscore30 +(1|ID) + (1|species/plot), data=d.all.complete)
m4.vp.RI <- lmer(RWC ~ P30.sml*VPD.zscore5 + P45.sml*VPD.zscore5 + P60.sml*VPD.zscore5 + (1|ID) + (1|species/plot), data=d.all.complete)
m5.vp.RI <- lmer(RWC ~ P30.sml*VPD.zscore15 + P45.sml*VPD.zscore15 + P60.sml*VPD.zscore15 + (1|ID) + (1|species/plot), data=d.all.complete)
m6.vp.RI <- lmer(RWC ~ P30.sml*VPD.zscore30 + P45.sml*VPD.zscore30 + P60.sml*VPD.zscore30 + (1|ID) + (1|species/plot), data=d.all.complete)

#Test precip*VPD models with precip random slope; set up global models for systematic AIC model comparison
m1.pv <- lmer(RWC~  precip.zscore15*VPD.zscore5 + precip.zscore15*VPD.zscore15 + precip.zscore15*VPD.zscore30 +(precip.zscore15|ID) + (1|species/plot), data=d.all.complete)
m2.pv <- lmer(RWC~  precip.zscore30*VPD.zscore5 + precip.zscore30*VPD.zscore15 + precip.zscore30*VPD.zscore30 +(precip.zscore30|ID) + (1|species/plot), data=d.all.complete)
m3.pv <- lmer(RWC~  precip.zscore45*VPD.zscore5 + precip.zscore45*VPD.zscore15 + precip.zscore45*VPD.zscore30 +(precip.zscore45|ID) + (1|species/plot), data=d.all.complete)
m4.pv <- lmer(RWC~  precip.zscore60*VPD.zscore5 + precip.zscore60*VPD.zscore15 + precip.zscore60*VPD.zscore30 +(precip.zscore60|ID) + (1|species/plot), data=d.all.complete)
m5.pv <- lmer(RWC~  P30.sml*VPD.zscore5 + P30.sml*VPD.zscore15 + P30.sml*VPD.zscore30 +(P30.sml|ID) + (1|species/plot), data=d.all.complete)
m6.pv <- lmer(RWC~  P45.sml*VPD.zscore5 + P45.sml*VPD.zscore15 + P45.sml*VPD.zscore30 +(P45.sml|ID) + (1|species/plot), data=d.all.complete)
m7.pv <- lmer(RWC~  P60.sml*VPD.zscore5 + P60.sml*VPD.zscore15 + P60.sml*VPD.zscore30 +(P60.sml|ID) + (1|species/plot), data=d.all.complete)

#Test precip*VPD models with VPD random slope; set up global models for systematic AIC model comparison
m1.vp <- lmer(RWC ~ precip.zscore15*VPD.zscore5 + precip.zscore30*VPD.zscore5 + precip.zscore45*VPD.zscore5 + precip.zscore60*VPD.zscore5 +(VPD.zscore5|ID) + (1|species/plot), data=d.all.complete)
m2.vp <- lmer(RWC ~ precip.zscore15*VPD.zscore15 + precip.zscore30*VPD.zscore15 + precip.zscore45*VPD.zscore15 + precip.zscore60*VPD.zscore15 +(VPD.zscore15|ID) + (1|species/plot), data=d.all.complete)
m3.vp <- lmer(RWC ~ precip.zscore15*VPD.zscore30 + precip.zscore30*VPD.zscore30 + precip.zscore45*VPD.zscore30 + precip.zscore60*VPD.zscore30 +(VPD.zscore30|ID) + (1|species/plot), data=d.all.complete)
m4.vp <- lmer(RWC ~ P30.sml*VPD.zscore5 + P45.sml*VPD.zscore5 + P60.sml*VPD.zscore5 + (VPD.zscore5|ID) + (1|species/plot), data=d.all.complete)
m5.vp <- lmer(RWC ~ P30.sml*VPD.zscore15 + P45.sml*VPD.zscore15 + P60.sml*VPD.zscore15 + (VPD.zscore15|ID) + (1|species/plot), data=d.all.complete)
m6.vp <- lmer(RWC ~ P30.sml*VPD.zscore30 + P45.sml*VPD.zscore30 + P60.sml*VPD.zscore30 + (VPD.zscore30|ID) + (1|species/plot), data=d.all.complete)

#Exhaustive AIC comparison of all global models m*.pv defined above
m1.pv.RI <- dredge(m1.pv.RI, m.lim=c(0,3))
m2.pv.RI <- dredge(m2.pv.RI, m.lim=c(0,3))
m3.pv.RI <- dredge(m3.pv.RI, m.lim=c(0,3))
m4.pv.RI <- dredge(m4.pv.RI, m.lim=c(0,3))
m5.pv.RI <- dredge(m5.pv.RI, m.lim=c(0,3))
m6.pv.RI <- dredge(m6.pv.RI, m.lim=c(0,3))
m7.pv.RI <- dredge(m7.pv.RI, m.lim=c(0,3))
m1.vp.RI <- dredge(m1.vp.RI, m.lim=c(0,3))
m2.vp.RI <- dredge(m2.vp.RI, m.lim=c(0,3))
m3.vp.RI <- dredge(m3.vp.RI, m.lim=c(0,3))
m4.vp.RI <- dredge(m4.vp.RI, m.lim=c(0,3))
m5.vp.RI <- dredge(m5.vp.RI, m.lim=c(0,3))
m6.vp.RI <- dredge(m6.vp.RI, m.lim=c(0,3))
m1.pv <- dredge(m1.pv, m.lim=c(0,3))
m2.pv <- dredge(m2.pv, m.lim=c(0,3))
m3.pv <- dredge(m3.pv, m.lim=c(0,3))
m4.pv <- dredge(m4.pv, m.lim=c(0,3))
m5.pv <- dredge(m5.pv, m.lim=c(0,3))
m6.pv <- dredge(m6.pv, m.lim=c(0,3))
m7.pv <- dredge(m7.pv, m.lim=c(0,3))
m1.vp <- dredge(m1.vp, m.lim=c(0,3))
m2.vp <- dredge(m2.vp, m.lim=c(0,3))
m3.vp <- dredge(m3.vp, m.lim=c(0,3))
m4.vp <- dredge(m4.vp, m.lim=c(0,3))
m5.vp <- dredge(m5.vp, m.lim=c(0,3))
m6.vp <- dredge(m6.vp, m.lim=c(0,3))

#Test precip * absolute VPD models with VPD random slope
m1.pva <- lmer(RWC~  precip.zscore15*VPD5 + precip.zscore15*VPD15 + precip.zscore15*VPD30 +(precip.zscore15|ID) + (1|species/plot), data=d.all.complete)
m2.pva <- lmer(RWC~  precip.zscore30*VPD5 + precip.zscore30*VPD15 + precip.zscore30*VPD30 +(precip.zscore30|ID) + (1|species/plot), data=d.all.complete)
m3.pva <- lmer(RWC~  precip.zscore45*VPD5 + precip.zscore45*VPD15 + precip.zscore45*VPD30 +(precip.zscore45|ID) + (1|species/plot), data=d.all.complete)
m4.pva <- lmer(RWC~  precip.zscore60*VPD5 + precip.zscore60*VPD15 + precip.zscore60*VPD30 +(precip.zscore60|ID) + (1|species/plot), data=d.all.complete)
m5.pva <- lmer(RWC~  P30.sml*VPD5 + P30.sml*VPD15 + P30.sml*VPD30 +(P30.sml|ID) + (1|species/plot), data=d.all.complete)
m6.pva <- lmer(RWC~  P45.sml*VPD5 + P45.sml*VPD15 + P45.sml*VPD30 +(P45.sml|ID) + (1|species/plot), data=d.all.complete)
m7.pva <- lmer(RWC~  P60.sml*VPD5 + P60.sml*VPD15 + P60.sml*VPD30 +(P60.sml|ID) + (1|species/plot), data=d.all.complete)

#Test precip * absolute VPD models with precip random slope
m1.vpa <- lmer(RWC~  VPD5*precip.zscore15 + VPD5*precip.zscore30 + VPD5*precip.zscore45 + VPD5*precip.zscore60 +(VPD5|ID) + (1|species/plot), data=d.all.complete)
m2.vpa <- lmer(RWC~  VPD15*precip.zscore15 + VPD15*precip.zscore30 + VPD15*precip.zscore45 + VPD15*precip.zscore60 +(VPD15|ID) + (1|species/plot), data=d.all.complete)
m3.vpa <- lmer(RWC~  VPD30*precip.zscore15 + VPD30*precip.zscore30 + VPD30*precip.zscore45 + VPD30*precip.zscore60 +(VPD30|ID) + (1|species/plot), data=d.all.complete)
m4.vpa <- lmer(RWC~  VPD5*P30.sml + VPD5*P45.sml + VPD5*P60.sml + (VPD5|ID) + (1|species/plot), data=d.all.complete)
m5.vpa <- lmer(RWC~  VPD15*P30.sml + VPD15*P45.sml + VPD15*P60.sml + (VPD15|ID) + (1|species/plot), data=d.all.complete)
m6.vpa <- lmer(RWC~  VPD30*P30.sml + VPD30*P45.sml + VPD30*P60.sml + (VPD30|ID) + (1|species/plot), data=d.all.complete)

#Test precip * absolute VPD models with VPD random slope
m1.pva.RI <- lmer(RWC~  precip.zscore15*VPD5 + precip.zscore15*VPD15 + precip.zscore15*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)
m2.pva.RI <- lmer(RWC~  precip.zscore30*VPD5 + precip.zscore30*VPD15 + precip.zscore30*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)
m3.pva.RI <- lmer(RWC~  precip.zscore45*VPD5 + precip.zscore45*VPD15 + precip.zscore45*VPD30 +(precip.zscore45|ID) + (1|species/plot), data=d.all.complete)
m4.pva.RI <- lmer(RWC~  precip.zscore60*VPD5 + precip.zscore60*VPD15 + precip.zscore60*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)
m5.pva.RI <- lmer(RWC~  P30.sml*VPD5 + P30.sml*VPD15 + P30.sml*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)
m6.pva.RI <- lmer(RWC~  P45.sml*VPD5 + P45.sml*VPD15 + P45.sml*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)
m7.pva.RI <- lmer(RWC~  P60.sml*VPD5 + P60.sml*VPD15 + P60.sml*VPD30 +(1|ID) + (1|species/plot), data=d.all.complete)

#Test precip * absolute VPD models with precip random slope
m1.vpa.RI <- lmer(RWC~  VPD5*precip.zscore15 + VPD5*precip.zscore30 + VPD5*precip.zscore45 + VPD5*precip.zscore60 +(1|ID) + (1|species/plot), data=d.all.complete)
m2.vpa.RI <- lmer(RWC~  VPD15*precip.zscore15 + VPD15*precip.zscore30 + VPD15*precip.zscore45 + VPD15*precip.zscore60 +(1|ID) + (1|species/plot), data=d.all.complete)
m3.vpa.RI <- lmer(RWC~  VPD30*precip.zscore15 + VPD30*precip.zscore30 + VPD30*precip.zscore45 + VPD30*precip.zscore60 +(1|ID) + (1|species/plot), data=d.all.complete)
m4.vpa.RI <- lmer(RWC~  VPD5*P30.sml + VPD5*P45.sml + VPD5*P60.sml + (1|ID) + (1|species/plot), data=d.all.complete)
m5.vpa.RI <- lmer(RWC~  VPD15*P30.sml + VPD15*P45.sml + VPD15*P60.sml + (1|ID) + (1|species/plot), data=d.all.complete)
m6.vpa.RI <- lmer(RWC~  VPD30*P30.sml + VPD30*P45.sml + VPD30*P60.sml + (1|ID) + (1|species/plot), data=d.all.complete)

#Exhaustive AIC comparison of all global models m*.pv defined above
dr.m1.pva <- dredge(m1.pva, m.lim=c(0,3))
dr.m2.pva <- dredge(m2.pva, m.lim=c(0,3))
dr.m3.pva <- dredge(m3.pva, m.lim=c(0,3))
dr.m4.pva <- dredge(m4.pva, m.lim=c(0,3))
dr.m5.pva <- dredge(m5.pva, m.lim=c(0,3))
dr.m6.pva <- dredge(m6.pva, m.lim=c(0,3))
dr.m7.pva <- dredge(m7.pva, m.lim=c(0,3))
dr.m1.vpa <- dredge(m1.vpa, m.lim=c(0,3))
dr.m2.vpa <- dredge(m2.vpa, m.lim=c(0,3))
dr.m3.vpa <- dredge(m3.vpa, m.lim=c(0,3))
dr.m4.vpa <- dredge(m4.vpa, m.lim=c(0,3))
dr.m5.vpa <- dredge(m5.vpa, m.lim=c(0,3))
dr.m6.vpa <- dredge(m6.vpa, m.lim=c(0,3))
dr.m1.pva.RI <- dredge(m1.pva.RI, m.lim=c(0,3))
dr.m2.pva.RI <- dredge(m2.pva.RI, m.lim=c(0,3))
dr.m3.pva.RI <- dredge(m3.pva.RI, m.lim=c(0,3))
dr.m4.pva.RI <- dredge(m4.pva.RI, m.lim=c(0,3))
dr.m5.pva.RI <- dredge(m5.pva.RI, m.lim=c(0,3))
dr.m6.pva.RI <- dredge(m6.pva.RI, m.lim=c(0,3))
dr.m7.pva.RI <- dredge(m7.pva.RI, m.lim=c(0,3))
dr.m1.vpa.RI <- dredge(m1.vpa.RI, m.lim=c(0,3))
dr.m2.vpa.RI <- dredge(m2.vpa.RI, m.lim=c(0,3))
dr.m3.vpa.RI <- dredge(m3.vpa.RI, m.lim=c(0,3))
dr.m4.vpa.RI <- dredge(m4.vpa.RI, m.lim=c(0,3))
dr.m5.vpa.RI <- dredge(m5.vpa.RI, m.lim=c(0,3))
dr.m6.vpa.RI <- dredge(m6.vpa.RI, m.lim=c(0,3))

#Extract best model from each dredge
#Global model m*.pv.RI (pv = holds precip constant; RI = random intercept only)
model.list <- list(m1.pv.RI,m2.pv.RI,m3.pv.RI,m4.pv.RI,m5.pv.RI,m6.pv.RI,m7.pv.RI)
model.AICc <- ldply(lapply(model.list, function(X) X[1,]$AICc))
AICc.bestmod <- cbind("m.pv.RI",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),])

#Global model m*.vp.RI (pv = holds VPD constant; RI = random intercept only)
model.list <- list(m1.vp.RI,m2.vp.RI,m3.vp.RI,m4.vp.RI,m5.vp.RI,m6.vp.RI)
model.AICc <- ldply(lapply(model.list, function(X) X[1,]$AICc))
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.vp.RI",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Global model m*.pv (pv = holds precip constant; random intercept + random slope for VPD anomaly)
model.list <- list(m1.pv,m2.pv,m3.pv,m4.pv,m5.pv,m6.pv,m7.pv)
model.AICc <- ldply(lapply(model.list, function(X) X[1,]$AICc))
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.pv",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Global model m*.vp (pv = holds VPD constant; random intercept + random slope for precipitation)
model.list <- list(m1.vp,m2.vp,m3.vp,m4.vp,m5.vp,m6.vp)
model.AICc <- ldply(lapply(model.list, function(X) X[1,]$AICc))
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.vp",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Global model m*.pva (pva = holds absolute VPD constant; random intercept + random slope for absolute VPD)
model.list <- list(dr.m1.vpa,dr.m2.vpa,dr.m3.vpa,dr.m4.vpa,dr.m5.vpa,dr.m6.vpa)
model.AICc <- ldply(lapply(model.list, function(X) X[2,]$AICc)) #first model contains two absolute VPD terms which is not allowed, so second best model selected
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.vpa",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Global model m*.pva (pva = holds absolute VPD constant; random intercept + random slope for precipitation)
model.list <- list(dr.m1.pva,dr.m2.pva,dr.m3.pva,dr.m4.pva,dr.m5.pva,dr.m6.pva,dr.m7.pva)
model.AICc <- ldply(lapply(model.list, function(X) X[2,]$AICc)) #first model contains two absolute VPD terms which is not allowed, so second best model selected
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.pva",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Global model m*.pva.RI (vpa = holds precip constant; random intercept + random slope for absolute VPD)
model.list <- list(dr.m1.vpa.RI,dr.m2.vpa.RI,dr.m3.vpa.RI,dr.m4.vpa.RI,dr.m5.vpa.RI,dr.m6.vpa.RI)
model.AICc <- ldply(lapply(model.list, function(X) X[2,]$AICc)) #first model contains two absolute VPD terms which is not allowed, so second best model selected
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.vpa.RI",which(model.AICc==min(model.AICc))[1],model.AICc[which(model.AICc==min(model.AICc))[1],]))

#Global model m*.pva.RI (pva = holds absolute VPD constant; random intercept + random slope for precipitation)
model.list <- list(dr.m1.pva.RI,dr.m2.pva.RI,dr.m3.pva.RI,dr.m4.pva.RI,dr.m5.pva.RI,dr.m6.pva.RI,dr.m7.pva.RI)
model.AICc <- ldply(lapply(model.list, function(X) X[2,]$AICc)) #first model contains two absolute VPD terms which is not allowed, so second best model selected
AICc.bestmod <- rbind(AICc.bestmod,cbind("m.pva.RI",which(model.AICc==min(model.AICc)),model.AICc[which(model.AICc==min(model.AICc)),]))

#Calculate the Akaike weights of each model
AIC.wts.globalmods <- (exp(-0.5*(as.numeric(AICc.bestmod[,3])-min(as.numeric(AICc.bestmod[,3])))))/sum((exp(-0.5*(as.numeric(AICc.bestmod[,3])-min(as.numeric(AICc.bestmod[,3]))))))
AICc.bestmod <- data.frame(AICc.bestmod[,1],as.numeric(AICc.bestmod[,2]),as.numeric(AICc.bestmod[,3]),AIC.wts.globalmods)
names(AICc.bestmod) <- c("GlobalModel","SubModel","AICc","AIC.wt")

#100% of AIC weight given to global model m.pv submodel 7. The model has the following structure: RWC ~ P60.sml + VPD.zscore5 + (P60.sml | ID) + (1 | species/plot)
getCall(m7.pv,1) #this function extracts the best model from the desired global model
summary(lmer(formula = RWC ~ P60.sml + VPD.zscore5 + (P60.sml | ID) + (1 | species/plot), data = d.all.complete))

#Testing different nesting structures for random effects
#Nesting structure of species/plot/ID results in perfect correlation of Species:Plot with P60.sml. To remove this we test different random effects structures. 
m7.pv.re0 <- lmer(formula = RWC ~ P60.sml + VPD.zscore5 + (P60.sml | species/plot/ID), data = d.all.complete) #perfect correlation of plot:species with P60.sml, thus not included
m7.pv.re1 <- lmer(RWC ~  P60.sml + VPD.zscore5 + (P60.sml | species/ID) + (P60.sml | plot), data=d.all.complete)
m7.pv.re2 <- lmer(RWC ~  P60.sml + VPD.zscore5 + (P60.sml | species/ID) + (1 | plot), data=d.all.complete)
m7.pv.re3 <- lmer(RWC ~  P60.sml + VPD.zscore5 + (P60.sml | ID) + (1 | species) + (1 | plot), data=d.all.complete)
m7.pv.re4 <- lmer(RWC ~  P60.sml + VPD.zscore5 + (P60.sml | ID) + (1 | species/plot), data=d.all.complete) #This model removes correlation of the random effects, minimizes AIC and results in most significant fixed effect precipitation response

#AIC comparison of best model with different random effect structures
AIC(m7.pv.re1,m7.pv.re2,m7.pv.re3,m7.pv.re4)