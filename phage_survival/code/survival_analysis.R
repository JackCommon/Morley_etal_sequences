#### Morley et al - Phage survival analysis ####
# Created: 8/5/18 by Jack Common

rm(list=ls())

#### Dependencies ####
#install.packages("survival")
#install.packages("rms")
#install.packages("car")
#install.packages("multcomp")
#install.packages("relaimpo")
#install.packages("dplyr")
#install.packages("magrittr")

library(survival)
library(rms)
library(car)
library(multcomp)
library(relaimpo)
library(dplyr)
library(magrittr)


#### Keplan-Meier ####

# Load data and attach the dataframe
phage<-read.csv("./phage_survival/analysis_data/phage_surv.csv", header=T)
phage$replicate %<>% as.factor()
phage$treatment %<>% as.factor()
attach(phage)
names(phage)

# KM ~ group
# Does the KM analysis and builds/saves the plot
summary(KM<-survfit(Surv(time_to_death,status)~treatment))

png("./figs/survplot.png", width=20, height=15, units="in", res=300)
par(mfrow=c(1,1), xpd=TRUE, oma=c(1,1,1,1), mai=c(1.02,.1,.82,0), bty="l", pty="s")

plot(survfit(Surv(phage$time_to_death,phage$status)~treatment), lty=c(1,3,5,6), 
     lwd=c(5,5,5,5), ylab="", xlab="", axes=FALSE, ylim=c(0,1), xlim=c(0,30))

axis(1, tcl=-0.1, pos=0, cex.axis=1, lwd=c(3), cex.axis=2)
axis(1, at=15, lab="Transfer", tcl=0, line=2, cex.axis=3)

axis(2, tcl=-0.1, pos=0, cex.axis=1, las=2, lwd=c(3), cex.axis = 2)
axis(2, at=0.5, lab="Proportion of replicates\nwith surviving phage", line=5, cex.axis=3, tcl=0)

legend(20,1, bty="o", title=c("Treatment"),
       legend=c(expression("10"*{}^{9}*""), 
                expression("10"*{}^{8}*""), 
                expression("10"*{}^{7}*""), 
                expression("10"*{}^{6}*"")),
       lty=c(6,5,3,1), lwd=c(5,5,5,5), cex=3, adj=0)

dev.off()


#### Cox PH model ####
cosph.mod<-coxph(Surv(time_to_death,status)~treatment)
summary(cosph.mod)
cosph.mod$loglik

anova(cosph.mod)
tapply(predict(cosph.mod),treatment,mean)

# Builds a Tukey comparison table that is copied to the clipboard (if you're on Mac OS)
# for easier copying to Excel
tukey <- summary(glht(cosph.mod, linfct = mcp(treatment = "Tukey")))
HRs <- exp(tukey$test$coefficients)
SEs <- exp(tukey$test$sigma)
Z <- tukey$test$tstat
P <- tukey$test$pvalues
HRs <- data.frame(HRs, SEs, Z, P)
clip = pipe('pbcopy', 'w')
write.table(HRs, file=clip, sep='\t', row.names = F, col.names = F)
close(clip)

# An alternative KM plot based on the Cox PH model
# Not entirely sure why it shows three lines?
plot(survfit(cosph.mod), xlim=c(0,30), lty=c(1,2,3,4))
