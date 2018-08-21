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


#### Functions ####
## Compare AIC values for model fit
# This function extracts the AIC for each GLM, and then compares the absolute relative 
# differences for each AIC to the model with the lowest AIC. This acts as a measure of 
# model fit. More can be found at: http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}
#### Keplan-Meier ####

# Load data and attach the dataframe
phage<-read.csv("./phage_survival/original_data/phage_surv.csv", header=T)
phage$replicate %<>% as.factor()
phage$treatment %<>% as.factor()
attach(phage)
names(phage)

# KM ~ group
# Does the KM analysis and builds/saves the plot
summary(KM<-survfit(Surv(time_to_death,status)~treatment))

png("./figs/phage/survplot.png", width=20, height=15, units="in", res=300)
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
m.null <- coxph(Surv(time_to_death, status)~1)
m1 <- coxph(Surv(time_to_death, status)~treatment)

AIC(m.null, m1) %>% compare_AICs()

summary(m1)


anova(m1)
tapply(predict(m1),treatment,mean)

# Builds a Tukey comparison table that is copied to the clipboard (if you're on Mac OS)
# for easier copying to Excel
tukey <- summary(glht(m1, linfct = mcp(treatment = "Tukey")))
HRs <- exp(tukey$test$coefficients)
SEs <- exp(tukey$test$sigma)
Z <- tukey$test$tstat
P <- tukey$test$pvalues
HRs <- data.frame(HRs, SEs, Z, P)
clip = pipe('pbcopy', 'w')
write.table(HRs, file=clip, sep='\t', row.names = T, col.names = T)
close(clip)

# An alternative KM plot based on the Cox PH model
# Not entirely sure why it shows three lines?
plot(survfit(m1), xlim=c(0,30), lty=c(1,2,3,4))


#### Binomial survival analysis ####
#install.packages("ggfortify")
#install.packages("survminer")
library(ggfortify)
library(survminer)

data <- read.csv("./phage_survival/original_data/binomial_survival_data.csv", header = TRUE)
data$Treatment %<>% as.factor

m.null <- glm(cbind(Alive, Dead)~1,family=binomial, data=data)
m1 <- glm(cbind(Alive, Dead)~Time,family=binomial, data=data)
m2 <- glm(cbind(Alive, Dead)~Treatment,family=binomial, data=data)
m3 <- glm(cbind(Alive, Dead)~Treatment+Time,family=binomial, data=data)
m.global <- glm(cbind(Alive, Dead)~Treatment*Time,family=binomial, data=data)

par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m3)
plot(m.global)

AIC(m.null, m1, m2, m3, m.global) %>% compare_AICs()

summary(m.global)
anova(m.global,test="Chisq")
summary(glht(m.global, linfct = mcp(Treatment = "Tukey")))
confint(m.global)

range(data$Time)
n<-seq(0,30,length.out=1000)
length(n)
intercept_TR1<-m.global$coef[1]
intercept_TR2<-m.global$coef[1]+m.global$coef[2]
intercept_TR3<-m.global$coef[1]+m.global$coef[3]
intercept_TR4<-m.global$coef[1]+m.global$coef[4]

slope_TR1 <- m.global$coef[5]
slope_TR2 <- m.global$coef[5]+m.global$coef[6]
slope_TR3 <- m.global$coef[5]+m.global$coef[7]
slope_TR4 <- m.global$coef[5]+m.global$coef[8]

fitted_TR1<-exp(intercept_TR1+slope_TR1*n)/(1+exp(intercept_TR1+slope_TR1*n))
fitted_TR2<-exp(intercept_TR2+slope_TR2*n)/(1+exp(intercept_TR2+slope_TR2*n))
fitted_TR3<-exp(intercept_TR3+slope_TR3*n)/(1+exp(intercept_TR3+slope_TR3*n))
fitted_TR4<-exp(intercept_TR4+slope_TR4*n)/(1+exp(intercept_TR4+slope_TR4*n))

par(mfrow=c(1,1))

par(las=1)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

binom_plot <- recordPlot()

plot(data$Time, data$Alive/(data$Alive+data$Dead),axes=FALSE, type="n",ylab = "",xlab = "", xlim=c(0,30), cex.axis=2.5, cex.lab=2.5)
axis(side = 1, lwd = 3, cex.axis=2.5, pos = c(0,1), tck= -.01)
axis(side = 2, lwd = 3, cex.axis=2.5, pos = c(-.15,1))
title(ylab="Mean proportion of replicates\nwith surviving phage",line=3.5, cex.lab=2.5, tcl=0)
title(xlab="Days post-infection (d.p.i.)",line=4, cex.lab=2.5, tcl=0)

binom_plot <- recordPlot()
points(fitted_TR1~n, type="l", lwd=2.5,lty=1, col="#718B2E")
points(fitted_TR2~n, type="l", lwd=2.5,lty=1, col= "#D39200")
points(fitted_TR3~n, type="l", lwd=2.5,lty=1, col= "red")
points(fitted_TR4~n, type="l", lwd=2.5,lty=1, col= "blue")

binom_plot <- recordPlot()
legend(20, 0.9, c("6"),lty = c(1), col=c('#D39200'),lwd=2.5, cex=1.5,bty='n')
legend(20, 0.8, c("7"),lty = c(1), col=c('#718B2E'),lwd=2.5, cex=1.5,bty='n')
legend(20, 0.7, c("8"),lty = c(1), col=c('red'),lwd=2.5, cex=1.5,bty='n')
legend(20, 0.6, c("9"),lty = c(1), col=c('blue'),lwd=2.5, cex=1.5,bty='n')

points(data$Time[data$Treatment=="6"], data$Alive[data$Treatment=="6"]/(data$Alive[data$Treatment=="6"]
                                                                                    +data$Dead[data$Treatment=="6"]), pch=19, cex=1,col=c('#718B2E'))
points(data$Time[data$Treatment=="7"], data$Alive[data$Treatment=="7"]/(data$Alive[data$Treatment=="7"]
                                                                                    +data$Dead[data$Treatment=="7"]), pch=19, cex=1,col=c('#D39200'))
points(data$Time[data$Treatment=="8"], data$Alive[data$Treatment=="8"]/(data$Alive[data$Treatment=="8"]
                                                                        +data$Dead[data$Treatment=="8"]), pch=19, cex=1,col=c('red'))
points(data$Time[data$Treatment=="9"], data$Alive[data$Treatment=="9"]/(data$Alive[data$Treatment=="9"]
                                                                        +data$Dead[data$Treatment=="9"]), pch=19, cex=1,col=c('blue'))
