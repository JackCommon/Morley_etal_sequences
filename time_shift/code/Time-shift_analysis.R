#### Morley et al - Time-shift analysis ####
# Created: 8/5/18 by Jack Common

rm(list=ls())

#### Dependencies ####
#install.packages("nlme")
#install.packages("lme4")
#install.packages("MuMIn")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("magrittr")

library(nlme)
library(lme4)
library(MuMIn)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)

#### Functions ####
# A function to quickly convert logit coefficients from a binomial GLM(M) 
# into more intuitive probability values
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

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


#### Data Wrangling ####
# Dan's original data is formatted as a matrix, which R will struggle with. 
# Therefore, that data needs to be converted into an R-readable long-format

# Step 1: Melt each wide-formatted matrix into a long-format dataframe and joing
# them based on host timepoint

# First make a vector for the phage timepoints
Phage.Timepoint <- c(rep("t1", 432), 
                     rep("t4", 432), 
                     rep("t9", 432)) %>% 
  rep(8)

# Replicate 2.1
two.one.pT1 <- read.csv("./time_shift/original_data/2.1/2.1_phageT1.csv", header=T) 
two.one.pT1 <- melt(two.one.pT1, id.vars = c("Replicate", "Phage", "Host.Environment", "Host.Timepoint"))

two.one.pT4 <- read.csv("./time_shift/original_data/2.1/2.1_phageT4.csv", header=T) 
two.one.pT4 <- melt(two.one.pT4, id.vars = c("Replicate", "Phage", "Host.Environment", "Host.Timepoint"))

two.one.pT9 <- read.csv("./time_shift/original_data/2.1/2.1_phageT9.csv", header=T) 
two.one.pT9 <- melt(two.one.pT9, id.vars = c("Replicate", "Phage","Host.Environment", "Host.Timepoint"))

two.one <- bind_rows(two.one.pT1, two.one.pT4, two.one.pT9)

Environment <- rep(two.one$Host.Environment, 8) 

# Replicate 2.2
two.two.pT1 <- read.csv("./time_shift/original_data/2.2/2.2_phageT1.csv", header=T) 
two.two.pT1 <- melt(two.two.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.two.pT4 <- read.csv("./time_shift/original_data/2.2/2.2_phageT4.csv", header=T) 
two.two.pT4 <- melt(two.two.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.two.pT9 <- read.csv("./time_shift/original_data/2.2/2.2_phageT9.csv", header=T) 
two.two.pT9 <- melt(two.two.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.two <- bind_rows(two.two.pT1, two.two.pT4, two.two.pT9)

# Replicate 2.3
two.three.pT1 <- read.csv("./time_shift/original_data/2.3/2.3_phageT1.csv", header=T) 
two.three.pT1 <- melt(two.three.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.three.pT4 <- read.csv("./time_shift/original_data/2.3/2.3_phageT4.csv", header=T) 
two.three.pT4 <- melt(two.three.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.three.pT9 <- read.csv("./time_shift/original_data/2.3/2.3_phageT9.csv", header=T) 
two.three.pT9 <- melt(two.three.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.three <- bind_rows(two.three.pT1, two.three.pT4, two.three.pT9)

# Replicate 2.4
two.four.pT1 <- read.csv("./time_shift/original_data/2.4/2.4_phageT1.csv", header=T) 
two.four.pT1 <- melt(two.four.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.four.pT4 <- read.csv("./time_shift/original_data/2.4/2.4_phageT4.csv", header=T) 
two.four.pT4 <- melt(two.four.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.four.pT9 <- read.csv("./time_shift/original_data/2.4/2.4_phageT9.csv", header=T) 
two.four.pT9 <- melt(two.four.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.four <- bind_rows(two.four.pT1, two.four.pT4, two.four.pT9)

# Replicate 2.5
two.five.pT1 <- read.csv("./time_shift/original_data/2.5/2.5_phageT1.csv", header=T) 
two.five.pT1 <- melt(two.five.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.five.pT4 <- read.csv("./time_shift/original_data/2.5/2.5_phageT4.csv", header=T) 
two.five.pT4 <- melt(two.five.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.five.pT9 <- read.csv("./time_shift/original_data/2.5/2.5_phageT9.csv", header=T) 
two.five.pT9 <- melt(two.five.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.five <- bind_rows(two.five.pT1, two.five.pT4, two.five.pT9)

# Replicate 2.6
two.six.pT1 <- read.csv("./time_shift/original_data/2.6/2.6_phageT1.csv", header=T) 
two.six.pT1 <- melt(two.six.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.six.pT4 <- read.csv("./time_shift/original_data/2.6/2.6_phageT4.csv", header=T) 
two.six.pT4 <- melt(two.six.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.six.pT9 <- read.csv("./time_shift/original_data/2.6/2.6_phageT9.csv", header=T) 
two.six.pT9 <- melt(two.six.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.six <- bind_rows(two.six.pT1, two.six.pT4, two.six.pT9)

# Replicate 2.7
two.seven.pT1 <- read.csv("./time_shift/original_data/2.7/2.7_phageT1.csv", header=T) 
two.seven.pT1 <- melt(two.seven.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.seven.pT4 <- read.csv("./time_shift/original_data/2.7/2.7_phageT4.csv", header=T) 
two.seven.pT4 <- melt(two.seven.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.seven.pT9 <- read.csv("./time_shift/original_data/2.7/2.7_phageT9.csv", header=T) 
two.seven.pT9 <- melt(two.seven.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.seven <- bind_rows(two.seven.pT1, two.seven.pT4, two.seven.pT9)

# Replicate 2.11
two.eleven.pT1 <- read.csv("./time_shift/original_data/2.11/2.11_phageT1.csv", header=T) 
two.eleven.pT1 <- melt(two.eleven.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.eleven.pT4 <- read.csv("./time_shift/original_data/2.11/2.11_phageT4.csv", header=T) 
two.eleven.pT4 <- melt(two.eleven.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.eleven.pT9 <- read.csv("./time_shift/original_data/2.11/2.11_phageT9.csv", header=T) 
two.eleven.pT9 <- melt(two.eleven.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

two.eleven <- bind_rows(two.eleven.pT1, two.eleven.pT4, two.eleven.pT9)

# Bind together each replicate dataframe
data <- bind_rows(two.one, two.two, two.three, two.four,
                  two.five, two.six, two.seven, two.eleven)

# Add in the new variables for downstream analysis
data$Phage.Timepoint <- as.factor(Phage.Timepoint)
data$Environment <- as.factor(Environment)

# Rename some variables so they make more visual sense
data <- plyr::rename(data, c("variable"="Host.Genotype", 
                             "value"="Infected",
                             "Phage"="Phage.Genotype"))
data$Replicate %<>% as.factor()
data$Phage.Genotype %<>% as.factor()

# Add the inverse of the infectivity data to get resistance data
data$Resistant <- ifelse(data$Infected=="1",0,1)

# Sort out the Host and Phage environment factors
data %<>% select(-Host.Environment)
data %<>% plyr::rename(c("Environment"="Host.Background"))
data$Phage.Background <- rep(NA, length(data$Host.Background))
for(i in seq(1,length(data$Host.Background),1)){
  if(data$Host.Background[i]=="Past"){
    data$Phage.Background[i]<-"Future"
  } 
}
for(i in seq(1,length(data$Host.Background),1)){
  if(data$Host.Background[i]=="Contemporary"){
    data$Phage.Background[i]<-"Contemporary"
  } 
}
for(i in seq(1,length(data$Host.Background),1)){
  if(data$Host.Background[i]=="Future"){
    data$Phage.Background[i]<-"Past"
  } 
}

# Host.Background -> If phage are challenged against hosts from their past, present or future
# Phage.Background -> If hosts are challenged against phage from their past, present or future

# Finally, save the data as a CSV to the working directory
write.csv(file="./time_shift/original_data/infectivity_resistance_long.csv", data, 
          row.names = F)

# Clear up intermediate dataframes
rm(two.one, two.one.pT1, two.one.pT4, two.one.pT9,
   two.two, two.two.pT1, two.two.pT4, two.two.pT9,
   two.three, two.three.pT1, two.three.pT4, two.three.pT9,
   two.four, two.four.pT1, two.four.pT4, two.four.pT9,
   two.five, two.five.pT1, two.five.pT4, two.five.pT9,
   two.six, two.six.pT1, two.six.pT4, two.six.pT9,
   two.seven, two.seven.pT1, two.seven.pT4, two.seven.pT9,
   two.eleven, two.eleven.pT1, two.eleven.pT4, two.eleven.pT9,
   Environment)


#### Analysis - GLM ####
# Reload the data 
data <- read.csv("./time_shift/original_data/infectivity_resistance_long.csv")
data$Replicate %<>% as.factor()
data$Phage.Genotype %<>% as.factor()
#data %<>% na.exclude # include to get marginal and conditional R2 from models
names(data)

# First build GLMs that tests the interaction between phage genotype and host genotype
# on average infectivity. Only the fixed effects are of interest for this 
# part of the analysis 

m.null <- glm(Infected~1,
            data=data,
            family=binomial(link="logit"))
par(mfrow=c(2,2))
plot(m1)

m1 <- glm(Infected~Host.Timepoint,
              data=data,
              family=binomial(link="logit"))
par(mfrow=c(2,2))
plot(m1)

m2 <- glm(Infected~Host.Timepoint*Phage.Timepoint,
          data=data,
          family=binomial(link="logit"))
par(mfrow=c(2,2))
plot(m2)

summary(m.null)
summary(m1)
summary(m2)

anova(m.null, m1, m2, test="Chisq")
logLik(m.null)
logLik(m1)
logLik(m2)

AIC(m.null, m1, m2) %>% compare_AICs()

anova(m2, test="Chisq")
model.tables(aov(m2), "mean")

data$Host.Timepoint %<>% relevel(ref="t9")
data$Host.Timepoint %<>% relevel(ref="t4")
data$Host.Timepoint %<>% relevel(ref="t1")

data$Phage.Timepoint %<>% relevel(ref="t9")
data$Phage.Timepoint %<>% relevel(ref="t4")
data$Phage.Timepoint %<>% relevel(ref="t1")

data$Host.Background %<>% relevel(ref="Past")
data$Host.Background %<>% relevel(ref="Contemporary")
data$Host.Background %<>% relevel(ref="Future")

m2 <- glm(Infected~Host.Timepoint*Phage.Timepoint,
          data=data,
          family=binomial(link="logit"))

logit2prob(m2$coefficients[1])
logit2prob(confint(m2))

#### Analysis - GLMMs of infectivity based on host background ####
# Host environment-only model
# Slope does not vary with respect to phage genotype
m1 <- glmer(Infected~Host.Background+(1|Host.Background),
            data=data,
            family=binomial())
summary(m1)
par(mfrow=c(2,2))
plot(m1)

# Genotype as a single random effect with no interaction
m2 <- glmer(Infected~Host.Background+(1|Phage.Genotype),
               data=data,
               family=binomial())
summary(m2)
plot(m2)

anova(m1,m2, test="Chisq")
anova(m2, test="Chisq")
R2 <- r.squaredGLMM(m2)
R2[1]/R2[2]*100

# Overall genotype x Environment model
# Slope varies for each phage genotype as a random effect

m3 <- glmer(Infected~Host.Background+(Host.Background|Phage.Genotype),
            data=data,
            family=binomial(link="logit"))
summary(m3)
anova(m3, test="Chisq")

plot(m3)

anova(m1, m2, m3, test="Chisq")
logLik(m1)
logLik(m2)
logLik(m3)
AIC(m1, m2, m3) %>% compare_AICs()
# Although the heteroskedacity can't seem to shift, I'll move ahead with model 2 based on the anova,
# log-likelihood and AIC comparisons
summary(m2)
logit2prob(fixef(m2)[[1]])
logit2prob(fixef(m2)[[1]]+fixef(m2)[[2]])
logit2prob(fixef(m2)[[1]]+fixef(m2)[[3]])

CIs <- confint(m2, parm="beta_")
CIs
logit2prob(CIs[1]+CIs[2])
logit2prob(CIs[1]+CIs[3])

logit2prob(CIs[4]+CIs[5])
logit2prob(CIs[4]+CIs[6])

#### Analysis - GLMMs of infectivity based on host background ####
# Host environment-only model
# Slope does not vary with respect to phage genotype
m1 <- glmer(Resistant~Phage.Background+(1|Phage.Background),
            data=data,
            family=binomial())
summary(m1)
par(mfrow=c(2,2))
plot(m1)

# Genotype as a single random effect with no interaction
m2 <- glmer(Resistant~Phage.Background+(1|Host.Genotype),
            data=data,
            family=binomial())
summary(m2)
plot(m2)

anova(m1,m2, test="Chisq")
r.squaredGLMM(m2)

# Overall genotype x Environment model
# Slope varies for each phage genotype as a random effect

m3 <- glmer(Resistant~Phage.Background+(Phage.Background|Host.Genotype),
            data=data,
            family=binomial(link="logit"))
summary(m3)
anova(m3, test="Chisq")

plot(m3)

anova(m1, m2, m3, test="Chisq")
logLik(m1)
logLik(m2)
logLik(m3)
AIC(m1, m2, m3) %>% compare_AICs()
# Although the heteroskedacity can't seem to shift, I'll move ahead with model 2 based on the anova,
# log-likelihood and AIC comparisons
summary(m2)
logit2prob(fixef(m2)[[1]])
logit2prob(fixef(m2)[[1]]+fixef(m2)[[2]])
logit2prob(fixef(m2)[[1]]+fixef(m2)[[3]])

CIs <- confint(m2, parm="beta_")
CIs
logit2prob(CIs[1]+CIs[2])
logit2prob(CIs[1]+CIs[3])

logit2prob(CIs[4]+CIs[5])
logit2prob(CIs[4]+CIs[6])

#### Analysis - Timepoint-specific E and GxE GLMMs ####
data.temp <- filter(data, Host.Timepoint=="t4")
m4 <- glmer(Resistant~Host.Background+(1|Phage.Genotype),
            data=data.temp,
            family=binomial())
summary(m4)
logit2prob(fixef(m4)[[1]])
logit2prob(fixef(m4)[[1]]+fixef(m4)[[2]])
logit2prob(fixef(m4)[[1]]+fixef(m4)[[3]])

CIs <- confint(m4, parm="beta_")
CIs
logit2prob(CIs[1]); logit2prob(CIs[4])
logit2prob(CIs[1]+CIs[2])
logit2prob(CIs[1]+CIs[3])

logit2prob(CIs[4]+CIs[5])
logit2prob(CIs[4]+CIs[6])

#### Figures ####
## Infectivity summary figure ####
infect_sum <- read.csv("./time_shift/summary_data/timepoint_contrasts.csv")
infect_plot <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=infect_sum)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity")+
  theme_bw()+
  labs(x="Host background", y="Proportion of hosts infected")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\nbackground",
                      breaks=c("t1", "t4", "t9"),
                      labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))

infect_plot

ggsave("infectivity_1.png", infect_plot, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

## Time-shift by timepoint figure ####
infect_sum_2 <- read.csv("./time_shift/summary_data/timepoint_timeshift_contrasts.csv")
infect_sum_2$Phage %<>% relevel(ref="Future")
infect_sum_2$Phage %<>% relevel(ref="Present")
infect_sum_2$Phage %<>% relevel(ref="Past")

infect_plot_2 <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=infect_sum_2)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  theme_bw()+
  
  labs(x="Timepoint", y="Infectivity")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Host\nbackground")+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))
infect_plot_2

ggsave("infectivity_2.png", infect_plot_2, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

### Time-shift infectivity contrast figure ####
coevo_infect <- read.csv("./time_shift/summary_data/timeshift_means.csv")

coevo_infect$Host.Background %<>% relevel(ref="Future")
coevo_infect$Host.Background %<>% relevel(ref="Present")
coevo_infect$Host.Background %<>% relevel(ref="Past")

coevo_infect_plot <- ggplot(aes(y=Mean.Infect, x=Host.Background, group=Group), data=coevo_infect)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_bw()+
  labs(x="Host background", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 0.8, 0.1)))
coevo_infect_plot

ggsave("coevo_infect.png", coevo_infect_plot, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

#### Resistance figures ####
## Timeshift resistance figure ####
resist_plot <- ggplot(aes(y=Mean.Resist, x=Host, Group=Phage), data=infect_sum)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity")+
  theme_bw()+
  labs(x="Host background", y="Resistance")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\nbackground",
                      breaks=c("t1", "t4", "t9"),
                      labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))

resist_plot

ggsave("resistance_1.png", resist_plot, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

## Resistance timeshift by timepoint figure ####

resist_plot_2 <- ggplot(aes(y=Mean.Resist, x=Host, Group=Phage), data=infect_sum_2)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  theme_bw()+
  
  labs(x="Transfer", y="Resistance")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Host\nbackground")+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))
resist_plot_2

ggsave("resistance_2.png", resist_plot_2, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

### Time-shift resistance figure ####
coevo_resist_plot <- ggplot(aes(y=Mean.Resist, x=Host.Background, group=Group), data=coevo_infect)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_bw()+
  labs(x="Host background", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,1))
coevo_resist_plot

ggsave("coevo_resist.png", coevo_resist_plot, path="./figs/coevo/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))


