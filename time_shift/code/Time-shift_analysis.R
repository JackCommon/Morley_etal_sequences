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

#### Data Wrangling ####
# Dan's original data is formatted as a matrix, which R will struggle with. 
# Therefore, that data needs to be converted into an R-readable long-format

# Step 1: Melt each wide-formatted matrix into a long-format dataframe and joing
# them based on host timepoint

# First make a vector for the phage timepoints
Phage.Timepoint <- c(rep("t1", 432), 
                     rep("t4", 432), 
                     rep("t9", 432)) %>% 
  rep(7)

# Replicate 2.1
two.one.pT1 <- read.csv("./time_shift/original_data/2.1/2.1_phageT1.csv", header=T) 
two.one.pT1 <- melt(two.one.pT1, id.vars = c("Replicate", "Phage", "Host.Environment", "Host.Timepoint"))

two.one.pT4 <- read.csv("./time_shift/original_data/2.1/2.1_phageT4.csv", header=T) 
two.one.pT4 <- melt(two.one.pT4, id.vars = c("Replicate", "Phage", "Host.Environment", "Host.Timepoint"))

two.one.pT9 <- read.csv("./time_shift/original_data/2.1/2.1_phageT9.csv", header=T) 
two.one.pT9 <- melt(two.one.pT9, id.vars = c("Replicate", "Phage","Host.Environment", "Host.Timepoint"))

two.one <- bind_rows(two.one.pT1, two.one.pT4, two.one.pT9)

Environment <- rep(two.one$Host.Environment, 7) # Change times to 8 when 2.5 data comes in

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
#two.five.pT1 <- read.csv("./time_shift/original_data/2.5/2.5_phageT1.csv", header=T) 
#two.five.pT1 <- melt(two.five.pT1, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

#two.five.pT4 <- read.csv("./time_shift/original_data/2.5/2.5_phageT4.csv", header=T) 
#two.five.pT4 <- melt(two.five.pT4, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

#two.five.pT9 <- read.csv("./time_shift/original_data/2.5/2.5_phageT9.csv", header=T) 
#two.five.pT9 <- melt(two.five.pT9, id.vars = c("Replicate", "Phage", "Host.Timepoint"))

#two.five <- bind_rows(two.five.pT1, two.five.pT4, two.five.pT9)

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
data <- bind_rows(two.one, two.two, two.three,
                  two.four, two.six, two.seven, two.eleven)

# Add the inverse of the infectivity data to get resistance data
data$Resistant <- ifelse(data$variable==1,0,1)

# Add in the new variables for downstream analysis
data$Phage.Timepoint <- as.factor(Phage.Timepoint)
data$Environment <- as.factor(Environment)

# Rename some variables so they make more visual sense
data <- plyr::rename(data, c("variable"="Host.Genotype", 
                             "value"="Infected",
                             "Phage"="Phage.Genotype"))
data$Replicate %<>% as.factor()
data$Phage.Genotype %<>% as.factor()

# Finally, save the data as a CSV to the working directory
write.csv(file="./time_shift/original_data/infectivity_resistance_long.csv", data, 
          row.names = F)

# Clear up intermediate dataframes
rm(two.one, two.one.pT1, two.one.pT4, two.one.pT9,
   two.two, two.two.pT1, two.two.pT4, two.two.pT9,
   two.three, two.three.pT1, two.three.pT4, two.three.pT9,
   two.four, two.four.pT1, two.four.pT4, two.four.pT9,
   #two.five, two.five.pT1, two.five.pT4, two.five.pT9,
   two.six, two.six.pT1, two.six.pT4, two.six.pT9,
   two.seven, two.seven.pT1, two.seven.pT4, two.seven.pT9,
   two.eleven, two.eleven.pT1, two.eleven.pT4, two.eleven.pT9,
   Environment)



#### Analysis - GLM ####

# A function to quickly convert logit coefficients from a binomial GLM(M) 
# into more intuitive probability values
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

names(data)

# First build a GLM that tests the interaction between phage genotype and host genotype
# on average infectivity. Because the only the fixed effects are of interest for this 
# part of the analysis, 

m1 <- glm(Infected~Phage.Timepoint*Host.Timepoint,
            data=data,
            family=binomial(link="logit"))
summary(m1)
anova(m1, test="Chisq")
model.tables(aov(m1), "mean")

data$Host.Timepoint %<>% relevel(ref="t9")
data$Host.Timepoint %<>% relevel(ref="t4")
data$Host.Timepoint %<>% relevel(ref="t1")

data$Phage.Timepoint %<>% relevel(ref="t9")
data$Phage.Timepoint %<>% relevel(ref="t4")
data$Phage.Timepoint %<>% relevel(ref="t1")

logit2prob(confint(m1))

#### Analysis - GLMMs of all data ####
# Host environment-only model
# Slope does not vary with respect to phage genotype
m2 <- glmer(Infected~Environment+(1|Environment),
            data=data,
            family=binomial())
summary(m2)
anova(m2, test="Chisq")

# Overall genotype x Environment model
# Slope varies for each phage genotype as a random effect

m3 <- glmer(Infected~Environment+(Environment|Phage.Genotype),
            data=data,
            family=binomial(link="logit"))
summary(m3)
anova(m3, test="Chisq")

# Calculate the relative importance of FSD to ARD by calculating the 
# ratio between the GxE mean square and the E mean square
CoEvoRatio <- function(model1, model2){
  f1 <- model1@call$formula
  f2 <- model2@call$formula
  Env.MS <- anova(model1)$`Mean Sq`
  GE.MS <- anova(model2)$`Mean Sq`
  Ratio <- GE.MS/Env.MS
  cat("Relative importance of FSD:ARD:", Ratio)
}

CoEvoRatio(m2,m3)

#### Analysis - Timepoint-specific E and GxE GLMMs ####
# As before, first just model the Environment as a fixed effect
m4 <- glmer(Infected~Environment+(1|Environment),
                  data=subset(data, Host.Timepoint="t9"),
                  family=binomial())
summary(m4)
anova(m4, test="Chisq")

# Then model the GxE interaction with phage genotype as a random effect
m5 <- glmer(Infected~Environment+(Environment|Phage.Genotype),
            data=subset(data, Host.Timepoint=="t9"),
            family=binomial())
summary(m5)
anova(m5, test="Chisq")

CoEvoRatio(m4, m5)
logit2prob(confint(m5))

#### Figures ####
## Infectivity summary figure
infect_sum <- read.csv("./time_shift/summary_data/infectivity_means.csv")
infect_plot <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=infect_sum)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0.2, size=1)+
  geom_path(stat="identity")+
  theme_bw()+
  labs(x="Host environment (days post-infection)", y="Mean phage infectivity")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\ngenotype",
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

ggsave("infectivity_1.png", infect_plot, path="./figs/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

## Timepoint contrast 2 summary figure
infect_sum_2 <- read.csv("./time_shift/summary_data/infectivity_means_2.csv")
infect_sum_2$Phage %<>% relevel(ref="Future")
infect_sum_2$Phage %<>% relevel(ref="Contemporary")
infect_sum_2$Phage %<>% relevel(ref="Past")

infect_plot_2 <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=infect_sum_2)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0.2, size=1)+
  theme_bw()+
  
  labs(x="Host environment (days post-infection)", y="Mean phage infectivity")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\ngenotype")+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=16))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))
infect_plot_2

ggsave("infectivity_2.png", infect_plot_2, path="./figs/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))

### Timepoint contrast figure
coevo_means <- read.csv("./time_shift/summary_data/comparison_means.csv")

coevo_means$Environment %<>% relevel(ref="Future")
coevo_means$Environment %<>% relevel(ref="Contemporary")
coevo_means$Environment %<>% relevel(ref="Past")

coevo_plot <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=coevo_means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0.1, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_bw()+
  labs(x="Host environment", y="Mean phage infectivity")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 0.7, 0.1)))
coevo_plot

ggsave("coevo_1.png", coevo_plot, path="./figs/",
       device="png", dpi=300, width=20, height=14, units = c("cm"))


