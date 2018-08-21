#### Morley et al - infectivity analysis ####
# Created: 10/5/18 by Jack Common
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

phage_stats = function(model){
  means = model.tables(aov(model), "means")
  genotype = c()
  coefs = c()
  for ( i in seq(1,12)){
    genotype[i] <- i
    coefs[i] <- means$tables$Phage.Genotype[i]
  }
  df <- data.frame(coefs)
  
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Phage infectivity coefficients copied to the clipboard')
  
}

host_stats = function(model){
  means = model.tables(aov(model), "means")
  genotype = c()
  coefs = c()
  for ( i in seq(1,12)){
    genotype[i] <- i
    coefs[i] <- means$tables$Host.Genotype[i]
  }
  df <- data.frame(coefs)
  
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Host resistance coefficients copied to the clipboard')
  
}

CopyCIs <- function(model){
  conf = logit2prob(confint(model, level=c(0.95)))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  CIs <- data.frame(conf[1,1], conf[1,2])
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(CIs, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Probability 95% CIs copied to the clipboard')
  
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

#### Data ####
# Import the full infection matrix dataset
data <- read.csv("./infectivity/original_data/infectivity_resistance_long.csv")
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor

#### Analysis - Proportion of hosts infected by each phage genotype ####
# Make subset dataframes of t1, t4 and t9
t1 <- filter(data, Host.Timepoint=="t1")
t4 <- filter(data, Host.Timepoint=="t4")
t9 <- filter(data, Host.Timepoint=="t9")

t1 <- filter(data, Phage.Timepoint=="t1")
t4 <- filter(data, Phage.Timepoint=="t4")
t9 <- filter(data, Phage.Timepoint=="t9")

# Then run a GLM to measure proportion of hosts infected by each phage genotype
# Change both the data and subset argument for timepoint and replicate, respecitively
m1 <- glm(Infected~Phage.Genotype,
          data=t1, subset=c(Replicate=="2.1"),
          family = binomial())
m2 <- glm(Resistant~Host.Genotype,
          data=t9, subset=c(Replicate=="2.7"),
          family = binomial())

# Copy the coefficients to the clipboard
phage_stats(m1)
host_stats(m2)
model.tables(aov(m2), "mean")

#### Load in replicate data parts and bind ####
one <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.1.csv")
two <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.2.csv")
three <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.3.csv")
four <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.4.csv")
five <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.5.csv")
six <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.6.csv")
seven <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.7.csv")
eleven <- read.csv("./infectivity/summary_data/infectivity_data_parts/2.11.csv")

# Bind and organise
data <- bind_rows(one, two, three, four, five, six, seven, eleven)

data$replicate %<>% as.factor()
data$phage %<>% as.factor()
data$host %<>% as.factor

# Rename some variables so they make more visual sense
data <- plyr::rename(data, c("replicate"="Replicate", 
                             "timepoint"="Timepoint",
                             "phage"="Phage.Genotype",
                             "host" = "Host.Genotype",
                             "infectivity"="Infected",
                             "resistance"="Resisted"))


# Save the bound data 
write.csv(file="./infectivity/summary_data/full_infect_resist_bygenotype.csv", data,
          row.names = F)

#### Analysis - Mean infectivity ####
data <- read.csv("./infectivity/summary_data/full_infect_resist_bygenotype.csv")
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor
data$Host.Genotype %<>% as.factor

# Nested models with Timepoint as a fixed effect
m1 <- glmer(Infected~Timepoint+(1|Timepoint),
            data=data_drop,
            family=binomial(link="identity"))

m2 <- glmer(Infected~Timepoint+(1|Replicate),
            data=data_drop,
            family=binomial(link="logit"))


# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)

# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
# Model 2 has the lowest AIC

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m2, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# On balance and despite the increased heteroskedacity in the second model,
# I'll use this one moving forwards due to the AIC and Chisq tests

# Get the model coefficients and confidence intervals
summary(m2)
fixed <- fixef(m2)
# Coefs
fixed[1]
fixed[1]+fixed[2]
fixed[1]+fixed[3]

data$Timepoint %<>% relevel(ref="t4")
m2 <- glmer(Infected~Timepoint+(1|Replicate),
            data=data,
            family=binomial(link="identity"))

CIs <- confint(m2, method="Wald", parm="beta_")
CIs[1]
CIs[1]+CIs[2]
CIs[1]+CIs[3]
CIs[4]
CIs[4]+CIs[5]
CIs[4]+CIs[6]


#### Analysis - Mean resistance ####
# Nested models with Timepoint as a fixed effect. Using a negative binomial family
# to account for zero inflation
m1 <- glmer(Resisted~Timepoint+(1|Timepoint),
               data=data,
               family=binomial(link="logit"))

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
               data=data_drop,
               family=binomial(link="logit"))


# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)

# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
# Model 2, with replicate as a random effect, has the lowest AIC

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m2, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# Model 2 is probably the best
# Heteroskedacity in Model 2 suggests it is most appropriate
# Get the model coefficients and confidence intervals
summary(m2)
fixed <- fixef(m2)
fixed
logit2prob(fixed[1])
logit2prob(fixed[1])+logit2prob(fixed[2])
logit2prob(fixed[1])+logit2prob(fixed[3])

data$Timepoint %<>% relevel(ref="t9")
data$Timepoint %<>% relevel(ref="t4")
data$Timepoint %<>% relevel(ref="t1")

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
            data=data,
            family=binomial(link="logit"))
CIs <- confint(m2, parm="beta_")

logit2prob(CIs[1])
logit2prob(CIs[1])+logit2prob(CIs[2])
logit2prob(CIs[1])+logit2prob(CIs[3])

logit2prob(CIs[4])
logit2prob(CIs[4])+logit2prob(CIs[5])
logit2prob(CIs[4])+logit2prob(CIs[6])
#### Figures ####

means <- read.csv("./infectivity/summary_data/infect_resist_means.csv")

infect_plot <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  #geom_path(stat="identity", size=.8, linetype=2)+
  theme_bw()+
  labs(x="Transfer", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))
infect_plot

resist_plot <- ggplot(aes(y=Mean.Resist, x=Timepoint, group=Group), data=means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  #geom_path(stat="identity", size=.8, linetype=2)+
  theme_bw()+
  labs(x="Transfer", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))
resist_plot

library(cowplot)
infect_resist_evo <- plot_grid(infect_plot+labs(x=""), resist_plot,
                               nrow=2, align = "hv", labels = c("A", "B"))

infect_resist_evo

detach(package:cowplot)

ggsave("infect_resist_evolution.png", infect_resist_evo, path="./figs/evo/",
       device="png", dpi=300, width=16, height = 24, units=c("cm"))
  