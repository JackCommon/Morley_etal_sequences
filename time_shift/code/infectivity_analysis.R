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
data <- read.csv("./time_shift/original_data/infectivity_resistance_long.csv")
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor

#### Analysis - Proportion of hosts infected by each phage genotype ####
# Make subset dataframes of t1, t4 and t9
t1 <- filter(data, Host.Timepoint=="t1")
t4 <- filter(data, Host.Timepoint=="t4")
t9 <- filter(data, Host.Timepoint=="t9")

# Then run a GLM to measure proportion of hosts infected by each phage genotype
# Change both the data and subset argument for timepoint and replicate, respecitively
m1 <- glm(Infected~Phage.Genotype,
          data=t9, subset=c(Replicate=="2.11"),
          family = binomial())
m2 <- glm(Resistant~Host.Genotype,
          data=t9, subset=c(Replicate=="2.2"),
          family = binomial())

# Copy the coefficients to the clipboard
phage_stats(m1)
host_stats(m2)

#### Load in replicate data parts and bind ####
one <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.1.csv")
two <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.2.csv")
three <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.3.csv")
four <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.4.csv")
five <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.5.csv")
six <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.6.csv")
seven <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.7.csv")
eleven <- read.csv("./time_shift/summary_data/infectivity_data_parts/2.11.csv")

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
write.csv(file="./time_shift/summary_data/full_infect_resist_bygenotype.csv", data,
          row.names = F)

#### Analysis - Mean infectivity ####
data <- read.csv("./time_shift/summary_data/full_infect_resist_bygenotype.csv")
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor
data$Host.Genotype %<>% as.factor

# Nested models with Timepoint as a fixed effect. Using a negative binomial family
# to account for zero inflation
m1 <- glmer.nb(Infected~Timepoint+(1|Timepoint),
            data=data,
            family=binomial())

m2 <- glmer.nb(Infected~Timepoint+(1|Replicate),
            data=data,
            family=binomial())

m3 <- glmer.nb(Infected~Timepoint+(Timepoint|Replicate),
            data=data,
            family=binomial())

# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)
# Maybe less fanning?
plot(m3)

# Compare AICs
AICs <- AIC(m1, m2, m3) %>% compare_AICs()
# Model 2, with replicate as a single random effect, has the lowest AIC

# ANOVA to compare effect contributions
anova(m1, m2, m3, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# On balance and despite the increased heteroskedacity in the second model,
# I'll use this one moving forwards due to the AIC and Chisq tests

# Get the model coefficients and confidence intervals
summary(m2)
# Coefs
logit2prob(-0.2551)
logit2prob(-0.2551-0.5364)
logit2prob(-0.2551-2.0177)

#### Analysis - Mean resistance ####

