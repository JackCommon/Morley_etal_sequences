### Sequence PWD analysis
## Created 25/4/18 by Jack Common

rm(list=ls())

# Packags
library(dplyr)
library(ggplot2)
library(scales)
library(magrittr)

# Functions
model_stats = function(model){
  sum = coef(model)
  conf = confint(model, level=c(0.95))
  print(c(sum[1]))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  
  stats = data.frame(conf[1,1], conf[1,2])
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Coefficients copied to the clipboard')
  
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

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

# Data
diversity <- read.csv("./sequences/summary_data/PWD.csv", header = T)
diversity$Replicate %<>% as.factor
diversity %<>% select(-PW_ID)

diversity$Replicate %<>% relevel(ref="2.11")
diversity$Replicate %<>% relevel(ref="2.7")
diversity$Replicate %<>% relevel(ref="2.6")
diversity$Replicate %<>% relevel(ref="2.5")
diversity$Replicate %<>% relevel(ref="2.4")
diversity$Replicate %<>% relevel(ref="2.3")
diversity$Replicate %<>% relevel(ref="2.2")
diversity$Replicate %<>% relevel(ref="2.1")

# Barplots

PWD_plot <- ggplot(aes(y=PWD, x=Replicate, group=Locus), data=diversity)+
  geom_bar(stat="identity", aes(fill=Locus))+
  facet_grid(~Timepoint)+
  theme_bw()
PWD_plot

mod.1 <- glm(PWD~Timepoint, data = diversity, family=binomial(link="sqrt"))
summary(mod.1)

