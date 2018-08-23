#### Phage protospacers sequence analysis 
## Created 20/8/18 by Jack Commo

rm(list=ls())

#### Dependencies ####
library(lme4)
library(MuMIn)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(broom)
library(cowplot)

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

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}
#### Data ####
data <- read.csv("./phage_sequences/summary_data/infectivity_protospacers_full.csv")

#data <- select(data, -Phage.Genotype, -Host.Timepoint, -Host.Genotype, -Phage.Timepoint, -Host.Background, -Phage.Background)

data$Replicate %<>% as.factor
#data$SpacerNumber %<>% as.factor
data$Phage.Genotype %<>% as.factor
data$Escape <- ifelse(is.na(data$Escape)==T, 0, data$Escape) %>% as.integer()
data$Other <- ifelse(is.na(data$Other)==T, 0, data$Other) %>% as.integer()
data$Protospacer <- ifelse(is.na(data$Protospacer)==T, 0, data$Protospacer) %>% as.integer()
data$SpacersTargetted <- ifelse(is.na(data$SpacersTargetted)==T, 0, data$SpacersTargetted) %>% as.integer()
#data$SpacersTargetted <- ifelse(is.na(data$SpacersTargetted)==T, 0, data$SpacersTargetted) %>% as.factor
data$NotTargetted <- data$SpacerNumber-data$SpacersTargetted
data$NotTargetted %<>% as.factor
data$Mutation <- ifelse(data$Escape==1,1,0) 
data$Mutation <- ifelse(data$Other==1,1,data$Mutation) %>% as.integer

infectious_success <- filter(data, NotTargetted==0) %>% 
  filter(Infected==1)
infectious_unsuccess <- filter(data, NotTargetted==0) %>% 
  filter(Infected==0)
not_success <- filter(data, NotTargetted!=0) %>% 
  filter(Infected==1)
not_unsuccess <- filter(data, NotTargetted!=0) %>% 
  filter(Infected==0)

infectious_success$Match <- c(rep("YY", length(infectious_success$Replicate))) %>% as.factor()
infectious_unsuccess$Match <- c(rep("YN", length(infectious_unsuccess$Replicate))) %>% as.factor
not_success$Match <- c(rep("NY", length(not_success$Replicate))) %>% as.factor
not_unsuccess$Match <- c(rep("NN", length(not_unsuccess$Replicate))) %>% as.factor

new_data <- bind_rows(infectious_success, infectious_unsuccess,
                      not_success, not_unsuccess)
new_data$Match %<>% as.factor

#### Effect of mutation on infection ####
m1 <- glmer(Infected~Mutation+(1|Replicate), data=new_data, family=binomial("logit"))
summary(multcomp::glht(m1))
summary(m1)

logit2prob(fixef(m1)[1])
logit2prob(fixef(m1)[1]+fixef(m1)[2])
CI <- confint(m1, parm="beta_")
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[3]+CI[4])

drop1(m1, test="Chisq")

#### Effect of escape mutations on infection ####

m2 <- glmer(Infected~Escape+(1|Replicate), data=data, family=binomial("logit"))
summary(multcomp::glht(m2))

logit2prob(fixef(m2)[1])
logit2prob(fixef(m2)[1]+fixef(m2)[2])
CI <- confint(m2, parm="beta_")
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[3]+CI[4])

model.tables(aov(m2), "mean")

anova(m2, test="Chisq")

#### Effect of spacer targetting on infection ####

m3 <- glmer(Infected~NotTargetted+(1|Replicate), data=data, family=binomial)
summary(multcomp::glht(m3))

logit2prob(fixef(m3)[1])
logit2prob(fixef(m3)[1]+fixef(m3)[2])
logit2prob(fixef(m3)[1]+fixef(m3)[3])
logit2prob(fixef(m3)[1]+fixef(m3)[4])
model.tables(aov(m3), "mean")

drop1(m3, test="Chisq")

new_data$NotTargetted %<>% relevel(ref="1")

CI <- confint.merMod(m3, parm="beta_")
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[1]+CI[3])
logit2prob(CI[4]+CI[5])
logit2prob(CI[4]+CI[6])


#### Effect of PWD on infection ####

data$PWD %<>% as.factor

m4 <- glmer(Infected~PWD+(1|Replicate), data=data, family=binomial("identity"))
summary(m4)
drop1(m4, test="Chisq")
wht <- glht(m4, linfct = mcp(PWD="Tukey"))
plot(print(confint(wht)))

## Signifcant negative relationship between probablity of infection and PWD
## Slope = -0.90272, Z = -9.134, p < 0.0001

logit2prob(fixef(m4)[2])
logit2prob(fixef(m4)[1]+fixef(m4)[2])
logit2prob(fixef(m4)[1]+fixef(m4)[3])
logit2prob(fixef(m4)[1]+fixef(m4)[4])
logit2prob(fixef(m4)[1]+fixef(m4)[5])
logit2prob(fixef(m4)[1]+fixef(m4)[6])

data$PWD %<>% relevel(ref="0.631")
CI <- confint(m4, parm="beta_")
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[1]+CI[3])
logit2prob(CI[1]+CI[4])
logit2prob(CI[1]+CI[5])
logit2prob(CI[1]+CI[6])

logit2prob(CI[7]+CI[8])
logit2prob(CI[7]+CI[9])
logit2prob(CI[7]+CI[10])
logit2prob(CI[7]+CI[11])
logit2prob(CI[7]+CI[12])
#### Effect of Match on Infectivity ####


#### Interaction b/w timeshift and escape mutations ####
data$Escape %<>% as.factor()
m5 <- glmer(Infected~NotTargetted*Phage.Background+(1|Phage.Background), data=filter(data, Mutation==1),
            family=binomial)
summary(multcomp::glht(m5))

logit2prob(fixef(m5)[1])
logit2prob(fixef(m5)[1]+fixef(m5)[2])
logit2prob(fixef(m5)[1]+fixef(m5)[3])
logit2prob(fixef(m5)[1]+fixef(m5)[4])

data$Escape %<>% relevel(ref="0")
data$Escape %<>% relevel(ref="1")

data$Phage.Background %<>% relevel(ref="Past")
data$Phage.Background %<>% relevel(ref="Contemporary")
data$Phage.Background %<>% relevel(ref="Future")

CI <- confint.merMod(m5, parm="beta_")
logit2prob(CI)

#### Test plots #### 
p0 <- ggplot(aes(x=NotTargetted, y=Infected), data=data)+
  geom_jitter(height=0.05)+
  geom_point()+
  binomial_smooth()+
  facet_wrap(~Phage.Background)+
  NULL
quartz()
p0

# pwd <- read.csv("./phage_sequences/summary_data/infections_PWD_sum.csv")
# pwd$PWD %<>% as.factor
# 
# p1 <- ggplot(aes(x=PWD, y=Mean), data=pwd)+
#   geom_point(size=3)+
#   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8)+
#   coord_cartesian(ylim=c(0,1))+
#   NULL
# p1

timeshift <- read.csv("./phage_sequences/summary_data/timeshift.csv")
timeshift$Mutation %<>% relevel(ref="Escape")
timeshift$Mutation %<>% relevel(ref="None")
timeshift$Phage.Background %<>% relevel(ref = "Future")
timeshift$Phage.Background %<>% relevel(ref = "Present")
timeshift$Phage.Background %<>% relevel(ref = "Past")

timeshift_fig <- ggplot(aes(x=Phage.Background, y=Mean, group=Mutation), data=timeshift)+
  geom_point(size=3, aes(colour=Mutation), position = position_dodge(.6))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, colour=Mutation), size=.8, width=0, position=position_dodge(.6))+
  coord_cartesian(ylim=c(0,1))+
  #ggtitle("Proportion of hosts that were\ninfected by or resisted a phage \nwith a mutation")+
  labs(x="Phage background", y="Proportion of hosts infected")+
  scale_colour_discrete(breaks=c("None", "Escape"),
                   labels=c("None", "Protospacer-\nassociated"))+
  theme_cowplot()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face="bold", size=16),
        legend.title.align = 0.5,
        legend.text = element_text(size=14),
        legend.text.align = 0.5,
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1.5, "cm"))+
  NULL
timeshift_fig

ggsave("Fig2.png", timeshift_fig, path="~/Desktop/", dpi=300, device="png",
       width=24, height=12, unit=c("cm"))


