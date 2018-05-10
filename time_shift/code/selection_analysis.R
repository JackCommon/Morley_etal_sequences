#### Morley et al - Selection analysis ####
# Created: 9/5/17 by Jack Common

#### Dependencies ####
#install.packages("lme4")
#install.packages("MuMIn")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("magrittr")

library(lme4)
library(MuMIn)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)

#### Functions ####
CoEvoRatio <- function(model1, model2){
  f1 <- model1@call$formula
  f2 <- model2@call$formula
  Env.MS <- anova(model1)$`Mean Sq`
  GE.MS <- anova(model2)$`Mean Sq`
  Ratio <- GE.MS/Env.MS
  cat("Relative importance of FSD:ARD:", Ratio)
}

#### 2.1 ####
one.t1 <- filter(data, Replicate == "2.1") %>% 
  filter(Host.Timepoint=="t1")

one.t4 <- filter(data, Replicate == "2.1") %>% 
  filter(Host.Timepoint=="t4")

one.t9 <- filter(data, Replicate == "2.1") %>% 
  filter(Host.Timepoint=="t9")
#### 2.2 ####
two.t1 <- filter(data, Replicate == "2.2") %>% 
  filter(Host.Timepoint=="t1")

two.t4 <- filter(data, Replicate == "2.2") %>% 
  filter(Host.Timepoint=="t4")

two.t9 <- filter(data, Replicate == "2.2") %>% 
  filter(Host.Timepoint=="t9")
#### 2.3 ####
three.t1 <- filter(data, Replicate == "2.3") %>% 
  filter(Host.Timepoint=="t1")

three.t4 <- filter(data, Replicate == "2.3") %>% 
  filter(Host.Timepoint=="t4")

three.t9 <- filter(data, Replicate == "2.3") %>% 
  filter(Host.Timepoint=="t9")
#### 2.4 ####
four.t1 <- filter(data, Replicate == "2.4") %>% 
  filter(Host.Timepoint=="t1")

four.t4 <- filter(data, Replicate == "2.4") %>% 
  filter(Host.Timepoint=="t4")

four.t9 <- filter(data, Replicate == "2.4") %>% 
  filter(Host.Timepoint=="t9")
#### 2.5 ####
#five.t1 <- filter(data, Replicate == "2.5") %>% 
#  filter(Host.Timepoint=="t1")

#five.t4 <- filter(data, Replicate == "2.5") %>% 
#  filter(Host.Timepoint=="t4")

#five.t9 <- filter(data, Replicate == "2.5") %>% 
#  filter(Host.Timepoint=="t9")
#### 2.6 ####
six.t1 <- filter(data, Replicate == "2.6") %>% 
  filter(Host.Timepoint=="t1")

six.t4 <- filter(data, Replicate == "2.6") %>% 
  filter(Host.Timepoint=="t4")

six.t9 <- filter(data, Replicate == "2.6") %>% 
  filter(Host.Timepoint=="t9")
#### 2.7 ####
seven.t1 <- filter(data, Replicate == "2.7") %>% 
  filter(Host.Timepoint=="t1")

seven.t4 <- filter(data, Replicate == "2.7") %>% 
  filter(Host.Timepoint=="t4")

seven.t9 <- filter(data, Replicate == "2.7") %>% 
  filter(Host.Timepoint=="t9")
#### 2.11 ####
eleven.t1 <- filter(data, Replicate == "2.11") %>% 
  filter(Host.Timepoint=="t1")

eleven.t4 <- filter(data, Replicate == "2.11") %>% 
  filter(Host.Timepoint=="t4")

eleven.t9 <- filter(data, Replicate == "2.11") %>% 
  filter(Host.Timepoint=="t9")

#### Analysis ####
# Environment-only effect
m1 <- glm(Infected~Phage.Timepoint,
            data=two.t4,
            family=binomial())
# Genotype * Environment effect
m2 <- glmer(Infected~Phage.Timepoint+(Phage.Timepoint|Phage.Genotype),
            data=two.t4,
            family=binomial())
summary(m1)
summary(m2)

# Get the FSD:ARD ratio
CoEvoRatio(m1, m2)

# If you get the annoying error on the env-only model then you can run it
# as a GLM and get the ratio as 1-deviance ratio of the two models
1-deviance(m2)/deviance(m1)

#### Models of FSD scores ####
scores <- read.csv("./time_shift/summary_data/FSD_scores.csv")
scores$Replicate %<>% as.factor()

m3 <- glm(Score~Timepoint,
          data = scores)
summary(m3)
model.tables(aov(m3), "mean")

#### Figures ####

FSD_means <- read.csv("./time_shift/summary_data/FSD_means.csv")

# With SE bars
FSD_SE <- ggplot(aes(x=Timepoint, y=Score), data=FSD_means)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Score-SE, ymax=Score+SE),
                size=.8, width=.05)+
  theme_bw()+
  labs(y="MS(GxE)/MS(E)")+
  ggtitle("Relative importance of fluctuating compared to\nescalating selection. Means and SEs (n=7)")+
  scale_x_discrete(name="Transfer",
                      breaks=c("t1", "t4", "t9"),
                      labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))
FSD_SE

ggsave("FSD_scores.png", FSD_SE, path="./figs/",
       device="png", dpi=300,
       height=15, width = 17, units = c("cm"))

# With 95% CIs
FSD_CI <- ggplot(aes(x=Timepoint, y=Score), data=FSD_means)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper),
                size=.8, width=.05)+
  theme_bw()+
  labs(y="MS(GxE)/MS(E)")+
  scale_x_discrete(name="Transfer",
                   breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0,0.5,0.1)))
FSD_CI
