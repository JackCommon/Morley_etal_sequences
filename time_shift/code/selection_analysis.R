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
  Env.MS <- anova(model1)$`Mean Sq`
  GE.MS <- anova(model2)$`Mean Sq`
  Ratio <- GE.MS/Env.MS
  print(Ratio)
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
five.t1 <- filter(data, Replicate == "2.5") %>% 
  filter(Host.Timepoint=="t1")

five.t4 <- filter(data, Replicate == "2.5") %>% 
  filter(Host.Timepoint=="t4")

five.t9 <- filter(data, Replicate == "2.5") %>% 
  filter(Host.Timepoint=="t9")
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
            data=five.t9,
            family=binomial())
# Genotype * Environment effect
m2 <- glmer(Infected~Phage.Timepoint+(Phage.Timepoint|Phage.Genotype),
            data=five.t9,
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

# Test to see if scores are normally-distributed
shapiro.test(scores$Score)
# Not normally distributed, so try a square-root transformation
shapiro.test(sqrt(scores$Score))
# Will move ahead with square-rooted residuals

m3 <- lmer(sqrt(Score)~Timepoint+(1|Timepoint),
          data = scores)

m4 <- lmer(sqrt(Score)~Timepoint+(1|Replicate),
            data=scores)
# Check for heteroskedacity
par(mfrow=c(2,2))
plot(m3)
plot(m4)

# Compare models with AIC, Log-Likelihood, and finally a chi-squared
AIC(m3, m4) %>% compare_AICs()
logLik(m3)
logLik(m4)
anova(m3, m4, test="Chisq")

# Looks like model 4, with Replicate as a random effect, is the best performing
summary(m4)
m4.coefs <- fixef(m4)
# Need to square the coefficients to put them back on the same scale as the data 
m4.coefs[[1]]^2
(m4.coefs[[1]]+m4.coefs[[2]])^2
(m4.coefs[[1]]+m4.coefs[[3]])^2
m4.CIs <- confint(m4, parm="beta_")

m4.CIs
(m4.CIs[1])^2; (m4.CIs[4])^2
(m4.CIs[1]+m4.CIs[2])^2
(m4.CIs[1]+m4.CIs[3])^2

(m4.CIs[4]+m4.CIs[5])^2
(m4.CIs[4]+m4.CIs[6])^2

#### Figures ####

FSD_means <- read.csv("./time_shift/summary_data/FSD_means.csv")

# With 95% CIs
FSD_CI <- ggplot(aes(x=Timepoint, y=Score), data=FSD_means)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Score-(1.96*SE), ymax=Score+(1.96*SE)),
                size=.8, width=.05)+
  theme_bw()+
  labs(y="MS(GxE)/MS(E)")+
  ggtitle("Relative importance of fluctuating compared to\nescalating selection. Means and 95% CIs (n=8)")+
  scale_x_discrete(name="Timepoint",
                   breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))
FSD_CI

ggsave("FSD_scores.png", FSD_CI, path="./figs/",
       device="png", dpi=300,
       height=15, width = 17, units = c("cm"))
