#### Phage protospacers sequence analysis ####
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
#### Data ####
data <- read.csv("./phage_sequences/summary_data/infectivity_protospacers_full.csv")

data <- select(data, -Phage.Genotype, -Host.Timepoint, -Host.Genotype, -Phage.Timepoint, -Host.Background, -Phage.Background)

data$Replicate %<>% as.factor
#data$SpacerNumber %<>% as.factor
data$Escape <- ifelse(is.na(data$Escape)==T, 0, data$Escape) %>% as.integer()
data$Other <- ifelse(is.na(data$Other)==T, 0, data$Other) %>% as.integer()
data$Protospacer <- ifelse(is.na(data$Protospacer)==T, 0, data$Protospacer) %>% as.integer()
#data$SpacersTargetted <- ifelse(is.na(data$SpacersTargetted)==T, 0, data$SpacersTargetted) %>% as.factor

data$Mutation <- ifelse(data$Escape==1,1,0) 
data$Mutation <- ifelse(data$Other==1,1,data$Mutation) %>% as.integer

data_copy <- data
data_copy$Escape <- ifelse(data_copy$Escape==1, "yes", "no") %>% as.factor
data_copy$Protospacer <- ifelse(data_copy$Protospacer==1, "yes", "no") %>% as.factor
data_copy$Mutation <- ifelse(data_copy$Mutation==1, "yes", "no") %>% as.factor

data_proto <- filter(data, Other==0)

data_infect <- filter(data, Infected==1)
data_resist <- filter(data, Resistant==1)

#### Analysis ####
# Very simple model to test things out

m1 <- glmer(Infected~Mutation+(1|Replicate), data=data, family=binomial("logit"))
summary(m1)

logit2prob(fixef(m1)[2])
logit2prob(fixef(m1)[1]+fixef(m1)[2])
CI <- confint(m1, parm="beta_")
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[3]+CI[4])

drop1(m1, test="Chisq")

m2 <- glm(Infected~Escape, data=data, family=binomial("logit"))
summary(m2)

logit2prob(coef(m2)[1])
logit2prob(coef(m2)[1]+coef(m2)[2])
CI <- confint(m2)
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[3]+CI[4])

model.tables(aov(m2), "mean")

anova(m2, test="Chisq")

data$Type <- ifelse(data$Escape==1,"Escape", NA) 
data$Type <- ifelse(data$Other==1, "Other", data$Type) %>% as.factor

m3 <- glm(Infected~as.factor(Mutation)*as.factor(Escape), data=data, family=binomial("logit"))
summary(m3)

model.tables(aov(m3), "mean")

logit2prob(coef(m3)[1])
logit2prob(coef(m3)[1]+coef(m3)[2])
logit2prob(fixef(m3)[1]+fixef(m3)[3])
CI <- confint(m3)
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[3]+CI[4])

anova(m3, test="Chisq")

m4 <- glm(Infected~NotTargetted, data=data, family=binomial())
summary(m4)

logit2prob(coef(m4)[1])
logit2prob(coef(m4)[1]+coef(m4)[2])
logit2prob(coef(m4)[1]+coef(m4)[3])
model.tables(aov(m4), "mean")

anova(m4, test="Chisq")

CI <- confint(m4)
logit2prob(CI)
logit2prob(CI[1]+CI[2])
logit2prob(CI[1]+CI[3])
logit2prob(CI[4]+CI[5])
logit2prob(CI[4]+CI[6])

#### Plots ####
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

p0 <- ggplot(aes(x=as.factor(Protospacer), y=Infected), fill=as.factor(Protospacer), data=filter(data, Escape==1))+
  #geom_point()+
  geom_bar(stat="identity", aes(fill=as.factor(Escape)))+
  #geom_text(aes(label=as.factor(Escape)))+
  #geom_boxplot(position="fill")+
  #geom_smooth(method="loess")+
  #geom_jitter(height = 0.05) +
  #facet_wrap(~Replicate)+
#  coord_flip()+
# binomial_smooth()+
  theme_cowplot()+
  labs(x="Location", y="Number infected")+
  scale_fill_discrete(name=c("Location"),
                      breaks=c(0,1),
                      labels=c("PAM", "Protospacer sequence"))+
  scale_x_discrete(breaks=c(0,1),
                   labels=c("PAM", "Protospacer sequence"))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold", size=14),
        legend.title.align = 0.5,
        legend.text = element_text(size=12),
        legend.position = "none")+

  NULL
p0

sumdat1 <- read.csv("./phage_sequences/summary_data/infections_targetted_sum.csv")
sumdat1$NotTargetted %<>% as.factor()

p1 <- ggplot(aes(x=NotTargetted, y=Mean), data=sumdat1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.8, width=0)+
  coord_cartesian(ylim=c(0,1))+
  #ggtitle("Proportion of hosts that were\ninfected by or resisted a phage \nwith a mutation")+
  labs(x="Number of host CRISPR spacers phage\nhad not evolved SNPs against", y="Proportion infected")+
  #facet_wrap(~Location,scales = "free_x")+
  theme_cowplot()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=14))+
  NULL
p1

all_sums <- read.csv("./phage_sequences/summary_data/all_summaries.csv")
all_sums$Mutation %<>% relevel(ref="Protospacer-associated")
all_sums$Mutation %<>% relevel(ref="Random")
all_sums$Mutation %<>% relevel(ref="None")


p2 <- ggplot(aes(x=Mutation, y=Mean), data=all_sums)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.8, width=0)+
  coord_cartesian(ylim=c(0,1))+
  #ggtitle("Proportion of hosts that were\ninfected by or resisted a phage \nwith a mutation")+
  labs(x="Number of host CRISPR spacers phage\nhad not evolved SNPs against", y="Proportion infected")+
  #facet_wrap(~Location,scales = "free_x")+
  theme_cowplot()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=14))+
  NULL
p2


p3 <- ggplot(aes(x=as.factor(Protospacer), fill=as.factor(Protospacer)), data=filter(data, Escape==1))+
  geom_bar()+
  theme_cowplot()+
  labs(x="Location", y="")+
  scale_fill_discrete(name=c("Location"),
                      breaks=c(0,1),
                      labels=c("PAM", "Protospacer sequence"))+
  scale_x_discrete(breaks=c(0,1),
                   labels=c("PAM", "Protospacer sequence"))+
  scale_y_continuous(breaks=c(seq(0,350,50)))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=14),
        legend.title = element_text(face="bold", size=14),
        legend.title.align = 0.5,
        legend.text = element_text(size=12),
        legend.position = "none")+
  
  NULL
p3

proto_counts <- read.csv("./phage_sequences/summary_data/proto_seq_summary_stats.csv")
proto_counts %<>% filter(Mutation!=c("Total"))
proto_counts$Mutation %<>% relevel(ref="PAM")
proto_counts$Mutation %<>% relevel(ref="Protospacer sequence")
proto_counts$Mutation %<>% relevel(ref="Protospacer-associated")
proto_counts$Mutation %<>% relevel(ref="Random")
proto_counts$Mutation %<>% relevel(ref="None")

just_proto <- select(proto_counts, -None, -Random, -Protospacer-associated)

m_data <- melt(proto_counts, id.vars = c("Mutation"))

p4 <- ggplot(aes(x=variable, y=value, fill=Mutation), data=m_data)+
  geom_col(position = position_dodge(1), colour="black")
p4

phage_plots <- plot_grid(p1, p3+ylab(""), 
                         p2, p4+ylab(""),
                         rel_widths = c(1,1.1,1,1.1), align = "hv",
                         labels = c("A", "B", "C", "D"), label_colour = "red")
phage_plots



ggsave("infection_sequence_plots.png", phage_plots, path="./figs/",
       device="png", width=25, height = 20, units=c("cm"), dpi=300)
ggsave("sequence_summary_plots.png", frequencies, path="./figs/",
       device="png", width=17, height = 25, units=c("cm"), dpi=300)


#### Notes ####

# Escape mutations (i.e. in the protospacers sequence or PAM) = 351
# Other mutations = 160
# No mutations detected = 257
# N = 768

# Protospacers SNPs = 357 (79.7%)
# PAM SNPs = 91 (20.3%)
