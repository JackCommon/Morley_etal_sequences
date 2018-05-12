#### Morley et al - spacer diversity analysis ####
# Created: 02/5/18 by Jack Common

rm(list=ls())

#### Dependencies ####

#install.packages("ggplot2")
#install.packages("scales")
#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("magrittr")
#install.packages("cowplot")

library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)

#### Data ####
data = read.csv("./sequences/summary_data/all_spacer_data.csv", header=T)
data$Replicate %<>% as.factor
data$Clone %<>% as.factor()
data$SpacerNumber %<>% as.factor()

# Calculate the median value between the spacer sequence blast hit start and end,
# and add it to the dataframe
medians <- c()
for(i in seq(1,212,1)){
  medians[i] <- median(c(data[i,7], data[i,8]))
}
data$SpacerMiddle <- medians %>% as.integer()

#### Unique spacers ####
## Summarise the data and get the counts of each unique median hit value for each factor
## combination
spacer.summary <- data %>% 
                  group_by(ProtospacerStart, ProtospacerEnd, Timepoint, Replicate, Locus) %>% 
                  count(SpacerMiddle)

# Take a look at the totals for each replicate
# Can filter by Timepoint as well 

spacer.summary %>% filter(Replicate=="2.1")
spacer.summary %>% filter(Replicate=="2.2")
spacer.summary %>% filter(Replicate=="2.3")
spacer.summary %>% filter(Replicate=="2.4")
spacer.summary %>% filter(Replicate=="2.5")
spacer.summary %>% filter(Replicate=="2.6")
spacer.summary %>% filter(Replicate=="2.7")
spacer.summary %>% filter(Replicate=="2.11")

#### Export the data (will be used later) ####
spacer.summary %>% 
  filter(Replicate=="2.1") %>%
  write.csv("./sequences/summary_data/2.1_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.2") %>%
  write.csv("./sequences/summary_data/2.2_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.3") %>%
  write.csv("./sequences/summary_data/2.3_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.4") %>%
  write.csv("./sequences/summary_data/2.4_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.5") %>%
  write.csv("./sequences/summary_data/2.5_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.6") %>%
  write.csv("./sequences/summary_data/2.6_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.7") %>%
  write.csv("./sequences/summary_data/2.7_spacers.csv", row.names = F)

spacer.summary %>% 
  filter(Replicate=="2.11") %>%
  write.csv("./sequences/summary_data/2.11_spacers.csv", row.names = F)

# And timepoint
spacer.summary %>% filter(Timepoint=="t1")
spacer.summary %>% filter(Timepoint=="t4")
spacer.summary %>% filter(Timepoint=="t9")

# And locus
spacer.summary %>% filter(Locus=="CR1")
spacer.summary %>% filter(Locus=="CR3")

#### Diversity calculations ####
## A function to calculate Simpson's Index of Diversity 
simpson  <- function(df){
  # Get the total "population" size (N) by totalling the
  # number of unique spacer sequences
  N <- sum(df$n)    
  # Initialize variable D at 0
  D = 0
  # Loop through the data and do the calculation for each genotype,
  # and add it to D
  for (i in seq(1, length(df$n),1)){
    D = D + ((df$n[i])/N)^2
  }
  # Print out the value of 1-D
  cat("Simpson's Index of Diversity (1-D) =", 1-D)
  # Copy D to the clipboard
  clip = pipe('pbcopy', 'w')
  write.table(1-D, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
}

## Diversity in each replicate
# 2.1
r1 <- filter(spacer.summary, Replicate=='2.1')
simpson(r1)
# 2.2
r2 <- filter(spacer.summary, Replicate=='2.2')
simpson(r2)
# 2.3
r3 <- filter(spacer.summary, Replicate=='2.3')
simpson(r3)
# 2.4
r4 <- filter(spacer.summary, Replicate=='2.4')
simpson(r4)
# 2.5
r5 <- filter(spacer.summary, Replicate=='2.5')
simpson(r5)
# 2.6
r6 <- filter(spacer.summary, Replicate=='2.6')
simpson(r6)
# 2.7
r7 <- filter(spacer.summary, Replicate=='2.7')
simpson(r7)
# 2.11
r11 <- filter(spacer.summary, Replicate=='2.11')
simpson(r11)

## Diversity in each timepoint
# t1
tp1 <- filter(spacer.summary, Timepoint=="t1")
simpson(tp1)
#t4
tp4 <- filter(spacer.summary, Timepoint=="t4", Locus=="CR3")
simpson(tp4)
#t9
tp9 <- filter(spacer.summary, Timepoint=="t9", Locus=="CR3")
simpson(tp9)

## Diversity in each locus
#CR1
cr1 <- filter(spacer.summary, Locus=="CR1")
simpson(cr1)
#CR3
cr3 <- filter(spacer.summary, Locus=="CR3")
simpson(cr3)

#### Figures ####
all_comps <- read.csv("./sequences/summary_data/diversity_all_combinations.csv")
all_comps$Replicate %<>% as.factor()

simpson_all <- ggplot(aes(x=Replicate, y=Simpson, group=Timepoint), 
                      data=all_comps)+
  geom_bar(stat="identity", aes(fill=Timepoint), 
           colour="black",
           position=position_dodge())+
  coord_cartesian(ylim=c(0,1))+
  labs(y="Spacer diversity\n(Simpson's Diversity Index)")
simpson_all

detach("package:cowplot")

ggsave("replicate_diversity.png", simpson_all, path="./figs/",
       device="png", dpi=300, width = 18, height=12, units = c("cm"))

timepoint_comps <- read.csv("./sequences/summary_data/diversity_timepoints.csv")
simpson_times <- ggplot(aes(x=Timepoint, y=Simpson, group=Locus), 
                        data=timepoint_comps)+
  geom_bar(stat="identity", aes(fill=Locus), 
           colour="black",
           position=position_dodge())+
  coord_cartesian(ylim=c(0,1))+
  labs(y="Spacer diversity\n(Simpson's Diversity Index)")
simpson_times

ggsave("timepoint_diversity.png", simpson_times, path="./figs/",
       device="png", dpi=300, width = 18, height=12, units = c("cm"))


#### Analysis ####
library(lme4)
str(all_comps)

m1 <- glm(Simpson~Timepoint,
            data=all_comps)
summary(m1)
anova(m1, test="Chisq")

sresid1 <- resid(m1, type="pearson")
hist(sresid1)
fitted.glm <- fitted(m1, level=1)
plot(sresid1~all_comps$Timepoint)
plot(m1)

#### Collated data ####
# Build a full dataset of the unique spacers in each replicate, timepoint and locus
one <- read.csv("./sequences/summary_data/2.1_spacers.csv")
two <- read.csv("./sequences/summary_data/2.2_spacers.csv")
three <- read.csv("./sequences/summary_data/2.3_spacers.csv")
four <- read.csv("./sequences/summary_data/2.4_spacers.csv")
five <- read.csv("./sequences/summary_data/2.5_spacers.csv")
six <- read.csv("./sequences/summary_data/2.6_spacers.csv")
seven <- read.csv("./sequences/summary_data/2.7_spacers.csv")
eleven <- read.csv("./sequences/summary_data/2.11_spacers.csv")

unique_spacers <- bind_rows(one, two, three, four,
                            five, six, seven, eleven)
unique_spacers$Replicate %<>% as.factor
unique_spacers$Timepoint %<>% as.factor
unique_spacers$Locus %<>% as.factor
unique_spacers %<>% select(-SpacerMiddle)

unique_spacers$Replicate %<>% relevel(ref="2.11")
unique_spacers$Replicate %<>% relevel(ref="2.7")
unique_spacers$Replicate %<>% relevel(ref="2.6")
unique_spacers$Replicate %<>% relevel(ref="2.5")
unique_spacers$Replicate %<>% relevel(ref="2.4")
unique_spacers$Replicate %<>% relevel(ref="2.3")
unique_spacers$Replicate %<>% relevel(ref="2.2")
unique_spacers$Replicate %<>% relevel(ref="2.1")

write.csv(file="./sequences/summary_data/unique_spacers_noseqs.csv", unique_spacers, row.names = F)
