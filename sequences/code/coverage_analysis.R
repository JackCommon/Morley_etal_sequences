#### Morley et al - Spacer coverage analysis ####
# Created: 30/4/18 by Jack Common

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

data$Replicate %<>% relevel(ref="2.11")
data$Replicate %<>% relevel(ref="2.7")
data$Replicate %<>% relevel(ref="2.6")
data$Replicate %<>% relevel(ref="2.5")
data$Replicate %<>% relevel(ref="2.4")
data$Replicate %<>% relevel(ref="2.3")
data$Replicate %<>% relevel(ref="2.2")
data$Replicate %<>% relevel(ref="2.1")
#### Functions ####

# Calculate the median value between the spacer sequence blast hit start and end,
# and add it to the dataframe
medians <- c()
for(i in seq(1,212,1)){
  medians[i] <- median(c(data[i,7], data[i,8]))
}
data$SpacerMiddle <- medians

# Some functions for better labels
timepoint_names_facet = list(
  't1' = '1 d.p.i.',
  't4' = '4 d.p.i.',
  't9' = '9 d.p.i.'
)

timepoint_names_legend = c("1","4","9")

timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

##### Figures ####
# Simple dot-plot of the centre-point of the blast hits. Number of points
# indicates how many times that hit occured.

# To get the grid arrangment, I need to make 13 individual plots and then arrange them 
t1 <- filter(data, Timepoint=="t1")
t4 <- filter(data, Timepoint=="t4")
t9 <- filter(data, Timepoint=="t9")

#### T1, Rep 2.2 - Plot1 ####
plot1 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=t1)+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.2 T1")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
 # theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  theme(legend.direction = "horizontal")+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())

#### T4, Rep 2.3 - Plot2 ####
plot2 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t4, Replicate=="2.3"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.3 T4")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())

#### T4, Rep 2.4 - Plot3 ####
plot3 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t4, Replicate=="2.4"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.4 T4")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T4, Rep 2.6 - Plot4 ####
plot4 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t4, Replicate=="2.6"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.6 T4")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T4, Rep 2.7 - Plot5 ####
# Need the blue fill colour for CR3 points
g <- ggplot_build(plot8)
unique(g$data[[1]]["fill"])

# Blue colour is #00BFC4

plot5 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t4, Replicate=="2.7"))+
  geom_dotplot(fill="#00BFC4", alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.7 T4")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.1 - Plot6 ####
plot6 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.1"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.1 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.2 - Plot7 ####
plot7 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.2"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.2 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.3 - Plot8 ####
plot8 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.3"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.3 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  theme(legend.direction = "horizontal")+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())



#### T9, Rep 2.4 - Plot9 ####
plot9 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.4"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.4 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())

#### T9, Rep 2.5 - Plot10 ####
plot10 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.5"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.5 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.6 - Plot11 ####
plot11 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.6"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.6 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.7 - Plot12 ####
plot12 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.7"))+
  geom_dotplot(fill="#00BFC4", alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.7 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
#### T9, Rep 2.11 - Plot13 ####
plot13 <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t9, Replicate=="2.11"))+
  geom_dotplot(aes(fill=Locus), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  #facet_wrap(Timepoint~Replicate,ncol = 3, nrow=5)+
  theme_cowplot()+
  ggtitle("2.11 T9")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
### Arrange plots ####
# Extract the legend from plot8 (with both CR1 and CR3 spacers) for later use
legend <- gtable_filter(ggplot_gtable(ggplot_build(plot8)), "guide-box")

blank <- ggplot(aes(x=SpacerMiddle, group=Locus), data=subset(t4, Replicate=="2.4"))+
  labs(x="")+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())
blank

t1_toprow <- plot_grid(plot1+
                         xlab("")+ylab("T1")+ggtitle("")+
                         theme(legend.position = "none",
                               axis.title.y = element_text(colour="red", size=40, angle = 0, vjust=0.5),
                               axis.text = element_text(size=30)),
                       blank, blank, blank,
                       labels = "A", label_size = 40,
                       nrow=1, scale=c(1,100,1,1))
t1_toprow

t4_midrow <- plot_grid(plot2+labs(x="", y="T4")+ggtitle("")+
                         theme(legend.position = "none",
                               axis.text = element_text(size=30),
                               axis.title.y = element_text(colour="red", size=40, angle = 0, vjust=0.5)),
                       plot3+labs(x="")+ggtitle("")+
                         theme(legend.position = "none",
                               axis.text = element_text(size=30)),
                       plot4+labs(x="")+ggtitle("")+
                         theme(legend.position = "none",
                               axis.text = element_text(size=30)),
                       plot5+labs(x="")+ggtitle("")+
                         theme(legend.position = "none",
                               axis.text = element_text(size=30)),
                       labels = c("B", "C", "D", "E"), label_size = 40,
                       ncol=4, nrow=1, align="v")
t4_midrow

t9_bottomrow <- plot_grid(plot6+labs(x="", y="T9")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30),
                                  axis.title.y = element_text(colour="red", size=40, angle = 0, vjust=.1)),
                          plot7+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          plot8+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          plot9+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          plot10+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          plot11+labs(x="Position on phage genome")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30),
                                  axis.title.x = element_text(size=50, hjust=-.8,
                                                              margin=margin(t=30))),
                          plot12+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          plot13+labs(x="")+ggtitle("")+
                            theme(legend.position = "none",
                                  axis.text = element_text(size=30)),
                          labels = c("F", "G", "H", "I", "J", "K", "L", "M"), label_size = 40,
                          ncol=4, nrow=4, align="hv")
t9_bottomrow

all_legend <- plot_grid(t1_toprow,
                  t4_midrow,
                  t9_bottomrow,
                  align = "v", nrow=3, ncol=1, 
                  rel_heights = c(1,1,4))
quartz()
all_legend


### Save figures ####
# Detach cowplot as it blocks ggsave
detach("package:cowplot")

ggsave("coverage_plot.png", plot1, path="./figs/",
       device="png", dpi=300,
       height=13, width=30, unit=c("cm"))

ggsave("all_legend.png", all_legend, path="./figs/",
       device="png", dpi=300,
       height=100, width=75, unit=c("cm"))

ggsave("legend.png", legend, path="./figs/",
       device="png", dpi=900,
       height=2, width=10, unit=c("cm"))
 



