### Dan's coevo stuff - spacer coverage analysis
# Created: 30/4/18 by Jack Common

rm(list=ls())

library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)

## Data
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
data$SpacerMiddle <- medians

## Coverage graphs
timepoint_names_facet = list(
  't1' = '1 d.p.i.',
  't4' = '4 d.p.i.',
  't9' = '9 d.p.i.'
)

timepoint_names_legend = c("1","4","9")

timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

plot1 <- ggplot(aes(x=SpacerMiddle, group=Replicate), data=data)+
  #geom_histogram(aes(fill=Locus), 
  #               position=position_dodge(),
  #               bins=50)+
  geom_dotplot(fill="white", alpha=0.5, colour="black",
               method="histodot", binwidth = 1,
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+
  facet_wrap(~Timepoint, labeller = timepoint_labeller)+
  theme_bw()+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704))+
  #scale_y_continuous(breaks=c(seq(0,1)),
  #                   labels=c(seq(0,20,5)))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(strip.text = element_text(face='bold', size=14))+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())

plot1

detach("package:cowplot")

ggsave("coverage_plot.png", plot1, path="./figs/",
       device="png", dpi=300,
       height=13, width=30, unit=c("cm"))



