#### Paper Figure Generation #####
## Morley et al
## Created 22/5/18 by Jack Common

rm(list=ls())

#### Dependencies ####
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(cowplot)
library(scales)
library(tidyr)

pd <- position_dodge(0.1)

### Figure 1 - Host & phage population dynamics ####
data <- read.csv("./phage_survival/original_data/phage_titre_Jack.csv")
data$Replicate %<>% as.factor()
#data$Timepoint %<>% as.factor()
data$Treatment %<>% as.factor()

# Convert negative values to zeros
data$OD600 <- ifelse(data$OD600<0,0,data$OD600)

# Melt data for PFU and OD600 dual plots
dataM <- melt(data, measure.vars = c("pfu", "OD600"))
dataM <- plyr::rename(dataM, c("variable"="measurement"))
dataM$Timepoint %<>% as.factor

phageIDs <- c(rep("p1",31) , rep("p2", 31), rep("p3", 31),
              rep("p4",31) , rep("p5", 31), rep("p6", 31),
              rep("p7",31) , rep("p8", 31), rep("p9", 31),
              rep("p10",31) , rep("p11", 31), rep("p12", 31)) %>% 
  rep(4)

hostIDs <- c(rep("h1",31) , rep("h2", 31), rep("h3", 31),
             rep("h4",31) , rep("h5", 31), rep("h6", 31),
             rep("h7",31) , rep("h8", 31), rep("h9", 31),
             rep("h10",31) , rep("h11", 31), rep("h12", 31)) %>% 
  rep(4)

ID2 <- c(phageIDs, hostIDs)

dataM$ID2 <- as.factor(ID2)


### 10^9 plot ##
nine_plot <- ggplot(aes(y=log10(pfu), x=Timepoint, group=Replicate), 
                    data=subset(data, Treatment == '9'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u./C.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{9}*" Treatment")))+
  facet_wrap(~Replicate)+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))+
  geom_line(aes(y=(OD600*1e+8)+1), colour="blue")+
  NULL
#quartz()
#nine_plot

### 10^8 plot ##
eight_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(data, Treatment == '8'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u./C.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{8}*" Treatment")))+
  facet_wrap(~Replicate)+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))+
  geom_line(aes(y=(OD600*1e+8)+1), colour="blue")+
  NULL
#quartz()
#eight_plot

### 10^7 plot ##
seven_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(data, Treatment == '7'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u./C.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{7}*" Treatment")))+
  facet_wrap(~Replicate)+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))+
  geom_line(aes(y=(OD600*1e+8)+1), colour="blue")+
  NULL
#quartz()
#seven_plot

### 10^6 plot ##
six_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                   data=subset(data, Treatment == '6'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u./C.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{6}*" Treatment")))+
  facet_wrap(~Replicate)+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))+
  geom_line(aes(y=(OD600*1e+8)+1), colour="blue")+
  NULL
#quartz()
#six_plot

### Arrange and save plots ##
Fig1 <- plot_grid(nine_plot+labs(x=""), 
                        eight_plot+labs(x=""),
                        seven_plot, 
                        six_plot,
                        labels = c("A", "B", "C", "D"), label_colour = "red")

ggsave("Fig1.png", Fig1, path="./figs/paper/",
       device="png", dpi=300,
       height=40, width=54, unit=c("cm"))

### Figure 2 - Evolution of resistance/infectivity ####
infect_resist_evo <- read.csv("./infectivity/summary_data/infect_resist_means.csv")

infect_plot <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=infect_resist_evo)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  #geom_path(stat="identity", size=.8, linetype=2)+
  theme_cowplot()+
  labs(x="Timepoint", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(axis.title.x = element_text(hjust=1.4))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,1))+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  NULL

resist_plot <- ggplot(aes(y=Mean.Resist, x=Timepoint, group=Group), data=infect_resist_evo)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  #geom_path(stat="identity", size=.8, linetype=2)+
  theme_cowplot()+
  labs(x="", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))

Fig2 <- plot_grid(infect_plot, resist_plot,
                  nrow=1, align = "hv", labels = c("A", "B"))

Fig2

ggsave("Fig2.png", Fig2, path="./figs/paper/",
       device="png", dpi=300, width=20, height=10, units = c("cm"))

### Figure 3 - Spacer diversity ####
#### Mean relative frequencies of different spacer nos 
#install.packages("wesanderson")
library(wesanderson)

pal <- wes_palette("Royal1",4, type = "discrete")

spacer_prop <- read.csv("./spacer_sequences/summary_data/spacer_prop.csv")
spacer_prop$SpacerNumber %<>% relevel(ref="Three")
spacer_prop$SpacerNumber %<>% relevel(ref="Two")
spacer_prop$SpacerNumber %<>% relevel(ref="One")
spacer_prop$SpacerNumber %<>% relevel(ref="Zero")
spacer_prop %<>% complete(Timepoint, SpacerNumber)

spacer_prop_plot <- ggplot(aes(x=Timepoint, y=Mean, group=SpacerNumber), data=spacer_prop)+
  geom_col(aes(fill=SpacerNumber), colour="black",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8, position = position_dodge(.9))+
  scale_y_continuous(breaks=seq(0,1,.25))+
  coord_cartesian(ylim=c(0,1))+
  
  theme_cowplot()+
  labs(x="Days post-infection (dpi)", y="Relative frequency")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  scale_fill_manual(name="Number of\nspacers",
                    values=pal)+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'right',
        legend.direction = "horizontal",
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  NULL

#spacer_prop_plot

#### Spcaer mean counts spacers 
spacer_counts <- read.csv("./spacer_sequences/summary_data/spacer_count.csv")
spacer_counts$Timepoint %<>% as.factor
spacer_counts$Group %<>% as.factor

spacer_count_plot <- ggplot(aes(x=Timepoint, y=Mean, group=Group), data=spacer_counts)+
  geom_point(size=3)+
  # geom_path(linetype=2, size=.8 )+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8, position = position_dodge(.9))+
  #scale_y_continuous(breaks=seq(0,1,.25))+
  #coord_cartesian(ylim=c(0,1))+
  
  theme_cowplot()+
  labs(x="Days post-infection (dpi)", y="Spacer number")+
  #scale_x_discrete(breaks=c("t1", "t4", "t9"),
  #                  labels=c("1", "4", "9"))+
  #scale_fill_discrete(name="Number of\nspacers")+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'right',
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  NULL

#spacer_count_plot

#### PWD among spacers by replicate and timepoint 
PWDs <- read.csv("./spacer_sequences/summary_data/diversity_all_combinations.csv")
PWDs$Replicate %<>% as.factor()
PWDs$Replicate %<>% relevel(ref="2.11")
PWDs$Replicate %<>% relevel(ref="2.7")
PWDs$Replicate %<>% relevel(ref="2.6")
PWDs$Replicate %<>% relevel(ref="2.5")
PWDs$Replicate %<>% relevel(ref="2.4")
PWDs$Replicate %<>% relevel(ref="2.3")
PWDs$Replicate %<>% relevel(ref="2.2")
PWDs$Replicate %<>% relevel(ref="2.1")

reps_og <- levels(PWDs$Replicate)
reps_labs <- c(seq(1,8,1)) %>% as.character()

PWDs %<>% complete(Replicate, Timepoint) #adds in NAs where data is missing

PWDs$PWD<- ifelse(is.na(PWDs$PWD)==T, 0, PWDs$PWD) #convert those NAs to zeros so that bar widths make sense

PWD_all <- ggplot(aes(x=Replicate, y=PWD), 
                  data=PWDs)+
  geom_col(aes(fill=Timepoint),
           position=position_dodge(), colour="black")+
  coord_cartesian(ylim=c(0,1))+
  labs(y="Pairwise difference")+
  scale_fill_discrete(breaks=c("t1", "t4", "t9"),
                      labels=c("1", "4", "9"))+
  scale_x_discrete(breaks=reps_og,
                   labels=reps_labs)+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        legend.position = 'right',
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  NULL
#PWD_all

#### Make and save Fig 3 

spacer_legend <- get_legend(spacer_prop_plot)
PWD_legend <- get_legend(PWD_all)

Fig3Legends <- plot_grid(spacer_legend, PWD_legend,
                     ncol=1, nrow=2)

Fig3 <- plot_grid(spacer_count_plot,
                          spacer_prop_plot+
                            theme(legend.position = "none"),
                          PWD_all+
                            theme(legend.position = "none"),
                          Fig3Legends,
                          ncol=2, nrow=2, align = "hv", axis="l", labels=c("A", "B", "C"),
                          rel_widths = c(1,1,1.4))
#Fig3

ggsave("Fig3.png", Fig3, path="./figs/paper/",
       device="png", dpi=300, width=30, height=25, units = c("cm"))


### Figure 4 - Spacer coverage ####
coverage <- read.csv("./spacer_sequences/summary_data/all_spacer_data.csv", header=T)
coverage$Replicate %<>% as.factor
coverage$Clone %<>% as.factor()
coverage$SpacerNumber %<>% as.factor()

medians <- c()
for(i in seq(1,212,1)){
  medians[i] <- median(c(coverage[i,7], coverage[i,8]))
}
coverage$SpacerMiddle <- medians %>% as.integer()


coverage$Replicate %<>% relevel(ref="2.11")
coverage$Replicate %<>% relevel(ref="2.7")
coverage$Replicate %<>% relevel(ref="2.6")
coverage$Replicate %<>% relevel(ref="2.5")
coverage$Replicate %<>% relevel(ref="2.4")
coverage$Replicate %<>% relevel(ref="2.3")
coverage$Replicate %<>% relevel(ref="2.2")
coverage$Replicate %<>% relevel(ref="2.1")

timepoint_names_legend = c("1","4","9")
timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

rep_og <- levels(coverage$Replicate)
rep_labs <- c(seq(1,8,1)) %>% as.character()
rep_labeller = function(variable, value) {
  return(rep_labs[value])
}


Fig4 <- ggplot(aes(x=SpacerMiddle, group=Timepoint), data=coverage)+
  geom_dotplot(aes(fill=Timepoint), alpha=0.5, colour="black",
               method="histodot", binwidth = 1,          # bindwidth = 1 to ensure no erroneous overlaps
               dotsize=2000)+
  coord_cartesian(xlim=c(seq(1,34704,1)))+               # y-axis is as long as the phage 2972 genome
  
  theme_cowplot()+
  ggtitle("")+
  labs(x="Position on phage genome", y="")+
  scale_x_continuous(breaks=c(0,10000,20000,30000,34704),
                     labels = c("0Kb", "10Kb", "20Kb", "30Kb", ""))+ # breaks every 10kb
  
  facet_wrap(~Replicate, ncol=4, nrow=4, labeller = rep_labeller)+
  theme(strip.text = element_text(face='bold', size=14))+
  
  theme(plot.title = element_text(face="bold", size=18))+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'right',
        legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14),
        legend.direction = "vertical")+
  scale_fill_discrete(name = c("Days\npost-infection"),
                      breaks=c("t1", "t4", "t9"),
                      labels=timepoint_names_legend)+
  # Remove the y-axis line as it's uninformative
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  
  theme(panel.grid.minor= element_blank(),
        panel.grid.major = element_blank())+
  NULL
Fig4

ggsave("Fig4.png", Fig4, path="./figs/paper/",
       device="png", dpi=300,
       height=15, width=33, unit=c("cm"))


### Figure 5 - Protospacer data ####
#### Figure 5A & B
proto_counts <- read.csv("./phage_sequences/summary_data/proto_seq_summary_stats.csv")
proto_counts %<>% filter(Mutation!=c("Total"))
proto_counts$Mutation %<>% relevel(ref="PAM")
proto_counts$Mutation %<>% relevel(ref="Sequence")
proto_counts$Mutation %<>% relevel(ref="Protospacer-associated")
proto_counts$Mutation %<>% relevel(ref="Random")
proto_counts$Mutation %<>% relevel(ref="None")

just_proto <- filter(proto_counts, Mutation%in%c("Sequence", "PAM"))
general_proto <- filter(proto_counts, Mutation%in%c("None", "Random", "Protospacer-associated"))

pal2 <- wesanderson::wes_palette("Royal1",2, type = "discrete")
pal3 <- wesanderson::wes_palette("Royal1",3, type = "discrete")

Fig5A <- ggplot(aes(x=Mutation, y=Count, fill=Mutation), data=general_proto)+
  geom_col(position = position_dodge(1), colour="black", width=.5)+
  theme_cowplot()+
  labs(x="Mutation", y="Count")+
  scale_x_discrete(breaks=c("None", "Random", "Protospacer-associated"),
                  labels=c("None", "Random", "Protospacer-\nassociated"))+
  scale_fill_manual(values=pal3)+
  coord_cartesian(ylim=c(0,50))+
  scale_y_continuous(breaks=c(seq(0,50,10)))+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.position = "none")+
  NULL
# Fig5A

Fig5B <- ggplot(aes(x=Mutation, y=Count, fill=Mutation), data=just_proto)+
  geom_col(position = position_dodge(1), colour="black", width=.5)+
  theme_cowplot()+
  labs(x="Location of SNP in protospacer", y="Count")+
  #scale_x_discrete(breaks=c("t1", "t4", "t9"),
  #                 labels=c("1", "4", "9"))+
  scale_fill_manual(values=pal2)+
  coord_cartesian(ylim=c(0,50))+
  scale_y_continuous(breaks=c(seq(0,50,10)))+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.position = "none")+
  NULL
# Fig5B

#### Fig 5C

all_sums <- read.csv("./phage_sequences/summary_data/all_summaries.csv")
all_sums$Mutation %<>% relevel(ref="Protospacer-associated")
all_sums$Mutation %<>% relevel(ref="None")


Fig5C <- ggplot(aes(x=Mutation, y=Mean), data=all_sums)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.8, width=0)+
  coord_cartesian(ylim=c(0,1))+
  #ggtitle("Proportion of hosts that were\ninfected by or resisted a phage \nwith a mutation")+
  labs(x="Mutation", y="Proportion of hosts infected")+
  scale_x_discrete(breaks=c("None", "Protospacer-associated"),
                   labels=c("None", "Protospacer-\nassociated"))+
  theme_cowplot()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(face="bold", size=16))+
  NULL
Fig5C

#### Fig 5D

ggsave("Fig1.png", Fig5C, path="~/Desktop/", dpi=300, device="png",
       width=15, height=12, unit=c("cm"))

proto_targets <- read.csv("./phage_sequences/summary_data/infections_targetted_sum.csv")
proto_targets$NotTargetted %<>% as.factor()

Fig5D <- ggplot(aes(x=NotTargetted, y=Mean), data=proto_targets)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=.8, width=0)+
  coord_cartesian(ylim=c(0,1))+
  #ggtitle("Proportion of hosts that were\ninfected by or resisted a phage \nwith a mutation")+
  labs(x="Number of host CRISPR spacers phage\nhad not evolved SNPs against", y="Probability of infection")+
  #facet_wrap(~Location,scales = "free_x")+
  theme_cowplot()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(face="bold", size=16))+
  NULL

#### Make and save Fig 5
Fig5 <- plot_grid(Fig5A+xlab(""), Fig5B+ylab(""), 
                  Fig5C, Fig5D+ylab(""),
                  ncol=2, nrow=2,
                  rel_widths = c(.9,1,.9,1), align = "hv",
                  labels = c("A", "B", "C", "D"))
#Fig5

ggsave("Fig5.png", Fig5, path="./figs/paper/",
       device="png", dpi=300, width=30, height=25, units = c("cm"))

### Figure 6 - Time-shift assay results ####
# Overall - phage ##
# timeshift_overall <- read.csv("./time_shift/summary_data/timeshift_means.csv")
# 
# timeshift_overall$Environment %<>% relevel(ref="Future")
# timeshift_overall$Environment %<>% relevel(ref="Present")
# timeshift_overall$Environment %<>% relevel(ref="Past")
# 
# coevo_infect_plot <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=timeshift_overall)+
#   geom_point(position = position_dodge(.5),
#              size=1.5)+
#   geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
#                 position = position_dodge(.5),
#                 width=0, size=1)+
#   geom_path(stat="identity", size=.8, linetype=2)+
#   theme_cowplot()+
#   labs(x="Host background", y="Proportion of\nhosts infected")+
#   theme(axis.title = element_text(face="bold", size=16))+
#   theme(axis.text = element_text(size=14))+
#   theme(legend.title = element_text(face='bold', size=14))+
#   theme(legend.title.align = 0.5)+
#   theme(legend.position = 'right')+
#   theme(legend.key.width = unit(2, 'cm'))+
#   theme(legend.key.height = unit(1, 'cm'))+
#   theme(legend.text = element_text(size=14))+
#   
#   coord_cartesian(ylim=c(0,1))+
#   scale_y_continuous(breaks=c(seq(0, 1, 0.25)))

# Overall - host ##
timeshift_host <- read.csv("./time_shift/summary_data/timeshift_means_phage_bkgrnd.csv")
timeshift_host$Environment %<>% relevel(ref="Future")
timeshift_host$Environment %<>% relevel(ref="Present")
timeshift_host$Environment %<>% relevel(ref="Past")

Fig6A <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=timeshift_host)+
  geom_point(position = position_dodge(.5),
             size=1.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_cowplot()+
  labs(x="Phage background", y="Proportion of\nhosts infected")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+
  coord_cartesian(ylim=c(0,1))

# By timepoint - phage ##
timepoint_TS_contrasts <- read.csv("./time_shift/summary_data/timepoint_contrasts.csv")
#infect_plot <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=timepoint_TS_contrasts)+
#  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
#           width=.5)+
#  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
#                position = position_dodge(.5),
#                width=0, size=1)+
#  geom_path(stat="identity")+
#  theme_cowplot()+
#  labs(x="Host background", y="Proportion of\nhosts infected")+
#  scale_x_discrete(breaks=c("t1", "t4", "t9"),
#                   labels=c("1", "4", "9"))+
#  coord_cartesian(ylim=c(0,1))+
#  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
#  scale_fill_discrete(name="Phage or host\nbackground",
  #                     breaks=c("t1", "t4", "t9"),
  #                     labels=c("1", "4", "9"))+
  # 
  # theme(axis.title = element_text(face="bold", size=16))+
  # theme(axis.text = element_text(size=14))+
  # theme(legend.title = element_text(face='bold', size=14))+
  # theme(legend.title.align = 0.5)+
  # theme(legend.position = 'right')+
  # theme(legend.key.width = unit(1, 'cm'))+
  # theme(legend.key.height = unit(1, 'cm'))+
  # theme(legend.text = element_text(size=14))


# By timepoint - host ##
Fig6B <- ggplot(aes(y=1-Mean.Resist, x=Phage, Group=Host), data=timepoint_TS_contrasts)+
  geom_bar(stat="identity",aes(fill=Host), position = position_dodge(.5),
           width=.5, colour="grey40")+
  geom_errorbar(aes(ymin=1-Resist.Lower, ymax=1-Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity")+
  theme_cowplot()+
  labs(x="Phage timepoint", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Host\ntimepoint",
                      breaks=c("t1", "t4", "t9"),
                      labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))

# Phage FSD scores ##
#FSD_means <- read.csv("./time_shift/summary_data/FSD_means.csv")
#FSD_phage <- ggplot(aes(x=Timepoint, y=Phage), data=FSD_means)+
#  geom_point(size=3)+
#  geom_errorbar(aes(ymin=Phage.Lower, ymax=Phage.Upper),
#                size=.8, width=0)+
#  theme_cowplot()+
#  labs(y="MS(GxE)/MS(E) phage")+
#  scale_x_discrete(name="Timepoint",
#                   breaks=c("t1", "t4", "t9"),
#                   labels=c("1", "4", "9"))+
  
#  theme(axis.title = element_text(face="bold", size=16))+
#  theme(axis.text = element_text(size=14))+
#  theme(plot.title = element_text(hjust=0.5, face="bold"))+
#  coord_cartesian(ylim=c(0,0.3727621))+
#  NULL
#FSD_phage


# Host FSD scores ##
#FSD_host <- ggplot(aes(x=Timepoint, y=Host), data=FSD_means)+
#  geom_point(size=3)+
#  geom_errorbar(aes(ymin=Host.Lower, ymax=Host.Upper),
#                size=.8, width=0)+
#  theme_cowplot()+
#  labs(y="MS(GxE)/MS(E) host")+
#  scale_x_discrete(name="Timepoint",
#                   breaks=c("t1", "t4", "t9"),
#                   labels=c("1", "4", "9"))+
  
#  theme(axis.title = element_text(face="bold", size=16))+
#  theme(axis.text = element_text(size=14))+
#  theme(plot.title = element_text(hjust=0.5, face="bold"))+
#  coord_cartesian(ylim=c(0, 0.3727621))+
#  NULL
#FSD_host

# Make Fig6 ##

#Fig6_left <- plot_grid(coevo_infect_plot, 
#                       infect_plot+
#                       theme(legend.position = c(0.65, -.25),
#                                legend.background = element_blank(),
#                                legend.direction = "horizontal",
#                                plot.margin = unit(c(0, 0, 2, 0), "cm")),
#                       ncol=1, align = "hv", axis = "l", labels=c("A", "C"))

#Fig6_right <- plot_grid(coevo_resist_plot, 
#                        resist_plot+
#                            theme(legend.position = "none",
#                                  plot.margin = unit(c(0, 0, 2, 0), "cm")), 
#                  ncol=1, align = "hv", axis = "l", labels = c("B", "D"))

#Fig6 <- plot_grid(Fig6_left, Fig6_right,
#                  ncol=2, align = "hv", axis="l")

Fig6 <- plot_grid(Fig6A,
                  Fig6B+ylab(""),
                  rel_widths = c(1,1.25), labels=c("A", "B"))

Fig6

ggsave("Fig6.png", Fig6, path="./figs/paper/",
       device="png", dpi=300, width=30, height=12, units = c("cm"))




