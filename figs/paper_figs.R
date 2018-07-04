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

pd <- position_dodge(0.1)

### Figure 1 - Host & phage population dynamics ####
# Phage data ####
raw_titres <- read.csv("./phage_survival/original_data/phage_titre_Jack.csv")
str(data)
raw_titres$Replicate %<>% as.factor()
raw_titres$Timepoint %<>% as.factor()
raw_titres$Treatment %<>% as.factor()

# 10^9 plot ####
nine_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                    data=subset(raw_titres, Treatment == '9'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Transfer', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{9}*" Treatment")))+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))
nine_plot

# 10^8 plot ####
eight_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(raw_titres, Treatment == '8'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Transfer', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{8}*" Treatment")))+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))
eight_plot

# 10^7 plot ####
seven_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(raw_titres, Treatment == '7'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Transfer', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{7}*" Treatment")))+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))
seven_plot

# 10^6 plot ####
six_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                   data=subset(raw_titres, Treatment == '6'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd, size=1)+
  
  labs(x='Transfer', y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  ggtitle(expression(bold("10"*{}^{6}*" Treatment")))+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  
  theme(axis.text = element_text(size=12))
six_plot

# Arrange phage plots ####
all_titres <- plot_grid(nine_plot+labs(x="", y="")+
                          theme(axis.title = element_text(size=21),
                                axis.text = element_text(size=16),
                                plot.title = element_text(size=28)),
                        eight_plot+labs(x="", y="")+
                          theme(axis.title = element_text(size=21),
                                axis.text = element_text(size=16),
                                plot.title = element_text(size=28)),
                        seven_plot+xlab("Transfer")+
                          theme(axis.title.x = element_text(size=25, hjust = 1.15, margin = margin(t=20)),
                                axis.title = element_text(size=21),
                                axis.title.y = element_text(margin=margin(r=20), hjust=1.5),
                                axis.text = element_text(size=16),
                                plot.title = element_text(size=28)),
                        six_plot+labs(x="", y="")+
                          theme(axis.title = element_text(size=21),
                                axis.text = element_text(size=16),
                                plot.title = element_text(size=28)),
                        nrow=2, ncol=2, align = "hv", labels = c("i.", "ii.", "iii.", "iv."), 
                        label_size = 28, label_colour = "red", label_x = .05)
# Host data ####
# Arrange and save Figure 1 ####
Fig1 <- plot_grid(all_titres, #all_cfu,
                  labels=c("A", "B"), label_size = 30,
                  nrow=2)

ggsave("Fig1.png", Fig1, path="./figs/paper/",
       dpi=300, device="png",
       width=54, height=75, units = c("cm"))

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
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))

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

### Figure 3 - Spacer and phage diversity ####
# The procedure to generate this figure is rather long, so see
# ./sequences/code/coverage_analysis.R

### Figure 4 - Time-shift assay results ####
# Overall - phage ####
timeshift_overall <- read.csv("./time_shift/summary_data/timeshift_means.csv")

timeshift_overall$Environment %<>% relevel(ref="Future")
timeshift_overall$Environment %<>% relevel(ref="Present")
timeshift_overall$Environment %<>% relevel(ref="Past")

coevo_infect_plot <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=timeshift_overall)+
  geom_point(position = position_dodge(.5),
             size=1.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_cowplot()+
  labs(x="Host background", y="Proportion of\nhosts infected")+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.25)))

# Overall - host ####
coevo_resist_plot <- ggplot(aes(y=Mean.Resist, x=Environment, group=Group), data=timeshift_overall)+
  geom_point(position = position_dodge(.5),
             size=1.5)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  theme_cowplot()+
  labs(x="Host background", y="Proportion of\nhosts resistant")+
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

# By timepoint - phage ####
timepoint_TS_contrasts <- read.csv("./time_shift/summary_data/timepoint_contrasts.csv")
infect_plot <- ggplot(aes(y=Mean.Infect, x=Host, Group=Phage), data=timepoint_TS_contrasts)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity")+
  theme_cowplot()+
  labs(x="Host background", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\nbackground",
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


# By timepoint - host ####
resist_plot <- ggplot(aes(y=Mean.Resist, x=Host, Group=Phage), data=timepoint_TS_contrasts)+
  geom_bar(stat="identity",aes(fill=Phage), position = position_dodge(.5),
           width=.5)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity")+
  theme_cowplot()+
  labs(x="Host background", y="Proportion of\nhosts resistant")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=c(seq(0,1,0.25)))+
  scale_fill_discrete(name="Phage\nbackground",
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

# Phage FSD scores ####
FSD_means <- read.csv("./time_shift/summary_data/FSD_means.csv")
FSD_phage <- ggplot(aes(x=Timepoint, y=Phage), data=FSD_means)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Phage-(1.96*Phage.SE), ymax=Phage+(1.96*Phage.SE)),
                size=.8, width=0)+
  theme_cowplot()+
  labs(y="MS(GxE)/MS(E) phage")+
  scale_x_discrete(name="Timepoint",
                   breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))+
  coord_cartesian(ylim=c(0,0.21166618))
FSD_phage


# Host FSD scores ####
FSD_host <- ggplot(aes(x=Timepoint, y=Host), data=FSD_means)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Host-(1.96*Host.SE), ymax=Host+(1.96*Host.SE)),
                size=.8, width=0)+
  theme_cowplot()+
  labs(y="MS(GxE)/MS(E) host")+
  scale_x_discrete(name="Timepoint",
                   breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))+
  coord_cartesian(ylim=c(0, 0.21166618))
FSD_host

# Make Fig4 ####

Fig4_left <- plot_grid(coevo_infect_plot, 
                       infect_plot+
                          theme(legend.position = c(0.5, -.35),
                                legend.background = element_blank(),
                                legend.direction = "horizontal",
                                plot.margin = unit(c(0, 0, 2, 0), "cm")), FSD_phage,
                       ncol=1, align = "hv", axis = "l", labels=c("A", "C", "E"))

Fig4_right <- plot_grid(coevo_resist_plot, 
                        resist_plot+
                            theme(legend.position = "none",
                                  plot.margin = unit(c(0, 0, 2, 0), "cm")), 
                        FSD_host,
                  ncol=1, align = "hv", axis = "l", labels = c("B", "D", "F"))

Fig4 <- plot_grid(Fig4_left, Fig4_right,
                  ncol=2, align = "hv", axis="l")
Fig4

ggsave("Fig4.png", Fig4, path="./figs/paper/",
       device="png", dpi=300, width=25, height=30, units = c("cm"))
#### Figure 5 - Phage and host genotype dynamics ####

#### Relative frequencies of spacers ####
unique_spacers <- read.csv("./sequences/summary_data/unique_spacers_noseqs.csv")

unique_spacers$RelativeFrequency <- unique_spacers$n/12

unique1 <- unique_spacers %>% 
  group_by(Timepoint, SpacerMiddle) %>% 
  slice(1L)

unique1$SpacerMiddle %<>% as.factor()

unique1$group[unique1$SpacerMiddle %in% c("883", "10732", "23585", "32567")] <- c("0.08", "",
                                                                                  "", "",
                                                                                  "0.83","1", 
                                                                                  "1", "0.08")

spacer_RelFreq <- ggplot(aes(y=RelativeFrequency, x=Timepoint, Group=SpacerMiddle), data=unique1)+
  geom_point(stat="identity", aes(colour=SpacerMiddle), position=position_dodge(.3), size=2.5)+
  geom_path(aes(group=SpacerMiddle), linetype=2, position = position_dodge(.3), size=.7)+
  geom_text(aes(label=group), position = position_identity(), size=5, face="bold")+
  #facet_wrap(~Replicate)+
  #geom_text(aes(label=group))+
  theme_cowplot()+
  labs(x="Timepoint", y="Frequency")+
  scale_x_discrete(breaks=c("t1", "t4", "t9"),
                   labels=c("1", "4", "9"))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.title.y = element_text(margin=margin(r=10)),
        axis.text = element_text(size=16))

spacer_RelFreq

Fig5 <- plot_grid(spacer_RelFreq, #phage_RelFreq
                  nrow=2, align = "hv", labels = c("A", "B"))
Fig5

ggsave("Fig5.png", Fig5, path="./figs/paper/",
       device="png", dpi=300, width=25, height=30, units = c("cm"))

