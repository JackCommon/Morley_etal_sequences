###############
### titre_figs.r
## Created: 18/5/18 by Jack Common

rm(list=ls())

### Dependencies ####
library(scales)
library(ggplot2)
library(cowplot)

pd <- position_dodge(0.1)

### Set up dataframe to add in raw values in long form ####
treat9 <- rep(9,31) %>% 
  rep(12)
treat8 <- rep(8,31) %>% 
  rep(12)
treat7 <- rep(7,31) %>% 
  rep(12)
treat6 <- rep(6,31) %>% 
  rep(12)

Treatment <- c(treat9, treat8, treat7, treat6)

rm(treat9, treat8, treat6, treat7)

Replicate <- NULL
for (i in seq(1,12)){
  Replicate <- c(Replicate, rep(i,31))
}
Replicate <- rep(Replicate,4)

Timepoint <- seq(0,30,1) %>% 
  rep(12*4)

pfu <- rep("", 1488)

data <- data.frame(Replicate, Timepoint, Treatment, pfu)
data

write.csv(file = "./phage_survival/original_data/phage_titre_Jack_nophage.csv", x=data, row.names = F)

### Read in formatted data ####
data <- read.csv("./phage_survival/original_data/phage_titre_Jack.csv")
str(data)
data$Replicate %<>% as.factor()
data$Timepoint %<>% as.factor()
data$Treatment %<>% as.factor()


### 10^9 plot ####
nine_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(data, Treatment == '9'))+
  
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

### 10^8 plot ####
eight_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                    data=subset(data, Treatment == '8'))+
  
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

### 10^7 plot ####
seven_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(data, Treatment == '7'))+
  
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

### 10^6 plot ####
six_plot <- ggplot(aes(y=pfu, x=Timepoint, group=Replicate), 
                     data=subset(data, Treatment == '6'))+
  
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

### Arrange and save plots ####
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
                        nrow=2, ncol=2, align = "hv")
quartz()
all_titres

detach(package:cowplot)
ggsave("all_titres.png", all_titres, path="./figs/phage/",
       device="png", dpi=300,
       height=30, width=54, unit=c("cm"))
