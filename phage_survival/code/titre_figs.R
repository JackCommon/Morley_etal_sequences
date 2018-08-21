###############
### titre_figs.r
## Created: 18/5/18 by Jack Common

rm(list=ls())

### Dependencies ####
library(scales)
library(ggplot2)
library(cowplot)
library(grid)
library(gridGraphics)
library(ggplotify)
library(reshape2)

pd <- position_dodge(0.1)

#### Functions ####
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

treatment_facet <- list(
  '6'  = expression(bold('10'^6*'')),
  '7' = expression(bold('10'^7*'')),
  '8' = expression(bold('10'^8*'')),
  '9' = expression(bold('10'^9*''))
)
treatment_labeller = function(variable, value) {
  return(treatment_facet[value])
}


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

write.csv(file = "./phage_survival/original_data/phage_titre_Jack_nophage.csv", x=data, row.names = F)

### Read in formatted data ####
data <- read.csv("./phage_survival/original_data/phage_titre_Jack.csv")
str(data)
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


### 10^9 plot ####
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
quartz()
nine_plot

### 10^8 plot ####
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
quartz()
eight_plot

### 10^7 plot ####
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
quartz()
seven_plot

### 10^6 plot ####
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
quartz()
six_plot


#### Summary model ####
library(multcomp)
library(relaimpo)
library(ggfortify)
library(survminer)

data <- read.csv("./phage_survival/original_data/binomial_survival_data.csv", header = TRUE)
data$Treatment %<>% as.factor

model <- glm(cbind(Alive, Dead)~Treatment*Time,family=binomial, data=data)

range(data$Time)
n<-seq(0,30,length.out=1000)
length(n)
intercept_TR1<-model$coef[1]
intercept_TR2<-model$coef[1]+model$coef[2]
intercept_TR3<-model$coef[1]+model$coef[3]
intercept_TR4<-model$coef[1]+model$coef[4]

slope_TR1 <- model$coef[5]
slope_TR2 <- model$coef[5]+model$coef[6]
slope_TR3 <- model$coef[5]+model$coef[7]
slope_TR4 <- model$coef[5]+model$coef[8]

fitted_TR1<-exp(intercept_TR1+slope_TR1*n)/(1+exp(intercept_TR1+slope_TR1*n))
fitted_TR2<-exp(intercept_TR2+slope_TR2*n)/(1+exp(intercept_TR2+slope_TR2*n))
fitted_TR3<-exp(intercept_TR3+slope_TR3*n)/(1+exp(intercept_TR3+slope_TR3*n))
fitted_TR4<-exp(intercept_TR4+slope_TR4*n)/(1+exp(intercept_TR4+slope_TR4*n))

#### Binomial summary fig ####
png("./figs/phage/binom_plot.png", width=20, height=15, units="in", res=300)
plot.new()
par(mfrow=c(1,1))
par(las=1)
mar.default <- c(8,8,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

plot(data$Time, data$Alive/(data$Alive+data$Dead),axes=F, type="n",ylab = "",xlab = "", xlim=c(0,30), cex.axis=2.5, cex.lab=2.5)
axis(side = 1, lwd = 3, cex.axis=2.5, pos = c(-.05,0), tck= -.02)
axis(side = 2, lwd = 3, cex.axis=2.5, pos = c(-.15,1))
title(ylab="Proportion of replicates\nwith surviving phage",line=5, cex.lab=4, tcl=0)
title(xlab="Days post-infection (d.p.i.)",line=4, cex.lab=4, tcl=0)

points(fitted_TR1~n, type="l", lwd=2.5,lty=1, col="#718B2E")
points(fitted_TR2~n, type="l", lwd=2.5,lty=1, col= "#D39200")
points(fitted_TR3~n, type="l", lwd=2.5,lty=1, col= "red")
points(fitted_TR4~n, type="l", lwd=2.5,lty=1, col= "blue")

legend(20, 0.9, c(expression("10"*{}^{6}*"")) ,lty = c(1), col=c('#D39200'),lwd=2.5, cex=4,bty='n')
legend(20, 0.8, c(expression("10"*{}^{7}*"")) ,lty = c(1), col=c('#718B2E'),lwd=2.5, cex=4,bty='n')
legend(20, 0.7, c(expression("10"*{}^{8}*"")) ,lty = c(1), col=c('red'),lwd=2.5, cex=4,bty='n')
legend(20, 0.6, c(expression("10"*{}^{9}*"")),lty = c(1), col=c('blue'),lwd=2.5, cex=4,bty='n')
dev.off()

binom_plot <- recordPlot()
binom_plot %<>% plot_to_gtable()
plot.new()
plot(binom_plot)
### Arrange and save plots ####
all_titres <- plot_grid(nine_plot+labs(x=""), 
                        eight_plot+labs(x=""),
                        seven_plot, 
                        six_plot,
                        labels = c("A", "B", "C", "D"), label_colour = "red")
quartz()
all_titres

full <- plot_grid(binom_plot, all_titres,
                  ncol=1, nrow=2)
quartz()
full

detach(package:cowplot)
ggsave("all_titres.png", all_titres, path="./figs/phage/",
       device="png", dpi=300,
       height=40, width=54, unit=c("cm"))

### OD600 and CFU models ####
comp <- read.csv("./phage_survival/original_data/OD_cfu.csv")
comp$Replicate %<>% as.factor()
#comp$Timepoint %<>% as.factor()
comp$Treatment %<>% as.factor()

comp %<>% filter(Timepoint %in% c(1,2,3,4))
comp$OD600 <- ifelse(comp$OD600<0,0,comp$OD600)

# Simple model of OD600 as a function of CFU
data_pos <- filter(data, pfu!=1)

m1 <- glm(log10(pfu)~sqrt(OD600), data=data, family=gaussian(link="identity"))
m2 <- glm(log10(pfu)~sqrt(OD600)+Treatment, data=data, family=gaussian(link="identity"))

summary(m1)
summary(m2)

par(mfrow=c(2,2))
plot(m1)
plot(m2)

anova(m1, test="Chisq")
anova(m2, test="F")

AIC(m1, m2) %>% compare_AICs()

quartz()
association_plot <- ggplot(aes(y=pfu, x=sqrt(OD600)), data=data)+
  geom_smooth(method="lm", formula = y~x, se=T, fullrange=T)+
  geom_point()+
  #stat_smooth(method="lm", formula = y~1, se=T, linetype=2, colour="grey2")+ # This gives a null slope
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  theme_cowplot()+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  labs(x=expression(bold(sqrt("OD"[600]))), y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  theme(strip.text.x = element_text(face="bold"))+
  NULL

ggsave("PFUxOD600.png", association_plot, path="./figs/phage/",
       device="png", dpi=300, width=25, height=15, units = c("cm"))
 