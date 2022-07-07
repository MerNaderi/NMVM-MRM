
rm(list=ls(all=TRUE))
library(gridExtra)
library(cowplot)
library(egg)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(grid)
library(gridtext)
#----------------------------------------------------------------
#BIC <- read.delim("D:/Dropbox/BIC.txt")
#table(apply(BIC[, c(5, 6)], 1, which.min))

mydata1<- c(0, 2, 2, 100, 95, 51,
            98, 0, 9, 100, 100, 97, 
            98, 91, 0, 100, 100, 97,
            0, 0, 0, 0, 0, 0,
            5, 0, 0, 100, 0, 3,
            49, 3, 3, 100, 97, 0)

#----------------------heatmap code--------------------------------------------------
cormat <- matrix(mydata1, ncol = 6, nrow = 6, byrow=T)
colnames(cormat) <- c('GHST-MRM','NIG-MRM','NMVBS-MRM', 'SCN-MRM', 'SSL-MRM', 'ST-MRM')
rownames(cormat) <- c('GHST-MRM','NIG-MRM','NMVBS-MRM', 'SCN-MRM', 'SSL-MRM', 'ST-MRM')
melted_cormat <- melt(cormat)

ggheatmap=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) + ggtitle("BIC") +
  geom_tile(color = "white")+ 
  scale_fill_gradient2(low = "white", high = "steelblue") + coord_fixed()+
  theme_minimal()+ theme(legend.title=element_blank())+
  theme(axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(angle =45, vjust = 1,size = 11, hjust = 1),
        axis.title.x=element_blank(),axis.title.y = element_blank(),
        legend.position="none")
#---------------------------heatmap and data--------------------------------
heatmap_1<-ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.ticks = element_blank())  
#==============================barplot=========================
d1=melted_cormat%>%
  group_by(Var1) %>% 
  summarize(value = sum(value)) 

bar_1 <- ggplot(data = d1, aes(x = Var1, y = value)) + 
  scale_fill_distiller(name = "Value", palette = "Blues", direction =1)+
  geom_bar(stat = "identity", aes(fill = value),width = 0.8)+ coord_flip() + theme_bw() +
  geom_text(aes(label = value), hjust = -0.2, size=5)+expand_limits(y = c(0,650))+
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        plot.margin = margin(0,0,0,-15), legend.position="none")
g1=ggarrange(
  heatmap_1,bar_1,widths = c(20,10), heights = c(20,20), byrow=FALSE)




mydata2 <- c(0, 60, 62, 67, 65, 59,
             40, 0, 77, 63, 62, 56, 
             38, 23, 0, 56, 58, 52,
             33, 37, 44, 0, 52, 38,
             35, 38, 42, 48, 0, 52, 
             41, 44, 48, 62, 48, 0)


#----------------------heatmap code--------------------------------------------------
cormat <- matrix(mydata2, ncol = 6, nrow = 6,byrow=T)
colnames(cormat) <- c('GHST-MRM','NIG-MRM','NMVBS-MRM', 'SCN-MRM', 'SSL-MRM', 'ST-MRM')
rownames(cormat) <- c('GHST-MRM','NIG-MRM','NMVBS-MRM', 'SCN-MRM', 'SSL-MRM', 'ST-MRM')
melted_cormat <- melt(cormat)

ggheatmap=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+ ggtitle("MCR") +
  geom_tile(color = "white")+ 
  scale_fill_gradient2(low = "white", high = "steelblue") + coord_fixed()+
  theme_minimal()+ theme(legend.title=element_blank())+
  theme(axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(angle =45, vjust = 1,size = 11, hjust = 1),
        axis.title.x=element_blank(),axis.title.y = element_blank(),
        legend.position="none")
#---------------------------heatmap and data--------------------------------
heatmap_1<-ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.ticks = element_blank())  
#==============================barplot=========================
d1=melted_cormat%>%
  group_by(Var1) %>% 
  summarize(value = sum(value)) 

bar_1 <- ggplot(data = d1, aes(x = Var1, y = value)) + 
  scale_fill_distiller(name = "Value", palette = "Blues", direction =1)+
  geom_bar(stat = "identity", aes(fill = value),width = 0.8)+ coord_flip() + theme_bw() +
  geom_text(aes(label = value), hjust = -0.2, size=5)+expand_limits(y = c(0,450))+
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        plot.margin = margin(0,0,0,-15), legend.position="none")
g1=ggarrange(
  heatmap_1,bar_1,widths = c(20,10), heights = c(20,20), byrow=FALSE)

