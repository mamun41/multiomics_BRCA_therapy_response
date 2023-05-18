library (ggplot2)
library (MASS)
library (readxl)
library (reshape2)
library (UpSetR)
library (vcd)
library (ggpubr)

alldata <- read.csv("D:/DOWNLOADBROWER/NAT-ML-main/inputs/all.csv", header = TRUE)

figure_font_size<-18


#Figure 2a
#wilcox
fig2a <- 
  ggplot(alldata,aes(x=RCBcategory,y=Swanton.PaclitaxelScore,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="Swanton.PaclitaxelScore")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+

  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2a

fig2b <- 
  ggplot(alldata,aes(x=RCBcategory,y=ESC.ssgsea.notnorm,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="ESC.ssgsea.notnorm")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+

  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2b

fig2c <- 
  ggplot(alldata,aes(x=RCBcategory,y=GGI.ssgsea.notnorm,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="GGI.ssgsea.notnorm")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2c

fig2d <- 
  ggplot(alldata,aes(x=RCBcategory,y=Danaher.Mast.cells,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="Danaher.Mast.cells")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2d

fig2e <- 
  ggplot(alldata,aes(x=RCBcategory,y=ESR1,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="ESR1.log2.tpm")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2e

fig2f <- 
  ggplot(alldata,aes(x=RCBcategory,y=PGR,fill=RCBcategory))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="PGR.log2.tpm")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2f


