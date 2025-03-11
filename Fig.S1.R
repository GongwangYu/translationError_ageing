
library(dplyr)
library(plyr)
library(readr)
library(parallel)
library(tidyr)
library(ggplot2)
library(data.table)

# Fig.S1a -------------------------------------------------------------------
load("rt_common.RData")
load("readcount_yeastPositiveSup45ts.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$wt<- composite(readcount[ names(readcount)[1:2]])
compRep$mut<- composite(readcount[ names(readcount)[3:4]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

wt<-RTscore$wt %>%dplyr::select(transcript,wt=rte_ext)
mut<-RTscore$mut %>% dplyr::select(transcript,mut=rte_ext)
dfmrg<-merge(wt,mut)%>%dplyr::select(2,3) 
dfmrg<-dfmrg%>%filter(wt<0.5 & mut<0.5 )

wilcox.test(dfmrg$mut,dfmrg$wt,paired = T, alternative = "g")$p.value
#n=53,p-value = 1.300441e-10

range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(wt > thres | mut > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(mut, wt, paired=T,alternative = "g")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

#dfplot$variable<-factor(dfplot$variable,levels = c("wt","mut"))

Fig.S1a<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("wt","mut"),
                    labels = c("WT","Mutant"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0.0000,0.2) )+
  scale_x_discrete(labels = c(0.001,0.01,0.02,0.05), 
                   breaks = c(0.001,0.01,0.02,0.05))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "1.3×" * 10**-10),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("1.3×" * 10**-10),size=2)+
  annotate("text",x=3,y=0.2,label=expression("1.3×" * 10**-10),size=2)+
  annotate("text",x=4,y=0.2,label=expression("1.9×" * 10**-10),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "53)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "53)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "53)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "52)" ),size=2)

Fig.S1a





# Fig.S1b -------------------------------------------------------------------
load("rt_common.RData")
load("readcount_yeastParo.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$paro<- composite(readcount[ names(readcount)[1:3]])
compRep$ypd<- composite(readcount[ names(readcount)[7:9]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

paro<-RTscore$paro %>%dplyr::select(transcript,paro=rte_ext)
ypd<-RTscore$ypd %>% dplyr::select(transcript,ypd=rte_ext)
dfmrg<-merge(paro,ypd)%>%dplyr::select(2,3)%>%filter(paro<0.5 & ypd<0.5 )

wilcox.test(dfmrg$paro,dfmrg$ypd,paired = T, alternative = "g")$p.value
#n=328,p-value = 1.245647e-44

range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(paro > thres | ypd > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(paro, ypd, paired=T,alternative = "g")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

dfplot$variable<-factor(dfplot$variable,levels = c("ypd","paro"))

Fig.S1b<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("ypd","paro"),
                    labels = c("Control","Paromomycin"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0.0000,0.2) )+
  scale_x_discrete(labels = c(0.001,0.01,0.02,0.05), 
                   breaks = c(0.001,0.01,0.02,0.05))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "1.2×" * 10**-44),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("7.5×" * 10**-44),size=2)+
  annotate("text",x=3,y=0.2,label=expression("7.1×" * 10**-43),size=2)+
  annotate("text",x=4,y=0.2,label=expression("1,2×" * 10**-39),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "314)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "301)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "290)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "261)" ),size=2)

Fig.S1b

# 合并图 ---------------------------------------------------------------------
library(cowplot)

Fig.S1 <- plot_grid(Fig.S1a,Fig.S1b,ncol = 2, labels = c("a","b"))  

ggsave("Fig.S1.pdf",Fig.S1,width =5, height =4)
