
library(dplyr)
library(plyr)
library(readr)
library(parallel)
library(tidyr)
library(ggplot2)
library(data.table)


# 最新 ----------------------------------------------------------------------

# fig2a -------------------------------------------------------------------
load("dfqsubs_flt_human.RData")#1127

range(dfqsubs_flt)
misincorp <- lapply(c(0.0001,0.0002,0.0005,0.001),function(thres){
  dfdat<-dfqsubs_flt %>%dplyr::filter(young > thres | old > thres)  %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()
misincorp$mw.p%>%unique()
misincorp$thres%>%table()

dfplot<-misincorp%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2a<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Misincorporation rates")+
  theme_test()+
  theme(plot.title=element_text(hjust=1),
        legend.position = "top",legend.key   = element_blank(),
        legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.005,0.01), labels = c("0","0.005","0.01"),
                     limits =c(0.0000,0.01) )+
  scale_x_discrete(breaks = c(10**-4,2*10**-4,5*10**-4,10**-3),
                   labels = c(expression( "1×"*10**-4),expression("2×"*10**-4),
                              expression("5×"*10**-4),expression("1×"*10**-3) ) )+
  annotate("text",x=1,y=0.008,label=expression(italic(P) * " = " * "1.1×" *10**-11 ),size=2)+
  annotate("text",x=2,y=0.008,label=expression("7.7×"*10**-11),size=2)+
  annotate("text",x=3,y=0.008,label=expression("4.6×"*10**-7),size=2)+
  annotate("text",x=4,y=0.008,label=expression("9.5×"*10**-4),size=2)+
  annotate("text",x=1,y=0.007,label=expression("(N" * " = " * "955)" ),size=2)+
  annotate("text",x=2,y=0.007,label=expression("(N" * " = " * "851)" ),size=2)+
  annotate("text",x=3,y=0.007,label=expression("(N" * " = " * "563)" ),size=2)+
  annotate("text",x=4,y=0.007,label=expression("(N" * " = " * "381)" ),size=2)

Fig2a


# fig2b -------------------------------------------------------------------
load("dfqsubs_flt_mouse.RData")#392

range(dfqsubs_flt)
misincorp <- lapply(c(0.0001,0.0002,0.0005,0.001),function(thres){
  dfdat<-dfqsubs_flt %>%dplyr::filter(young > thres | old > thres)  %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()
misincorp$mw.p%>%unique()
misincorp$thres%>%table()
dfplot<-misincorp%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2b<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Misincorporation rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.005,0.01), labels = c("0","0.005","0.01"),
                     limits =c(0.0000,0.01) )+
  scale_x_discrete(breaks = c(10**-4,2*10**-4,5*10**-4,10**-3),
                   labels = c(expression( "1×"*10**-4),expression("2×"*10**-4),
                              expression("5×"*10**-4),expression("1×"*10**-3) ) )+
  annotate("text",x=1,y=0.0095,label=expression(italic(P) * " = " * "2.04×" *10**-3 ),size=2)+
  annotate("text",x=2,y=0.0095,label=expression("2.02×"*10**-3),size=2)+
  annotate("text",x=3,y=0.0095,label=expression("1.86×"*10**-3),size=2)+
  annotate("text",x=4,y=0.0095,label=expression("1.09×"*10**-3),size=2)+
  annotate("text",x=1,y=0.0085,label=expression("(N" * " = " * "385)" ),size=2)+
  annotate("text",x=2,y=0.0085,label=expression("(N" * " = " * "380)" ),size=2)+
  annotate("text",x=3,y=0.0085,label=expression("(N" * " = " * "365)" ),size=2)+
  annotate("text",x=4,y=0.0085,label=expression("(N" * " = " * "322)" ),size=2)

Fig2b


# fig2c -------------------------------------------------------------------
load("dfqsubs_flt_yeast.RData")#57

range(dfqsubs_flt)
misincorp <- lapply(c(0.0001,0.0002,0.0005,0.001),function(thres){
  dfdat<-dfqsubs_flt %>%dplyr::filter(young > thres | old > thres)  %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()
misincorp$mw.p%>%unique()
misincorp$thres%>%table()
dfplot<-misincorp%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2c<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Misincorporation rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.005,0.01), labels = c("0","0.005","0.01"),
                     limits =c(0.0000,0.01) )+
  scale_x_discrete(breaks = c(10**-4,2*10**-4,5*10**-4,10**-3),
                   labels = c(expression( "1×"*10**-4),expression("2×"*10**-4),
                              expression("5×"*10**-4),expression("1×"*10**-3) ) )+
  annotate("text",x=1,y=0.008,label=expression(italic(P) * " = " * "0.015" ),size=2)+
  annotate("text",x=2,y=0.008,label=expression("0.018"),size=2)+
  annotate("text",x=3,y=0.008,label=expression("0.023"),size=2)+
  annotate("text",x=4,y=0.008,label=expression("0.025"),size=2)+
  annotate("text",x=1,y=0.007,label=expression("(N" * " = " * "53)" ),size=2)+
  annotate("text",x=2,y=0.007,label=expression("(N" * " = " * "47)" ),size=2)+
  annotate("text",x=3,y=0.007,label=expression("(N" * " = " * "43)" ),size=2)+
  annotate("text",x=4,y=0.007,label=expression("(N" * " = " * "37)" ),size=2)

Fig2c

# Fig2.d human ------------------------------------------------------------------

load("rt_common.RData")
load("readcount_humanCell.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0;
compRep <- list() # Composite replicates read count
names(readcount)
compRep$old<- composite(readcount[ names(readcount)[1:2]])
compRep$young<- composite(readcount[ names(readcount)[3:4]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3, exten_m5 )

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg<-merge(young,old)%>%dplyr::select(2,3)%>%
  filter(young<0.5 & old <0.5) 

wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#n=866,p-value = 2.333434e-53


range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2d<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0,0.2) )+
  scale_x_discrete(labels = c(0.001,0.01,0.02,0.05), 
                   breaks = c(0.001,0.01,0.02,0.05))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "2.4×" * 10**-53),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("3.7×" * 10**-53),size=2)+
  annotate("text",x=3,y=0.2,label=expression("3.5×" * 10**-52),size=2)+
  annotate("text",x=4,y=0.2,label=expression("2.9×" * 10**-49),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "856)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "843)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "789)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "664)" ),size=2)

Fig2d


# Fig2.e mouse ------------------------------------------------------------------

load("rt_common.RData")
load("readcount_mouseAge.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[26:28]])
compRep$old<- composite(readcount[ names(readcount)[19:21]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)
young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)

dfmrg<-merge(young,old)%>%dplyr::select(2,3) %>% filter(young<0.5 & old <0.5)

wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#n=557,p-value = 1.056362e-12

range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2e<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0.0000,0.2) )+
  scale_x_discrete(labels = c(0.001,0.01,0.02,0.05), 
                   breaks = c(0.001,0.01,0.02,0.05))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "1.1×" * 10**-12),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("4.7×" * 10**-11),size=2)+
  annotate("text",x=3,y=0.2,label=expression("1.3×" * 10**-8),size=2)+
  annotate("text",x=4,y=0.2,label=expression("2.3×" * 10**-5),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "532)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "458)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "341)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "174)" ),size=2)

Fig2e


# Fig2f-yeast -------------------------------------------------------------------

load("rt_common.RData")
load("readcount_yeastChrAge.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[1:4]])
compRep$old<- composite(readcount[ names(readcount)[9:12]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)

dfmrg<-merge(young,old)%>%dplyr::select(2,3) %>% filter(young<0.5 & old <0.5)
wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#n=786,p-value = 1.300069e-16


range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig2f<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("young","old"),
                    labels = c("Young","Aged"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0.0000,0.2) )+
  scale_x_discrete(labels = c(0.001,0.01,0.02,0.05), 
                   breaks = c(0.001,0.01,0.02,0.05))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "1.3×" * 10**-16),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("9.5×" * 10**-17),size=2)+
  annotate("text",x=3,y=0.2,label=expression("1.3×" * 10**-17),size=2)+
  annotate("text",x=4,y=0.2,label=expression("3.4×" * 10**-18),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "786)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "767)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "688)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "412)" ),size=2)
Fig2f


# Fig2g -------------------------------------------------------------------
load("dfsubs_rt.cor.RData")

cor.test(dfsubs_rt$young_misinRate,dfsubs_rt$young_RTrate,method = "s")
#n=28 ,rho=0.14,p=0.4755
Fig2g<-ggplot(dfsubs_rt,aes(x=young_misinRate,y=young_RTrate))+geom_point()+
  xlab("Misincorporation rate\nin young samples")+
  ylab("Readthrough rate\nin young samples")+
  theme_test()+
  scale_y_continuous(breaks=c(0.0,0.2,0.4), labels = c(0.0,0.2,0.4),
                     limits =c(0,0.4) )+
  scale_x_continuous(breaks=c(0,0.0005,0.001), labels = c(0.0,"0.0005",0.001),
                     limits =c(0,0.001) )+
  annotate("text",x=0.0005,y=0.39,label=expression(italic(rho) * " = " * 0.14))+
  annotate("text",x=0.0005,y=0.35,label=expression(italic(P) * " = " *0.48))

Fig2g
# Fig2h -------------------------------------------------------------------

load("dfsubs_rt.cor.RData")

cor.test(dfsubs_rt$old_misinRate,dfsubs_rt$old_RTrate,method = "s")
#rho=0.04107339,p=0.8356
Fig2h<-ggplot(dfsubs_rt,aes(x=old_misinRate,y=old_RTrate))+geom_point()+
  xlab("Misincorporation rate\nin aged samples")+
  ylab("Readthrough rate\nin aged samples")+
  theme_test()+
  scale_y_continuous(breaks=c(0.0,0.3,0.6), labels = c(0.0,0.3,0.6),
                     limits =c(0,0.6) )+
  scale_x_continuous(breaks=c(0,0.0005,0.001), labels = c(0.0,"0.0005",0.001),
                     limits =c(0,0.001) )+
  annotate("text",x=0.0005,y=0.59,label=expression(italic(rho) * " = " * 0.04))+
  annotate("text",x=0.0005,y=0.50,label=expression(italic(P) * " = " *0.83))

Fig2h

# Fig2i -------------------------------------------------------------------

load("dfsubs_rt.cor.RData")

cor.test(dfsubs_rt$sub_O_Yratio,dfsubs_rt$rt_O_Yratio,method = "s")
#rho=-0.06609742,p=0.7382
Fig2i<-ggplot(dfsubs_rt,aes(x=sub_O_Yratio,y=rt_O_Yratio))+geom_point()+
  xlab("Aged-to-young relative increase\nfor misincorporation rate")+
  ylab("Aged-to-young relative increase\nfor readthrough rate")+
  theme_test()+
  scale_y_continuous(breaks=c(-2,-1,0,1,2), labels = c(-2,-1,0,1,2),
                     limits =c(-2,2) )+
  scale_x_continuous(breaks=c(-2,-1,0,1,2), labels = c(-2,-1,0,1,2),
                     limits =c(-2,2) )+
  annotate("text",x=-1,y=2,label=expression(italic(rho) * " = " * -0.07))+
  annotate("text",x=-1,y=1.5,label=expression(italic(P) * " = " *0.74))
Fig2i
# 合并图 ---------------------------------------------------------------------
library(cowplot)

Fig2 <- plot_grid(Fig2a,Fig2b,Fig2c,Fig2d,Fig2e,Fig2f,Fig2g,Fig2h,Fig2i,ncol = 3,
                  labels = c("a","b","c","d","e","f","g","h","i"))  

ggsave("Fig2.pdf",Fig2,width =8, height =8)


# 结束 ----------------------------------------------------------------------

