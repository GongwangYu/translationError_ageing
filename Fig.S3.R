


# fig.s3a -----------------------------------------------------------------
#human skeletal muscle
load("readcount_humanTissueAge.RData")
load("rt_common.RData")

cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[1:3]])
compRep$old<- composite(readcount[ names(readcount)[4:5]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg<-merge(young,old)%>%dplyr::select(2,3) %>% filter(young<0.5 & old <0.5)
wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#n=114,p-value = 0.0002524633


range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig.S3a<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
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
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "2.5×" * 10**-4),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("2.2×" * 10**-4),size=2)+
  annotate("text",x=3,y=0.2,label=expression("2.1×" * 10**-4),size=2)+
  annotate("text",x=4,y=0.2,label=expression("1.7×" * 10**-4),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "104)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "93)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "86)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "70)" ),size=2)

Fig.S3a



#yeast
load("readcount_yeastReplicAge.RData")
load("rt_common.RData")

cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[1:3]])
compRep$old<- composite(readcount[ names(readcount)[4:6]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg<-merge(young,old)%>%dplyr::select(2,3) %>% filter(young<0.5 & old <0.5)
wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#n=967,p-value = 6.303939e-23


range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig.S3b<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
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
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "6.3×" * 10**-23),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("8.1×" * 10**-23),size=2)+
  annotate("text",x=3,y=0.2,label=expression("5.9×" * 10**-23),size=2)+
  annotate("text",x=4,y=0.2,label=expression("4.4×" * 10**-19),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "962)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "957)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "927)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "743)" ),size=2)

Fig.S3b

# 合并图 ---------------------------------------------------------------------
library(cowplot)

Fig.S3 <- plot_grid(Fig.S3a,Fig.S3b,ncol = 2, labels = c("a","b"))  

ggsave("Fig.S3.pdf",Fig.S3,width =5, height =4)

