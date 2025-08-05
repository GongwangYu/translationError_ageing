library(dplyr)
library(plyr)
library(readr)
library(parallel)
library(tidyr)
library(ggplot2)
library(data.table)

# Fig1.b ------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_str.csv")

dfqsubs_flt<-lapply(1:nrow(dfqsubs), function(i){
  dfBP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"
  
  dfDP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Modification.intensity")))%>%t()%>%
    as.data.frame()%>%tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%dplyr::select(2,3)
  colnames(dfDP)[2]<-"dp"
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T))%>%
    filter(dfratio<0.01)
  
  wt<-dfmerg%>%filter(grepl("wt",sample,ignore.case=T))%>%getElement("dfratio")%>%mean(na.rm=T)
  ep<-dfmerg%>%filter(grepl("ep",sample,ignore.case=T))%>%getElement("dfratio")%>%mean(na.rm=T)
  
  return(data.frame(wt,ep))
  
})%>%rbind.fill()%>%filter(wt >0 | ep >0)%>%na.omit()

wilcox.test(dfqsubs_flt$wt,dfqsubs_flt$ep,paired = T, alternative = "l")
#N=58,p-value = 3.922e-11, cutoff=0.01

###
range(dfqsubs_flt)
misincorp <- lapply(c(0.0001,0.001,0.002,0.005),function(thres){
  dfdat<-dfqsubs_flt %>%dplyr::filter(wt > thres | ep > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(wt, ep, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()
misincorp$mw.p%>%unique()
misincorp$thres%>%table()
dfplot<-misincorp%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig1b<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("wt","ep"),
                    labels = c("Control","Streptomycin"))+
  labs(x="Threshold for error rates (>x)",y="Misincorporation rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.005,0.01), labels = c("0","0.005","0.01"),
                     limits =c(0.0000,0.01) )+
  scale_x_discrete(breaks = c(10**-4,10**-3,2*10**-3,5*10**-3),
                   labels = c(expression( "1×"*10**-4),expression("1×"*10**-3),
                              expression("2×"*10**-3),expression("5×"*10**-3) ) )+
  annotate("text",x=1,y=0.01,label=expression(italic(P) * " = " * "3.9×" * 10**-11 ),size=2)+
  annotate("text",x=2,y=0.01,label=expression("1.8×" * 10**-10 ),size=2)+
  annotate("text",x=3,y=0.01,label=expression("4.3×" * 10**-12),size=2)+
  annotate("text",x=4,y=0.01,label=expression("3.8×" * 10**-5),size=2)+
  annotate("text",x=1,y=0.0095,label=expression("(N" * " = " * "58)" ),size=2)+
  annotate("text",x=2,y=0.0095,label=expression("(N" * " = " * "53)" ),size=2)+
  annotate("text",x=3,y=0.0095,label=expression("(N" * " = " * "42)" ),size=2)+
  annotate("text",x=4,y=0.0095,label=expression("(N" * " = " * "17)" ),size=2)

Fig1b



# fig1d -------------------------------------------------------------------
load("rt_common.RData")
load("readcount_yeastPositive.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$Mut<- composite(readcount[ names(readcount)[1:2]])
compRep$WT<- composite(readcount[ names(readcount)[3:4]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5 )

WT<-RTscore$WT %>%dplyr::select(transcript,WT=rte_ext)
Mut<-RTscore$Mut %>% dplyr::select(transcript,Mut=rte_ext)

dfmrg<-merge(WT,Mut)%>%dplyr::select(2,3)%>%filter(WT<0.5 & Mut<0.5 )

median(dfmrg$WT) #0.003987241
median(dfmrg$Mut) #0.0441141
wilcox.test(dfmrg$WT,dfmrg$Mut,paired = T, alternative = "l")$p.value
#n=221,p-value = 4.387029e-38
range(dfmrg)

readthr <- lapply( c(0.001,0.002,0.005,0.01),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(WT > thres | Mut > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(WT, Mut, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()

dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig1d<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF", "#E64B35FF"),
                    breaks = c("WT","Mut"),
                    labels = c("WT","Mutant"))+
  labs(x="Threshold for error rates (>x)",y="Readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),
                     limits =c(0.0000,0.2) )+
  scale_x_discrete(labels = c(expression("1×"*10**-3),expression("2×"*10**-3),
                              expression("5×"*10**-3),expression("1×"*10**-2)), 
                   breaks = c(10**-3,2*10**-3,5*10**-3,10**-2))+
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "4.4×" * 10**-38),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("6.2×" * 10**-38),size=2)+
  annotate("text",x=3,y=0.2,label=expression("1.8×" * 10**-37),size=2)+
  annotate("text",x=4,y=0.2,label=expression("1.6×" * 10**-36),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "221)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "220)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "217)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "211)" ),size=2)

Fig1d

#figure ---------------------------------------------------------------------
library(cowplot)
Fig1<-plot_grid("","",Fig1b,"","",Fig1d,ncol = 3)  
ggsave("Fig1.pdf",Fig1,width =8, height =6)

