
# Fig.S2 -------------------------------------------------------------------

# 最终的代码 -------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_humanbone.csv")

dfqsubs_flt<-mclapply(1:nrow(dfqsubs), function(i){
  dfBP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"
  
  dfDP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Modification.intensity")))%>%t()%>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfDP)[2]<-"dp"
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  return(data.frame(young=young,old=old))
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()

wilcox.test(dfqsubs_flt$young,dfqsubs_flt$old, alternative = "l")
#N=96,p-value = 0.04618  cutoff=0.01


range(dfqsubs_flt)
misincorp <- lapply(c(0.0001,0.0002,0.0005,0.001),function(thres){
  dfdat<-dfqsubs_flt %>%dplyr::filter(young > thres | old > thres)  %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old,alternative = "l")$p.value)
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
