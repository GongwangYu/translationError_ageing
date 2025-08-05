library(data.table)
# common codes -------------------------------------------------------------------
load("rt_common.RData")

lvls_AA<-c("A","C","D","E","F","G","H","K","I/L",
           "M","N","P","Q","R","S","T","V","W","Y","allAA")

#Classification by AA and codon preference
lvls_codon<-c("allCodon","GCC","GCU","GCA","GCG","UGC","UGU","GAC","GAU","GAG","GAA","UUC","UUU",
              "GGC","GGA","GGG","GGU","CAC","CAU","AAG","AAA","AUC","AUU","AUA","CUG","CUC","CUU",
              "UUG","CUA","UUA","AUG","AAC","AAU","CCC","CCU", "CCA","CCG","CAG","CAA","CGG","AGA",
              "AGG","CGC","CGA","CGU","AGC","UCC","UCU","UCA","AGU","UCG","ACC","ACA","ACU","ACG",
              "GUG","GUC","GUU","GUA","UGG","UAC","UAU")

# Fig3a -------------------------------------------------------------------

load("misin.human_disc.heatmap.RData")

dfallMtr$destination<-factor(dfallMtr$destination,levels =lvls_AA)
dfallMtr$codon<-factor(dfallMtr$codon,levels = lvls_codon)

Fig3a<-ggplot(dfallMtr,aes(x=destination,y=codon,fill=effeSize))+
  geom_tile(aes(fill=effeSize))+
  geom_text(aes(label=stars),color="black",size=4,
            hjust="middle",vjust="bootom")+
  scale_fill_gradient2("Aged-to-young relative increase\nin translation error rates",
                       na.value="grey",low = "#4DBBD5FF",high = "#E64B35FF",mid = "white",
                       guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                            barwidth =4,barheight = 0.5, ticks.colour="black",
                                            ticks.linewidth = 1 ))+
  labs(x="Misincorporated amino acids",y="Original codon")+
  scale_y_discrete(breaks=lvls_codon,labels=c("",lvls_codon[-1]))+  
  scale_x_discrete(breaks=c("A","C","D","E","F","G","H","K","I/L",
                            "M","N","P","Q","R","S","T","V","W","Y","allAA"),
                   labels=c("A","C","D","E","F","G","H","K","I/L",
                            "M","N","P","Q","R","S","T","V","W","Y",""))+  
  theme( axis.ticks.y = element_blank(),
         panel.background=element_blank())+
  theme_test()+
  theme(legend.position = "top",legend.title  =element_text(size = 8))

Fig3a
# Fig3b -------------------------------------------------------------------

load("misin.mouse.heatmap.RData")
dfallMtr$destination<-factor(dfallMtr$destination,levels =lvls_AA)
dfallMtr$codon<-factor(dfallMtr$codon,levels = lvls_codon)

Fig3b<-ggplot(dfallMtr,aes(x=destination,y=codon,fill=effeSize))+
  geom_tile(aes(fill=effeSize))+
  geom_text(aes(label=stars),color="black",size=4,
            hjust="middle",vjust="bootom")+
  scale_fill_gradient2("Aged-to-young relative increase\nin translation error rates",
                       na.value="grey",low = "#4DBBD5FF",high = "#E64B35FF",mid = "white",
                       guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                            barwidth =4,barheight = 0.5, ticks.colour="black",ticks.linewidth = 1 ))+
  labs(x="Misincorporated amino acids",y="Original codon")+
  scale_y_discrete(breaks=lvls_codon,labels=c("",lvls_codon[-1]))+  
  scale_x_discrete(breaks=c("A","C","D","E","F","G","H","K","I/L",
                            "M","N","P","Q","R","S","T","V","W","Y","allAA"),
                   labels=c("A","C","D","E","F","G","H","K","I/L",
                            "M","N","P","Q","R","S","T","V","W","Y",""))+  
  theme( axis.ticks.y = element_blank(),
         panel.background=element_blank())+
  theme_test()+
  theme(legend.position = "top",legend.title  =element_text(size = 8))

Fig3b

# Fig3c/d/e -------------------------------------------------------------------
#human
load("rt_common.RData")
load("readcount_humanCell.RData")
load("stopCdnSeq_human.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$old<- composite(readcount[ names(readcount)[1:2]])
compRep$young<- composite(readcount[ names(readcount)[3:4]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3, exten_m5 )

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg_human<-merge(young,old)%>%merge(stopCdnSeq)%>%
  dplyr::mutate(relRatio=(old-young)/((old+young)/2),spec="human")%>%na.omit() %>%
  dplyr::filter(young < 0.5 & old <0.5)

summary(dfmrg_human$relRatio)

table(dfmrg_human$stopcdn)
#TAA TAG TGA 
#246 167 444 
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_human%>%filter(stopcdn=="TAG"))$young)
# p-value = 0.07285 
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$young)
# p-value = 0.3031
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAG"))$young,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$young)
#p-value = 0.002993 

wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_human%>%filter(stopcdn=="TAG"))$old)
#p-value = 0.0147
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$old)
#p-value = 0.501
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAG"))$old,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$old)
#p-value = 0.04252 

wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_human%>%filter(stopcdn=="TAG"))$relRatio)
#p-value = 0.3343
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.03778
wilcox.test((dfmrg_human%>%filter(stopcdn=="TAG"))$relRatio,
            (dfmrg_human%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.4089

#mouse
load("readcount_mouseAge.RData")
load("stopCdnSeq_mouse.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[26:28]])
compRep$old<- composite(readcount[ names(readcount)[19:21]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg_mouse<-merge(young,old)%>%merge(stopCdnSeq)%>%
  dplyr::mutate(relRatio=(old-young)/((old+young)/2),spec="mouse")%>%na.omit() %>%
  dplyr::filter(young <0.5 & old < 0.5)

summary(dfmrg_mouse$relRatio)

table(dfmrg_mouse$stopcdn)
#TAA TAG TGA 
#140 119 274 
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_mouse%>%filter(stopcdn=="TAG"))$young)
# p-value = 0.3305
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$young)
# p-value = 0.01076
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAG"))$young,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$young)
#p-value = 0.1569

wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_mouse%>%filter(stopcdn=="TAG"))$old)
#p-value = 0.01766
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$old)
#p-value = 0.003535
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAG"))$old,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$old)
#p-value = 0.9649

wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_mouse%>%filter(stopcdn=="TAG"))$relRatio)
#p-value = 0.3542
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.6163
wilcox.test((dfmrg_mouse%>%filter(stopcdn=="TAG"))$relRatio,
            (dfmrg_mouse%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.1188



#yeast
load("readcount_yeastChrAge.RData")
load("stopCdnSeq_yeast.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[1:4]])
compRep$old<- composite(readcount[ names(readcount)[9:12]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young%>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)
dfmrg_yeast<-merge(young,old)%>%separate(col = "transcript",into = c("transcript"))%>%
  merge(stopCdnSeq)%>%
  dplyr::mutate(relRatio=(old-young)/((old+young)/2),spec="yeast")%>%na.omit() %>%
  dplyr::filter(young < 0.5 & old <0.5)

summary(dfmrg_yeast$relRatio)
table(dfmrg_yeast$stopcdn)
#TAA TAG TGA 
#364 175 221 

wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_yeast%>%filter(stopcdn=="TAG"))$young)
#p-value = 0.001553
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$young,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$young)
# p-value = 9.641e-06
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAG"))$young,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$young)
#p-value = 0.4361

wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_yeast%>%filter(stopcdn=="TAG"))$old)
#p-value = 0.1909
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$old,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$old)
#p-value = 6.012e-05
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAG"))$old,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$old)
#p-value = 0.01654

wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_yeast%>%filter(stopcdn=="TAG"))$relRatio)
#p-value = 0.1459
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAA"))$relRatio,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.6416
wilcox.test((dfmrg_yeast%>%filter(stopcdn=="TAG"))$relRatio,
            (dfmrg_yeast%>%filter(stopcdn=="TGA"))$relRatio)
#p-value = 0.0878

##combined data

alldat<-rbind(dfmrg_human,dfmrg_mouse)%>%rbind(dfmrg_yeast) %>%
  dplyr::mutate(stopcdn=(case_when(stopcdn=="TAA" ~ "UAA",stopcdn=="TAG" ~ "UAG",
                                   stopcdn=="TGA" ~ "UGA")),
                context3=case_when(context3=="T" ~"U",context3!="T"~ context3),
                context5=case_when(context5=="T" ~"U",context5!="T"~ context5)
  )
alldat$spec<-factor(alldat$spec,levels = c("human","mouse","yeast"),
                    labels =c("Human","Mouse","Yeast") )

Fig3c<-ggplot(alldat,aes(x=factor(spec),y=young,fill=stopcdn))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF","#4DBBD5FF","#E64B35FF"))+
  labs(x="Species",y="Readthrough rates in young cells")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),limits = c(0,0.2))+
  scale_x_discrete()

Fig3c


Fig3d<-ggplot(alldat,aes(x=factor(spec),y=old,fill=stopcdn))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF","#4DBBD5FF","#E64B35FF"))+
  labs(x="Species",y="Readthrough rates in aged cells")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(breaks=c(0,0.1,0.2), labels = c(0,0.1,0.2),limits = c(0,0.2))+
  scale_x_discrete()

Fig3d

Fig3e<-ggplot(alldat,aes(x=factor(spec),y=relRatio,fill=stopcdn))+
  geom_boxplot(outlier.colour=NA,position=position_dodge(0.8),
               key_glyph="rect",width=0.5,lwd=0.2,fatten =0.8)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.2,lwd=0.2,outlier.shape = NA,
               position =position_dodge(0.8))+
  scale_fill_manual(name="",values = c("#00A087FF","#4DBBD5FF","#E64B35FF"))+
  labs(x="Species",y="Aged-to-young relative increase\nin readthrough rates")+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))+
  scale_y_continuous(limits = c(-2,2.7))+
  scale_x_discrete()+
  theme_test()+
  theme(legend.position = "top",legend.key   = element_blank(),legend.key.size = unit(3,"mm"))

Fig3e
# combinded figures ---------------------------------------------------------------------
library(cowplot)
Fig3ab <- plot_grid(Fig3a,Fig3b,ncol = 2,
                    labels = c("a","b"))  
Fig3cde <- plot_grid(Fig3c,Fig3d,Fig3e,ncol = 1,labels = c("c","d","e"))  
Fig3 <- plot_grid("",Fig3cde,ncol = 2,rel_widths  = c(2.5,1))
ggsave("Fig3.pdf",Fig3,width =8, height =8)
