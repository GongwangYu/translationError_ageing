library(data.table)
library(dplyr)
library(readr)
library(tidyverse)
library(riboWaltz)
library(plyr)
library(factoextra)
library(circlize)
library(ComplexHeatmap)
library(ggplotify)
library(cowplot)

load("rt_common.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5=0;frame_ = "0"

# 1 Fig5a  ----------------------------------------------------------------
# 1.1 遗传操作改变寿命 ----------------------------------------------------

load("readcount_lifespan.RData")
compRep <- list() 
names(readcount)
compRep$S1_long_lifespan<- composite(readcount[ names(readcount)[1:3] ])
compRep$S2_long_lifespan<- composite(readcount[ names(readcount)[7:9] ])
compRep$S1_short_lifespan<- composite(readcount[names(readcount)[4:6]  ]) 
compRep$S2_short_lifespan<- composite(readcount[ names(readcount)[10:11] ])
compRep$S3_short_lifespan<- composite(readcount[ names(readcount)[12:14] ])
compRep$wt<- composite(readcount[ names(readcount)[15:17] ])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3, exten_m5 )
S1_long<-RTscore$S1_long_lifespan %>%dplyr::select(transcript,S1_long=rte_ext)
S2_long<-RTscore$S2_long_lifespan %>%dplyr::select(transcript,S2_long=rte_ext)
#rpl6a/rpl15b/rpl27是否导致寿命减短还不确定，因此去掉。
# S1_short<-RTscore$S1_short_lifespan %>%dplyr::select(transcript,S1_short=rte_ext)
# S2_short<-RTscore$S2_short_lifespan %>%dplyr::select(transcript,S2_short=rte_ext)
# S3_short<-RTscore$S3_short_lifespan %>%dplyr::select(transcript,S3_short=rte_ext)
wt<-RTscore$wt %>%dplyr::select(transcript,wt=rte_ext)

###条件摸索
{
  # dfmrg<-merge(S1_long,S2_long,all = T)%>%merge(S1_short,all = T)%>%
  #   merge(S2_short,all = T)%>%merge(S3_short,all = T)%>%merge(wt,all = T)
  # dfmrg[is.na(dfmrg)]<-0
  
  # dfmrg<-merge(S1_long,S2_long)%>%merge(S1_short,)%>%
  #   merge(S2_short)%>%merge(S3_short)%>%merge(wt)
  
  dfmrg<-merge(S1_long,S2_long)%>%merge(wt)
}

dfmrg<-dfmrg%>%dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

dfdat<-dfmrg%>%column_to_rownames(var = "gene_name")
dfdat[dfdat==0]<-0.000001
dfdat<-(dfdat-dfdat$wt)/((dfdat+dfdat$wt)/2)

wilcox.test(dfdat$S1_long,mu = 0,alternative = "l")$p.value #1.27894e-07
wilcox.test(dfdat$S2_long,mu = 0,alternative = "l")$p.value #86.015491e-08
# wilcox.test(dfdat$S1_short,mu = 0,alternative = "g")$p.value 
# wilcox.test(dfdat$S2_short,mu = 0,alternative = "g")$p.value 
# wilcox.test(dfdat$S3_short,mu = 0,alternative = "g")$p.value 


dfplot1<-data.frame(treated=c("S1_long","S2_long"),
                    mean_error=c(mean(dfdat$S1_long),mean(dfdat$S2_long)),
                    sd_error=c(sd(dfdat$S1_long),sd(dfdat$S2_long)),
                    se_error=c(sd(dfdat$S1_long)/sqrt(nrow(dfdat)),
                               sd(dfdat$S2_long)/sqrt(nrow(dfdat)) ),
                    pvalue=c(1.27894e-07,6.015491e-08))

plot2dat<-dfdat%>%rownames_to_column(var = "gene")%>%dplyr::select(-wt)
# GR限制 ----------------------------------------------------------------------

load("readcount_yeastGlucoseRestriction.RData")
compRep <- list() 
names(readcount)
#将protocol A和B合并
compRep$SD<- composite(readcount[ names(readcount)[4:6]])
compRep$GR<- composite(readcount[ names(readcount)[1:3]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3, exten_m5 )
SD<-RTscore$SD %>% dplyr::select(transcript,SD=rte_ext)
GR<-RTscore$GR %>% dplyr::select(transcript,GR=rte_ext)


dfmrg<-merge(GR,SD)%>%dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

dfdat<-dfmrg%>%column_to_rownames(var = "gene_name")
dfdat[dfdat==0]<-0.000001
dfdat<-(dfdat-dfdat$SD)/((dfdat+dfdat$SD)/2)

wilcox.test(dfdat$GR,mu = 0,alternative = "l")$p.value #8.225053e-19

dfplot1<-rbind(dfplot1,
               data.frame(treated=c("GR"),
                          mean_error=mean(dfdat$GR),
                          sd_error=sd(dfdat$GR),
                          se_error=sd(dfdat$GR)/sqrt(nrow(dfdat)),
                          pvalue=8.225053e-19
               ) )

plot2dat2<-dfdat%>%rownames_to_column(var = "gene")%>%dplyr::select(-SD)%>%merge(plot2dat,all = T)
# Met 限制饮食 ----------------------------------------------------------------
load("readcount_yeastMetRestriction.RData")
compRep <- list() 
names(readcount)
compRep$MetR<- composite(readcount[ names(readcount)[1]])
compRep$SD<- composite(readcount[ names(readcount)[2]])

RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3 ,exten_m5)
SD<-RTscore$SD %>%dplyr::select(transcript,SD=rte_ext)
MetR<-RTscore$MetR %>%dplyr::select(transcript,MetR=rte_ext)

dfmrg<-merge(SD,MetR)%>%dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

dfdat<-dfmrg%>%column_to_rownames(var = "gene_name")
dfdat[dfdat==0]<-0.000001
dfdat<-(dfdat-dfdat$SD)/((dfdat+dfdat$SD)/2)


wilcox.test(dfdat$MetR,mu = 0,alternative = "l")$p.value #0.03900637

dfplot1<-rbind(dfplot1,
               data.frame(treated=c("MetR"),
                          mean_error=mean(dfdat$MetR),
                          sd_error=sd(dfdat$MetR),
                          se_error=sd(dfdat$MetR)/sqrt(nrow(dfdat)),
                          pvalue=0.03900637
               ) )

plot2dat3<-dfdat%>%rownames_to_column(var = "gene")%>%dplyr::select(-SD)%>%merge(plot2dat2,all = T)

# 雷帕霉素处理 ------------------------------------------------------------------
load("readcount_yeastRap1h3h.RData")
compRep <- list() 
names(readcount)
compRep$ypd<- composite(readcount[ names(readcount)[1:2]])
compRep$rapa<- composite(readcount[ names(readcount)[5]])

RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3 ,exten_m5)

ypd<-RTscore$ypd %>%dplyr::select(transcript,ypd=rte_ext)
rapa<-RTscore$rapa %>%dplyr::select(transcript,rapa=rte_ext)

dfmrg<-merge(ypd,rapa)%>%dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

dfdat<-dfmrg%>%column_to_rownames(var = "gene_name")
dfdat[dfdat==0]<-0.000001
dfdat<-(dfdat-dfdat$ypd)/((dfdat+dfdat$ypd)/2)

wilcox.test(dfdat$rapa,mu = 0,alternative = "l")$p.value #1.171764e-13

dfplot1<-rbind(dfplot1,
               data.frame(treated=c("rapa"),
                          mean_error=mean(dfdat$rapa),
                          sd_error=sd(dfdat$rapa),
                          se_error=sd(dfdat$rapa)/sqrt(nrow(dfdat)),
                          pvalue=1.171764e-13) )
plot2dat4<-dfdat%>%rownames_to_column(var = "gene")%>%dplyr::select(-ypd)%>%merge(plot2dat3,all = T)
# 巴龙霉素处理 ------------------------------------------------------------------
load("readcount_yeastParo.RData")
compRep <- list() 
names(readcount)
compRep$ypd<- composite(readcount[ names(readcount)[7:9] ])
compRep$paro<- composite(readcount[ names(readcount)[1:3] ])
compRep$rapa<- composite(readcount[ names(readcount)[4:6] ])

RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3 ,exten_m5)
ypd<-RTscore$ypd %>%dplyr::select(transcript,ypd=rte_ext)
paro<-RTscore$paro %>%dplyr::select(transcript,paro=rte_ext)
rapa<-RTscore$rapa %>%dplyr::select(transcript,rapa=rte_ext)

dfmrg<-merge(ypd,rapa)%>%merge(paro)%>%
  dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

dfdat<-dfmrg%>%column_to_rownames(var = "gene_name")
dfdat[dfdat==0]<-0.000001
dfdat<-(dfdat-dfdat$ypd)/((dfdat+dfdat$ypd)/2)


wilcox.test(dfdat$rapa,mu = 0,alternative = "l")$p.value #0.007410955
wilcox.test(dfdat$paro,mu = 0,alternative = "g")$p.value #7.279755e-23

dfplot1<-rbind(dfplot1,
               data.frame(treated=c("paro"),
                          mean_error=mean(dfdat$paro),
                          sd_error=sd(dfdat$paro),
                          se_error=sd(dfdat$paro)/sqrt(nrow(dfdat)),
                          pvalue=7.279755e-23
               ) )

plot2dat5<-dfdat%>%rownames_to_column(var = "gene")%>%dplyr::select(-c(ypd,rapa))%>%merge(plot2dat4,all = T)
# 绘图 ----------------------------------------------------------------------

Fig5a<-ggplot(dfplot1%>%filter(treated!="MetR"),
              aes(x=treated,y=mean_error,fill=treated))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin =mean_error-se_error,ymax = mean_error+se_error ),
                width=0.1)+
  #geom_hline(yintercept = 0,linetype="dashed")+
  scale_x_discrete(limits=c("paro","GR","rapa",
                            "S1_long","S2_long"),
                   labels=c("Paromomycin","Glucose restriction","Rapamycin",
                            "Long lived S1","Long lived S2"))+
  scale_fill_manual(values = c( "#BF1D2D","#293890","#BF1D2D","#BF1D2D","#BF1D2D"))+
  annotate("text",x=1.2,y=1.5,label=expression(italic(P) * "<" *10^{-22}))+
  annotate("text",x=2,y=-0.7,label=expression(10^{-18} ))+
  annotate("text",x=3,y=-0.7,label=expression( 10^{-12}))+
  annotate("text",x=4,y=-0.8,label=expression(10^{-6}))+
  annotate("text",x=5,y=0-0.9,label=expression(10^{-7}))+
  theme_test()+
  #labs(title = "Yeast")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Treatment for lifespan changes")+
  ylab("Treatment-to-WT relative\nincrease in readthrough rate") +
  theme(axis.text.x = element_text(size=10,hjust = 1,
                                   vjust = 1,angle = 45),
        legend.position = "none")


Fig5a


# fig5b -------------------------------------------------------------------

plot5bDat<-plot2dat5 %>%rowwise()%>%
  dplyr::mutate(long=mean(c(S1_long,S2_long),na.rm=T))%>%
  column_to_rownames(var = "gene")%>%
  dplyr::select(-c(paro,MetR,S1_long,S2_long))
plot5bDat[is.na(plot5bDat)]<-0

cutoff=quantile(as.matrix(plot5bDat),probs = 0.25)

commu=((plot5bDat < cutoff) %>%rowSums(na.rm=T))==3
dfcommu=plot5bDat[commu,]

spec=((plot5bDat < cutoff) %>%rowSums(na.rm=T))==1
dfspec=plot5bDat[spec,]

rapa=dfspec[(dfspec$rapa< cutoff),] 
gr=dfspec[(dfspec$GR< cutoff),] 
long=dfspec[(dfspec$long< cutoff),]


dfclustered<-rbind(dfcommu,gr)%>%rbind(rapa)%>%rbind(long)
colnames(dfclustered)<-c("Rapamycin","Glucose restriction","Long lived")
mat<-dfclustered%>%as.matrix()
col_fun = colorRamp2(c(-1, 0, 1), c( "#293890","white", "#BF1D2D"))

group <- factor(c(rep("Comm",nrow(dfcommu)),rep("GR",nrow(gr)),rep("Rapa",nrow(rapa)),
                  rep("S",nrow(long))))

Fig5b<-Heatmap(mat, col = col_fun,cluster_rows = F,cluster_columns = F,
               show_row_names = F,
               row_title_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               #row_split = group, # 根据基因功能分组信息分割热图
               heatmap_legend_param = list( 
                 title = "Treatment-to-WT relative\nincrease in readthrough rate",
                 title_position ="leftcenter-rot",  
                 direction="vertical",
                 labels_gp = gpar(fontsize = 6),
                 title_gp = gpar(fontsize = 8),
                 grid_width = unit(2, "mm"),grid_height = unit(6, "mm")
                 
               )
)


Fig5b<-as.ggplot(Fig5b)


# plot5c ------------------------------------------------------------------
dfdat_go=list(dfcommu,rapa,gr,long)
dfname=c("Common","Rapamycin","Glucose restriction","Long lived")
go_all<-lapply(1:length(dfdat_go), function(i){
  gene=rownames(dfdat_go[[i]]) 
  gene=bitr(gene,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Sc.sgd.db")
  gene=dplyr::distinct(gene,ENSEMBL,.keep_all=T)
  go=enrichGO(gene = gene$ENTREZID,OrgDb ="org.Sc.sgd.db",keyType = "ENTREZID",ont = "BP")
  go_result<-go@result%>%filter(pvalue<0.05 &Count >=1 )
  if(nrow(go_result)==0){return(NULL)}
  go_result$type=dfname[i]
  return(go_result)
  
})%>%rbind.fill()


go_all_top5<-go_all%>%dplyr::select(c(Description,pvalue,type))%>%
  filter(!grepl("vacuolar acidification",
                Description,ignore.case=T))%>%
  dplyr::mutate(pvalue=-log10(pvalue))%>%
  dplyr::group_by(type)%>%arrange(type, pvalue)%>%
  top_n(n=5,wt=pvalue)%>%
  dplyr::group_by(type)%>%
  sample_n(5)%>%arrange(type, pvalue)%>%
  rownames_to_column()


go_all_top5$Description <- str_c(str_to_title(str_sub(go_all_top5$Description, 1, 1)), 
                                 str_sub(go_all_top5$Description, 2))
go_all_top5$Description<- factor(as.integer(go_all_top5$rowname),labels = go_all_top5$Description)
go_all_top5$type<-factor(go_all_top5$type,
                         levels = c("Common","Glucose restriction","Long lived","Rapamycin"))




Fig5c<-ggplot(go_all_top5,aes(x=Description,y=pvalue))+
  geom_bar(stat="identity",width=0.8,fill="#a6b3d7",alpha=0.9)+#绘制条形图
  geom_text(aes(y=0,#控制文本标签起始位置
                label=Description),
            size=3,hjust=0)+#hjust=0左对齐
  scale_y_continuous(expand=c(0,0.05))+
  facet_grid(type~. ,scales = "free_y",space = "free_y")+
  labs(x=NULL,y="-log10 (p.adjust)")+
  theme_bw()+coord_flip()+
  theme(
    axis.text.y=element_blank(),##删去y轴label
    axis.ticks.y=element_blank(),##删去y轴刻度线
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank()
  )+
  theme(strip.background=element_rect(fill=c("#2b3990")),##分面背景更改为红色
        strip.text=element_text(size=10,colour="white"))

Fig5c
# 合并------------------------------------------------------------------

plot5ab<-plot_grid(Fig5a,Fig5b,ncol = 2,nrow = 1,rel_widths = c(1.5,1))
Fig5<-plot_grid(plot5ab,Fig5c,ncol = 1,nrow = 2,rel_heights = c(1,1.5))
ggsave("Fig5.pdf",Fig5,width =6, height =10)
