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
library(clusterProfiler)
# Fig.4a ------------------------------------------------------------------
load("dfsubs_compile_desc.RData")
dfsubs_compile_desc<-dfsubs_compile
load("dfsubs_compile_bone.RData")
dfsubs_compile_bone<-dfsubs_compile

young_desc<-dfsubs_compile_desc%>%dplyr::select(substype,young_desc=young)%>%group_by(substype)%>%
  dplyr::summarise(young_desc=mean(young_desc))
young_bone<-dfsubs_compile_bone%>%dplyr::select(substype,young_bone=young)%>%group_by(substype)%>%
  dplyr::summarise(young_bone=mean(young_bone))

dfmrg<-merge(young_desc,young_bone,all = T)  
dfmrg[is.na(dfmrg)]<-0
dfmrg<-dfmrg%>%filter(young_desc>0 |young_bone>0 )
cor.test(dfmrg$young_desc,dfmrg$young_bone,method = "p")
#rho=-0.02670859  ,p-value = 0.7335,method="p"

Fig4a<-ggplot(dfmrg,aes(x=young_desc,y=young_bone))+geom_point()+
  xlab("Misincorporation rate\nin young intervertebral disc")+
  ylab("Misincorporation rate\nin young bone")+
  theme_test()+
  scale_y_continuous(breaks=c(0,0.0025,0.005), labels = c(0,0.0025,0.005),
                     limits =c(0,0.005) )+
  scale_x_continuous(breaks=c(0,0.0025,0.005), labels = c(0,0.0025,0.005),
                     limits =c(0,0.005))+
  annotate("text",x=0.004,y=0.0045,label=expression(italic(rho) * " = " * -0.03))+
  annotate("text",x=0.004,y=0.004,label=expression(italic(P) * " = " *0.73))

Fig4a

# Fig.4b ------------------------------------------------------------------

old_desc<-dfsubs_compile_desc%>%dplyr::select(substype,old_desc=old)%>%
  group_by(substype)%>%dplyr::summarise(old_desc=mean(old_desc))
old_bone<-dfsubs_compile_bone%>%dplyr::select(substype,old_bone=old)%>%
  group_by(substype)%>%dplyr::summarise(old_bone=mean(old_bone))

dfmrg<-merge(old_desc,old_bone,all=T)
dfmrg[is.na(dfmrg)]<-0
dfmrg<-dfmrg%>%filter(old_desc>0 |old_bone>0 )
cor.test(dfmrg$old_desc,dfmrg$old_bone,method = "p")
#rho=-0.1223429   p-value = 0.09265,method = "p", ALL

Fig4b<-ggplot(dfmrg,aes(x=old_desc,y=old_bone))+geom_point()+
  xlab("Misincorporation rate\nin aged intervertebral disc")+
  ylab("Misincorporation rate\nin aged bone")+
  theme_test()+
  scale_y_continuous(breaks=c(0,0.004,0.008), labels = c(0,0.004,0.008),
                     limits =c(0,0.008) )+
  scale_x_continuous(breaks=c(0,0.004,0.008), labels = c(0,0.004,0.008),
                     limits =c(0,0.008) )+
  annotate("text",x=0.006,y=0.007,label=expression(italic(rho) * " = " * -0.12))+
  annotate("text",x=0.006,y=0.006,label=expression(italic(P) * " = " *0.09))

Fig4b

# Fig.4c ------------------------------------------------------------------

rel_incr_desc<-dfsubs_compile_desc%>%dplyr::mutate(
  rel_incr_desc=(old-young) / ( (old+young)/2 ) )%>%
  group_by(substype) %>%
  dplyr::summarise(desc_mean=mean(rel_incr_desc,na.rm = T))

rel_incr_bone<-dfsubs_compile_bone%>%dplyr::mutate(
  rel_incr_bone=(old-young) / ( (old+young)/2 ) )%>%
  group_by(substype) %>%
  dplyr::summarise(bone_mean=mean(rel_incr_bone,na.rm = T))


dfmrg<-merge(rel_incr_desc,rel_incr_bone,all=T)
dfmrg[is.na(dfmrg)]<-0
cor.test(dfmrg$desc_mean,dfmrg$bone_mean,method="p")
#rho=0.005748063 ,p-value = 0.9324,method="s"

Fig4c<-ggplot(dfmrg,aes(x=desc_mean,y=bone_mean))+geom_point()+
  xlab("Aged−to−young relative\nincrease in intervertebral disc")+
  ylab("Aged−to−young relative\nincrease in bone")+
  theme_test()+
  annotate("text",x=-1.2,y=1.6,label=expression(italic(rho) * " = " * 0.006))+
  annotate("text",x=-1.2,y=1.2,label=expression(italic(P) * " = " *0.93))

Fig4c

# Fig.4d ------------------------------------------------------------------
load("rt_common.RData")
load("readcount_tissuecompar_xiezhi.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5=0; 

compRep <- list() # Composite replicates read count
names(readcount)
compRep$Brain_old<- composite(readcount[ names(readcount)[1:2]])
compRep$Heart_old<- composite(readcount[ names(readcount)[3:4]])
compRep$Kidney_old<- composite(readcount[ names(readcount)[5:6]])
compRep$Liver_old<- composite(readcount[ names(readcount)[7:8]])
compRep$Lung_old<- composite(readcount[ names(readcount)[9:10]])
compRep$Brain_young<- composite(readcount[ names(readcount)[11:12]])
compRep$Heart_young<- composite(readcount[ names(readcount)[13:14]])
compRep$Kidney_young<- composite(readcount[ names(readcount)[15:16]])
compRep$Liver_young<- composite(readcount[ names(readcount)[17:18]])
compRep$Lung_young<- composite(readcount[ names(readcount)[19:20]])

RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5 )

dfmrg<-lapply(names(RTscore), function(thistissue){
  eachtissue<-RTscore[[thistissue]]%>%dplyr::select(transcript,rte_ext)%>%
    filter(rte_ext<=0.1)
  eachtissue$tissue<-thistissue
  return(eachtissue)
})%>%rbind.fill() %>% 
  pivot_wider(id_cols = "transcript",names_from = "tissue",values_from = "rte_ext")
dfmrg[is.na(dfmrg)]<-0

dfmrg_flt<-dfmrg%>%rowwise()%>%
  dplyr::mutate(sum_young=sum(c(Brain_young,Heart_young,Kidney_young,
                                Liver_young,Lung_young),na.rm=T))%>%
  dplyr::mutate(sum_old=sum(c(Brain_old,Heart_old,Kidney_old,
                              Liver_old,Lung_old),na.rm=T))%>%
  dplyr::mutate(Brain_rel=(Brain_old - Brain_young)/( (Brain_old + Brain_young)/2 ),
                Heart_rel=(Heart_old - Heart_young)/( (Heart_old + Heart_young)/2 ),
                Kidney_rel=(Kidney_old - Kidney_young)/( (Kidney_old + Kidney_young)/2 ),
                Liver_rel=(Liver_old - Liver_young)/( (Liver_old + Liver_young)/2 ),
                Lung_rel=(Lung_old - Lung_young)/( (Lung_old + Lung_young)/2 ) ) %>%
  dplyr::mutate(dfvar=var(c(Brain_rel,Heart_rel,Kidney_rel,Liver_rel,Lung_rel) ,na.rm = T))%>%
  filter(sum_young>0 & sum_old>0,dfvar!=0)%>%dplyr::select(-c(sum_young,sum_old,dfvar)) %>%
  column_to_rownames(var="transcript")
dfmrg_flt[is.na(dfmrg_flt)]<-0



###rel
{
  dfdat_rel<-dfmrg_flt%>% dplyr::select(contains("rel"))
  colnames(dfdat_rel)<-c("Brain","Heart","Kidney","Liver","Lung")
  
  #cutoff=1
  cutoff=quantile(as.matrix(dfdat_rel),probs = 0.75)
  
  
  commu_up=((dfdat_rel > cutoff) %>%rowSums(na.rm=T))==5
  dfcommu_up=dfdat_rel[commu_up,]
  
  
  spec=((dfdat_rel > cutoff) %>%rowSums(na.rm=T))==1
  dfspec=dfdat_rel[spec,]
  
  brain=dfspec[(dfspec$Brain > cutoff),] 
  heart=dfspec[(dfspec$Heart > cutoff),] 
  kidney=dfspec[(dfspec$Kidney > cutoff),] 
  liver=dfspec[(dfspec$Liver > cutoff),] 
  lung=dfspec[(dfspec$Lung > cutoff),] 
  
  
  dfclustered<-rbind(dfcommu_up)%>%rbind(brain)%>%
    rbind(heart)%>%rbind(kidney)%>%rbind(liver)%>%rbind(lung)
  
  
  mat<-dfclustered%>%as.matrix()
  col_fun = colorRamp2(c(-1, 0, 1), c("#293890", "white", "#BF1D2D"))
  
  
  rel<-Heatmap(mat, col = col_fun,
               #column_title = "Relative increase",
               cluster_rows = F,cluster_columns = F,  show_row_names = F,
               heatmap_legend_param = list( 
                 title = "Aged-to-young relative\nincrease in readthrough rate",
                 title_position ="leftcenter-rot", 
                 direction="vertical",
                 labels_gp = gpar(fontsize = 10),
                 title_gp = gpar(fontsize = 10))
               )
  
  Fig4d_rel<-grid.grabExpr( draw(rel, heatmap_legend_side = "right"))
  
  
}
##young
{
dfdat_young<-dfmrg_flt%>% dplyr::select(contains("young"))
dfdat_young<-dfdat_young[rownames(dfclustered),]
colnames(dfdat_young)<-c("Brain","Heart","Kidney","Liver","Lung")

mat<-dfdat_young%>%as.matrix()
col_fun_yong = colorRamp2(c(0, 0.05), c("white", "#4DBBD5FF"))

young<-Heatmap(mat, col=col_fun_yong,
               #column_title = "Young",
               cluster_rows = F,cluster_columns = F,show_row_names = F,
               heatmap_legend_param = list( 
                 title = "Readthrough rate\nin young samples",
                 title_position ="leftcenter-rot", 
                 direction="vertical",
                 labels_gp = gpar(fontsize = 10),
                 title_gp = gpar(fontsize = 10)
               ))
Fig4d_young<-grid.grabExpr( draw(young, heatmap_legend_side = "right"))
}
##old
{
dfdat_old<-dfmrg_flt%>% dplyr::select(contains("old"))
dfdat_old<-dfdat_old[rownames(dfclustered),]
colnames(dfdat_old)<-c("Brain","Heart","Kidney","Liver","Lung")

mat<-dfdat_old%>%as.matrix()

col_fun_old = colorRamp2(c(0, 0.05), c("white", "#00A087FF"))


old<-Heatmap(mat, col=col_fun_old,
             cluster_rows = F,cluster_columns = F, show_row_names = F,
             heatmap_legend_param = list( 
               title = "Readthrough rate\nin aged samples",
               title_position ="leftcenter-rot", 
               direction="vertical",
               labels_gp = gpar(fontsize = 10),
               title_gp = gpar(fontsize = 10)
             ))
Fig4d_old<-grid.grabExpr( draw(old, heatmap_legend_side = "right") )
}


# fig4e -------------------------------------------------------------------

dfdat_go=list(dfcommu_up,brain,heart,kidney,liver,lung)
dfname=c("Common","Brain","Heart","Kidney","Liver","Lung")
go_all<-lapply(1:length(dfdat_go), function(i){
  gene=rownames(dfdat_go[[i]]) 
  gene=bitr(gene,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
  gene=dplyr::distinct(gene,ENSEMBL,.keep_all=T)
  go=enrichGO(gene = gene$ENTREZID,OrgDb ="org.Mm.eg.db",keyType = "ENTREZID",ont = "ALL")
  go_result<-go@result%>%filter(p.adjust<0.05 &Count >=2 )
  if(nrow(go_result)==0){return(NULL)}
  go_result$type=dfname[i]
  return(go_result)
  
})%>%rbind.fill() 

go_all_top5<-go_all%>%dplyr::select(c(Description,p.adjust,type))%>%
  filter(!grepl("gliogenesis|tRNA|mitochondrial protein|RNA polymerase|FNIP-folliculin RagC/D GAP",
                Description,ignore.case=T))%>%
  dplyr::mutate(pvalue=-log10(p.adjust))%>%
  dplyr::group_by(type)%>%arrange(type, pvalue)%>%
  top_n(n=5,wt=pvalue)%>%
  dplyr::group_by(type)%>%
  sample_n(5)%>%arrange(type, pvalue)%>%
  rownames_to_column()


go_all_top5$Description <- str_c(str_to_title(str_sub(go_all_top5$Description, 1, 1)), 
                                 str_sub(go_all_top5$Description, 2))
go_all_top5$Description<- factor(as.integer(go_all_top5$rowname),labels = go_all_top5$Description)
go_all_top5$type<-factor(go_all_top5$type,levels = c("Common","Brain","Heart","Kidney","Liver","Lung"))
Fig4e<-ggplot(go_all_top5,aes(x=Description,y=pvalue))+
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
        strip.text=element_text(size=10,colour="white")) ##分面字体更改为白色

Fig4e

# Fig.4f ------------------------------------------------------------------
load("rt_common.RData")
load("readcount_tissuecompar_xiezhi.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5=0; 
RTscore <- lapply(readcount, rt_efficiency, cds_m5, cds_m3,exten_m5 )
names(RTscore)
dfmrg<-lapply(names(RTscore)[1:20], function(thistissue){
  eachtissue<-RTscore[[thistissue]]%>%dplyr::select(transcript,rte_ext)
  eachtissue$tissue<-thistissue
  return(eachtissue)
})%>%rbind.fill() %>% 
  pivot_wider(id_cols = "transcript",names_from = "tissue",values_from = "rte_ext")

dfmrg[is.na(dfmrg)]<-0

dfyoung<-dfmrg%>%select(contains("young"))%>%t()
dfyoung_dist=dist(dfyoung,method = "euclidean")%>%as.numeric()

dfold<-dfmrg%>%select(contains("old"))%>%t()
dfold_dist=dist(dfold,method = "euclidean")%>%as.numeric()

wilcox.test(dfyoung_dist,dfold_dist)
#p-value = 2.502e-05
mean(dfold_dist)/mean(dfyoung_dist) #1.28


dfdis_sameTissue<-lapply(c("liver","kidney","brain","lung","heart"),
                         function(thisTissue){
                           dfdat<-dfmrg%>%dplyr::select(contains(thisTissue)) %>%t()
                           dfdist=dist(dfdat,method = "euclidean") %>%as.matrix()%>%
                             as.data.frame()%>%
                             select(contains("old"))%>%rownames_to_column()%>%
                             filter(grepl("young",rowname)) %>%
                             column_to_rownames("rowname")%>%
                             as.matrix()%>%as.numeric()
                           
                           return(data.frame(tissue=thisTissue, dfdist)) 
                         })%>%rbind.fill()%>%getElement("dfdist")

t.test(dfyoung_dist,dfdis_sameTissue)
#p-value = 0.02315
t.test(dfold_dist,dfdis_sameTissue)
#p-value = 0.2321
dfdatplot<-data.frame(dfvalue=c(dfyoung_dist,dfold_dist,dfdis_sameTissue),
                      dftype=c(rep("Young",length(dfyoung_dist)),
                               rep("Old",length(dfyoung_dist)),
                               rep("Young vs Old",length(dfdis_sameTissue)) ))
dfdatplot$dftype<-factor(dfdatplot$dftype,levels = c("Young","Old","Young vs Old"))
dfdatplot$group<-factor(c(rep("Between different tissues at the same age",length(c(dfyoung_dist,dfyoung_dist))),
                          rep("Between different ages for the same tissue",length(c(dfdis_sameTissue)))),
                        levels = c("Between different tissues at the same age",
                                   "Between different ages for the same tissue"))
# Fig4f<-ggplot(dfdatplot,aes(x=dftype,y=dfvalue,fill = dftype))+
#   stat_boxplot(geom = "errorbar",width=0.15)+
#   geom_boxplot(fill="white",outlier.fill = "white",outlier.color = "white")+
#   geom_jitter(aes(fill=dftype),width =0.2,shape = 21,size=2.5)+
#   scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442"))+  
#   ylab("Readthrough rate difference of a tissue pair\n(Euclidean distance of shared readthrough genes)")+
#   xlab(NULL)+
#   scale_x_discrete(breaks=c("Young","Old","Young vs Old"),
#                    labels=c("Between different tissues\nin young samples",
#                             "Between different tissues\nin aged samples",
#                             "Between young and aged\nsamples for the same tissues"))+
#   theme_test()+
#   coord_flip()+ 
#   annotate("segment",x=1,xend = 1.95,y=7.8,yend = 7.8)+
#   annotate("text",x=1.5,y=8.1,label=expression(italic(P) * " = " * "3×" * 10**-5),angle = -90)+
#   annotate("segment",x=1,xend = 3,y=8.5,yend = 8.5)+
#   annotate("text",x=2,y=8.8,label=expression(italic(P) * " = " * "0.02"),angle = -90)+
#   annotate("segment",x=2.05,xend = 3,y=7.8,yend = 7.8)+
#   annotate("text",x=2.5,y=8.1,label=expression(italic(P) * " = " * "0.23"),angle = -90)+
#   theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
#         axis.text.y = element_text(size = 10,angle = 30,
#                                    hjust = 1,vjust = 0))
# Fig4f

Fig4f<-ggplot(dfdatplot,aes(x=dftype,y=dfvalue,fill = dftype))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_boxplot(fill="white",outlier.fill = "white",outlier.color = "white")+
  geom_jitter(aes(fill=dftype),width =0.2,shape = 21,size=2.5)+
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF","#00A087FF"))+  
  ylab("Readthrough rate difference of\na tissue pair(Euclidean distance of\nshared readthrough genes)")+
  xlab(NULL)+
  scale_x_discrete(breaks=c("Young","Old","Young vs Old"),
                   labels=c("Between different tissues\nin young samples",
                            "Between different tissues\nin aged samples",
                            "Between young and aged\nsamples for the same tissues"))+
  theme_test()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
        axis.text.x = element_text(size = 10,angle = 40,
                                   hjust = 0.9,vjust = 1))+
  annotate("segment",x=1,xend = 1.95,y=7.8,yend = 7.8)+
  annotate("text",x=1.5,y=8.1,label=expression(italic(P) * " = " * "3×" * 10**-5))+
  annotate("segment",x=1,xend = 3,y=8.5,yend = 8.5)+
  annotate("text",x=2,y=8.8,label=expression(italic(P) * " = " * "0.02"))+
  annotate("segment",x=2.05,xend = 3,y=7.8,yend = 7.8)+
  annotate("text",x=2.5,y=8.1,label=expression(italic(P) * " = " * "0.23"))
 
Fig4f

# Fig4g 线性模型 --------------------------------------------------------------

#save(dfmrg_flt,file="mouse_tissue_rt.RData")
library(tidyverse)
load("stopCdnSeq_mouse.RData")
names(stopCdnSeq)[1]<-"gene"
load("mouse_tissue_rt.RData")
dfdat_rel<-dfmrg_flt%>% dplyr::select(contains("rel"))
colnames(dfdat_rel)<-c("Brain","Heart","Kidney","Liver","Lung")
dfdat<-dfdat_rel%>%rownames_to_column(var="gene")

dfmrg<-merge(dfdat,stopCdnSeq)%>%dplyr::select(1:7)%>%
  pivot_longer(cols = c("Brain","Heart","Kidney","Liver","Lung"),
               names_to = "tissue",values_to ="rel_rate"  )

lmRes <- lm(rel_rate ~ stopcdn * tissue, data=dfmrg)
summary(lmRes)

dfanv<-anova(lmRes)
dfanv$`Pr(>F)`

dfmrg_tissue<-dfmrg%>%group_by(tissue)%>%
  dplyr::summarise(dfmean=mean(rel_rate), 
                   se=sd(rel_rate)/sqrt(length(rel_rate)))
colnames(dfmrg_tissue)[1]<-"type"
dfmrg_stopcdn<-dfmrg%>%group_by(stopcdn)%>%
  dplyr::summarise(dfmean=mean(rel_rate), 
                   se=sd(rel_rate)/sqrt(length(rel_rate)))
colnames(dfmrg_stopcdn)[1]<-"type"

dfpool<-rbind(dfmrg_tissue,dfmrg_stopcdn)
dfpool$type<-factor(dfpool$type,levels = c("Brain","Lung","Liver",
                                           "Heart","Kidney","TGA","TAG","TAA"))
dfpool$dffactor<-c(rep("tisuue",5),rep("ycodon",3))


Fig4g<-ggplot(dfpool,aes(x=type,y=dfmean))+
  geom_errorbar(aes(ymax=dfmean+se,ymin=dfmean-se,width=0.2))+
  geom_point(size=0.8)+
  geom_vline(xintercept = 5.5, color = "red", linetype = "dashed", size = 1)+
  annotate("segment",x=1,xend = 5,y=0.5,yend = 0.5)+
  annotate("segment",x=6,xend = 8,y=0.5,yend = 0.5)+
  annotate("text",x=3,y=0.6,label=expression("ANOVA " *italic(P) * " = " * "3.2×" * 10**-63),
           size=3)+
  annotate("text",x=7,y=0.6,label=expression("0.03"),size=3)+
  labs(x=NULL,y="Average aged−to−young relative\nincrease for readthrough rate" )+
  theme_test()+
  theme(axis.text.x = element_text(size = 10,angle = 40,
                                   hjust = 0.9,vjust = 1))

Fig4g

# 合并图 ---------------------------------------------------------------------
library(cowplot)
layout <- matrix(c(1, 2, 3, 
                   4,5,6), nrow = 2, byrow = TRUE)

Fig4a_d=grid.arrange(ggplotGrob(Fig4a),ggplotGrob(Fig4b),ggplotGrob(Fig4c),
             Fig4d_rel,Fig4d_young,Fig4d_old,
             layout_matrix = layout)

Fig4fg <- plot_grid(Fig4f,Fig4g,nrow = 2,rel_heights = c(2,1.2)) 
Fig4efg<- plot_grid(Fig4e,Fig4fg,nrow = 1,rel_widths  = c(1,1)) 
Fig4 <- plot_grid(Fig4a_d,Fig4efg,ncol = 1 ,rel_heights = c(2,2))
ggsave("Fig4.pdf",Fig4,width =20, height =29,units = "cm")

