# human-rt ----------------------------------------------------------------

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

dfmrg<-merge(young,old)%>%dplyr::rename(ensg=transcript)

ensg_name<-read_tsv("ensg_name_human.txt")
names(ensg_name)<-c("ensg","gene_name")

expr<-read_tsv("human_rna_seq_genes_tpm.txt") %>%
  dplyr::select("Gene_Id","Young1","Young2","Senescent1","Senescent2")%>%
  separate(col = "Gene_Id",into = c("a","gene_name"),sep = "_")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(expr_young=mean(c(Young1,Young2)),
                expr_old=mean(c(Senescent1,Senescent2)))%>%
  dplyr::select(gene_name,expr_young,expr_old)%>%merge(ensg_name)


dfmrg<-dplyr::left_join(dfmrg,expr,by="ensg") %>%
  filter(young<0.5 & old <0.5) %>%
  filter(!((expr_old/expr_young)>1)) %>%
  dplyr::select(young,old)

wilcox.test(dfmrg$young,dfmrg$old,paired = T,alternative = "l")$p.value
#n=426,4.723705e-30

range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig.S4a<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
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
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "4.7×" * 10**-30),size=2)+
  annotate("text",x=2.1,y=0.2,label=expression("5.5×" * 10**-30),size=2)+
  annotate("text",x=3,y=0.2,label=expression("1.4×" * 10**-29),size=2)+
  annotate("text",x=4,y=0.2,label=expression("2.4×" * 10**-28),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "422)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "416)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "392)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "336)" ),size=2)

Fig.S4a

# mouse-rt ----------------------------------------------------------------

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

dfmrg<-merge(young,old)%>%dplyr::rename(ensg=transcript)

wilcox.test(dfmrg$young,dfmrg$old,paired = T, alternative = "l")$p.value
#p-value = 3.199819e-13

# 去除表达量的影响 ----------------------------------------------------------------
if(FALSE){
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF("gencode.vM35.annotation.gtf",format="gtf")
  exons.list.per.gene <- exonsBy(txdb, by = "gene")
  
  #通过reduce函数避免重复计算重叠区
  exonic.gene.sizes <- lapply(exons.list.per.gene,
                              function(x){sum(width(reduce(x)))})
  #生成的geneID为ensemble编号
  eff_length <- do.call(rbind,lapply(exonic.gene.sizes, data.frame))
  eff_length <- data.frame(gene_id = rownames(eff_length),effLen = eff_length[,1])
  rownames(eff_length)<-eff_length$gene_id 
  rownames(eff_length) <- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]
  eff_length$gene_id<- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]
  head(eff_length)
  
  count<-read_tsv("mouse_rna_seq_counts.txt")%>%
    dplyr::select("gene_id","liver_3m_1_rna","liver_3m_2_rna","liver_3m_3_rna",
                  "liver_32m_1_rna","liver_32m_2_rna","liver_32m_3_rna")
  
  count$gene_id<-do.call(rbind,strsplit(as.character(count$gene_id),'\\.'))[,1]
  rownames(count)<-count$gene_id
  
  dfmrg<-merge(count,eff_length)
  
  count<-dfmrg%>%dplyr::select(2:7)
  effLen<-dfmrg$effLen
  # counts:转录组的count矩阵，行为基因，列为样本
  # effLen：一个数值型向量，值是基因长度，顺序应该与count的列一致对应。
  Counts2TPM <- function(count, effLen){
    rate <- log(count) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }
  
  trans_tpm <- apply(count, 2, Counts2TPM, effLen = effLen) %>%as.data.frame()
  trans_tpm$gene_id<-dfmrg$gene_id
  
  expr<-trans_tpm%>%dplyr::group_by(ensg=gene_id)%>%
    dplyr::summarise(expr_young=log2(mean(c(liver_3m_1_rna,liver_3m_2_rna,liver_3m_3_rna),na.rm=T)+1),
                     expr_old=log2(mean(c(liver_32m_1_rna,liver_32m_2_rna,liver_32m_3_rna),na.rm=T)+1))
  save(expr,file = "mouse_liver_expr.RData")
}
load("mouse_liver_expr.RData")
dfmrg<-dplyr::left_join(dfmrg,expr,by="ensg") %>%
  filter(young<0.5 & old <0.5) %>%
  filter(!((expr_old/expr_young)>1)) %>%
  dplyr::select(young,old)

wilcox.test(dfmrg$young,dfmrg$old,paired = T,alternative = "l")$p.value
#n=218,0.0001617667

range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig.S4b<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
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
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "1.6×" * 10**-4),size=2)+
  annotate("text",x=2,y=0.2,label=expression("7.2×" * 10**-4),size=2)+
  annotate("text",x=3,y=0.2,label=expression("6.8×" * 10**-3),size=2)+
  annotate("text",x=4,y=0.2,label=expression("9.6×" * 10**-2),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "203)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "169)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "120)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "54)" ),size=2)

Fig.S4b

# yeast-rt ----------------------------------------------------------------
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

dfmrg<-merge(young,old)%>%dplyr::rename(gene_name=transcript)%>%
  dplyr::mutate(gene_name=gsub("_mRNA","",gene_name))

# 去除表达量的影响 -------------------------------------------------------------
expr<-read_xlsx("yeast_chronoAge_rpkm.xlsx")
colnames(expr)[1]<-"gene_name"

expr<-expr%>%dplyr::group_by(gene_name)%>%
  dplyr::summarise(expr_young=mean(c(BY_D1_1,BY_D1_2,BY_D1_3)),
                   expr_old=mean(c(BY_D4_1,BY_D4_2,BY_D4_3)))

dfmrg<-dplyr::left_join(dfmrg,expr,by="gene_name") %>%
  dplyr::filter(young<0.5 & old <0.5) %>%
  dplyr::filter(!((expr_old/expr_young)>1)) %>%
  dplyr::select(young,old) 

wilcox.test(dfmrg$young,dfmrg$old,paired = T,alternative = "l")$p.value
#n=389,4.960404e-10


range(dfmrg)
readthr <- lapply( c(0.001,0.01,0.02,0.05),function(thres){
  dfdat<-dfmrg %>%dplyr::filter(young > thres | old > thres) %>% 
    dplyr::mutate(thres = thres) %>%
    dplyr::mutate(mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value)
}) %>% rbind.fill()

readthr$mw.p%>%unique()
readthr$thres%>%table()
dfplot<-readthr%>%reshape2::melt(id.vars =c("thres","mw.p") )

Fig.S4c<-ggplot(dfplot,aes(x=factor(thres),y=value,fill=variable))+
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
  annotate("text",x=1.1,y=0.2,label=expression(italic(P) * " = " * "4.9×" * 10**-10),size=2)+
  annotate("text",x=2,y=0.2,label=expression("4.6×" * 10**-10),size=2)+
  annotate("text",x=3,y=0.2,label=expression("2.6×" * 10**-10),size=2)+
  annotate("text",x=4,y=0.2,label=expression("4.6×" * 10**-10),size=2)+
  annotate("text",x=1,y=0.19,label=expression("(N" * " = " * "389)" ),size=2)+
  annotate("text",x=2,y=0.19,label=expression("(N" * " = " * "384)" ),size=2)+
  annotate("text",x=3,y=0.19,label=expression("(N" * " = " * "344)" ),size=2)+
  annotate("text",x=4,y=0.19,label=expression("(N" * " = " * "217)" ),size=2)

Fig.S4c



# 合并图 ---------------------------------------------------------------------

library(cowplot)

Fig.S4 <- plot_grid(Fig.S4a,Fig.S4b,Fig.S4c,ncol = 3,
                    labels = c("a","b","c"))  

ggsave("Fig.S4.pdf",Fig.S4,width =8, height =4)
