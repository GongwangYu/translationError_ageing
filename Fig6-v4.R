
# human minsincorp --------------------------------------------------------
dfqsubs<-read.csv("detct_quant_human.csv")

dfqsubs_flt<-lapply(1:nrow(dfqsubs), function(i){
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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%
    dplyr::mutate(error=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(error<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("error")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("error")%>%mean(na.rm=T)
  
  protein=dfqsubs[i,]$proteins %>%str_split( " ") %>%unlist() %>%unique()
  if(length(protein)>1){return(NULL)}else{
    protein=gsub("\\.\\d+","",protein)
    return(data.frame(protein,young=young,old=old)) 
  }
  
})%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%
  dplyr::group_by(ensg=protein)%>%
  dplyr::summarise(young=mean(young,na.rm = T),old=mean(old,na.rm = T),
                   site_n=length(ensg))%>%
  filter(site_n>=1)%>%rowwise()%>%
  mutate(rel_incr=(old-young)/((old+young)/2)) 

#椎间盘组织中基因的表达量
human_enst_ensg<-read_tsv("human_biomart_anno.txt") %>%dplyr::select(1,3)
colnames(human_enst_ensg)<-c("ensg","name")

human_desc_expr<-read_tsv("GSE245147_Degenerated_NO_Degenerated_RPKM.txt")%>%
  rowwise()%>%dplyr::mutate(expr=log2(mean(c(NO_Degenerated_1,NO_Degenerated_2,NO_Degenerated_3))))%>%
  dplyr::select(name=Geneid,len=Length,expr) %>%left_join(human_enst_ensg) %>%
  dplyr::select(ensg,expr)%>%unique()%>%na.omit()

dfmrg<-left_join(dfqsubs_flt,human_desc_expr,by="ensg")%>%na.omit() #N=140
cor.test(dfmrg$young,dfmrg$expr,method = "s")    #rho=0.3288222   ,p=7.288e-05
cor.test(dfmrg$old,dfmrg$expr,method = "s")      #rho=-0.1706945   ,p=0.04376
cor.test(dfmrg$rel_incr,dfmrg$expr,method = "s") #rho=-0.3972936   ,p=1.172e-06

Fig6a<- ggplot(dfmrg,aes(x=expr,y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin misincorporation rate")+
  theme_test()+
  annotate("text",x=-4,y=2.5,label=expression(italic(rho) * " = " * -0.40))+
  annotate("text",x=3,y=2.5,label=expression(italic(P) * " < " *10**-5))



# human readthr -----------------------------------------------------------
human_enst_ensg<-read_tsv("human_biomart_anno.txt") %>%dplyr::select(1,3)%>%unique()
colnames(human_enst_ensg)<-c("ensg","name")

human_expr<-read_tsv("human_rna_seq_genes_tpm.txt") %>%
  dplyr::select("Gene_Id","Young1","Young2")%>%
  separate(col = "Gene_Id",into = c("a","name"),sep = "_")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(expr=log2(mean(c(Young1,Young2))))%>%
  dplyr::select(name,expr)%>%merge(human_enst_ensg)%>%
  dplyr::select(ensg,expr)


load("rt_common.RData")
load("readcount_humanCell.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0; 
compRep <- list() # Composite replicates read count
names(readcount)
compRep$old<- composite(readcount[ names(readcount)[1:2]])
compRep$young<- composite(readcount[ names(readcount)[3:4]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3, exten_m5 )

young<-RTscore$young %>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)

dfmrg<-merge(young,old)%>%rowwise()%>%dplyr::rename(ensg=transcript)%>%
  dplyr::mutate(rel_incr=(old-young)/((old+young)/2))%>%left_join(human_expr)%>%na.omit()#N=903

cor.test(dfmrg$young,dfmrg$expr,method = "s")    #rho=-0.3803647  ,p < 2.2e-16
cor.test(dfmrg$old,dfmrg$expr,method = "s")      #rho=-0.3134226  ,p < 2.2e-16
cor.test(dfmrg$rel_incr,dfmrg$expr,method = "s") #rho=0.04228733  ,p=0.2042

Fig6b<- ggplot(dfmrg,aes(x=expr,y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin readthrough rate")+
  theme_test()+
  annotate("text",x=5,y=2.5,label=expression(italic(rho) * " = " * 0.04))+
  annotate("text",x=10,y=2.5,label=expression(italic(P) * " = " *0.20))


# mouse -------------------------------------------------------------------

{
  mouse_lung_epxr<-read_tsv("GSE124872_raw_counts_whole_lung_bulk.txt")%>%dplyr::select(1,6:12)
  mouse_lung_epxr_flt<-(mouse_lung_epxr[3:8]/(mouse_lung_epxr$Length/1000))
  mouse_lung_epxr_flt$old_wholelung_1<-mouse_lung_epxr_flt$old_wholelung_1/
    (sum(mouse_lung_epxr$old_wholelung_1)/10^6)
  mouse_lung_epxr_flt$old_wholelung_2<-mouse_lung_epxr_flt$old_wholelung_2/
    (sum(mouse_lung_epxr$old_wholelung_2)/10^6)
  mouse_lung_epxr_flt$old_wholelung_3<-mouse_lung_epxr_flt$old_wholelung_3/
    (sum(mouse_lung_epxr$old_wholelung_3)/10^6)
  mouse_lung_epxr_flt$young_wholelung_1<-mouse_lung_epxr_flt$young_wholelung_1/
    (sum(mouse_lung_epxr$young_wholelung_1)/10^6)
  mouse_lung_epxr_flt$young_wholelung_2<-mouse_lung_epxr_flt$young_wholelung_2/
    (sum(mouse_lung_epxr$young_wholelung_2)/10^6)
  mouse_lung_epxr_flt$young_wholelung_3<-mouse_lung_epxr_flt$young_wholelung_3/
    (sum(mouse_lung_epxr$young_wholelung_3)/10^6)
  mouse_lung_epxr<-cbind(mouse_lung_epxr[,1],mouse_lung_epxr_flt)%>%
    dplyr::rename(ensg = Geneid)%>%rowwise()%>%
    dplyr::mutate(expr_young=mean(young_wholelung_1,young_wholelung_2,young_wholelung_3),
                  expr_old=mean(old_wholelung_1,old_wholelung_2,old_wholelung_3))%>%dplyr::select(1,8,9)
  
}

dfqsubs<-read.csv("detct_quant_mouse.csv")

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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%
    dplyr::mutate(error=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(error<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("error")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("error")%>%mean(na.rm=T)
  
  protein=dfqsubs[i,]$proteins %>%str_split( " ") %>%unlist() %>%unique()
  if(length(protein)>1){return(NULL)}else{
    protein=gsub("\\.\\d+","",protein)
    return(data.frame(protein,young=young,old=old)) 
  }
  
},mc.cores = 3)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%
  dplyr::group_by(ensg=protein)%>%
  dplyr::summarise(young=mean(young,na.rm = T),old=mean(old,na.rm = T),
                   site_n=length(ensg))%>%
  filter(site_n>=1)%>%
  rowwise()%>%dplyr::mutate(rel_incr=(old-young)/((old+young)/2))

dfmrg<-left_join(dfqsubs_flt,mouse_lung_epxr,by="ensg")%>%na.omit() #N=131

cor.test(dfmrg$young,dfmrg$expr_young,method = "s")    #rho=0.2018316     ,p=0.02079
cor.test(dfmrg$old,dfmrg$expr_young,method = "s")      #rho=-0.06887791     ,p=0.4344
cor.test(dfmrg$rel_incr,dfmrg$expr_young,method = "s") #rho=-0.1760575     ,p=0.04427

Fig6c<- ggplot(dfmrg,aes(x=log2(expr_young),y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin misincorporation rate")+
  theme_test()+
  annotate("text",x=0,y=2.5,label=expression(italic(rho) * " = " * -0.18))+
  annotate("text",x=6,y=2.5,label=expression(italic(P) * " = " *0.04))


# mouse readthr -----------------------------------------------------------
load("mouse_liver_expr.RData")
load("rt_common.RData")
load("readcount_mouseAge.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[26:28]])
compRep$old<- composite(readcount[ names(readcount)[19:21]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young %>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)

dfmrg<-merge(young,old)%>%rowwise()%>% dplyr::rename(ensg=transcript)%>%
  dplyr::mutate(rel_incr=(old-young)/((old+young)/2))%>%left_join(expr)%>%na.omit() #N=535

cor.test(dfmrg$young,dfmrg$expr_young,method = "s")      #rho=-0.3286979     , p-value = 6.051e-15
cor.test(dfmrg$old,dfmrg$expr_young,method = "s")        #rho= -0.3152487      , p-value = 8.275e-14
cor.test(dfmrg$rel_incr,dfmrg$expr_young,method = "s")   #rho= 0.05603377        , p-value = 0.1956

Fig6d<- ggplot(dfmrg,aes(x=expr_young,y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin readthrough rate")+
  theme_test()+
  annotate("text",x=3,y=2.5,label=expression(italic(rho) * " = " * 0.06))+
  annotate("text",x=10,y=2.5,label=expression(italic(P) * " = " *0.19))


# yeast -------------------------------------------------------------------

# yeast_misincor ----------------------------------------------------------

dfqsubs<-read.csv("detct_quant_yeast.csv")
#来自../substitution_AA/yeast_ageing/combined_all_FDR0.02/txt/detct_quant.csv"

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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(error_ratio=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(error_ratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("Young_rep",sample,ignore.case=T))%>%
    getElement("error_ratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("Old_rep",sample,ignore.case=T))%>%
    getElement("error_ratio")%>%mean(na.rm=T)
  
  protein=dfqsubs[i,]$proteins %>%str_split( " ") %>%unlist() %>%unique()
  if(length(protein)>1){return(NULL)}else{
    protein=gsub("_mRNA","",protein)
    return(data.frame(protein,young=young,old=old))}
  
},mc.cores = 3)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%
  dplyr::group_by(ensg=protein)%>%
  dplyr::summarise(young=mean(young,na.rm = T),old=mean(old,na.rm = T),
                   site_n=length(ensg))%>% 
  filter(site_n>=1)%>%rowwise()%>%
  dplyr::mutate(rel_incr=(old-young)/((old+young)/2))


load("01.yeast.RData")
expr<-  data.frame(ensg=toPlot.yeast$orf,
                   expr=toPlot.yeast$logexpr)

dfmrg<-left_join(dfqsubs_flt,expr)%>%na.omit() #N=26
cor.test(dfmrg$young,dfmrg$expr,method = "s")       # rho=0.457808   ,p-value = 0.01868
cor.test(dfmrg$old,dfmrg$expr,method = "s")         # rho=-0.4311939  ,p-value = 0.02786
cor.test(dfmrg$rel_incr,dfmrg$expr,method = "s")    # rho=-0.5509538    ,p-value = 0.003535

Fig6e<- ggplot(dfmrg,aes(x=expr,y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin misincorporation rate")+
  theme_test()+
  annotate("text",x=4,y=2.5,label=expression(italic(rho) * " = " * -0.55))+
  annotate("text",x=6.5,y=2.5,label=expression(italic(P) * " = " *0.003))


# yeast_readthr -----------------------------------------------------------

load("rt_common.RData")
load("readcount_yeastChrAge.RData")
cds_m5 = 15;cds_m3 = 33;exten_m5 = 0
compRep <- list() # Composite replicates read count
names(readcount)
compRep$young<- composite(readcount[ names(readcount)[1:4]])
compRep$old<- composite(readcount[ names(readcount)[9:12]])
RTscore <- lapply(compRep, rt_efficiency, cds_m5, cds_m3,exten_m5)

young<-RTscore$young %>%dplyr::select(transcript,young=rte_ext)
old<-RTscore$old %>%dplyr::select(transcript,old=rte_ext)

dfmrg<-merge(young,old)%>% separate(transcript,into = c("ensg"),sep = "_")%>%rowwise()%>%
  dplyr::mutate(rel_incr=(old-young)/((old+young)/2))%>%left_join(expr) #792


cor.test(dfmrg$young,dfmrg$expr,method = "s")      #rho=-0.224274  , p-value = 3.123e-08
cor.test(dfmrg$old,dfmrg$expr,method = "s")        #rho= -0.2063598   , p-value = 3.733e-07
cor.test(dfmrg$rel_incr,dfmrg$expr,method = "s")   #rho= 0.005991014     , p-value = 0.884

Fig6f<- ggplot(dfmrg,aes(x=expr,y=rel_incr))+geom_point()+
  xlab("Expression level")+
  ylab("Aged-to-young relative increase\nin readthrough rate")+
  theme_test()+
  annotate("text",x=2,y=2.5,label=expression(italic(rho) * " = " * 0.006))+
  annotate("text",x=6,y=2.5,label=expression(italic(P) * " = " *0.88))

# 合并图 ---------------------------------------------------------------------

library(cowplot)
load("fig6g.RData")
Fig6af <- plot_grid(Fig6a,Fig6b,Fig6c,Fig6d,Fig6e,Fig6f
                    ,ncol = 2,byrow = T,
                  labels = c("a","d","b","e","c","f")) 

Fig6g<-plot_grid(fig6g,"",ncol=2,rel_widths = c(1.5,1),
                labels = c("g",""))
Fig6<-plot_grid(Fig6af,fig6g,ncol=1,labels = c("","g"),rel_heights  = c(4,1))
ggsave("Fig6.pdf",Fig6,width =5, height =10)
