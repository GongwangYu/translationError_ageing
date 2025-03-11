
# 最终的代码 -------------------------------------------------------------------
dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant.csv")

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

wilcox.test(dfqsubs_flt$young,dfqsubs_flt$old,paired=T, alternative = "l")
#N=1127,p-value = 1.689e-12  cutoff=0.01

save(dfqsubs_flt,file = "../RData.file/dfqsubs_flt_human.RData")


misincorp <- dfqsubs_flt %>% dplyr::rowwise()%>%
  dplyr::mutate(old_increased = (old - young) ) %>%
  getElement("old_increased")%>%median()
#8.562229e-05




# 去除RNA Editing ------------------------------------------------------------
rm(list=ls())
dfrawdat<-read_tsv("../../downdata/TABLE1_hg38.txt")%>%
  filter(Gene.refGene == "nonsynonymous SNV")%>%
  getElement("ExonicFunc.wgEncodeGencodeBasicV34")

dfRNAediting<-mclapply(1:length(dfrawdat), function(i){
  eachrow<-dfrawdat[i]
  enst<-str_split(eachrow,",")%>%unlist() %>%str_extract_all(":ENST\\d+") %>%unlist()
  edit<-str_split(eachrow,",")%>%unlist() %>%str_extract_all("p\\.\\w+") %>%unlist()
  dfdat<-data.frame(enst,edit)
  dfdat$enst<-gsub(":","",dfdat$enst)
  dfdat$edit<-gsub("p\\.","",dfdat$edit)
  
  return(dfdat)
},mc.cores = 50)%>%rbind.fill()


dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant_enst.csv")
dfqsubs_flt_mouse<-mclapply(1:nrow(dfqsubs), function(i){
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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T))%>%
    filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  origin<-dfqsubs[i,]$origin
  destination<-dfqsubs[i,]$destination
  positions<-dfqsubs[i,]$positions %>% str_split(" ") %>%unlist()
  ensts<-gsub("\\.+\\d+","",dfqsubs[i,]$proteins) %>% str_split(" ") %>%unlist()
  
  expected_RNAEditing<-paste(origin,positions,destination,sep = "")
  
  return(data.frame(enst=ensts,expected_RNAEditing,young=young,old=old))
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()

dfmrg<-left_join(dfqsubs_flt_mouse,dfRNAediting,by="enst")%>%
  filter(expected_RNAEditing ==edit) #没有RNA editing



####----去除SNP

dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant.csv") 
ensg<-mclapply(1:nrow(dfqsubs), function(i){
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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T))%>%
    filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  
  ensgs<-gsub("\\.+\\d+","",dfqsubs[i,]$proteins) %>% str_split(" ") %>%unlist() %>%unique()
  
  return(data.frame(enst=ensgs,young=young,old=old))
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%dplyr::select(1)
write.table(ensg,file = "ensg",row.names = F,col.names = F,quote = F)


snp<-read_tsv("../../downdata/human_snp_anon.txt")
colnames(snp)<-c("ensg","enst","enst_pos","protein_pos","variant","var_cons","protein_all")

snp<-snp%>%dplyr::filter(var_cons=="missense_variant")%>%
  separate(col = "protein_all",into = c("orign","dest"),sep = "\\/")%>%dplyr::rowwise()%>%
  dplyr::mutate(snp_obs= paste0(c(orign,protein_pos,dest),collapse = ""))%>%
  dplyr::select(1,9)%>%unique()



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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T))%>%
    filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  origin<-dfqsubs[i,]$origin
  destination<-dfqsubs[i,]$destination
  positions<-dfqsubs[i,]$positions %>% str_split(" ") %>%unlist()
  ensgs<-gsub("\\.+\\d+","",dfqsubs[i,]$proteins) %>% str_split(" ") %>%unlist()
  
  expected_snp<-paste(origin,positions,destination,sep = "")
  
  
  return(data.frame(ensg=ensgs,expected_snp,young=young,old=old))
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%unique()

ensg_snp_remov<-left_join(dfqsubs_flt,snp,by="ensg")%>%
  filter((expected_snp ==snp_obs)) %>%getElement("ensg")#4

#去除snp相关基因
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
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T))%>%
    filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  ensgs<-gsub("\\.+\\d+","",dfqsubs[i,]$proteins) %>% str_split(" ") %>%unlist()%>%unique()
  

  
  if(length(intersect(ensgs,ensg_snp_remov))==0){ 
    return(data.frame(young=young,old=old))}else{return(NULL)}
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()%>%unique()

wilcox.test(dfqsubs_flt$young,dfqsubs_flt$old,paired=T, alternative = "l")

#N=1012,p-value = 1.455e-06  cutoff=0.01

# 结束 ----------------------------------------------------------------------


load("~/disk3/Nascent_RNA/mystyle.print.RData")
library(ggthemes)
############################
##############
##统计质量
dfqsubs<-read.csv("raw/combined/txt/detct_quant.csv")%>%
  filter(protein=="ENST00000356839")%>%
  dplyr::select(c(starts_with("Base.intensity"),starts_with("Modification.intensity")))%>%t()%>%
  as.data.frame()
dfqsubs[dfqsubs==0]<-NA

dfname<-c(paste("base_old",1:11,sep = "-"),paste("base_young",1:11,sep = "-"),
          paste("base_old",12:22,sep = "-"),paste("base_young",12:22,sep = "-"),
          paste("base_old",23:33,sep = "-"),paste("base_young",23:33,sep = "-"),
          paste("mod_old",1:11,sep = "-"),paste("mod_young",1:11,sep = "-"),
          paste("mod_old",12:22,sep = "-"),paste("mod_young",12:22,sep = "-"),
          paste("mod_old",23:33,sep = "-"),paste("mod_young",23:33,sep = "-"))

colnames(dfqsubs)<-dfname
dfqsubs_young<-dfqsubs%>%dplyr::select(contains("young"))
dfqsubs_old<-dfqsubs%>%dplyr::select(contains("old"))
dfall<-cbind(dfqsubs_young,dfqsubs_old)
missmap(dfall, col=c("black", "grey"), main = 'Missing Map',rank.order = F,
        margins = c(5, 5))



dfsum<-read_tsv("txt_AllWindow/txt/summary.txt")%>%
  dplyr::select(1:3,23,27,31,35)

dfsum_stat<-fread("txt_AllWindow/txt/summary.txt")%>%filter(`Raw file` !="Total") %>% rowwise()%>%
  mutate(Sample=stringr::str_match(Experiment,pattern = "oung|old")%>%as.character())%>%
  group_by(Sample)%>%dplyr::summarise(dfcover=sum(`MS/MS submitted`))


AllCodons<-c()
allnucle<-c("A","T","C","G")
for (a in allnucle){
  a=a
  for (b in allnucle) {
    b=b
    for (c in allnucle) {
      c=c
      df=paste(a,b,c,sep = "")
      AllCodons<-c(AllCodons,df)
    }
  }
}

dfsubs_type<-read.csv("txt_AllWindow/txt/detct_quant.csv")%>%
  mutate(mysubs=paste(codon,destination,sep = " to "))%>%
  dplyr::select(c(!matches("intensity|Cluster|Amount|features")))%>%
  group_by(mysubs)%>%mutate(dfnum=length(mysubs))%>%
  filter(grepl("TGC",mysubs))

setdiff(AllCodons, dfsubs_type$codon)

###########################
dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant.csv")

dfqsubs_flt_human<-mclapply(1:nrow(dfqsubs), function(i){
  dfBP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"
  #dfBP$bp[dfBP$bp==0]<-sample(dfBP$bp[dfBP$bp>0],length( dfBP$bp[dfBP$bp==0]),replace = T)#填补缺失值
  # dfBP_young<-dfBP%>%dplyr::filter(grepl("young",sample,ignore.case=T))
  # if(length(dfBP_young$bp[dfBP_young$bp>0])>0){
  #   dfBP_young$bp[dfBP_young$bp==0]<-sample(dfBP_young$bp[dfBP_young$bp>0],
  #                                           length( dfBP_young$bp[dfBP_young$bp==0]),replace = T)  
  # }
  # 
  # dfBP_old<-dfBP%>%dplyr::filter(grepl("old",sample,ignore.case=T))
  # if(length(dfBP_old$bp[dfBP_old$bp>0])>0){
  #   dfBP_old$bp[dfBP_old$bp==0]<-sample(dfBP_old$bp[dfBP_old$bp>0],
  #                                       length( dfBP_old$bp[dfBP_old$bp==0]),replace = T)  
  # }
  # 
  # dfBP<-rbind(dfBP_young,dfBP_old)
  
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

wilcox.test(dfqsubs_flt_human$young,dfqsubs_flt_human$old,paired=T, alternative = "l")
save(dfqsubs_flt_human,file = "../RData.file/dfqsubs_flt_human.RData")
#N=1127,p-value = 1.689e-12  cutoff=0.01.
#N=1135,p-value = 1.689e-12  cutoff=0.01,不去除NAN.
#N=1135,p-value = 6.545162e-20  cutoff=0.01,不去除NAN，不paired=T
human_misincorp<-dfqsubs_flt_human
save(human_misincorp,file = "/mnt/data3/disk/yugongwang/translation_error/substitution_AA/human_aging_elife_2020/human_misincorp.RData")

#####################
load("/mnt/data3/disk/yugongwang/translation_error/substitution_AA/human_aging_elife_2020/human_misincorp.RData")

misincorp <- human_misincorp %>% 
  dplyr::rename(young=young_rate,old=old_rate) %>% 
  mutate(old_increased = (old - young)/(old + young) ) %>%
  mutate(avg = (old + young)/2)

range(misincorp$avg)
dfdat=lapply(10 ** seq(-7,-3,by=1),function(thres){
  misincorp %>% 
    filter(avg > thres) %>%
    summarise(
      n = length(old_increased),
      n_increased = sum(old_increased > 0), 
      n_decreased = sum(old_increased < 0),
      frac_increased = n_increased / n,
      frac_decreased = n_decreased / n,
      mean_increase = mean(old_increased),
      se=sd(old_increased)/sqrt(n),
      mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value) %>%
    mutate(thres = thres)
}) %>% rbind.fill()

dfplot<-dfdat%>%select(thres,frac_increased,frac_decreased,mean_increase,se,mw.p)%>%
  melt(id.vars =c("thres","mean_increase","se","mw.p") )

plot1f<-ggplot(dfplot,aes(x=factor(thres),y=ifelse(variable=="frac_increased",value,-value)) )+
  geom_bar( aes(fill=variable),stat = "identity")+
  geom_point(aes(x=factor(thres), y=ifelse(variable=="frac_increased",mean_increase,NA) ) )+
  geom_errorbar(aes(ymin=ifelse(variable=="frac_increased",mean_increase-se,NA) ,
                    ymax=ifelse(variable=="frac_increased",mean_increase+se,NA)),width=0.3)+
  geom_line(data=dfplot%>%filter(variable=="frac_increased"),
             aes(x=factor(thres), y=mean_increase),group=1 )+
  # geom_text(aes(x=factor(thres), y=ifelse(variable=="frac_increased",value,NA),
  #               label=ifelse(variable=="frac_increased" & mw.p<0.05,"*",NA)) ,size=8)+
  scale_fill_manual(name="",values = c("#F39B7FFF", "#8491B4FF"),
                    breaks = c("frac_increased","frac_decreased"),
                    labels = c("Increased","Decreased"))+
  scale_y_continuous( breaks = c(-0.5,-0.25,0,0.25,0.5), labels = abs, 
                     sec.axis = sec_axis(name = "Aged-to-young relative increase\nin translation error rates",
                                         trans=~.*1,breaks =c(-0.5,-0.25,0,0.25,0.5),  
                                         labels = c(-0.5,-0.25,0,0.25,0.5)),
                     expand = expansion(mult = c(0.1,0.1)))+
  geom_hline(aes(yintercept = 0),linetype=5,color="#4DBBD5FF")+
  labs(x="Threshold for translation error rates (>x)",y="Fraction of translation error sites")+
  scale_x_discrete(labels = c(expression(10**-7),expression(10**-6),expression(10**-5),
                                expression(10**-4),expression(10**-3)),
                                breaks = c(10**-7,10**-6,10**-5,10**-4,10**-3))+
  theme_test()+
  labs(title = "Human")+theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.y.left  = element_text(color =c(rep("#8491B4FF",2),"#4DBBD5FF",rep("#F39B7FFF",2))),
        axis.text.y.right = element_text(color =c(rep("black",2),"#4DBBD5FF",rep("black",2))),
        legend.position = "none",legend.key = element_blank())

fig1<-list()
fig1[[6]]<-plot1f

########################

p + theme(axis.text.y = element_text(colour=x_cols))
dfqsubs_flt_human_plot<-dfqsubs_flt_human%>%mutate(diff=old_rate-young_rate)

dfqsubs_flt_human_plot<-dfqsubs_flt_human_plot[order(dfqsubs_flt_human_plot$diff),]
dfqsubs_flt_human_plot$range<-1:nrow(dfqsubs_flt_human_plot)
dfqsubs_flt_human_plot$dftype<-dfqsubs_flt_human_plot$diff>0

dfqsubs_flt_human_plot$dftype%>%table()
ggplot(dfqsubs_flt_human_plot, aes(range, diff,fill = dftype)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("grey","red"))+
  annotate("text",x=900,y=0.002,label="684",color="red")+
  annotate("text",x=300,y=-0.002,label="443",color="grey")+
  ylab("Difference of misincorporation rate\nbwteen aged and young cells)")+
  xlab("Rank of genes")

#QQ plot
x_human <- dfqsubs_flt_human$young_rate
y_human <- dfqsubs_flt_human$old_rate
sx_human <- sort(x_human); sy_human <- sort(y_human)

plot1D<-ggplot() + geom_abline(a=0,b=1,lwd=1,lty=2,col="gray")+
  geom_point(aes(x=log2(x_human), y=log2(y_human)),cex=0.8)+
  # scale_x_continuous(limits = c(0,0.01),breaks = c(0,0.005,0.01),labels =c(0,0.005,0.01))+
  # scale_y_continuous(limits = c(0,0.01),breaks = c(0,0.005,0.01),labels =c(0,0.005,0.01))+
  # annotate("text",x=0.007,y=0.002,label=expression(italic(P) * " < " * 10^{-11}))+
  xlab("log2 (Misincorporation rate\nin young cells)")+
  ylab("log2 (Misincorporation rate\nin aged cells)")+
  mystyle.print()+
  labs(title = "Human")+
  theme(plot.title=element_text(hjust=0.5))





{
  
  dfdat_plot<-(dfqsubs_flt+10^-5)%>%pivot_longer(cols = contains("rate"))
  
  ggplot(dfdat_plot,aes(x=name,y=value,fill=name))+
    geom_boxplot(outlier.shape = NA,linetype="dashed")+ #去除极值
    stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                     fill=name),
                 outlier.shape = NA)+
    stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),
                 width=0.2,outlier.shape = NA)+
    stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),
                 width=0.2,outlier.shape = NA)+
    scale_fill_manual(values=c( "red","grey"))+
    scale_y_continuous(breaks=trans_breaks("log10",function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)),
                       trans=log10_trans(),limits = c(-6,0.1))+
    scale_x_discrete(limits=c("young_rate","old_rate"),
                     labels=c("young_rate"="Young","old_rate"="Old"))+
    geom_segment(aes(x=1,y=0.02,xend=2,yend=0.02))+
    annotate("text",x=1.5,y=0.05,label=expression(italic(P) * " < " * 10^{-11}))+
    theme_bw(base_size = 12)+
    labs(title = "Human")+
    theme(legend.position = "",plot.title=element_text(hjust=0.5))+
    xlab("")+ylab("Misincorporation rate") +
    theme(axis.text.x = element_text(size=12))
}

##图的比例：200x300



{
# {dfdat_plot<-(dfqsubs_flt+10^-5)%>%pivot_longer(cols = contains("rate"))
#   ggplot(dfdat_plot,aes(x=name,y=log10(value),fill=name))+
#     geom_boxplot()+
#     geom_signif(comparisons = list(c("young_rate","old_rate")),
#                 annotations = "P = 10e-15",
#                 test=wilcox.test,
#                 textsize = 5)+
#     
#     scale_fill_manual(values=c("red", "grey"))+
#     scale_y_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
#                        limits=c(-5,-1.5))+
#     scale_x_discrete(limits=c("young_rate","old_rate"),
#                      labels=c("young_rate"="Young","old_rate"="Old"))+
#     theme_bw()+
#     labs(title = "Human")+
#     
#     theme(legend.position = "",plot.title=element_text(hjust=0.5,size=16))+
#     xlab("")+ylab("Misincorporation rate") +
#     theme(axis.title = element_text(size=16),
#           axis.text = element_text(size = 16),
#           legend.text =  element_text(size = 16)
#           
#     )
# }

}
###########
#error分析
#per site
#dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant.csv")%>%
dfqsubs_flt<-dfqsubs%>%  
  dplyr::select(c(starts_with("Ratio")))%>%as.data.frame()
dfqsubs_flt[dfqsubs_flt=="Inf"]<-NA
dfqsubs_flt[dfqsubs_flt>0.001]<-NA

young<-dfqsubs_flt%>% dplyr::select(c(contains("young")))%>%
  dplyr::summarise(young_rate=rowMeans(.,na.rm = T))
old<-dfqsubs_flt%>%dplyr::select(c(contains("old")))%>%
  dplyr::summarise(old_rate=rowMeans(.,na.rm = T))

dfdat<-cbind(young,old)%>%filter(young_rate >0 & old_rate >0)%>%na.omit()

wilcox.test(dfdat$young_rate,dfdat$old_rate,paired = T,alternative = "l")
#N=1014,p-value = 2.290547e-58, cutoff=0.1
#N=1125,p-value = 1.149091e-75  cutoff=0.01

t.test(dfdat$young_rate,dfdat$old_rate,paired = T,alternative = "l")
#N=1482,p-value = 1 cutoff=0.1
#N=1125,p-value = 4.065e-07 cutoff=0.01


##绘图
{
#散点图
{
  dfdat_plot<-dfdat+10^-5
  dfmerg<-lapply(1:nrow(dfdat_plot), function(i){
    eachgene<-dfdat_plot[i,]
    if(nrow(eachgene)==0){return(NULL)}
    if(eachgene$old_rate > eachgene$young_rate){
      eachgene$dftype<-"UP"} 
    else if(eachgene$old_rate < eachgene$young_rate){
      eachgene$dftype<-"down"} 
    else{eachgene$dftype<-"NO"} 
    return(eachgene)
  }) %>%rbind.fill
  
  table(dfmerg$dftype)
  binom.test(685,1125,alternative = "greater")
  
  ggplot(dfmerg,aes(x=log10(old_rate ),y=log10(young_rate) ,color=dftype) )+
    geom_point(shape=1,alpha=1)+
    scale_color_manual(values = c("black","red"))+
    scale_x_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
                       limits=c(-5,-2))+
    scale_y_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
                       limits=c(-5,-2))+
    geom_abline(slope = 1,intercept = 0,lty = 'dashed',size = 0.8) +
    annotate("text",x=-4,y=-2.2,label="440",color="black")+
    annotate("text",x=-2.2,y=-4,label="685",color="red")+
    annotate("text",x=-3,y=-4.8,label=expression(italic(P) * " = " * 10^-11),color="black")+
    xlab("Misincorporation rate (Old)")+
    ylab("Misincorporation rate (Young)")+
    theme_test()+
    #mystyle.print()+
    labs(title = "Human")+
    theme(legend.position = "",
          plot.title=element_text(hjust=0.5))
  
}

{
  dfdat_plot<-(dfdat+10^-5)%>%pivot_longer(cols = contains("rate"))
  ggplot(dfdat_plot,aes(x=name,y=log10(value),fill=name))+
    geom_boxplot()+
    geom_signif(comparisons = list(c("young_rate","old_rate")),
                annotations = "P < 10e-11",
                test=wilcox.test)+
    scale_fill_manual(values=c("red", "grey"))+
    scale_y_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
                       limits=c(-5,-1.5))+
    scale_x_discrete(limits=c("young_rate","old_rate"),
                     labels=c("young_rate"="Young","old_rate"="Old"))+
    theme_test()+
    labs(title = "Human")+
    
    theme(legend.position = "",plot.title=element_text(hjust=0.5))+
    xlab("")+ylab("Misincorporation rate") 
}













##绘图
dfdat_plot<-(dfdat+10^-5)%>%pivot_longer(cols = contains("rate"))
ggplot(dfdat_plot,aes(x=name,y=log10(value),fill=name))+
  geom_boxplot()+
  scale_fill_manual(values=c("red", "grey"))+
  scale_x_discrete(limits=c("young_rate","old_rate"),
                   labels=c("young_rate"="Young","old_rate"="Old"))+
  theme_test()+
  labs(title = "Human")+
  
  theme(legend.position = "")+
  xlab("")+ylab("Log10 (misincorporation rate)") 




dfdat<-(dfdat+10^-5)%>%pivot_longer(cols = contains("rate"))

ggplot(dfdat,aes(x=name,y=log10(value),fill=name))+
  geom_point(position = "jitter",shape=21,alpha=0.3)+
  stat_summary(fun.y=median, geom="point", shape=18,
               size=8, color="red")+
  scale_fill_manual(values=c("blue", "grey"))+
  scale_x_discrete(limits=c("young_rate","old_rate"),
                   labels=c("young_rate"="Young","old_rate"="Old"))+
  theme_test()+

  theme(legend.position = "")+
  xlab("")+ylab("Log10 (misincorporation rate)")
 
  
ggplot(dfdat,aes(x=name,y=log10(value),fill=name))+
  geom_boxplot()+
  scale_fill_manual(values=c("red", "grey"))+
  scale_x_discrete(limits=c("young_rate","old_rate"),
                   labels=c("young_rate"="Young","old_rate"="Old"))+
  theme_test()+
  
  theme(legend.position = "")+
  xlab("")+ylab("Log10 (misincorporation rate)")  




  
  ggplot(dfdat, aes(x=name,y=log10(value))) + 
  geom_dotplot(binaxis='y', stackdir='center',stackratio=0.5, dotsize=0.5)+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")

  ggplot(ToothGrowth, aes(x=dose, y=len, fill=dose)) +
    geom_dotplot(binaxis='y', stackdir='center')

dfmerg<-lapply(1:nrow(dfdat), function(i){
  eachgene<-dfdat[i,]
  if(nrow(eachgene)==0){return(NULL)}
  if(eachgene$old_rate > eachgene$young_rate){
    eachgene$dftype<-"UP"} 
  else if(eachgene$old_rate < eachgene$young_rate){
    eachgene$dftype<-"down"} 
  else{eachgene$dftype<-"NO"} 
  return(eachgene)
}) %>%rbind.fill

table(dfmerg$dftype)
binom.test(636,1044,alternative = "greater")

ggplot(dfmerg,aes(x=log10(old_rate ),y=log10(young_rate) ,color=dftype) )+
  geom_point(shape=1,alpha=1)+
  scale_color_manual(values = c("blue","red"))+
  scale_x_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
                     limits=c(-5,-2))+
  scale_y_continuous(breaks=c(-5,-4,-3,-2),label=(c("10e-5","10e-4","10e-3","10e-2")),
                     limits=c(-5,-2))+
  geom_abline(slope = 1,intercept = 0,lty = 'dashed',size = 0.8) +
  annotate("text",x=-4.5,y=-2.5,label="408",color="blue")+
  annotate("text",x=-2.5,y=-4.5,label="636",color="red")+
  annotate("text",x=-3,y=-4.8,label=expression(italic(P) * " < " * 10^{-6}),color="black")+
  xlab("Misincorporation rate (Old)")+
  ylab("Misincorporation rate (Young)")+
  mystyle.print()+
  labs(title = "Human")+
  theme(legend.position = "",
        plot.title=element_text(hjust=0.5))

}

##########################
#per substitution level


dfsubs_compile<-mclapply(1:nrow(dfqsubs),function(i){
  eachrow<-dfqsubs[i,]
  protein_index<-sample(1:length(str_split(eachrow$proteins," ")%>%unlist()),1 )
  codon<-((str_split(eachrow$codons," ")%>%unlist())[protein_index])%>%gsub("T","U",.)
  destination<-eachrow$destination
  
  dfBP<-eachrow%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"
  
  dfDP<-eachrow%>%dplyr::select(c(starts_with("Modification.intensity")))%>%t()%>%
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
  
  return(data.frame(codon,destination,young,old))
},mc.cores = 50)%>%rbind.fill()%>%filter(young>0 | old >0)%>%na.omit()%>%
  mutate(destination=ifelse((destination=="I" | destination=="L") ,"I/L",destination))%>%
  mutate(substype=paste(codon,destination,sep = " to "))

save(dfsubs_compile,file ="../RData.file/dfsubs_compile_desc.RData" )

###############################














dfsubs_compile<-mclapply(1:nrow(dfqsubs),function(i){
  eachrow<-dfqsubs[i,]
  protein_index<-sample(1:length(str_split(eachrow$proteins," ")%>%unlist()),1 )
  eachrow$protein<-(str_split(eachrow$proteins," ")%>%unlist())[protein_index]
  eachrow$codon<-((str_split(eachrow$codons," ")%>%unlist())[protein_index])%>%gsub("T","U",.)
  eachrow$position<-(str_split(eachrow$positions," ")%>%unlist())[protein_index]
  return(eachrow)
},mc.cores = 60)%>%rbind.fill()

dfsubs_compile[dfsubs_compile=="Inf"]<-NA

cutoff=0.01

young<-dfsubs_compile%>%
  mutate(mysubs=paste(codon,destination,sep = " to "))%>%
  dplyr::select(c(mysubs,starts_with("Ratio")))%>%
  dplyr::select(c(mysubs,contains("young")))%>%
  pivot_longer(cols = starts_with("Ratio"),names_to = "sample",
               values_to = "young_rate")%>%
  filter(young_rate<=cutoff)%>%
  group_by(mysubs)%>%dplyr::summarise(young_rate=mean(young_rate,na.rm=T))
  

old<-dfsubs_compile%>%
  mutate(mysubs=paste(codon,destination,sep = " to "))%>%
  dplyr::select(c(mysubs,starts_with("Ratio")))%>%
  dplyr::select(c(mysubs,contains("old")))%>%
  pivot_longer(cols = starts_with("Ratio"),names_to = "sample",
               values_to = "old_rate")%>%
  filter(old_rate<=cutoff)%>%
  group_by(mysubs)%>%dplyr::summarise(old_rate=mean(old_rate,na.rm=T))


dfdat<-merge(young,old,all=T)%>%na.omit()%>%filter(young_rate >0 | old_rate >0)

wilcox.test(dfdat$young_rate,dfdat$old_rate,paired = T,alternative = "l")
#N=186,p-value = 1.613e-15 cutoff=0.1
#N=142,p-value = 5.136e-16 cutoff=0.01


t.test(dfdat$young_rate,dfdat$old_rate,paired = T,alternative = "l")
#N=235,p-value = 0.2438 cutoff=0.1
#N=172,p-value = 0.08568 cutoff=0.01

###############################################end

#########################################
#画热图
dfqsubs_type<-dfsubs_compile%>%
  mutate(mysubs=paste(codon,destination,sep = " to "))

AA=dfqsubs_type%>%dplyr::select(mysubs)

amino_acids = c("A","C","D","E","F","G","H","K","I/L",
                "M","N","P","Q","R","S","T","V","W","Y")


AllCodons<-c()
allnucle<-c("U","C","A","G")
for (a in allnucle){
  a=a
  for (b in allnucle) {
    b=b
    for (c in allnucle) {
      c=c
      df=paste(a,b,c,sep = "")
      AllCodons<-c(AllCodons,df)
    }
  }
}
amino_acids_64 = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' %>%
  str_split("")%>%unlist()

dfcodon_AA_coord<-data.frame(codon=AllCodons,AA=amino_acids_64)

allSubsType<-c()
for (a in AllCodons){
  a=a
  for (b in amino_acids) {
    b=b
      df=paste(a,b,sep = " to ")
      allSubsType<-c(allSubsType,df)
    
  }
}

thistype<-allSubsType[10]
thistype<-"GCC to S"

dfmtr<-mclapply(allSubsType, function(thistype){
  eachtype<-dfqsubs_type%>%filter(mysubs==thistype)
  if(nrow(eachtype)==0){return(data.frame(thistype,myratio_O_Y=NA,mypvalue=NA))}
  
  young<-eachtype%>% dplyr::select(c(mysubs,starts_with("Ratio")))%>%
    dplyr::select(c(mysubs,contains("young")))%>%
    pivot_longer(cols = starts_with("Ratio"),names_to = "sample",
                 values_to = "young_rate")%>%filter(!is.infinite(young_rate))%>%
    filter(young_rate<0.01)
  
  if(nrow(young)==0){return(data.frame(thistype,myratio_O_Y=NA,mypvalue=NA))}
  
  old<-eachtype%>% dplyr::select(c(mysubs,starts_with("Ratio")))%>%
    dplyr::select(c(mysubs,contains("old")))%>%
    pivot_longer(cols = starts_with("Ratio"),names_to = "sample",
                 values_to = "old_rate")%>%filter(!is.infinite(old_rate))%>%
    filter(old_rate<0.01)
  
  if(nrow(old)==0){return(data.frame(thistype,myratio_O_Y=NA,mypvalue=NA))}
  
 mypvalue<- wilcox.test(young$young_rate,old$old_rate)
 mypvalue<-mypvalue$p.value  
 myratio_O_Y<-mean(old$old_rate,na.rm=T)/mean(young$young_rate,na.rm=T)
 return(data.frame(thistype,myratio_O_Y,mypvalue))
},mc.cores = 60)%>%rbind.fill()



dfmtr_fil<-dfmtr%>%separate(col=thistype,into = c("codon","destination"),sep = " to ")

dfmtr_fil$myratio_O_Y <-log2(dfmtr_fil$myratio_O_Y)

dfmtr_fil$myratio_O_Y[dfmtr_fil$myratio_O_Y=="Inf"]<- 
  max(dfmtr_fil$myratio_O_Y[dfmtr_fil$myratio_O_Y !="Inf"],na.rm = T)
dfmtr_fil$myratio_O_Y[dfmtr_fil$myratio_O_Y=="-Inf"]<- 
  min(dfmtr_fil$myratio_O_Y[dfmtr_fil$myratio_O_Y !="-Inf"],na.rm = T)

dfmtr_fil$myratio_O_Y[dfmtr_fil$myratio_O_Y=="NaN"]<- NA

dfmtr_fil$stars<-ifelse(dfmtr_fil$mypvalue<0.05,"*","")


ggplot(dfmtr_fil,aes(x=destination,y=codon))+
  geom_tile(aes(fill=myratio_O_Y))+
  scale_fill_gradient2("Ratio(Old/young)",
                       na.value=NA,low = "blue",high = "red",mid = "white",
                       guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                            barwidth =8,barheight = 1, ticks.colour="black",ticks.linewidth = 1 ))+
  labs(x="Destination amino acdis",y="Original codon")+
  theme(
        axis.ticks.y = element_blank(),
        panel.background=element_blank())+
  theme_test()+
  theme(legend.position = "top")+
  geom_text(aes(label=stars),color="black",size=4,
            hjust="middle",vjust="bootom")


###################################














aa=dfmtr%>%filter(myratio_O_Y>1 & mypvalue<0.05)
bb=dfmtr%>%filter(myratio_O_Y<1 & mypvalue<0.05)


dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")
dfqsubs$Base[is.na(dfqsubs$Base)]<-0
dfqsubs<-dfqsubs%>%mutate(errorFre=DP/(Base+DP))


original_aa<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
destination_aa<-c("A","C","D","E","F","G","H","I","L","K","M","N","P","Q","R","S","T","V","W","Y")

dfmatrix<-lapply(original_aa, function(thisAA){
  each_original_aa<-dfqsubs%>%filter(origin==thisAA)
  
  young<-each_original_aa%>%filter(grepl("young",Sample))%>%
    group_by(destination)%>%dplyr::summarise(rate_Y=mean(errorFre,na.rm = T))
  old<-each_original_aa%>%filter(grepl("old",Sample))%>%
    group_by(destination)%>%dplyr::summarise(rate_O=mean(errorFre,na.rm = T))
  
  df_old_young_ratio<-merge(young,old)%>%mutate( ratio=log10( (rate_O+ 10^-4)/(rate_Y+10^-4) ) )
  dfMean=data.frame(des_aa="Mean", ratio=mean(df_old_young_ratio$ratio,na.rm=T))
  
  dfdat<-data.frame(des_aa=destination_aa,ratio=0)
  dfdat$ratio<-df_old_young_ratio$ratio[match(dfdat$des_aa,df_old_young_ratio$destination)]
  dfdat<-rbind(dfdat,dfMean)
  dfdat= t(dfdat[,2])%>%as.data.frame()
  colnames(dfdat)<-c(destination_aa,"Origin_Mean")
  return(dfdat)
})%>%rbind.fill()
rownames(dfmatrix)<-original_aa

dfcolMean=data.frame(Destination_Mean=colMeans(dfmatrix,na.rm = T))  
dfcolMean=dfcolMean%>%t()%>%as.data.frame()

dfplot<-rbind(dfmatrix,dfcolMean)
dfplot[21,21]<-NA
dt<-as.matrix(dfplot)







#############################

dfmerg<-lapply(1:nrow(dfdat), function(i){
  eachgene<-dfdat[i,]
  if(nrow(eachgene)==0){return(NULL)}
  if(eachgene$rate_old >= eachgene$rate_young){
    eachgene$dftype<-"UP"} 
  else if(eachgene$rate_old < eachgene$rate_young){
    eachgene$dftype<-"down"} 
  return(eachgene)
}) %>%rbind.fill

table(dfmerg$dftype)
binom.test(90,154)

ggplot(dfmerg,aes(x=( rate_old ),y=(rate_young) ,color=dftype) )+
  geom_point(shape=1,alpha=0.3)+
  scale_color_manual(values = c("blue","red"))+
  scale_x_log10(limits=c(10^-5,0.01))+
  scale_y_log10(limits=c(10^-5,0.01))+
  geom_abline(slope = 1,intercept = 0,lty = 'dashed',size = 0.8) +
  annotate("text",x=10^-4.5,y=10^-2.5,label="90",color="blue")+
  annotate("text",x=10^-2.5,y=10^-4.5,label="64",color="red")+
  annotate("text",x=10^-3,y=10^-5,label="P = 0.02",color="black")+
  xlab("Misincorporation rate (Old)")+
  ylab("Misincorporation rate (Young)")+
  mystyle.print()+
  labs(title = "Human intervertebral discs")+
  theme(legend.position = "",
        plot.title=element_text(hjust=0.5))


###########################end



















dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")
#dfqsubs=read.csv("txt_AllWindow/txt/qSubs.csv")

young=dfqsubs%>%filter(grepl("young",Sample))%>%filter(!is.na(Base))%>%
  filter(Base>0 & DP >0)%>%select(c(14,15))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_Y=mean(Ratio,na.rm = T))

old=dfqsubs%>%filter(grepl("old",Sample))%>%filter(!is.na(Base))%>%
  filter(Base>0 & DP >0)%>%select(c(14,15))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_O=mean(Ratio,na.rm = T))

dfdat<-merge(young,old,all = T)
dfdat[is.na(dfdat)]<-0
wilcox.test(dfdat$rate_Y,dfdat$rate_O,paired = T,alternative = "less")
#p-value = 1.071e-06

dfmerg_human<-lapply(1:nrow(dfdat), function(i){
  eachgene<-dfdat[i,]
  if(nrow(eachgene)==0){return(NULL)}
  if(eachgene$rate_O > eachgene$rate_Y){
    eachgene$dftype<-"UP"} 
  else if(eachgene$rate_O < eachgene$rate_Y){
    eachgene$dftype<-"down"} 
  else{eachgene$dftype<-"NO"} 
  return(eachgene)
}) %>%rbind.fill
dfmerg_human$spec<-"Human"
table(dfmerg$dftype)
binom.test(330,479,alternative = "greater")

plot2a<-ggplot(dfmerg,aes(x=log2( rate_O+10^-4 ),y=log2(rate_Y+10^-4) ,color=dftype) )+geom_point(alpha=0.3)+
  scale_color_manual(values = c("blue","red"))+
  xlim(c(-15,10))+
  ylim(c(-15,10))+
  geom_abline(slope = 1,intercept = 0,lty = 'dashed',size = 0.8) +
  annotate("text",x=4,y=-9,label="330",color="red")+
  annotate("text",x=-9,y=4,label="149",color="blue")+
  annotate("text",x=-5,y=8,label=expression(italic(P) * " < " * 10^{-5}),color="black")+
  xlab("")+
  ylab("Misincorporation rate (Young)")+
  mystyle.print()+
  labs(title = "Human")+
  theme(legend.position = "",
        plot.title=element_text(hjust=0.5))


#########################################
dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")
dfqsubs$Base[is.na(dfqsubs$Base)]<-0
dfqsubs<-dfqsubs%>%mutate(errorFre=DP/(Base+DP))


original_aa<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
destination_aa<-c("A","C","D","E","F","G","H","I","L","K","M","N","P","Q","R","S","T","V","W","Y")

dfmatrix<-lapply(original_aa, function(thisAA){
  each_original_aa<-dfqsubs%>%filter(origin==thisAA)
  
  young<-each_original_aa%>%filter(grepl("young",Sample))%>%
    group_by(destination)%>%dplyr::summarise(rate_Y=mean(errorFre,na.rm = T))
  old<-each_original_aa%>%filter(grepl("old",Sample))%>%
    group_by(destination)%>%dplyr::summarise(rate_O=mean(errorFre,na.rm = T))
  
  df_old_young_ratio<-merge(young,old)%>%mutate( ratio=log10( (rate_O+ 10^-4)/(rate_Y+10^-4) ) )
  dfMean=data.frame(des_aa="Mean", ratio=mean(df_old_young_ratio$ratio,na.rm=T))
  
  dfdat<-data.frame(des_aa=destination_aa,ratio=0)
  dfdat$ratio<-df_old_young_ratio$ratio[match(dfdat$des_aa,df_old_young_ratio$destination)]
  dfdat<-rbind(dfdat,dfMean)
  dfdat= t(dfdat[,2])%>%as.data.frame()
  colnames(dfdat)<-c(destination_aa,"Origin_Mean")
  return(dfdat)
})%>%rbind.fill()
rownames(dfmatrix)<-original_aa

dfcolMean=data.frame(Destination_Mean=colMeans(dfmatrix,na.rm = T))  
dfcolMean=dfcolMean%>%t()%>%as.data.frame()

dfplot<-rbind(dfmatrix,dfcolMean)
dfplot[21,21]<-NA
dt<-as.matrix(dfplot)


##########
thisAA<-original_aa[1]
dfpvalue<-lapply(original_aa, function(thisAA){
  each_original_aa<-dfqsubs%>%filter(origin==thisAA)
  
  young<-each_original_aa%>%filter(grepl("young",Sample)) %>%select(7,17)%>%
    dplyr::rename(ratio_Y=errorFre)
  old<-each_original_aa%>%filter(grepl("old",Sample))%>%select(7,17)%>%
    dplyr::rename(ratio_O=errorFre)
  
  dfaa<-intersect(young$destination%>%unique,old$destination%>%unique)
  eachpvalue<-lapply(dfaa, function(aa){
    eachAA_Y<-young%>%filter(destination==aa)%>%na.omit()
    eachAA_O<-old%>%filter(destination==aa)%>%na.omit()
    if(nrow(eachAA_O)==0 | nrow(eachAA_Y)==0 ){return(data.frame(dfdes=aa,pvalue=NA))}
    dftest=wilcox.test(eachAA_Y$ratio_Y,eachAA_O$ratio_O)
    return(data.frame(dfdes=aa,pvalue=dftest$p.value))
  })%>%rbind.fill()
  
  dfdat<-data.frame(des_aa=destination_aa,pvalue=NA)
  if(is.null(eachpvalue)){
    dfMean=data.frame(des_aa="Origin_Mean", pvalue=NA)
    dfdat<-rbind(dfdat,dfMean)
    dfdat= t(dfdat[,2])%>%as.data.frame()
    colnames(dfdat)<-c(destination_aa,"Origin_Mean") 
    return(dfdat)
  }
  
  dfalltest<-wilcox.test(young$ratio_Y,old$ratio_O)
  dfMean=data.frame(des_aa="Origin_Mean", pvalue=dfalltest$p.value)
  dfdat$pvalue<-eachpvalue$pvalue[match(dfdat$des_aa,eachpvalue$dfdes)]
  dfdat<-rbind(dfdat,dfMean)
  dfdat= t(dfdat[,2])%>%as.data.frame()
  
  colnames(dfdat)<-c(destination_aa,"Origin_Mean")
  return(dfdat)
  
})%>%rbind.fill()



thisAA<-destination_aa[2]
dfpvalue_des<-lapply(c(destination_aa,"Mean"), function(thisAA){
  each_destination_aa<-dfqsubs%>%filter(destination==thisAA)
  if(nrow(each_destination_aa)==0){return( data.frame(des_aa=thisAA, pvalue=NA) )}
  young<-each_destination_aa%>%filter(grepl("young",Sample))%>%select(7,15)%>%
    dplyr::rename(ratio_Y=Ratio)
  old<-each_destination_aa%>%filter(grepl("old",Sample))%>%select(7,15)%>%
    dplyr::rename(ratio_O=Ratio)
  dfalltest<-wilcox.test(young$ratio_Y,old$ratio_O)
  dfMean=data.frame(des_aa=thisAA, pvalue=dfalltest$p.value)
  return(dfMean)
})%>%rbind.fill()

dfpvalue_des= t(dfpvalue_des[,2])%>%as.data.frame()
colnames(dfpvalue_des)<-c(destination_aa,"Origin_Mean")
dfpvalue_all<-rbind(dfpvalue,dfpvalue_des)
rownames(dfpvalue_all)<-c(original_aa,"Destination_Mean")

dfpvalue_all_mt<- as.matrix(dfpvalue_all)

getSig <- function(dc) {
  if(is.na(dc)) sc <- ''
  else if (dc <= 0.01) sc <- '***'
  else if (dc <= 0.05) sc <- '**'
  else if (dc > 0.05) sc <- ''
  
}
sig_mat <- matrix( sapply(dfpvalue_all_mt, getSig), nrow=nrow(dfpvalue_all_mt) )


pheatmap(dt,cluster_rows=F,cluster_cols=F,
         display_numbers=sig_mat,
         border="white",
         fontsize_number=12,
         fontsize_col = 8,
         fontsize_row = 8,
         angle_col = 90,
         legend_breaks=c(-5,-2.5,0,2.5,5),
         legend_labels=c("-5","-2.5","0","2.5","5"),
         na_col="white",
         
         fontsize=6
         
         
)















































#########
#dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")%>%filter(!grepl("Inf",Ratio))
dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")
dfqsubs$Base[is.na(dfqsubs$Base)]<-0
dfqsubs<-dfqsubs%>%mutate(Ratio=DP/Base)

young=dfqsubs%>%filter(grepl("young",Sample))%>%select(c(14,15))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_Y=mean(Ratio,na.rm = T))

old=dfqsubs%>%filter(grepl("old",Sample))%>%select(c(14,15))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_O=mean(Ratio,na.rm = T))

dfdat<-merge(young,old)

dfdat$rate_Y[is.infinite(dfdat$rate_Y)]<-median(dfdat$rate_Y)
dfdat$rate_O[is.infinite(dfdat$rate_Y)]<-median(dfdat$rate_O)

wilcox.test(dfdat$rate_Y,dfdat$rate_O,paired = T,alternative = "less")
#n=458,p-value = 8.339e-09


##############
dfqsubs=read.csv("raw/combined_oneSample/txt/qSubs.csv")
dfqsubs$Base[is.na(dfqsubs$Base)]<-0
dfqsubs<-dfqsubs%>%mutate(errorFre=DP/(Base+DP))

young=dfqsubs%>%filter(grepl("young",Sample))%>%select(c(14,17))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_Y=mean(errorFre,na.rm = T))

old=dfqsubs%>%filter(grepl("old",Sample))%>%select(c(14,17))%>%
  group_by(substitution_index)%>%dplyr::summarise(rate_O=mean(errorFre,na.rm = T))

dfdat<-merge(young,old)

wilcox.test(dfdat$rate_Y,dfdat$rate_O,paired = T,alternative = "less")
#p-value = 7.912e-07


