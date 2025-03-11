

# common ------------------------------------------------------------------
{
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
  amino_acids_64 = 'FFLLSSSSYY--CC-WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' %>%
    str_split("")%>%unlist()
  dfcodon_AA_coord<-data.frame(codon=AllCodons,AA=amino_acids_64)  
}


# human -------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_human.csv")

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
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("\\.\\d+","",dfqsubs[i,]$protein)
  rel_incr=(old-young)/( (old+young) / 2 )
  return(data.frame(ensg,codon,young,old,rel_incr)) 

},mc.cores = 6)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()


###表达量
#椎间盘组织中基因的表达量
human_enst_ensg<-read_tsv("human_biomart_anno.txt") %>%dplyr::select(1,3)
colnames(human_enst_ensg)<-c("ensg","name")

human_desc_expr<-read_tsv("GSE245147_Degenerated_NO_Degenerated_RPKM.txt")%>%
  rowwise()%>%dplyr::mutate(expr=log2(mean(c(NO_Degenerated_1,NO_Degenerated_2,NO_Degenerated_3))+1))%>%
  dplyr::select(name=Geneid,len=Length,expr) %>%left_join(human_enst_ensg) %>%
  dplyr::select(ensg,expr)%>%unique()%>%na.omit()

dfmrg<-left_join(dfqsubs_flt,human_desc_expr,by="ensg")%>%na.omit() #N=436

lmRes <- lm(rel_incr ~ codon * expr, data=dfmrg)
summary(lmRes)

###rscu计算
# library(parazitCUB)
 theHost <- read.host("Homo_sapiens.GRCh38.cds.all.fa.gz", sep = "|")
# rscu <- RSCU.values(theHost)
# save(rscu,file = "rscu_human.RData")
load("rscu_human.RData")
rscu_flt<-rscu$Homo_sapiens.GRCh38.cds.all.fa.gz %>%rownames_to_column(var="gene")%>%
  filter(grepl("transcript_biotype:protein_coding",gene)) %>%
  dplyr::mutate(gene= str_match(gene, "ENSG[0-9]+"))%>%
  filter((atg==1 & taa %in% c(0,3) & tag %in% c(0,3) & tga %in% c(0,3)))%>%
  dplyr::select(-1) %>%colMeans()%>%
  as.data.frame()%>%rownames_to_column(var="codon") %>%dplyr::mutate(codon=toupper(codon))%>%
  dplyr::mutate(codon=gsub("T","U",codon))%>%merge(dfcodon_AA_coord,by="codon")
colnames(rscu_flt)[2]<-"rscu"
###统计分析
{
  # rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=ifelse(rscu>1,"optimal","nonoptimal"))%>%
  #   filter(!grepl("M|W",AA))
  cutoff_opt=1.2
  cutoff_nonopt=0.8
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ))%>%filter(!grepl("M|W",AA)) %>%filter(prefer!="Unknown")
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()
  
  opt<-dfqsubs_flt_pefer%>%filter(prefer=="optimal")
  nonopt<-dfqsubs_flt_pefer%>%filter(prefer=="nonoptimal")
  wilcox.test(opt$rel_incr,nonopt$rel_incr,alternative = "l")
  #p-value = 0.004625
  mean(opt$rel_incr)
  #0.2282027
  mean(nonopt$rel_incr)
  #0.7545854
  ###OR分析
  table(opt$rel_incr>0)
  table(nonopt$rel_incr>0)
  data <- matrix(c(table(opt$rel_incr>0), table(nonopt$rel_incr>0)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("decrease", "increase")
  rownames(data) <- c("opt", "nonopt")
  fisher.test(data, alternative = "g")
  #OR=1.85203,p-value = 0.01933
  
  bottom25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.25)
  top25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.75)
  bottom<-dfqsubs_flt_pefer%>%filter(rel_incr<= bottom25)
  table(bottom$prefer)
  top<-dfqsubs_flt_pefer%>%filter(rel_incr>= top25)
  table(top$prefer)
  data <- matrix(c(table(top$prefer),table(bottom$prefer)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt", "opt")
  rownames(data) <- c( "higherError","lowerError")
  fisher.test(data, alternative = "g")
  #OR=2.117665,p-value = 0.02292
  
  mdn<-mean(dfqsubs_flt_pefer$rel_incr)
  dftable<-dfqsubs_flt_pefer%>%dplyr::mutate(ctgy=case_when(
    rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
    rel_incr >=mdn & prefer=="optimal" ~ "b",
    rel_incr <mdn & prefer=="nonoptimal" ~ "c",
    rel_incr <mdn & prefer=="optimal" ~ "d"
  ))
  data<-matrix( table(dftable$ctgy), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt","opt")
  rownames(data) <- c("higherError", "lowerError") 
  fisher.test(data, alternative = "g")
  #1.920149 ,p-value = 0.01171  
}

# mouse -------------------------------------------------------------------
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
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("\\.\\d+","",dfqsubs[i,]$protein)
  rel_incr=(old-young)/((old+young)/2)
  return(data.frame(ensg,codon,young,old,rel_incr)) 
  
  
},mc.cores = 6)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()


###表达量
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
    dplyr::mutate(expr_young=log2(mean(young_wholelung_1,young_wholelung_2,young_wholelung_3)+1),
                  expr_old=log2(mean(old_wholelung_1,old_wholelung_2,old_wholelung_3)+1))%>%
    dplyr::select(1,8,9)
  
  }

dfmrg<-left_join(dfqsubs_flt,mouse_lung_epxr,by="ensg")%>%na.omit() #N=436

lmRes <- lm(rel_incr ~ codon * expr_young, data=dfmrg)
summary(lmRes)

###rscu计算
# library(parazitCUB)
# theHost <- read.host("Mus_musculus.GRCm39.cds.all.fa.gz", sep = "|")
# rscu <- RSCU.values(theHost)
# save(rscu,file = "rscu_mouse.RData")
load("rscu_mouse.RData")
rscu_flt<-rscu$Mus_musculus.GRCm39.cds.all.fa.gz %>%rownames_to_column(var="gene")%>%
  filter(grepl("transcript_biotype:protein_coding",gene)) %>%
  dplyr::mutate(gene= str_match(gene, "ENSMUSG[0-9]+"))%>%
  filter((atg==1 & taa %in% c(0,3) & tag %in% c(0,3) & tga %in% c(0,3)))%>%
  dplyr::select(-1) %>%colMeans()%>%
  as.data.frame()%>%rownames_to_column(var="codon") %>%dplyr::mutate(codon=toupper(codon))%>%
  dplyr::mutate(codon=gsub("T","U",codon))%>%merge(dfcodon_AA_coord,by="codon")
colnames(rscu_flt)[2]<-"rscu"

###统计分析
{
  # rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=ifelse(rscu>1,"optimal","nonoptimal"))%>%
  #   filter(!grepl("M|W",AA))
  cutoff_opt=1
  cutoff_nonopt=1
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ))%>%filter(!grepl("M|W",AA)) %>%filter(prefer!="Unknown")
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()
  
  opt<-dfqsubs_flt_pefer%>%filter(prefer=="optimal")
  nonopt<-dfqsubs_flt_pefer%>%filter(prefer=="nonoptimal")
  wilcox.test(opt$rel_incr,nonopt$rel_incr,alternative = "l")
  #p-value = 0.89
  mean(opt$rel_incr)
  #0.3138528
  mean(nonopt$rel_incr)
  #-0.04557991
  ###OR分析
  table(opt$rel_incr>0)
  table(nonopt$rel_incr>0)
  data <- matrix(c(table(opt$rel_incr>0), table(nonopt$rel_incr>0)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("decrease", "increase")
  rownames(data) <- c("opt", "nonopt")
  fisher.test(data, alternative = "g")
  #OR=0.5917908,p-value = 0.9393
  
  bottom25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.25)
  top25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.75)
  bottom<-dfqsubs_flt_pefer%>%filter(rel_incr<= bottom25)
  table(bottom$prefer)
  top<-dfqsubs_flt_pefer%>%filter(rel_incr>= top25)
  table(top$prefer)
  data <- matrix(c(table(top$prefer),table(bottom$prefer)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt", "opt")
  rownames(data) <- c( "higherError","lowerError")
  fisher.test(data, alternative = "g")
  #OR=0.5607127,p-value = 0.9116
  
  mdn<-median(dfqsubs_flt_pefer$rel_incr)
  dftable<-dfqsubs_flt_pefer%>%dplyr::mutate(ctgy=case_when(
    rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
    rel_incr >=mdn & prefer=="optimal" ~ "b",
    rel_incr <mdn & prefer=="nonoptimal" ~ "c",
    rel_incr <mdn & prefer=="optimal" ~ "d"
  ))
  data<-matrix( table(dftable$ctgy), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt","opt")
  rownames(data) <- c("higherError", "lowerError") 
  fisher.test(data, alternative = "g")
  #0.6813633 ,p-value = 0.8828  
}

# yeast -------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_yeast.csv")

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
  
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("_mRNA","",dfqsubs[i,]$protein)
  rel_incr= (old-young) /( (old+young) /2 )
  return(data.frame(ensg,codon,young,old,rel_incr)) 

},mc.cores = 6)%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()


###表达量
load("01.yeast.RData")
expr<-  data.frame(ensg=toPlot.yeast$orf,
                   expr=toPlot.yeast$logexpr)

dfmrg<-left_join(dfqsubs_flt,expr,by="ensg")%>%na.omit() #N=436

lmRes <- lm(rel_incr ~ codon * expr, data=dfmrg)
summary(lmRes)

###rscu计算
# library(parazitCUB)
# theHost <- read.host("Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz", sep = "|")
# rscu <- RSCU.values(theHost)
# save(rscu,file = "rscu_yeast.RData")
load("rscu_yeast.RData")
rscu_flt<-rscu$`Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz` %>%rownames_to_column(var="gene")%>%
  filter(grepl("transcript_biotype:protein_coding",gene)) %>%
  dplyr::mutate(gene= str_match(gene, "^Y.+_mRNA"))%>%na.omit()%>%
  filter((atg==1 & taa %in% c(0,3) & tag %in% c(0,3) & tga %in% c(0,3)))%>%
  dplyr::mutate(gene=gsub("_mRNA","",gene))%>%dplyr::select(-1) %>%colMeans()%>%
  as.data.frame()%>%rownames_to_column(var="codon") %>%dplyr::mutate(codon=toupper(codon))%>%
  dplyr::mutate(codon=gsub("T","U",codon))%>%merge(dfcodon_AA_coord,by="codon")
colnames(rscu_flt)[2]<-"rscu"

###统计分析
{
  # rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=ifelse(rscu>1,"optimal","nonoptimal"))%>%
  #   filter(!grepl("M|W",AA))
  cutoff_opt=1.2
  cutoff_nonopt=0.8
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ))%>%filter(!grepl("M|W",AA)) %>%filter(prefer!="Unknown")
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()
  
  opt<-dfqsubs_flt_pefer%>%filter(prefer=="optimal")
  nonopt<-dfqsubs_flt_pefer%>%filter(prefer=="nonoptimal")
  wilcox.test(opt$rel_incr,nonopt$rel_incr,alternative = "l")
  #p-value = 0.7274
  mean(opt$rel_incr)
  #1.007232
  mean(nonopt$rel_incr)
  #0.5639471
  ###OR分析
  table(opt$rel_incr>0)
  table(nonopt$rel_incr>0)
  data <- matrix(c(table(opt$rel_incr>0), table(nonopt$rel_incr>0)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("decrease", "increase")
  rownames(data) <- c("opt", "nonopt")
  fisher.test(data, alternative = "g")
  #OR=1.048361,p-value = 0.7317
  
  bottom25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.25)
  top25<-dfqsubs_flt_pefer$rel_incr%>%quantile(.,probs=0.75)
  bottom<-dfqsubs_flt_pefer%>%filter(rel_incr<= bottom25)
  table(bottom$prefer)
  top<-dfqsubs_flt_pefer%>%filter(rel_incr>= top25)
  table(top$prefer)
  data <- matrix(c(table(top$prefer),table(bottom$prefer)), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt", "opt")
  rownames(data) <- c( "higherError","lowerError")
  fisher.test(data, alternative = "g")
  #OR=0.8296783,p-value = 0.799
  
  mdn<-mean(dfqsubs_flt_pefer$rel_incr)
  dftable<-dfqsubs_flt_pefer%>%dplyr::mutate(ctgy=case_when(
    rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
    rel_incr >=mdn & prefer=="optimal" ~ "b",
    rel_incr <mdn & prefer=="nonoptimal" ~ "c",
    rel_incr <mdn & prefer=="optimal" ~ "d"
  ))
  data<-matrix( table(dftable$ctgy), ncol = 2, byrow = TRUE)
  colnames(data) <- c("nonopt","opt")
  rownames(data) <- c("higherError", "lowerError") 
  fisher.test(data, alternative = "g")
  #0.5988078 ,p-value = 0.8514  
}

