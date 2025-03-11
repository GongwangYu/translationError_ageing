


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

combind_species<-data.frame()
# human -------------------------------------------------------------------
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
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("\\.\\d+","",dfqsubs[i,]$protein)
  rel_incr=(old-young)/( (old+young) / 2 )
  return(data.frame(ensg,codon,young,old,rel_incr)) 
  
})%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()

###rscu计算
# library(parazitCUB)
# theHost <- read.host("Homo_sapiens.GRCh38.cds.all.fa.gz", sep = "|")
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
  cutoff_opt=1
  cutoff_nonopt=1
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ))%>%
    filter(!grepl("M|W",AA)) %>%
    filter(prefer!="Unknown")
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()%>%
    filter(!grepl("UGC|UGU",codon))
  AA<-dfqsubs_flt_pefer$AA%>%unique()
  
  #方法1
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-median(eachAA$rel_incr)
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.203188,p-value = 0.1618
  }
  #方法2
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-0
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.379872,p-value = 0.0316
  }
  
  #方法3
  { 
    cutoff=0.25
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      bottom<-eachAA$rel_incr%>%quantile(.,probs=cutoff)
      top<-eachAA$rel_incr%>%quantile(.,probs=1-cutoff)
      
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=top & prefer=="nonoptimal" ~ "a",
        rel_incr >=top & prefer=="optimal" ~ "b",
        rel_incr <=bottom & prefer=="nonoptimal" ~ "c",
        rel_incr <=bottom & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    mydata.frame <- dfMH
    #mydata.frame <- mydata.frame+1
    my.matrix <- as.matrix(mydata.frame)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.502112,p-value = 0.03452
   
    df_1000_sd<- mclapply(1:1000, function(x){
      sample_data<- sample_n(mydata.frame,nrow(mydata.frame),replace = T);
      my.matrix <- as.matrix(sample_data)
      myvector <- as.vector(t(my.matrix))
      dims=c(2,2,nrow(mydata.frame));
      myarray <- array(myvector,dims)
      mytest<- mantelhaen.test(myarray)
      return(data.frame(commonOR=mytest$estimate,pvalue=mytest$p.value))
    },mc.cores = 40)%>% rbind.fill()
    my.se<-sd(df_1000_sd$commonOR)
    #my.se=0.4771917
    
    combind_species<-rbind(dfMH,combind_species)
  }
  
}  


# mouse -------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_mouse.csv")

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
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("\\.\\d+","",dfqsubs[i,]$protein)
  rel_incr=(old-young)/((old+young)/2)
  return(data.frame(ensg,codon,young,old,rel_incr)) 
  
  
})%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()

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
  cutoff_opt=1
  cutoff_nonopt=1
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ) )%>%
    filter( !grepl("M|W",AA) ) %>%
    filter( prefer != "Unknown" )
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()%>%
    filter( !grepl("GCC|GAG",codon) )
    #filter( !grepl("GCC",codon) )
  AA<-dfqsubs_flt_pefer$AA%>%unique()
  
  #方法1
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-median(eachAA$rel_incr)
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.161985,p-value = 0.3267
  }
  #方法2
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-0
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.060998,p-value = 0.464
  }
  
  #方法3
  { 
    cutoff=0.25
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      bottom<-eachAA$rel_incr%>%quantile(.,probs=cutoff)
      top<-eachAA$rel_incr%>%quantile(.,probs=(1-cutoff))
      
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=top & prefer=="nonoptimal" ~ "a",
        rel_incr >=top & prefer=="optimal" ~ "b",
        rel_incr <=bottom & prefer=="nonoptimal" ~ "c",
        rel_incr <=bottom & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    mydata.frame <- dfMH
    #mydata.frame <- mydata.frame+1
    my.matrix <- as.matrix(mydata.frame)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.253497,p-value = 0.2966
    
    df_1000_sd<- mclapply(1:1000, function(x){
      sample_data<- sample_n(mydata.frame,nrow(mydata.frame),replace = T);
      my.matrix <- as.matrix(sample_data)
      myvector <- as.vector(t(my.matrix))
      dims=c(2,2,nrow(mydata.frame));
      myarray <- array(myvector,dims)
      mytest<- mantelhaen.test(myarray)
      return(data.frame(commonOR=mytest$estimate))
    },mc.cores = 40)%>% rbind.fill()
    my.se<-sd(df_1000_sd$commonOR)
    #my.se=0.2202951
    
    combind_species<-rbind(dfMH,combind_species)
  }
}

# yeast -------------------------------------------------------------------
dfqsubs<-read.csv("detct_quant_yeast.csv")

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
  
  
  codon<-gsub("T","U",dfqsubs[i,]$codon) 
  ensg= gsub("_mRNA","",dfqsubs[i,]$protein)
  rel_incr= (old-young) /( (old+young) /2 )
  return(data.frame(ensg,codon,young,old,rel_incr)) 
  
})%>%rbind.fill()%>%
  filter(young >0 | old >0)%>%na.omit()


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
  cutoff_opt=1
  cutoff_nonopt=1
  rscu_pref<-rscu_flt%>%dplyr::mutate(prefer=case_when(
    rscu>cutoff_opt ~"optimal",
    rscu<cutoff_nonopt ~"nonoptimal",
    TRUE ~ "Unknown"
  ))%>%
    filter(!grepl("M|W",AA)) %>%
    filter(prefer!="Unknown")
  
  dfqsubs_flt_pefer<-dfqsubs_flt%>%left_join(rscu_pref,by="codon") %>%na.omit()%>%
    filter(!grepl("GCU",codon))
  AA<-dfqsubs_flt_pefer$AA%>%unique()
  
  #方法1
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-median(eachAA$rel_incr)
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.291228,p-value = 0.3281
  }
  #方法2
  {
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      mdn<-0
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=mdn & prefer=="nonoptimal" ~ "a",
        rel_incr >=mdn & prefer=="optimal" ~ "b",
        rel_incr <mdn & prefer=="nonoptimal" ~ "c",
        rel_incr <mdn & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    my.matrix <- as.matrix(dfMH)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.353319,p-value = 0.2869
  }
  
  #方法3
  { 
    cutoff=0.25
    dfMH<-mclapply(AA,function(thisAA){
      eachAA<-dfqsubs_flt_pefer%>%filter(AA==thisAA)
      bottom<-eachAA$rel_incr%>%quantile(.,probs=cutoff)
      top<-eachAA$rel_incr%>%quantile(.,probs=1-cutoff)
      
      dftable<-eachAA%>%dplyr::mutate(ctgy=case_when(
        rel_incr >=top & prefer=="nonoptimal" ~ "a",
        rel_incr >=top & prefer=="optimal" ~ "b",
        rel_incr <=bottom & prefer=="nonoptimal" ~ "c",
        rel_incr <=bottom & prefer=="optimal" ~ "d"
      ))
      data<-matrix( table(dftable$ctgy), ncol = 4, byrow = TRUE)%>%as.data.frame()
      colnames(data) <- c("a","b","c","d")
      return(data)
      
    },mc.cores = 6)%>%rbind.fill() 
    mydata.frame <- dfMH
    #mydata.frame <- mydata.frame+1
    my.matrix <- as.matrix(mydata.frame)
    myvector <- as.vector(t(my.matrix))
    dims=c(2,2,nrow(dfMH));
    myarray <- array(myvector,dims)
    mantelhaen.test(myarray,alternative = "g") 
    #or=1.299639,p-value = 0.3267
    
    df_1000_sd<- mclapply(1:1000, function(x){
      sample_data<- sample_n(mydata.frame,nrow(mydata.frame),replace = T);
      my.matrix <- as.matrix(sample_data)
      myvector <- as.vector(t(my.matrix))
      dims=c(2,2,nrow(mydata.frame));
      myarray <- array(myvector,dims)
      mytest<- mantelhaen.test(myarray)
      return(data.frame(commonOR=mytest$estimate))
    },mc.cores = 40)%>% rbind.fill()
    my.se<-sd(df_1000_sd$commonOR)
    #my.se=0.3954973
    
    combind_species<-rbind(dfMH,combind_species)
  }
  
}  


# 三个物种合并 ------------------------------------------------------------------


mydata.frame <- combind_species
my.matrix <- as.matrix(mydata.frame)
myvector <- as.vector(t(my.matrix))
dims=c(2,2,nrow(mydata.frame));
myarray <- array(myvector,dims)
mantelhaen.test(myarray,alternative = "g") 
#or=1.39737,p=0.02481

df_1000_sd<- mclapply(1:1000, function(x){
  sample_data<- sample_n(mydata.frame,nrow(mydata.frame),replace = T);
  my.matrix <- as.matrix(sample_data)
  myvector <- as.vector(t(my.matrix))
  dims=c(2,2,nrow(mydata.frame));
  myarray <- array(myvector,dims)
  mytest<- mantelhaen.test(myarray)
  return(data.frame(commonOR=mytest$estimate))
},mc.cores = 40)%>% rbind.fill()
my.se<-sd(df_1000_sd$commonOR)
#my.se=0.2696257

# 绘图 ----------------------------------------------------------------------
OR_data<-data.frame(dftype=c("Human","Mouse","Yeast","ALL"),
                    or=c(1.502112,1.253497,1.299639,1.39737),
                   se=c(0.4771917,0.2202951,0.3954973,0.2696257),
                   pvalue=c(0.03452,0.2966,0.3267,0.02481))

OR_data$dftype<-factor(OR_data$dftype,levels =c("Human","Mouse","Yeast","ALL") )
fig6g<-ggplot(OR_data,aes(x=dftype,y=or))+
  geom_bar(stat = "identity",fill = "white", color = "black") +
  geom_errorbar(aes(ymax=or+se,ymin=or-se,width=0.2))+
  coord_flip(ylim = c(0.8,2.3))+
  labs(y="Combined odds ratio",x="Species")+
  geom_hline(aes(yintercept=1),lty=2)+
  theme_test()+
  annotate("text",x=1,y=2.2,label=expression(italic(P)* " = " * 0.03),size=3)+
  annotate("text",x=2,y=2.2,label=expression( italic(P)* " = " * 0.29),size=3)+
  annotate("text",x=3,y=2.2,label=expression( italic(P)* " = " * 0.32),size=3)+
  annotate("text",x=4,y=2.2,label=expression( italic(P)* " = " * 0.02),size=3)+
  scale_y_continuous(breaks = c(1.0,1.5,2,2.5))

fig6g

save.image("fig6g.RData")
load("fig6g.RData")
