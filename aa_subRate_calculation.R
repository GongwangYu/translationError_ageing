library(dplyr)
library(plyr)
library(readr)
library(parallel)
library(tidyr)
library(ggplot2)
library(data.table)


dfqsubs<-read.csv("detct_quant.csv")

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


save(dfqsubs_flt,file = "dfqsubs_flt.RData")

