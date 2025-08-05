

library(stringr)

dfqsubs<-read.csv("txt_AllWindow/txt/detct_quant.csv")

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
dfcodon_AA_coord<-data.frame(codon=AllCodons,AA=amino_acids_64)%>%
  mutate(AA_poll=ifelse((AA=="I" | AA=="L") ,"I/L",AA))%>%
  mutate(AAcor=paste(codon,AA_poll,sep = " to "))

allSubsType<-c()
for (a in AllCodons){
  a=a
  for (b in amino_acids) {
    b=b
    df=paste(a,b,sep = " to ")
    allSubsType<-c(allSubsType,df)
    
  }
}

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

dfmtr_onecodonToOneAA<-mclapply(allSubsType, function(thistype){
  eachtype<-dfsubs_compile%>%filter(substype==thistype)
  codon=strsplit(thistype," to ")[[1]][1]
  destination=strsplit(thistype," to ")[[1]][2]
  if(nrow(eachtype)==0){return(data.frame(codon,destination,effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  codon = eachtype$codon%>%unique()
  destination = eachtype$destination%>%unique()
  return(data.frame(codon,destination,effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill() 


dfmtr_allCodnToOneAA<-mclapply(amino_acids, function(thistype){
  eachtype<-dfsubs_compile%>%filter(destination==thistype)
  if(nrow(eachtype)==0){return(data.frame(codon="allCodon",destination=thistype,effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  return(data.frame(codon="allCodon", destination=thistype,effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill()


dfmtr_oneCodonToallAA<-mclapply(AllCodons, function(thistype){
  eachtype<-dfsubs_compile%>%filter(codon==thistype)
  if(nrow(eachtype)==0){return(data.frame(codon=thistype,destination="allAA",effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  return(data.frame(codon=thistype,destination="allAA",effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill()


dfallMtr<-rbind(dfmtr_onecodonToOneAA,dfmtr_allCodnToOneAA)%>%rbind(dfmtr_oneCodonToallAA)%>%
  mutate(effeSize=ifelse(is.na(effeSize),0,effeSize))%>%
  mutate(substype=paste(codon,destination,sep = " to "))%>%
  mutate(nobio=ifelse((substype %in% dfcodon_AA_coord$AAcor)|(codon %in%c("UAA","UAG","UGA")) ,T,F))%>%
  mutate(effeSize=ifelse((substype %in% dfcodon_AA_coord$AAcor)|(codon %in%c("UAA","UAG","UGA")) ,NA,effeSize))%>%
  mutate(stars=ifelse(pvalue<0.05,"*",""))%>%
  filter(!(codon %in%c("UAA","UAG","UGA") ))


save(dfallMtr,file = "../RData.file/misin.human_disc.heatmap.RData")

ggplot(dfallMtr,aes(x=destination,y=codon,fill=effeSize))+
  geom_tile(aes(fill=effeSize))+
  geom_text(aes(label=stars),color="black",size=4,
            hjust="middle",vjust="bootom")+
  scale_fill_gradient2("Aged-to-young relative increase\nin translation error rates",
                       na.value="grey",low = "#4DBBD5FF",high = "#E64B35FF",mid = "white",
                       guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                            barwidth =8,barheight = 1, ticks.colour="black",ticks.linewidth = 1 ))+
  labs(x="Misincorporation of amino acids",y="Original codon")+
  scale_y_discrete(breaks=dfallMtr$codon%>%unique()%>%sort(),labels=c("",(dfallMtr$codon%>%unique()%>%sort())[-1]))+               
  theme(
    axis.ticks.y = element_blank(),
    panel.background=element_blank())+
  theme_test()+
  theme(legend.position = "top")
