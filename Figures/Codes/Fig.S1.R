library(ggtext)
library(ggplot2)
library(data.table)

# Fig.S1 ------------------------------------------------------------------
load("plot.RData")
OR_codon_plot$Species<-c("Human","Yeast","Mouse")
# plot
Fig.S1A <- ggplot(OR_codon_plot, aes(x = OR - 1, y = Species)) +  # 将 OR 平移使 1 为起点
  geom_bar(stat = "identity", fill = "white", color = "grey35", 
           width = 0.6, position = position_dodge(0.8),linewidth =1) +
  geom_errorbar(aes(xmin = CI_low-1, xmax = CI_high-1), width = 0.2, color = "black",
                position = position_dodge(0.8), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # 将虚线保持在0的位置
  scale_x_continuous(
    limits = c(-0.6, 2.45),  # 
    breaks = seq(-0.4, 1.4, 0.4),  
    labels = c("0.6", "1.0", "1.4","1.8","2.2")  
  ) +
  scale_y_discrete(position = "left") +  # 
  labs(title = "", x = "Combined odds ratio", y = "Species") +
  #theme_bw(base_size = 16) +
  theme_classic() +
  theme(
    axis.text.y = element_markdown()
  )+
  annotate("text",x=2,y=3,label=expression(italic(P) * " = 5x" *10**-239),size=3)+
  annotate("text",x=2,y=2,label=expression(" N.S. "),size=3)+
  annotate("text",x=2,y=1,label=expression(italic(P) * " = 7x" * 10**-220),size=3)
Fig.S1A

# Fig.S1B -----------------------------------------------------------------
library(stringr)
library(tidyr)

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

dfdat<-data.frame(substype=allSubsType) %>%
  filter(!grepl("UAA|UAG|UGA",substype))%>%
  filter(!(substype %in% (dfcodon_AA_coord$AAcor%>%unique()) )) %>%
  separate(col = substype,into = c("codon","dest"),sep=" to ")

dfcognate<-lapply(1:nrow(dfdat), function(i){
  eachrow<-dfdat[i,]
  codon<-eachrow$codon
  dest<-eachrow$dest
  dest_codon<-dfcodon_AA_coord%>%filter(AA==dest) %>%getElement("codon")
  
  {
    check_match <- function(codon, dest_codon) {
     
      for (target in dest_codon) {
       
        match_count <- 0
        
        
        for (i in 1:3) {
          if (substr(codon, i, i) == substr(target, i, i)) {
            match_count <- match_count + 1
          }
        }
        
        # 如果有 2 个字符匹配，返回 1
        if (match_count >= 2) {
          return(1)
        }
      }
      
      
      return(0)
    }
    
   
    cognate <- check_match(codon, dest_codon)
  }
  
  eachrow$cognate<-cognate
  return(eachrow)
})%>%rbind.fill()
dfcognate$cognate%>%table() %>%as.numeric()

load("dfmisincorpration_info.RData")

dfhuman<-dfmisincorpration_info%>%filter(species=="human")

dfhuman$cognate%>%table() %>%as.numeric()

dffish<-matrix(c(dfcognate$cognate%>%table() %>%as.numeric(),
                 dfhuman$cognate%>%table() %>%as.numeric()),nrow = 2)
colnames(dffish)<-c("expected","observed")
rownames(dffish)<-c("NoCE","NeCE")
a=fisher.test(dffish)
a$p.value
#p=4.53856e-25

dfmouse<-dfmisincorpration_info%>%filter(species=="mouse")
dfmouse$cognate%>%table() %>%as.numeric()
dffish<-matrix(c(dfcognate$cognate%>%table() %>%as.numeric(),
                 dfmouse$cognate%>%table() %>%as.numeric()),nrow = 2)
colnames(dffish)<-c("expected","observed")
rownames(dffish)<-c("NoCE","NeCE")
a=fisher.test(dffish)
a$p.value
#p=1.2011e-16

dfyeast<-dfmisincorpration_info%>%filter(species=="yeast")
dfyeast$cognate%>%table() %>%as.numeric()
dffish<-matrix(c(dfcognate$cognate%>%table() %>%as.numeric(),
                 dfyeast$cognate%>%table() %>%as.numeric()),nrow = 2)
colnames(dffish)<-c("expected","observed")
rownames(dffish)<-c("NoCE","NeCE")
a=fisher.test(dffish)
a$p.value
#p=2.327526e-09


data <- data.frame(
  species = rep(c("Expected","Human", "Mouse", "Yeast"), each = 2),
  category = rep(c("NoCE", "NeCE"), times = 4),
  value = c(788, 310, 416, 437, 119, 149, 9, 29)
)
# 
data$total <- ave(data$value, data$species, FUN = sum)

# 
data$percent <- data$value / data$total * 100
# 
Fig.S1B<-ggplot(data, aes(x = species, y = percent, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_y_continuous(labels = scales::percent) +  # 将y轴标签转换为百分比格式
  labs(y = "Fraction of detected misincorporations") +
  theme_classic() +
  theme(legend.position = "top") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), 
            position = position_stack(vjust = 0.5),  # 
            size = 3.5, color = "white")+
  scale_y_continuous(breaks = c(0,25,50,75,100),labels =c(0,25,50,75,100) )+
  theme(axis.title.x = element_blank(),legend.title = element_blank())+
  annotate("segment", x = 1, xend = 2, y = 102, yend = 102)+
  annotate("text",x=1.5,y=107,label=expression(italic(P) * " = " * "4×" *10**-25 ),size=3)+
  annotate("segment", x = 1, xend = 3, y = 112, yend = 112)+
  annotate("text",x=2,y=117,label=expression(italic(P) * " = " * "1×" *10**-16 ),size=3)+
  annotate("segment", x = 1, xend = 4, y = 122, yend = 122)+
  annotate("text",x=2.5,y=127,label=expression(italic(P) * " = " * "2×" *10**-9 ),size=3)
Fig.S1B


#  ----------------------------------------------------------------------
library(cowplot)
Fig.S1 <- plot_grid(Fig.S1A,Fig.S1B,ncol = 2, labels = c("a","b"))  

ggsave("Fig.S1.pdf",Fig.S1,width =8, height =4)
