library(devtools)
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
library(riboWaltz)
library(readr)
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)
library(pheatmap)

load("nis_yeast.RData")

read_count <- function(dt, annotation, cds_m5, cds_m3, exten_m5, frame_ ) {
  dt <- data.table(dt)
  # Subset data
  message("Subsetting data")
  if (frame_ == "all") {
    message("\tFrame: 0, 1, 2")
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3, ]
    utr3_sub <- dt[psite_region == "3utr", ]
  } else if (frame_ < 0 | frame_ > 2) {
    message("\tInvalid frame number. Please put either 0, 1, 2, or all")
  } else {
    message(paste0("\tFrame: ", frame_))
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3 & frame %in% frame_, ]
    utr3_sub <- dt[psite_region == "3utr" & psite_from_stop > exten_m5 & frame %in% frame_, ]
  }
  utr3_sub <- merge(utr3_sub, annotation[, c("transcript", "nis_pos")], by = "transcript", all.x = TRUE) # Add nis information
  ext_sub <-  utr3_sub[psite_from_stop < nis_pos, ]   # extension
  distal_sub <- utr3_sub[psite_from_stop >= nis_pos, ]  # distal 3'-UTR
  
  
  # Count reads
  message("Counting reads")
  cds_tab <- cds_sub[, list(c_cds = .N), by = list(transcript)]
  utr3_tab <- utr3_sub[, list(c_utr3 = .N), by = list(transcript)]
  ext_tab <- ext_sub[, list(c_ext = .N), by = list(transcript)]
  distal_tab <- distal_sub[, list(c_dist = .N), by = list(transcript)]
  
  # Combine data
  message("Combining read count tables of different regions")
  cu_tab <- Reduce(function(df1, df2) merge(df1, df2, by = "transcript", all.x = TRUE), 
                   list(annotation, cds_tab, utr3_tab, ext_tab, distal_tab))
  message("Replacing NA count with 0")
  cu_tab[is.na(cu_tab)] <- 0
  
  # Calculate length of extension and distal 3'UTR
  cu_tab$l_ext <- pmin(cu_tab$nis_pos-1, cu_tab$l_utr3)
  cu_tab$l_dist <- cu_tab$l_utr3 - cu_tab$l_ext
  cu_tab$nis_pos <- NULL
  
  cu_tab <- data.table(cu_tab %>% mutate_if(is.character, as.factor))
  message("Done\n")
  return(cu_tab)
}

composite <- function(data_list) {
  data_comp <- rbindlist(data_list, use.names = TRUE)[,lapply(.SD, sum),
                                                      by = list(transcript, l_tr, l_utr5, l_cds, l_utr3, l_ext, l_dist)]
  return(data_comp)
}

rt_efficiency <- function(dt, cds_m5, cds_m3,exten_m5) {
  # Calculate RPKM
  lib_size <- (sum(dt$c_cds) + sum(dt$c_utr3))/10^6
  dt$rpkm_cds <- dt$c_cds/(lib_size * ((dt$l_cds - cds_m5 - cds_m3)/10^3))
  dt$rpkm_utr3 <- dt$c_utr3/(lib_size * ((dt$l_utr3)/10^3))
  dt$rpkm_ext <- dt$c_ext/(lib_size * ((dt$l_ext - exten_m5)/10^3))
  dt$rpkm_dist <- dt$c_dist/(lib_size * ((dt$l_dist)/10^3))
  
  # Calculate readthrough efficiency
  dt$rte_utr3 <- dt$rpkm_utr3/dt$rpkm_cds
  dt$rte_ext <- dt$rpkm_ext/dt$rpkm_cds
  
  dt<-dt%>%
    filter(l_ext>30)%>%
    filter(rpkm_cds>5 )%>%
    filter(rpkm_utr3>0.5 )%>%
    filter(rpkm_cds>rpkm_ext)
  return(dt)
}


rt_common<-function(dir){
  
  ###Load bam files
  reads_list <- bamtolist(bamfolder = dir, annotation = annotation_file)
  #生成"transcript" "end5" "end3" "length" "cds_start" "cds_stop"的数据  
  
  ###Selection of read lengths
  filtered_list <- length_filter(data = reads_list,
                                 length_filter_mode = "custom",
                                 length_range = c(20:23,27:32))
  
  ###Determine P-site offset for each read length
  psite_offset <- psite(filtered_list, flanking = 6, extremity = "3end")
  
  #Calculate P-site position for each read
  reads_psite_list <- psite_info(filtered_list, psite_offset)
  
  # Calculate reading frame based on P-site position from start codon
  for (i in 1:length(reads_psite_list)) {
    reads_psite_list[[i]][["frame"]] <- reads_psite_list[[i]][["psite_from_start"]] %% 3
  }
  
  
  #  -------------------------------------------------------------------
  # step1 -----------------------------------------------------------
  ###Count read in the CDS and extension region for each transcript
  ##stop codon是计算到extension区域内的长度
  # cds_m5: the number of nucleotides from the start codon (5' end of CDS) to EXCLUDE from the CDS count
  # cds_m3: the number of nucleotides from the stop codon (3' end of CDS) to EXCLUDE from the CDS count
  # exten_m5 : the number of nucleotides following the stop codon(NTC) to EXCLUDE
  read_count <- function(dt, annotation, cds_m5, cds_m3, exten_m5, frame_ ) {
    dt <- data.table(dt)
    # Subset data
    message("Subsetting data")
    if (frame_ == "all") {
      message("\tFrame: 0, 1, 2")
      cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3, ]
      utr3_sub <- dt[psite_region == "3utr", ]
    } else if (frame_ < 0 | frame_ > 2) {
      message("\tInvalid frame number. Please put either 0, 1, 2, or all")
    } else {
      message(paste0("\tFrame: ", frame_))
      cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3 & frame %in% frame_, ]
      utr3_sub <- dt[psite_region == "3utr" & psite_from_stop > exten_m5 & frame %in% frame_, ]
    }
    utr3_sub <- merge(utr3_sub, annotation[, c("transcript", "nis_pos")], by = "transcript", all.x = TRUE) # Add nis information
    ext_sub <-  utr3_sub[psite_from_stop < nis_pos, ]   # extension
    distal_sub <- utr3_sub[psite_from_stop >= nis_pos, ]  # distal 3'-UTR
    
    
    # Count reads
    message("Counting reads")
    cds_tab <- cds_sub[, list(c_cds = .N), by = list(transcript)]
    utr3_tab <- utr3_sub[, list(c_utr3 = .N), by = list(transcript)]
    ext_tab <- ext_sub[, list(c_ext = .N), by = list(transcript)]
    distal_tab <- distal_sub[, list(c_dist = .N), by = list(transcript)]
    
    # Combine data
    message("Combining read count tables of different regions")
    cu_tab <- Reduce(function(df1, df2) merge(df1, df2, by = "transcript", all.x = TRUE), 
                     list(annotation, cds_tab, utr3_tab, ext_tab, distal_tab))
    message("Replacing NA count with 0")
    cu_tab[is.na(cu_tab)] <- 0
    
    # Calculate length of extension and distal 3'UTR
    cu_tab$l_ext <- pmin(cu_tab$nis_pos-1, cu_tab$l_utr3)
    cu_tab$l_dist <- cu_tab$l_utr3 - cu_tab$l_ext
    cu_tab$nis_pos <- NULL
    
    cu_tab <- data.table(cu_tab %>% mutate_if(is.character, as.factor))
    message("Done\n")
    return(cu_tab)
  }
  
  ## Pool replicates read count
  # data_list: list of data tables that are replicates to combine
  composite <- function(data_list) {
    data_comp <- rbindlist(data_list, use.names = TRUE)[,lapply(.SD, sum),
                                                        by = list(transcript, l_tr, l_utr5, l_cds, l_utr3, l_ext, l_dist)]
    return(data_comp)
  }
  
  ## Calculate readthrough efficiency and RPKM for each region
  # dt: read count data table
  rt_efficiency <- function(dt, cds_m5, cds_m3,exten_m5) {
    # Calculate RPKM
    lib_size <- (sum(dt$c_cds) + sum(dt$c_utr3))/10^6
    dt$rpkm_cds <- dt$c_cds/(lib_size * ((dt$l_cds - cds_m5 - cds_m3)/10^3))
    dt$rpkm_utr3 <- dt$c_utr3/(lib_size * ((dt$l_utr3)/10^3))
    dt$rpkm_ext <- dt$c_ext/(lib_size * ((dt$l_ext - exten_m5)/10^3))
    dt$rpkm_dist <- dt$c_dist/(lib_size * ((dt$l_dist)/10^3))
    
    # Calculate readthrough efficiency
    dt$rte_utr3 <- dt$rpkm_utr3/dt$rpkm_cds
    dt$rte_ext <- dt$rpkm_ext/dt$rpkm_cds
    
    dt<-dt%>%
      filter(l_ext>30)%>%
      filter(rpkm_cds>5 )%>%
      filter(rpkm_utr3>0.5 )%>%
      filter(rpkm_cds>rpkm_ext)
    
    return(dt)
  }
  
  
  # step 3.2 ----------------------------------------------------------
  
  
  readcount<-lapply(reads_psite_list, read_count, annotation = nis, 
                    cds_m5 , cds_m3 , exten_m5, frame_ )
  
}


cds_m5 = 15;cds_m3 = 33;exten_m5=0; frame_ = "0"
readcount<-rt_common(dir = "./align/")
save(readcount,file="readcount.RData")

