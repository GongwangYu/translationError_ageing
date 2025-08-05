#!/bin/bash
rawdataDir="/mnt/data/home/yugongwang/disk3/translation_error/readthrough/rawdata/nature_yeast_2022"
Index_ncRNA="/home/yugongwang/disk3/translation_error/readthrough/reference_readthrough/yeast/yeast_ncRNA_bt1_index/yeast_NCBI_ncRNA"
Index_mRNA="/home/yugongwang/disk3/translation_error/readthrough/reference_readthrough/yeast/yeast_mRNA_bt1_index/gffread_transcripts_from_v5_transcriptome_13AUG20218_removed_premRNA"
run="4"
adapter="CTGTAGGCACCATCAAT"

mkdir -p align${run}
mkdir -p log${run}


cut -f1,2 id_name_list_WT|while read id name
#cut -f1,2 test|while read id name
do
echo $id	$name

cutadapt -a ${adapter} ${rawdataDir}/${id}.fastq.gz -o align${run}/${name}_trimed.fq -j 80 \
-u 1 \
-q 10 \
--trim-n \
-m 15 \
--max-n=2 \
--discard-untrimmed  \
1> log${run}/${name}.trimed.log 


##mapped to ncRNA
bowtie $Index_ncRNA align${run}/${name}_trimed.fq  align${run}/${name}_MappedTo_ncRNA.sam \
--un align${run}/${name}_unMappedTo_ncRNA.fq \
--norc \
-v 2 \
-p 80  \
2> log${run}/${name}.mappedTo_ncRNA.log 

##mapped to mRNA
bowtie $Index_mRNA align${run}/${name}_unMappedTo_ncRNA.fq align${run}/${name}_MappedTo_mRNA.sam \
--norc \
-v 2 -m 1 -a --best --strata \
-S \
-p 80 \
2> log${run}/${name}.mappedTo_mRNA.log 


samtools view -@ 80 -bF 4 -q 30 align${run}/${name}_MappedTo_mRNA.sam -o align${run}/${name}.bam 


done

