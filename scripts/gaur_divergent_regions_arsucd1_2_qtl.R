#------------------------------------------------------
# Program name: gaur_divergent_regions_arsucd1_2_qtl.R
# Objective: check which region in qtl the divergent
#           region hits
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ape)
library(GenomicRanges)
library(stringr)

##### Divergent region overlaps with ARS-UCD1.2 QTL #####
#read in qtl ARS-UCD1.2
qtl_gff3 <- read.gff("/Users/lloyd/Documents/lloyd_2021/cattle_ARS_UCD1.2_QTL/arsucd1_2_qtl.gff")

#rename chr name
qtl_gff3$seqid <- as.character(qtl_gff3$seqid)

qtl_gff3$seqid <- gsub("Chr.","",qtl_gff3$seqid)

names(qtl_gff3)[1] <- "chr"

#work with complete cases in start and end
compl_cases <- complete.cases(qtl_gff3$start,qtl_gff3$end)
qtl_gff3 <- qtl_gff3[compl_cases,]

#need a strand
qtl_gff3$strand <- rep("*",length(qtl_gff3$strand))
qtl_gff3 <- as_tibble(qtl_gff3)

# divergent region
divergent_region_hereford <- data.frame(chr = "15",start = 77900000, end = 80300000, strand = "*")

# Make GRanges object
divergent_region_hereford_GR <- makeGRangesFromDataFrame(divergent_region_hereford, keep.extra.columns = TRUE,
                                                         seqnames.field="chr", start.field="start",
                                                         end.field="end", strand.field="strand")

#QTL
QTL_GR <- makeGRangesFromDataFrame(qtl_gff3, keep.extra.columns = TRUE,
                                   seqnames.field="chr", start.field="start",
                                   end.field="end",strand.field="strand")

#overlap divergent region with QTL
overlap_divergent_QTL_df <- mergeByOverlaps(divergent_region_hereford_GR,QTL_GR)

overlap_divergent_QTL_as_df <- as.data.frame(overlap_divergent_QTL_df)

# found 2 reproduction related traits in the region between 77900000 and 80300000
# which is when OR4B1 starts and ends with OR5AK30
overlap_divergent_QTL_as_df[1,20]
#[1] "QTL_ID=176438;Name=Conception rate;Abbrev=CONCRATE;PUBMED_ID=31139206;trait_ID=1075;trait=Conception rate;FlankMarker=rs41626653;breed=holstein;VTO_name=fertility trait;Map_Type=Genome;Model=Mendelian;Test_Base=Genome-wise;Significance=Significant;P-value=<0.05;Dominance_Effect=0.76"

overlap_divergent_QTL_as_df[15,20]
#[1] "QTL_ID=30531;Name=Stillbirth;Abbrev=SB;PUBMED_ID=22888914;trait_ID=1153;trait=Stillbirth;FlankMarker=rs109187331,rs43197822,rs109503638;breed=holstein;VTO_name=stillborn offspring quantity;Map_Type=Genome;Model=Mendelian;Test_Base=Genome-wise;peak_cM=110;Significance=Significant;P-value=<0.05"

#the two relevant papers for the region
#stillbirth https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-13-71
#conception rate https://www.frontiersin.org/articles/10.3389/fgene.2019.00412/full 

#there is a Table S1 in the stillbirth paper above that showed all the significant SNPs,
#using what I think is the COMBINED dataset method
#The stillbirth QTL is chr15:78742837-78742841, which is the nearest reproduction related trait to the
#divergent region

#search in cattle QTLdb
#https://www.animalgenome.org/cgi-bin/QTLdb/BT/qdetails?QTL_ID=30531&submit=go 
