#------------------------------------------------------
# Program name: gaur_annotation_gtf_flipping_coordinates.R
# Objective: create reverse comp gtf coordinates for ribbon plot
#           
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

# library(ggplot2)
library(refGenome)
library(dplyr)
library(readr)

# names(brahmanGTF) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
#using refGenome to read in genome gtf

##### Brahman #####
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190320/")

#Brahman EBI first draft v1.96
read.gtf(ens,"Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.96.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Brahman_gene <- getGenePositions(ens)
dim(Brahman_gene)
# 

#get chr15 and rev comp annotation
#length of chr15 84270218
#region of interest is 3748952 to 5140465
left <- 3748952
right <- 5140465
length <- 5140465 - 3748952 + 1
selectedRegionName <- "UOA_Brahman_1_chr15_3748952_5140465_rc"

Brahman_gene_chr15 <- Brahman_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

#rev comp
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = length - newstart +1) %>% 
  mutate(newend = length - newend +1)
Brahman_gene_chr15$newstrand <- ifelse(Brahman_gene_chr15$strand == "+","-","+")

Brahman_gene_chr15$selectedRegionName <- selectedRegionName

newcoorDF <- Brahman_gene_chr15 %>% 
  select(selectedRegionName,newend,newstart,gene_id,score,newstrand,gene_biotype)

#write out the results as bed for overlap in ribbon
# write_tsv(newcoorDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/UOA_Brahman_1_chr15_3748952_5140465_rc.bed",
#           col_names = FALSE)

##### Gaur #####
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2019/Research/gaur/EBI_annotation/gaur/")

#Bos_gaurus.ARS_UOA_Gaur_1.101.gtf
read.gtf(ens,"Bos_gaurus.ARS_UOA_Gaur_1.101.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Gaur_gene <- getGenePositions(ens)
dim(Gaur_gene)

#get chr15 and rev comp annotation
#length of chr15 91621889
#region of interest is 45417056 to 46744541
left <- 45417056
right <- 46744541
length <- 46744541 - 45417056 + 1
selectedRegionName <- "gaur_15_44877349_46204834"

Gaur_gene_chr15 <- Gaur_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Gaur_gene_chr15 <- Gaur_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

#rev comp
Gaur_gene_chr15 <- Gaur_gene_chr15 %>% mutate(newstart = length - newstart +1) %>% 
  mutate(newend = length - newend +1)
Gaur_gene_chr15$newstrand <- ifelse(Gaur_gene_chr15$strand == "+","-","+")

Gaur_gene_chr15$selectedRegionName <- selectedRegionName

newcoorDF <- Gaur_gene_chr15 %>% 
  dplyr::select(selectedRegionName,newend,newstart,gene_id,score,newstrand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorDF,"/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/ribbon_plots/UOA_Gaur_1_chr15_46744541_45417056_rc.bed",
          col_names = FALSE)

# read the gap coor file and insert coordinates into above annotation dataframe
gaur_gap <- read_tsv(file = "/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/ribbon_plots/nucmer_divergent_chr15/gaur_revcomp_15_44877349_46204834.coor",
                     col_names = FALSE)

colnames(gaur_gap) <- c("selectedRegionName","newstart","newend")

gaur_gap$gene_id <- "gap"
gaur_gap$score <- "."
gaur_gap$newstrand <- "+"
gaur_gap$gene_biotype <- "gap"

write_tsv(gaur_gap,"/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/ribbon_plots/UOA_Gaur_1_chr15_46744541_45417056_rc_gapOnly.bed",
          col_names = FALSE)

# #Angus
# # create ensemblGenome object for storing Ensembl genomic annotation data
# ens <- ensemblGenome()
# 
# setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Paternal_20190320/")
# 
# #Brahman EBI first draft v1.96
# read.gtf(ens,"Bos_indicus_x_bos_taurus_Pat.UOA_Angus_1.96.gtf")
# 
# # counts all annotations on each seqname
# tableSeqids(ens)
# 
# # # create table of genes
# Angus_gene <- getGenePositions(ens)
# dim(Angus_gene)
# 
# #get chr15 annotation
# #length of chr15 84012196
# #region of interest is 78799177 to 80168904
# left <- 78799177
# right <- 80168904
# length <- 80168904 - 78799177 + 1
# selectedRegionName <- "UOA_Angus_1_chr15_78799177_80168904"
# 
# Angus_gene_chr15 <- Angus_gene %>% filter(seqid == 15) %>% 
#   filter(start >= left) %>% filter(end <= right)
# 
# #new zoom in coordinates
# Angus_gene_chr15 <- Angus_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
#   mutate(newend = end -left + 1)
# 
# Angus_gene_chr15$selectedRegionName <- selectedRegionName
# 
# newcoorAngusDF <- Angus_gene_chr15 %>% 
#   select(selectedRegionName,newstart,newend,gene_id,score,strand,gene_biotype)
# 
# #write out the results as bed for overlap in ribbon
# write_tsv(newcoorAngusDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/UOA_Angus_1_chr15_78799177_80168904.bed",
#           col_names = FALSE)

#Hereford
# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/arsucd1.2/")

#ARS-UCD1.2 
read.gtf(ens,"Bos_taurus.ARS-UCD1.2.95.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Hereford_gene <- getGenePositions(ens)
dim(Hereford_gene)

#get chr15 annotation
#length of chr15 85007780
#region of interest is 78791037 to 80120961
left <- 78791037
right <- 80120961
length <- 80120961 - 78791037 + 1
selectedRegionName <- "ARSUCD1_2_chr15_78791037_80120961"

Hereford_gene_chr15 <- Hereford_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Hereford_gene_chr15 <- Hereford_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

Hereford_gene_chr15$selectedRegionName <- selectedRegionName

newcoorHerefordDF <- Hereford_gene_chr15 %>% 
  select(selectedRegionName,newstart,newend,gene_id,score,strand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorHerefordDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/ARSUCD1_2_chr15_78791037_80120961.bed",
          col_names = FALSE)

##### Find all genes in hereford chr15 77.4 MB – 80.6 MB #####

#get chr15 annotation
#length of chr15 85007780
left <- 77400000
right <- 80600000
length <- 80600000 - 77400000 + 1
selectedRegionName <- "ARSUCD1_2_chr15_77400000_80600000"

Hereford_gene_chr15 <- Hereford_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Hereford_gene_chr15 <- Hereford_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

Hereford_gene_chr15$selectedRegionName <- selectedRegionName

newcoorHerefordDF <- Hereford_gene_chr15 %>% 
  select(selectedRegionName,newstart,newend,gene_id,score,strand,gene_biotype)

#write out the results to mine cattle gene atlas
write_tsv(newcoorHerefordDF,"/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/misassembly_issues/interesting_region/brahman_chr15_FADS2P1/synteny_draw/ARSUCD1_2_chr15_77400000_80600000_geneatlas.bed",
          col_names = FALSE)

##### Find all genes in brahman chr15 3.2 MB to 6.6 MB #####

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

setwd("/Users/lloyd/Documents/lloyd_2017/Research/Brahman_Angus/Assembly_version/final_to_correct_20180905/EBI_annotation/ensembl_annotation_Maternal_20190320/")

#Brahman EBI first draft v1.96
read.gtf(ens,"Bos_indicus_x_bos_taurus_Mom.UOA_Brahman_1.96.gtf")

# counts all annotations on each seqname
tableSeqids(ens)

# # create table of genes
Brahman_gene <- getGenePositions(ens)
dim(Brahman_gene)
# 

#get chr15 and rev comp annotation
#length of chr15 84270218
#region of interest is 3200000 to 6600000
left <- 3200000
right <- 6600000
length <- 6600000 - 3200000 + 1
selectedRegionName <- "UOA_Brahman_1_chr15_3200000_6600000_rc"

Brahman_gene_chr15 <- Brahman_gene %>% filter(seqid == 15) %>% 
  filter(start >= left) %>% filter(end <= right)

#new zoom in coordinates
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = start -left + 1) %>%
  mutate(newend = end -left + 1)

#rev comp
Brahman_gene_chr15 <- Brahman_gene_chr15 %>% mutate(newstart = length - newstart +1) %>% 
  mutate(newend = length - newend +1)
Brahman_gene_chr15$newstrand <- ifelse(Brahman_gene_chr15$strand == "+","-","+")

Brahman_gene_chr15$selectedRegionName <- selectedRegionName

newcoorDF <- Brahman_gene_chr15 %>% tibble() %>%
  dplyr::select(selectedRegionName,newend,newstart,gene_id,score,newstrand,gene_biotype)

#write out the results as bed for overlap in ribbon
write_tsv(newcoorDF,"/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/zoom_in_OR_region_brahman/UOA_Brahman_1_chr15_3200000_6600000_Equiv_to_Hereford_geneatlas.tsv",
          col_names = TRUE)
