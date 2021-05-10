#------------------------------------------------------
# Program name: gaur_cafe_match_gene_famlies.R
# Objective: group the orthongroups into PANTHER gene families
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

library(magrittr)
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(gtools)
library(tibble)
library(data.table)
library(biomaRt)
'%!in%' <- function(x,y){!('%in%'(x,y))}


# Get PTHR_UniProtKB match


#PANTHER15 <- read.table("PANTHER15.txt", header = F)

#PANTHER15 <- lapply(c(1:length(PANTHER15$V1)), function(x){gsub(".*=","",PANTHER15$V1[x])})

#PANTHER15 <- unlist(PANTHER15)%>%
# data.frame()%>%
# set_colnames("X1")



# PTHR_rank <- PANTHER15$X1%>%grep("PTHR",.)
# 
# ID_match <- lapply(c(1:length(PTHR_rank)), function(x){
#   if (x <= 15701){
#     GAP <- (as.numeric(PTHR_rank[x+1]-PTHR_rank[x])-2)
#   }else{
#     GAP <- "-1"
#   }
#   if (GAP < 0){
#     OUTPUT <- "NULL"
#     OUTPUT
#   }else{
#     PANTHER15$X1[(c(as.numeric(PTHR_rank[x])+1)):((c(as.numeric(PTHR_rank[x])+1))+GAP)] 
#   }
# })%>%
#   set_names(PANTHER15$X1[PTHR_rank])
# 
# 
# PANTHER15 <- lapply(1:length(ID_match), FUN = function(x, object = ID_match){
#   object[[x]]%>%
#     as.data.frame()%>%
#     set_colnames("UniProtKB")%>%
#     mutate(PANTHER_ID = rep(names(object[x]), length(object[[x]])))
# })%>%do.call("rbind",.)
# 
# PTHR_UniProtKB <- subset(PANTHER15, UniProtKB %!in% "NULL")
# PTHR_UniProtKB$PANTHER_ID%>%unique%>%length() #15701

PTHR_UniProtKB <- read_rds("Input/PTHR_15/PTHR_UniProtKB.rds") 


# Match by the ID (No blast method)


# detlete other columns before reading
OGtbl <- read.table("input/From_OrthoFinder/OGtbl.tsv",header = T)%>%
  set_colnames(c("orthogroup_ID", "geneID", "spc"))

OGtbl$orthogroup_ID <- gsub("OG1v", "OG",OGtbl$orthogroup_ID)

# Check number of seqs for each species
table(OGtbl$spc)
# How many orthogroup ID
OGtbl$orthogroup_ID%>%unique()%>%length()
OGtbl$geneID <- as.character(OGtbl$geneID)


## Human

#Human = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

### GO terms match their entrezgene_id and ensembl_gene_id
#ALL_Human_genes <- getBM(
#  attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name","gene_biotype", "description","entrezgene_id"),mart= Human)

#saveRDS(ALL_Human_genes, "ALL_Human_genes.rds")
ALL_Human_genes <- read_rds("/Users/kellydren/Documents/Kelly_annotation/ALL_Human_genes.rds")%>%
  subset(gene_biotype %in% "protein_coding")

ALL_Human_genes$entrezgene_id <- as.character(ALL_Human_genes$entrezgene_id)

OGtbl_Human <- OGtbl%>%subset(spc %in% c("Hsap"))%>%
  left_join(ALL_Human_genes[,c(5,6,9)],by = c("geneID" = "ensembl_gene_id"))

dim(OGtbl_Human)

# 19127 entrezgene_id found in Human
(!is.na(OGtbl_Human$entrezgene_id))%>%table()


## Hbta


#Cattle = useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")

### GO terms match their entrezgene_id and ensembl_gene_id
# ALL_Hbta_genes <- getBM(
#   attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name","gene_biotype", "description","entrezgene_id"),mart= Cattle)

#saveRDS(ALL_Hbta_genes, "ALL_Hbta_genes.rds")
ALL_Hbta_genes <- read_rds("/Users/kellydren/Documents/Kelly_annotation/ALL_Hbta_genes.rds")%>%subset(gene_biotype %in% "protein_coding")

ALL_Hbta_genes$entrezgene_id <- as.character(ALL_Hbta_genes$entrezgene_id)

OGtbl_Cattle <- OGtbl%>%subset(spc %in% c("Hbta"))%>%
  left_join(ALL_Hbta_genes[,c(5,6,9)],by = c("geneID" = "ensembl_gene_id"))

dim(OGtbl_Cattle)

# 18493 entrezgene_id found in Cattle
(!is.na(OGtbl_Cattle$entrezgene_id))%>%table()


## Chir


#goat = useMart(biomart="ensembl", dataset="chircus_gene_ensembl")

### GO terms match their entrezgene_id and ensembl_gene_id
# ALL_goat_genes <- getBM(
#   attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name","gene_biotype", "description","entrezgene_id"),mart= goat)

#saveRDS(ALL_goat_genes, "ALL_goat_genes.rds")
ALL_goat_genes <- read_rds("/Users/kellydren/Documents/Kelly_annotation/ALL_goat_genes.rds")%>%subset(gene_biotype %in% "protein_coding")

ALL_goat_genes$entrezgene_id <- as.character(ALL_goat_genes$entrezgene_id)

OGtbl_Goat <- OGtbl%>%subset(spc %in% c("Chir"))%>%
  left_join(ALL_goat_genes[,c(5,6,9)],by = c("geneID" = "ensembl_gene_id"))

dim(OGtbl_Goat)

# 17401 entrezgene_id found in Goat
(!is.na(OGtbl_Goat$entrezgene_id))%>%table()


## Sscr


#Pig = useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")

### GO terms match their entrezgene_id and ensembl_gene_id
#ALL_pig_genes <- getBM(
#  attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name","gene_biotype", "description","entrezgene_id"),mart= Pig)

#saveRDS(ALL_pig_genes, "ALL_pig_genes.rds")
ALL_pig_genes <- read_rds("/Users/kellydren/Documents/Kelly_annotation/ALL_pig_genes.rds")%>%subset(gene_biotype %in% "protein_coding")

ALL_pig_genes$entrezgene_id <- as.character(ALL_pig_genes$entrezgene_id)

OGtbl_Pig <- OGtbl%>%subset(spc %in% c("Sscr"))%>%
  left_join(ALL_pig_genes[,c(5,6,9)],by = c("geneID" = "ensembl_gene_id"))

dim(OGtbl_Pig)

# 17936 entrezgene_id found in Pig
(!is.na(OGtbl_Pig$entrezgene_id))%>%table()


# Join_matching

# Check format for rbind
head(OGtbl_Human)
head(OGtbl_Cattle)
head(OGtbl_Goat)
head(OGtbl_Pig)

# match the colnames
OGtbl_Cattle <- OGtbl_Cattle[,c(1,2,3,4,5)]%>%
  set_colnames(colnames(OGtbl_Human))
dim(OGtbl_Cattle)

OGtbl_Goat <- OGtbl_Goat[,c(1,2,3,4,5)]%>%
  set_colnames(colnames(OGtbl_Human))
dim(OGtbl_Goat)

OGtbl_Pig <- OGtbl_Pig[,c(1,2,3,4,5)]%>%
  set_colnames(colnames(OGtbl_Human))
dim(OGtbl_Pig)

# unique each OG for one species
OGtbl_select <- rbind(OGtbl_Human,subset(OGtbl_Cattle,orthogroup_ID %!in% c(OGtbl_Human$orthogroup_ID%>%unique())))

OGtbl_select <- rbind(OGtbl_select,subset(OGtbl_Pig,orthogroup_ID %!in% c(OGtbl_select$orthogroup_ID%>%unique())))

OGtbl_select <- rbind(OGtbl_select,subset(OGtbl_Goat,orthogroup_ID %!in% c(OGtbl_select$orthogroup_ID%>%unique())))

# How many for each species
OGtbl_select$spc%>%table()
# Chir  Hbta  Hsap  Sscr 
#  638  2706 20027  1587  
OGtbl_select$orthogroup_ID%>%unique%>%length()

OGtbl_select%>%
  write_csv("output/CAFE/OGtbl_select.csv")



Entrez_ToUniProtKB <- read.table("input/ToUniProtKB/Entrez_ToUniProtKB.txt",header = T)%>%
  set_colnames(c("Entrez_ID", "UniProtKB_ID"))

# Fix mutiple Entrez ID to one UniProtKB ID
multi_Entrez_ToUniProtKB <- Entrez_ToUniProtKB[Entrez_ToUniProtKB$Entrez_ID%>%grep(",",.),]


multi_Entrez_ToUniProtKB <- lapply(c(1:dim(multi_Entrez_ToUniProtKB)[1]), function(x){
  str_split(multi_Entrez_ToUniProtKB$Entrez_ID[x], pattern = ",")%>%
    as.data.frame()%>%
    mutate(UniProtKB=rep(multi_Entrez_ToUniProtKB$UniProtKB_ID[x]))%>%
    set_colnames(c("Entrez_ID","UniProtKB_ID"))
})%>%do.call("rbind",.)

Entrez_ToUniProtKB <- Entrez_ToUniProtKB[setdiff(c(1:dim(Entrez_ToUniProtKB)[1]), Entrez_ToUniProtKB$Entrez_ID%>%grep(",",.)),]%>%
  rbind(multi_Entrez_ToUniProtKB)

Entrez_ToUniProtKB$Entrez_ID%>%unique()%>%length() #19368
Entrez_ToUniProtKB$UniProtKB_ID%>%unique()%>%length() #32763


match_table_species <- OGtbl_select%>%
  drop_na()%>%
  left_join(Entrez_ToUniProtKB, by = c( "entrezgene_id"="Entrez_ID"))%>%
  left_join(PTHR_UniProtKB, by = c( "UniProtKB_ID"="UniProtKB"))

ALL_orthogroup_PTHRfam <- match_table_species[,c(1,7)]



## Find the genes in orthogroups for all species


# the match_table_species contain the matching information between orthogroup_ID, geneID, UniProtKB_ID and PANTHER_ID for selected sequences in script CAFE_Match_gene_family.Rmd

# make a overall table for orthogroup_ID and PANTHER_ID and their genes for each species
match_table_species[,c("orthogroup_ID","PANTHER_ID")]%>%
  drop_na()%>%
  unique()%>%
  left_join(OGtbl,by = "orthogroup_ID")%>%
  write_csv("output/CAFE/genes_in_orthogroup_all_species.csv")


# fix duplicated PTHRfam for the same orthogroup

ALL_orthogroup_PTHRfam[ALL_orthogroup_PTHRfam$orthogroup_ID%>%duplicated(),]%>%extract2("orthogroup_ID")%>%unique()%>%length() #8613

dup_ALL_orthogroup_PTHRfam <- ALL_orthogroup_PTHRfam[ALL_orthogroup_PTHRfam$orthogroup_ID%>%duplicated(),]

dup_list <- subset(ALL_orthogroup_PTHRfam, orthogroup_ID %in% dup_ALL_orthogroup_PTHRfam$orthogroup_ID)%>%extract2("orthogroup_ID")%>%unique()

Themost_pickedPTHRID <- lapply(c(1:length(dup_list)), function(x){subset(ALL_orthogroup_PTHRfam, orthogroup_ID %in% dup_list[x])%>%
    extract2("PANTHER_ID")%>%
    table()%>%
    which.max()%>%
    names()})%>%
  set_names(dup_list)

Themost_pickedPTHRID <- Themost_pickedPTHRID%>%
  do.call("rbind",.)%>%
  as.data.frame()%>%
  rownames_to_column("orthogroup_ID")%>%
  set_colnames(c("orthogroup_ID","PANTHER_ID"))

head(Themost_pickedPTHRID)

# Unique the duplicated OG
ALL_orthogroup_PTHRfam <- subset(ALL_orthogroup_PTHRfam, orthogroup_ID %!in% unique(dup_ALL_orthogroup_PTHRfam$orthogroup_ID))
ALL_orthogroup_PTHRfam <- ALL_orthogroup_PTHRfam%>%
  rbind(Themost_pickedPTHRID)

OGtbl$orthogroup_ID%>%
  unique()%>%
  length()

ALL_orthogroup_PTHRfam$orthogroup_ID%>%
  duplicated()%>%
  table()

ALL_orthogroup_PTHRfam$PANTHER_ID%>%
  unique()%>%
  length()

# No duplicated match 17536 out of 19920 matched


# Combine for cafe


# The Orthogroups.GeneCount.tsv table is from OrthoFinder
cafe_input_data <- read.table("input/From_OrthoFinder/Orthogroups.GeneCount.tsv", header = T)
cafe_input_data <- cafe_input_data%>%left_join(ALL_orthogroup_PTHRfam, by=c("Orthogroup" = "orthogroup_ID"))%>%unique()

cafe_input_data$PANTHER_ID%>%is.na()%>%table()
# FALSE  TRUE 
# 16926  2994 
cafe_input_data$PANTHER_ID%>%unique()%>%length() #7414

cafe_input_data <- cafe_input_data[!cafe_input_data$PANTHER_ID%>%is.na(),]

cafe_input_data <- cafe_input_data[,c(1,11,2:10)]%>%
  set_colnames(colnames(cafe_input_data[,c(1, 11,2:10)]))

colnames(cafe_input_data) <- colnames(cafe_input_data)%>%gsub("PANTHER_ID","FAMILY",.)
head(cafe_input_data)




com_cafe_input_data <- cafe_input_data%>%
  split(., f = as.factor(cafe_input_data$FAMILY))

com_cafe_input_data <- lapply(c(1:length(com_cafe_input_data)), function(x){
  com_cafe_input_data[[x]][,c(1:2)][1,]%>%cbind(com_cafe_input_data[[x]][,c(3:11)]%>%
                                                  colSums()%>%
                                                  as.data.frame()%>%
                                                  set_rownames(names(colSums(com_cafe_input_data[[x]][,c(3:11)])))%>%t())})%>%
  do.call("rbind",.)

com_cafe_input_data$FAMILY%>%is.na()%>%table()
com_cafe_input_data$FAMILY%>%unique()%>%length() # match into 7413 gene families
write_csv(com_cafe_input_data, "output/CAFE/com_cafe_input_data.csv") #7395 at last

# nsites

# Check number of sites used for r8s
nsite <- read.table("nsite.txt",header = F )
sum(nsite$V2)


