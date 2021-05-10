#------------------------------------------------------
# Program name: gaur_find_genes_in_GO_terms.R
# Objective: The GO and KEGG pathway analysis for families on gaur branch 
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

# Loading packages

library(ggplot2)
library(magrittr)
library(dplyr)
library(tibble)
library(biomaRt)
# Analysis

# Cattle database

Cattle = useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")

### GO terms match their entrezgene_id and ensembl_gene_id
cattle_genes_GO <- getBM(
  attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name", "description","go_id","name_1006", "definition_1006"),mart= Cattle)

saveRDS(cattle_genes_GO, "output/CAFE/ALL_cattle_genes.rds")
head(cattle_genes_GO)

# Read GO terms
cattle_genes_GO <- read_rds("output/CAFE/ALL_cattle_genes.rds")
head(cattle_genes_GO)


Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO <- read_csv("output/CAFE/Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO.csv")%>%
  as.data.frame()

# detlete other columns before reading
OGtbl <- read.table("input/From_OrthoFinder/OGtbl.tsv",header = T)%>%
  set_colnames(c("orthogroup_ID", "geneID", "spc"))
OGtbl$orthogroup_ID <- gsub("OG1v", "OG",OGtbl$orthogroup_ID)

Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch <- read.csv("output/CAFE/Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch.csv",header = T)

OGtbl <- subset(OGtbl, orthogroup_ID %in% Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch$orthogroup_ID)
cattle_genes_GO$InData <- cattle_genes_GO$ensembl_gene_id %in% OGtbl$geneID

#subset(cattle_genes_GO, go_id %in% Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO$GO_ID)%>%
#  left_join(Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch, by = c("ensembl_gene_id" = "geneID"))%>%
#  write_csv("output/CAFE/Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO_genes.csv")

Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO_genes <- read_csv("output/CAFE/Hbtaref_genes_in_orthogroup_of_sigFamilies_CAFE_GO_genes.csv")%>%
  as.data.frame()

ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep <- read_csv("input/From_blast/ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep.csv")%>%as.data.frame()

