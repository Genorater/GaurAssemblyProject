#------------------------------------------------------
# Program name: gaur_paml_go_kegg_pathway_analysis.R
# Objective: GO and KEGG pathway analysis for positive selection results
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

# Packages

library(magrittr)
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(gtools)
library(tibble)
library(limma)

'%!in%' <- function(x,y)!('%in%'(x,y))


# Analysis

# Read data

OGtbl <- read.table("input/From_OrthoFinder/OGtbl.tsv",header = T)%>%
  set_colnames(c("orthogroup_ID", "geneID", "spc"))

OGtbl$orthogroup_ID <- gsub("OG1v", "OG",OGtbl$orthogroup_ID)

# Check species
table(OGtbl$spc)



results_Positive_select <- read.csv("output/PAML/branch_results_Positive_select.csv")
results_Positive_select%<>%subset(InnateDB %in% "TRUE")



colnames(results_Positive_select) <- gsub("spc", "species",colnames(results_Positive_select))



ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond <- readRDS("input/From_blast/ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond.rds")
head(ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond)



# How many orthogroup ID
OGtbl$orthogroup_ID%>%unique()%>%length()
OGtbl$geneID <- as.character(OGtbl$geneID)



results_Positive_select_OGtbl_Bgau <- subset(OGtbl, orthogroup_ID %in% results_Positive_select$ID)%>%
  subset(spc %in% "Bgau")%>%
  left_join(ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond,by = c("geneID"= "Bosg_gene_id_nover"))%>%
  unique()

head(results_Positive_select_OGtbl_Bgau)


## Gaur (Cattle database)


results_Positive_select_OGtbl_Hbta <- subset(OGtbl, orthogroup_ID %in% results_Positive_select$ID)%>%
  subset(spc %in% "Hbta")%>%
  unique()

head(results_Positive_select_OGtbl_Hbta)



#ah <- AnnotationHub()
#saveRDS(ah,"All_annotation.rds")
ah <- read_rds("/Users/kellydren/Documents/Kelly_annotation/All_annotation.rds") 
#subset(ah, rdataclass == "EnsDb" & species == "Bos taurus")
ensDb <- ah[["AH83145"]]
ensDb



genesGR <- GenomicFeatures::genes(ensDb)
genesGR


{r }
cols2Keep <- c("gene_id","gene_name", "gene_biotype", "description", "entrezid")
mcols(genesGR) <- mcols(genesGR)[, cols2Keep]

# get gene annotation
Genes <- genesGR%>%
  as.data.frame()
Genes$entrezid <- Genes$entrezid%>%as.character()



cattle_ALL_entrezID <- genesGR %>% 
  subset(!is.na(entrezid)) %>%
  mcols() %>%
  .[["entrezid"]] %>%
  unlist() %>%
  unique() 


# trans to entrezid
results_Positive_select_OGtbl_Bgau$Bost_gene_id %<>% gsub("\\..","",. )

results_Positive_select_OGtbl_Bgau <- results_Positive_select_OGtbl_Bgau%>%left_join(Genes[,c("gene_id","entrezid")], by = c("Bost_gene_id" = "gene_id")) 

results_Positive_select_OGtbl_Bgau_entrezid <- results_Positive_select_OGtbl_Bgau$entrezid%>%
  gsub("c","",.)%>%
  gsub("\\(","",.)%>%
  gsub("\\)","",.)%>%
  str_split( ", ")%>%
  unlist()


### GO
goRes <- goana(results_Positive_select_OGtbl_Bgau_entrezid, cattle_ALL_entrezID, species = "Bt")
Gaur_goRes <- goRes%>%
  rownames_to_column("GO_ID")%>%
  mutate(fdr = p.adjust(P.DE, "fdr"))%>%
  subset(fdr < 0.05)%>%
  arrange(fdr)

Gaur_goRes%>%
  dim()

head(Gaur_goRes)


### KEGG
keggRes <- kegga(results_Positive_select_OGtbl_Bgau_entrezid, cattle_ALL_entrezID, species = "Bt")

Gaur_keggRes <- keggRes%>% 
  rownames_to_column("KEGG_ID")%>%
  mutate(fdr = p.adjust(P.DE, method = "fdr"))%>%
  subset(fdr < 0.05)%>%
  arrange(fdr)

Gaur_keggRes%>%
  dim()

head(Gaur_keggRes)

# ECM-receptor interaction pathways were the most upregulated gene-enriched signaling pathways. They play an important role in the process of tumor shedding, adhesion, degradation, movement and hyperplasia. The role of ECM in other cancers has been proved.

