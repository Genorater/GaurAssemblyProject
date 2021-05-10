#------------------------------------------------------
# Program name: gaur_cafe_go_kegg_pathway_analysis.R
# Objective: The GO and KEGG pathway analysis for families on gaur branch 
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

# Loading packages

library(magrittr)
library(readr)
library(dplyr)
library(gtools)
library(tibble)
library(biomaRt)
library(limma)
'%!in%' <- function(x,y)!('%in%'(x,y))

# Analysis

## Cattle database


#ah <- AnnotationHub()
#saveRDS(ah,"All_annotation.rds")
ah <- read_rds("/Users/kellydren/Documents/Kelly_annotation/All_annotation.rds") 
#subset(ah, rdataclass == "EnsDb" & species == "Bos taurus")
ensDb <- ah[["AH83145"]]
ensDb

genesGR <- GenomicFeatures::genes(ensDb)
genesGR


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

CAFE_sig_families_Bgaurbranch <- rbind(read_csv("output/CAFE/Bgaur_sig_families_gains.csv")%>%
                                         as.data.frame(),
                                       read_csv("output/CAFE/Bgaur_sig_families_losses.csv")%>%
                                         as.data.frame())

genes_in_orthogroup_all_species <- read_csv("output/CAFE/genes_in_orthogroup_all_species.csv")

Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch <- subset(genes_in_orthogroup_all_species, PANTHER_ID %in% CAFE_sig_families_Bgaurbranch$FAMILY)%>%
  subset(spc %in% "Hbta")

# There are 16 families without cattle gene
write.csv(Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch,"output/CAFE/Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch.csv")


## Bgaur CAFE PANTHER gene family 
### Check all 373 families


Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch <- read_csv("output/CAFE/Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch.csv")%>%
  as.data.frame()

# There are 16 families without cattle gene
# result in 357 families
Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch$PANTHER_ID%>%unique()%>%length()
Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch$orthogroup_ID%>%unique()%>%length()
Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch$geneID%>%unique()%>%length()



# trans to entrezid
Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch <- Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch%>%left_join(Genes[,c("gene_id","entrezid")], by = c("geneID" = "gene_id")) 

Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch_entrezid <- Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch$entrezid%>%
  gsub("c","",.)%>%
  gsub("\\(","",.)%>%
  gsub("\\)","",.)%>%
  str_split( ", ")%>%
  unlist()


#### GO

goRes <- goana(Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch_entrezid, cattle_ALL_entrezID, species = "Bt")

cattle_goRes <- goRes%>%
  rownames_to_column("GO_ID")%>%
  mutate(fdr = p.adjust(P.DE, "fdr"))%>%
  subset(fdr < 0.05)%>%
  arrange(fdr)

cattle_goRes%>%
  dim()

cattle_goRes


keggRes <- kegga(Hbtaref_genes_in_orthogroup_CAFE_Bgaurbranch_entrezid, cattle_ALL_entrezID, species = "Bt") 

cattle_keggRes <- keggRes%>% 
  rownames_to_column("KEGG_ID")%>%
  mutate(fdr = p.adjust(P.DE, method = "fdr"))%>%
  subset(fdr < 0.05)%>%
  arrange(fdr)

cattle_keggRes%>%
  dim()

cattle_keggRes

