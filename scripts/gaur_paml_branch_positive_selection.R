#------------------------------------------------------
# Program name: gaur_paml_branch_positive_selection.R
# Objective: Read the positive selection results from PAML and find the immune related sites.
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

library(ggplot2)
library(magrittr)
library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(gtools)
library(limma)
'%!in%' <- function(x,y)!('%in%'(x,y))


## branch selection

Likelihood_branch_A <- read.delim("input/From_paml/all_results_PAML_Likelihood_branch_A_gaur.txt", header = FALSE)%>%
  as.data.frame()%>%
  set_colnames(c("ID", "logLikelihood_branch_A"))

Likelihood_branch_A$ID <- gsub("all_results/","",Likelihood_branch_A$ID)

Likelihood_branch_A$ID%>%length()



# match ID
OG_PAML_list <- read.delim("input/From_paml/OG_PAML_list.txt", header = FALSE)%>%
  as.data.frame()

OG_PAML_list$V1%>%length()

stopifnot(setdiff(Likelihood_branch_A$ID,OG_PAML_list$V1) == 0)

left_join(OG_PAML_list,Likelihood_branch_A, by = c("V1" = "ID"))



Likelihood_branch_B <- read.delim("input/From_paml/all_results_PAML_Likelihood_branch_B_gaur.txt", header = FALSE)%>%
  as.data.frame()%>%
  set_colnames(c("ID", "logLikelihood_branch_B"))

Likelihood_branch_B$ID <- gsub("all_results/","",Likelihood_branch_B$ID)

Likelihood_branch_B$ID%>%length()



# The Alternative model - Null model
Sorted_Likelihood_branch_A <- Likelihood_branch_A
Sorted_Likelihood_branch_B <- Likelihood_branch_B
# Join into one table
Positive_select <- Sorted_Likelihood_branch_B%>%
  left_join(Sorted_Likelihood_branch_A, by = "ID")

Positive_select$logLikelihood_branch_B <- as.numeric(Positive_select$logLikelihood_branch_B)
Positive_select$logLikelihood_branch_A <- as.numeric(Positive_select$logLikelihood_branch_A)

# Calculate the p values
Positive_select$log_diff <- 2*(Positive_select$logLikelihood_branch_A - Positive_select$logLikelihood_branch_B)
Positive_select$p.values <- pchisq(2*(Positive_select$logLikelihood_branch_A - Positive_select$logLikelihood_branch_B),df=1, lower.tail=FALSE)/2



branch_results_Positive_select <- Positive_select%>%
  mutate(fdr = p.adjust(p.values, "fdr"),
         bonferroni = p.adjust(p.values, "bonferroni"))%>%
  dplyr::filter(fdr < 0.05)%>%
  arrange(fdr)

match_table_species <- read_csv("output/CAFE/match_table_species.csv")
length(branch_results_Positive_select$ID%>%unique()) #665 Positive selected genes

results_Positive_select <- branch_results_Positive_select%>%
  left_join(match_table_species,by = c("ID" = "orthogroup_ID"))



colnames(results_Positive_select) <- gsub("spc", "species",colnames(results_Positive_select))
results_Positive_select <- unique(results_Positive_select[,c(1,2,3,4,5,6,7,8,9,10,11,13)])

results_Positive_select <- results_Positive_select[,1:11]%>%
  unique()%>%
  left_join(drop_na(results_Positive_select[,c(1,12)]))




# detlete other columns before reading
# OGtbl <- read.table("input/From_OrthoFinder/OGtbl.tsv",header = T)%>%
#   set_colnames(c("orthogroup_ID", "geneID", "spc"))
# 
# OGtbl$orthogroup_ID <- gsub("OG1v", "OG",OGtbl$orthogroup_ID)
# 
# # Check number of seqs for each species
# table(OGtbl$spc)
# # How many orthogroup ID
# OGtbl$orthogroup_ID%>%unique()%>%length()
# OGtbl$geneID <- as.character(OGtbl$geneID)
# 
# results_Positive_select%<>%left_join(subset(OGtbl,spc %in% "Bgau"), by = c("ID" = "orthogroup_ID"))
# write_csv(branch_results_Positive_select,"output/PAML/branch_results_Positive_select.csv")




InnateDB_immune <- read_csv("input/Immune_gene_ref/InnateDB_genes_all.csv")%>%
  as.data.frame()
InnateDB_immune$entrez <- as.character(InnateDB_immune$entrez)

results_Positive_select$InnateDB <- results_Positive_select$geneID %in% InnateDB_immune$ensembl
results_Positive_select%>%write_csv("output/PAML/branch_results_Positive_select.csv")




ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond <- readRDS("input/From_blast/ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond.rds")
head(ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond)



# detlete other columns before reading
OGtbl <- read.table("input/From_OrthoFinder/OGtbl.tsv",header = T)%>%
  set_colnames(c("orthogroup_ID", "geneID", "spc"))

OGtbl$orthogroup_ID <- gsub("OG1v", "OG",OGtbl$orthogroup_ID)



results_Positive_select_OGtbl_Bgau <- subset(OGtbl, orthogroup_ID %in% results_Positive_select$ID)%>%
  subset(spc %in% "Bgau")%>%
  left_join(ARS_UOA_Gaur_1pep_vs_ARS_UCD1_2pep_correspond,by = c("geneID"= "Bosg_gene_id_nover"))%>%
  unique()

head(results_Positive_select_OGtbl_Bgau)
write_csv(results_Positive_select_OGtbl_Bgau, "input/From_paml/results_Positive_select_OGtbl_Bgau.csv")

