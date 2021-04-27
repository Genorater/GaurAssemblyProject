#------------------------------------------------------
# Program name: gaur_species_specific_genes.R
# Objective: to find any gaur specific genes using the
#           rationale the orthogroup has only gaur genes
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(dplyr)
library(readr)

setwd("/Users/lloyd/Documents/lloyd_2019/Research/gaur/species_specific_genes/")

# read in Orthogroups.tsv
Orthogroups <- read_tsv("Orthogroups.tsv", col_names = TRUE)

# hereford
# check that hereford can reproduce the result ENSBTAG00000037878, which is Hbta specific
# loop through and ask the question whether is.na for all except the wanted species

selc <- c()

for (i in 1:nrow(Orthogroups)){
  logic <- is.na(Orthogroups[i,c(2:6,8:10)])
  hereford_specific_logic <- sum(logic) == 8
  selc <- c(selc, hereford_specific_logic)
}

hereford_specific_OG <- Orthogroups$Orthogroup[selc]
hereford_specific_OG

# gaur

selc <- c()

for (i in 1:nrow(Orthogroups)){
  logic <- is.na(Orthogroups[i,c(2,4:10)])
  gaur_specific_logic <- sum(logic) == 8
  selc <- c(selc, gaur_specific_logic)
}

gaur_specific_OG <- Orthogroups$Orthogroup[selc]
gaur_specific_OG

gaur_specific_OG_gp1 <- Orthogroups %>% filter(Orthogroup == gaur_specific_OG[1]) 
gaur_specific_OG_gp1$Bgau
# [1] "ENSBGUG00000009969, ENSBGUG00000013779, ENSBGUG00000021745"

gaur_specific_OG_gp2 <- Orthogroups %>% filter(Orthogroup == gaur_specific_OG[2]) 
gaur_specific_OG_gp2$Bgau
# [1] "ENSBGUG00000013971, ENSBGUG00000018258"

gaur_specific_OG_gp3 <- Orthogroups %>% filter(Orthogroup == gaur_specific_OG[3]) 
gaur_specific_OG_gp3$Bgau
# [1] "ENSBGUG00000009775, ENSBGUG00000009810"
