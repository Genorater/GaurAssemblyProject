#------------------------------------------------------
# Program name: gaur_reciprocal_best_hit_cattle.R
# Objective: find the longest isoform from gaur that
#           best blastp hit the hereford and vice versa
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(magrittr)
library(biomaRt)

setwd("/Users/lloyd/Documents/lloyd_2019/Research/gaur/gaur_BlastDB/reciprocal_best_hit_gaur_vs_hereford")

# read in data, gaur blast ref hereford
gaur_blast_ref_hereford <- read_tsv("ARS_UOA_Gaur_1longpep_vs_ARS_UCD1_2longpep.blpm6", col_names = FALSE)

colnames(gaur_blast_ref_hereford) <- 
  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

# read in data, hereford blast ref gaur
hereford_blast_ref_gaur <- read_tsv("ARS_UCD1_2longpep_vs_ARS_UOA_Gaur_1longpep.blpm6", col_names = FALSE)

colnames(hereford_blast_ref_gaur) <- 
  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

# loop through gaur blast ref hereford table
# for each gaur longest isoform protein, what is the hereford protein?
# using the matching hereford protein, does it match back to the exact gaur protein
# in other words, reciprocal best hit for longest isoform, both ways

# Initialize rbh vector to keep all reciprocal best hit matches
rbh <- c()

for (i in 1:nrow(gaur_blast_ref_hereford)){
  gaur_blast_ref_hereford_tbl_gaur_id <- gaur_blast_ref_hereford$qseqid[i]
  gaur_blast_ref_hereford_tbl_hereford_id <- gaur_blast_ref_hereford$sseqid[i]
  
  hereford_selc_logic <- hereford_blast_ref_gaur$qseqid %in% gaur_blast_ref_hereford_tbl_hereford_id
  hereford_blast_ref_gaur_tbl_gaur_id <- hereford_blast_ref_gaur$sseqid[hereford_selc_logic]
  
  if (identical(hereford_blast_ref_gaur_tbl_gaur_id, character(0))){
    hereford_blast_ref_gaur_tbl_gaur_id <- "notavailable"
  } 
  
  if (gaur_blast_ref_hereford_tbl_gaur_id == hereford_blast_ref_gaur_tbl_gaur_id){
    dummy <- "yes"
    rbh <- c(rbh,dummy)
  } else {
    dummy <- "no"
    rbh <- c(rbh,dummy)
  }
}

gaur_blast_ref_hereford_rbh <- gaur_blast_ref_hereford %>% mutate(rbh = rbh)

# read in biomart to get descriptions
Cattle = useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")
# mart = useMart('ensembl')
# listDatasets(mart)

### cattle genes match their ensembl_gene_id
cattle_genes <- getBM(
  attributes= c("chromosome_name","start_position","end_position","strand","ensembl_gene_id","external_gene_name", "description"),mart= Cattle)

head(cattle_genes)

gaur_blast_ref_hereford_rbh_desc <- left_join(gaur_blast_ref_hereford_rbh,cattle_genes, by = c("sseqid"= "ensembl_gene_id")) # %>% as.data.frame()

chr_name_int <- as.integer(gaur_blast_ref_hereford_rbh_desc$chromosome_name)

gaur_blast_ref_hereford_rbh_desc_sort <- gaur_blast_ref_hereford_rbh_desc %>% mutate(chr_name_int = chr_name_int) %>%
  arrange(chr_name_int,start_position)

gaur_blast_ref_hereford_rbh_desc_sort_chr15 <- gaur_blast_ref_hereford_rbh_desc_sort %>% filter(chromosome_name == "15")

gaur_blast_ref_hereford_rbh_desc_sort_chr15_divergent_expanded <- gaur_blast_ref_hereford_rbh_desc_sort_chr15 %>%
  filter(start_position >= 77409614) %>% filter(end_position <= 80584707)

write_tsv(gaur_blast_ref_hereford_rbh_desc_sort_chr15_divergent_expanded, 
          "gaur_blast_ref_hereford_rbh_desc_sort_chr15_divergent_expanded.tsv", col_names = TRUE)
