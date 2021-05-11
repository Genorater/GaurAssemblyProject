#------------------------------------------------------
# Program name: gaur_find_Btau_gene_atlas.R
# Objective: load and display cattle gene atlas
#           
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(magrittr)

#SLC
# read in all gene atlas csv files
path_to_atlas <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/gene_gains_losses_postCAFE_specific_orthogroups/ruminant_specific/SLC/"

#create a list to read all *csv file
file_list <- list.files(path = path_to_atlas, pattern = "*.csv")

final_DF <- c()

for (i in 1:length(file_list)){
  file <- assign(file_list[i], read_tsv(paste0(path_to_atlas,file_list[i]),col_names = TRUE))
  prefix <- gsub(".csv","",file_list[i])
  DF <- file
  DF$prefix <- prefix
  
  final_DF <- rbind(final_DF, DF)
}

# mean FPKM after group by genes and then tissue
FPKM_mean <- final_DF %>% group_by(prefix,Tissue) %>% summarise(mean_FPKM = mean(FPKM))

#heatmap
# set factor so the gene will be plot in order
FPKM_mean$prefix <- factor(FPKM_mean$prefix,levels = unique(FPKM_mean$prefix))

# plot in two colors
tiff(filename = "heatmap_ruminant_specific_SLC.tiff", width = 700, height = 700)
ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= mean_FPKM)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissue")+
  labs(fill = "FPKM") +
  scale_fill_gradientn(colours = c("white", "blue"),
                       values = scales::rescale(c(0,100,250))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))
dev.off()

rm(list=ls())

#lysozyme
# read in all gene atlas csv files
path_to_atlas <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/gene_gains_losses_postCAFE_specific_orthogroups/ruminant_specific/lysozyme/"

#create a list to read all *csv file
file_list <- list.files(path = path_to_atlas, pattern = "*.csv")

final_DF <- c()

for (i in 1:length(file_list)){
  file <- assign(file_list[i], read_tsv(paste0(path_to_atlas,file_list[i]),col_names = TRUE))
  prefix <- gsub(".csv","",file_list[i])
  DF <- file
  DF$prefix <- prefix
  
  final_DF <- rbind(final_DF, DF)
}

# mean FPKM after group by genes and then tissue
FPKM_mean <- final_DF %>% group_by(prefix,Tissue) %>% summarise(mean_FPKM = mean(FPKM))

#heatmap
# set factor so the gene will be plot in order
FPKM_mean$prefix <- factor(FPKM_mean$prefix,levels = unique(FPKM_mean$prefix))
FPKM_mean$Tissue %<>% gsub("_"," ",.) 
# plot in two colors
#with log2 of values
# tiff(filename = "heatmap_ruminant_specific_SLC.tiff", width = 700, height = 700)
# ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= log2(mean_FPKM+0.1))) + 
#   geom_tile() + 
#   xlab("Gene") + 
#   ylab("Tissue")+
#   labs(fill = "FPKM") +
#   scale_fill_gradientn(colours = c("white", "blue"),
#                        values = scales::rescale(c(0,100,250))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))
# dev.off()

#straight plotting with FPKM values
ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= mean_FPKM)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissue")+
  labs(fill = "FPKM") +
  scale_fill_gradientn(colours = c("white", "blue"),
                       values = scales::rescale(c(0,100,250))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))
