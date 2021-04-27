#------------------------------------------------------
# Program name: gaur_heatmap.R
# Objective: to compare the gaur missing genes in divergent
#           region chr15 against hereford and brahman
# Author: Kelly Ren
# Email add: KellyDREN@gmail.com
#------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(viridis)  
# Read the file
setwd("~/Documents/Kelly_2021/Research/Lloyd_Gaur/Cooperation")
ARSUCD1_2_chr15_77400000_80600000_geneatlas <- read_tsv("ARSUCD1_2_chr15_77400000_80600000_geneatlas.bed", col_names = FALSE)
FPKM_mean <- read_tsv("FPKM_mean.tsv", col_names = TRUE)

head(FPKM_mean)
# prefix             Tissue             mean_FPKM    no
# <chr>              <chr>                  <dbl> <dbl>
#   1 ENSBTAG00000018742 Abomasum                54.4     1
# 2 ENSBTAG00000018742 Adipose                 22.4     1
# 3 ENSBTAG00000018742 Adrenal                 29.2     1
# 4 ENSBTAG00000018742 Ampulla                 24.1     1
# 5 ENSBTAG00000018742 Anterior_pituitary      23.6     1
# 6 ENSBTAG00000018742 Aorta                   19.9     1

# set factor so the gene will be plot in order
FPKM_mean$prefix <- factor(FPKM_mean$prefix,levels = unique(FPKM_mean$prefix))

# plot in three colors
ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= mean_FPKM)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissues")+
  labs(fill = "mean FPKM") +
  scale_fill_gradientn(colours = c("white", "blue", "red"),
                      values = scales::rescale(c(0,60,150,250))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))

# plot in two colors
ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= mean_FPKM)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissues")+
  labs(fill = "mean FPKM") +
  scale_fill_gradientn(colours = c("white", "blue"),
                       values = scales::rescale(c(0,100,250))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))

# calculate the log2 values
FPKM_mean$log_mean <- log2(FPKM_mean$mean_FPKM)
# set infinite values to -14 == white
FPKM_mean$log_mean[is.infinite(FPKM_mean$log_mean)] <- -14

ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= log_mean)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissues")+
  labs(fill = "log2(mean FPKM)") +
  scale_fill_gradientn(colours = c("white", "blue"),
                       values = scales::rescale(c(-14,0,7))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))

