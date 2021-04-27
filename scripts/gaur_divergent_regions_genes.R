#------------------------------------------------------
# Program name: gaur_divergent_regions_genes.R
# Objective: to compare the gaur missing genes in divergent
#           region chr15 against hereford and brahman
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(dplyr)
library(readr)

setwd("/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/ribbon_plots/nucmer_divergent_chr15")

#hereford
hereford_15_78791037_80120961 <- read_tsv("hereford_15_78791037_80120961.bed", col_names = FALSE)

colnames(hereford_15_78791037_80120961) <- 
  c("selectedRegionName","newstart","newend","gene_id","score","strand","gene_biotype")

HEREFORD_VS_GAUR_divergeChr15.delta.coords <- read_tsv("HEREFORD_VS_GAUR_divergeChr15.delta.coords",
                                                       col_names = FALSE)

names(HEREFORD_VS_GAUR_divergeChr15.delta.coords) <- c("Ref_start",	"Ref_end", "Query_start", "Query_end", 
                                                       "Ref_align_len",	"Query_len","Perc_id", "LEN_R", 
                                                       "LEN_Q",	"Ref", "Query")

# from 453091 to 459847, there is 95.62% match
# therefore, below gene is matching well with gaur
# hereford_15_78791037_80120961 454616 455575 ENSBTAG00000035985 . + protein_coding

# from inspection of the ribbon plot, the first 18 genes/pseudo only have two genes matching the hereford
# one is the ENSBTAG00000030920, another is ENSBTAG00000035985. I'm referring to cattle genes here
# see Word doc annotating_genes_divergent_chr15_on_ribbonplot for details

# divergent missing genes in gaur but expanded in hereford
hereford_15_78791037_80120961_diverge <- hereford_15_78791037_80120961[c(2:14,16:18),]

hereford_diverge_genes <- hereford_15_78791037_80120961_diverge$gene_id

hereford_diverge_genes
# [1] "ENSBTAG00000052005" "ENSBTAG00000014513" "ENSBTAG00000046527" "ENSBTAG00000046285" "ENSBTAG00000001291"
# [6] "ENSBTAG00000054791" "ENSBTAG00000053432" "ENSBTAG00000037878" "ENSBTAG00000048469" "ENSBTAG00000051546"
# [11] "ENSBTAG00000048987" "ENSBTAG00000054804" "ENSBTAG00000035986" "ENSBTAG00000050889" "ENSBTAG00000053685"
# [16] "ENSBTAG00000035059"

hereford_diverge_genes_symbol <- c("Uncharacterized protein1","pseudogene1","OR5BE5","OR8I2","OR8H13",
                                   "OR5T21","OR5T1C","OR8J2E","OR8J12","OR8K67",
                                   "OR8K3B","OR8K64","OR8K5B","OR8J19","OR8J16",
                                   "OR8J3F")

#brahman
brahman_15_79129754_80521267 <- read_tsv("brahman_15_79129754_80521267.bed", col_names = FALSE)

colnames(brahman_15_79129754_80521267) <- 
  c("selectedRegionName","newstart","newend","gene_id","score","strand","gene_biotype")

brahman_15_79129754_80521267_arranged <- brahman_15_79129754_80521267 %>% arrange(newstart)

BRAHMAN_VS_GAUR_divergeChr15.delta.coords <- read_tsv("BRAHMAN_VS_GAUR_divergeChr15.delta.coords",
                                                       col_names = FALSE)

names(BRAHMAN_VS_GAUR_divergeChr15.delta.coords) <- c("Ref_start",	"Ref_end", "Query_start", "Query_end", 
                                                       "Ref_align_len",	"Query_len","Perc_id", "LEN_R", 
                                                       "LEN_Q",	"Ref", "Query")

# from 249736 to 256467, from 259029 to 265383, these two regions have more than 96% match
# therefore, below gene is matching well with gaur
# brahman_15_79129754_80521267 254846 256388 ENSBIXG00005002614 . + protein_coding

# skipping ENSBIXG00005021668, which is the FADS2P1 leftmost in ribbon plot, then ENSBIXG00005002614
# brahman divergent region is below

brahman_15_79129754_80521267_arranged_diverge <- brahman_15_79129754_80521267_arranged[c(2:5,7:9),]

brahman_diverge_genes <- brahman_15_79129754_80521267_arranged_diverge$gene_id

brahman_diverge_genes
# [1] "ENSBIXG00005021709" "ENSBIXG00005007373" "ENSBIXG00005021772" "ENSBIXG00005021792" "ENSBIXG00005007839"
# [6] "ENSBIXG00005003679" "ENSBIXG00005007613"

brahman_diverge_genes_symbol <- c("Uncharacterized protein1","pseudogene1","OR8H1-like","pseudogene2","Uncharacterized protein2",
                                  "Uncharacterized protein3","FADS2P1")

##### Getting cattle gene atlas info #####
all_hereford_genes_in_divergent <- hereford_15_78791037_80120961$gene_id
all_hereford_genes_in_divergent_atlas <- c("y","na","pseudo","y","y",
                                            "y","na","na","y","na",
                                            "na","na","na","y","y",
                                            "na","na","y","y","na",
                                            "pseudo","na","na","na","na",
                                            "na","y","y","na","na",
                                            "na","y","y","pseudo","na",
                                            "y","na","y","y","y",
                                            "pseudo","y","y","na","na",
                                            "na","y","y","y")

all_hereford_genes_in_divergent_DF <- 
  as.data.frame(cbind(all_hereford_genes_in_divergent,all_hereford_genes_in_divergent_atlas),
                stringsAsFactors = FALSE)

#make the 49 genes in this divergent region full of gene symbols
hereford_15_78791037_80120961$gene_id %in% hereford_diverge_genes

hereford_diverge_genes_symbol_all <- c("FADS2P1_a","Uncharacterized protein1","OR8H16P","OR5BE5","OR8I2",
                                       "OR8H13","OR5T21","OR5T1C","OR8J2E","OR8J12",
                                       "OR8K67","OR8K3B","OR8K64","OR8K5B","OR8K1",
                                       "OR8J19","OR8J16","OR8J3F","OR5AL8","OR8K62",
                                       "OR8K65P","OR8K66","OR5T2","Uncharacterized protein2","OR9G4",
                                       "OR8K60","OR8K5","Uncharacterized protein3","OR8J17","Uncharacterized protein4",
                                       "OR8J3","OR8U1","OR8U9","OR5AL9P","OR5AL2",
                                       "OR8U3","OR5M3","OR5M11","OR5M13D","OR5M10",
                                       "OR5AP2","OR5AR1","OR2AH1","OR9G3","OR9G4D",
                                       "Uncharacterized protein5","Uncharacterized protein6","OR5G34","FADS2P1_b")

hereford_15_78791037_80120961$gene_id

all_hereford_genes_in_divergent_DF$symbol <- hereford_diverge_genes_symbol_all

##### Find all genes in hereford chr15 77.4 MB – 80.6 MB #####

#hereford
ARSUCD1_2_chr15_77400000_80600000_geneatlas <- 
  read_tsv("/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/gene_atlas_search/ARSUCD1_2_chr15_77400000_80600000_geneatlas.bed", 
           col_names = FALSE)

colnames(ARSUCD1_2_chr15_77400000_80600000_geneatlas) <- 
  c("selectedRegionName","newstart","newend","gene_id","score","strand","gene_biotype")

ARSUCD1_2_chr15_77400000_80600000_geneatlas$no <- 1:nrow(ARSUCD1_2_chr15_77400000_80600000_geneatlas)

#need to add all gene symbols to ARSUCD1_2_chr15_77400000_80600000_geneatlas
write_tsv(ARSUCD1_2_chr15_77400000_80600000_geneatlas,"/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/zoom_in_OR_region/ARSUCD1_2_chr15_77400000_80600000_geneatlas.tsv")

# read in all gene atlas in the region flanking the divergent region
path_to_atlas <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/gene_atlas_search/"

# test <- read_tsv(paste0(path_to_atlas,"ENSBTAG00000048141.csv"))

#create a list to read all *cov file
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

FPKM_mean$no <- rep(0,nrow(FPKM_mean))

for (j in 1:nrow(FPKM_mean)){
  selc <- FPKM_mean$prefix[j] == ARSUCD1_2_chr15_77400000_80600000_geneatlas$gene_id
  the_number <- ARSUCD1_2_chr15_77400000_80600000_geneatlas$no[selc]
  FPKM_mean$no[j] <- the_number
}

FPKM_mean <- FPKM_mean %>% arrange(no)

sperm_only <- FPKM_mean %>% filter(Tissue == "Sperm")

plot(sperm_only$no,sperm_only$mean_FPKM, type = "l")

#write out the results for heat map
write_tsv(FPKM_mean,"/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/gene_atlas_search/FPKM_mean.tsv",
          col_names = TRUE)

##### heatmap #####
library(ggplot2)
library(viridis)  

# Read the file
# setwd("~/Documents/Kelly_2021/Research/Lloyd_Gaur/Cooperation")
ARSUCD1_2_chr15_77400000_80600000_geneatlas <- 
  read_tsv("/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/gene_atlas_search/ARSUCD1_2_chr15_77400000_80600000_geneatlas.bed", col_names = FALSE)
FPKM_mean <- 
  read_tsv("/Users/lloyd/Documents/lloyd_2019/Research/gaur/divergent_region/chr15_cattle_species/gene_atlas_search/FPKM_mean.tsv", col_names = TRUE)

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
tiff(filename = "heatmap_expanded_divergent_region.tiff", width = 700, height = 700)
ggplot(FPKM_mean, aes(x = prefix, y = gsub("_","",Tissue), fill= mean_FPKM)) + 
  geom_tile() + 
  xlab("Gene") + 
  ylab("Tissue")+
  labs(fill = "FPKM") +
  scale_fill_gradientn(colours = c("white", "blue"),
                       values = scales::rescale(c(0,100,250))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y=element_text(size=7))
dev.off()

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
