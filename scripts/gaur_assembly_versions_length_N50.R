#------------------------------------------------------
# Program name: gaur_assembly_versions_length_N50.R
# Objective: Analysis of contig/scaffold lengths and N50 
#          of various gaur assemblies
# Author: Lloyd
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(ggplot2)
library(dplyr)

#My own N50 function
N50 <- function(assembly_df){
  assembly_df <- assembly_df %>% arrange(desc(length))
  vect <- c()
  for (i in 1:nrow(assembly_df)){
    N50length <- 0.5*sum(as.numeric(assembly_df$length))
    vect <- c(vect,assembly_df$length[i])
    cumsumlength <- sum(as.numeric(vect))
    if (cumsumlength >= N50length){
      print(assembly_df$length[i])
      break
    } 
  }
}

###############################################################################

# reading tab files of read lengths
# directory containing read lengths
dir1 <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/assembly_versions/"

# path to filename
path1 <- paste0(dir1,"gaur_contigs.arrow2.rl")
path2 <- paste0(dir1,"PGA_assembly-002.rl")

# PB gaur assembly
gaur_PB <- read_tsv(path1,col_names = FALSE)
gaur_PB <- gaur_PB[c(1,3,4,5)]
names(gaur_PB) <- c("name","gap","length","perc_gap")

#sort for length, longest at the top row
gaur_PB <- gaur_PB %>% arrange(desc(length))

# create var for gap /100 to count gaps
gaur_PB <- gaur_PB %>% mutate(no_of_gap = gap/100)

sum(gaur_PB$length)
#2717269665

N50(gaur_PB)
#13786623/1e6 = 13.8 Mb 

# PB + HIC gaur assembly
gaur_PB_HIC <- read_tsv(path2,col_names = FALSE)
gaur_PB_HIC <- gaur_PB_HIC[c(1,3,4,5)]
names(gaur_PB_HIC) <- c("name","gap","length","perc_gap")

#sort for length, longest at the top row
gaur_PB_HIC <- gaur_PB_HIC %>% arrange(desc(length))

# create var for gap /100 to count gaps
gaur_PB_HIC <- gaur_PB_HIC %>% mutate(no_of_gap = gap/100)

#all scaffolds
sum(gaur_PB_HIC$length)
#2717496547

sum(gaur_PB_HIC$no_of_gap)
#2269

#chromosome-level or more than 1 Mb scaffolds
sum(gaur_PB_HIC$length[1:30])
#2687673489

sum(gaur_PB_HIC$no_of_gap[1:30])
#2269

#only ~1.1% of bases not included in scaffolds
(2717496547 - 2687673489)/2717496547
#[1] 0.01097446

N50(gaur_PB_HIC)
#104474243/1e6 = 104.5 Mb 

########################################################################
# 20200624
# After Ben and Derek polished assembly

# path to filename
path3 <- paste0(dir1,"ARS-UOA_Gaur_PB_HiC_Arrow_p2.rl")

# PB gaur polished assembly
gaur_PB_polished <- read_tsv(path3,col_names = FALSE)
gaur_PB_polished <- gaur_PB_polished[c(1,3,4,5)]
names(gaur_PB_polished) <- c("name","gap","length","perc_gap")

#sort for length, longest at the top row
gaur_PB_polished <- gaur_PB_polished %>% arrange(desc(length))

sum(gaur_PB_polished$length)
#2717364261

N50(gaur_PB_polished)
#104472007/1e6 = 104.5 Mb

#percentage of bases in unplaced
unplaced_length <- sum(gaur_PB_polished$length) - sum(gaur_PB_polished$length[1:30])

(unplaced_length / sum(gaur_PB_polished$length))*100
#1.1%
#accounting for mito, still ~1.1%
#(29829050-16346)/2717364261
#[1] 0.01097118

# path to filename
path4 <- paste0(dir1,"ARS_UOA_Gaur_1.rl")

# PB gaur polished assembly, with mito fixed and chr 19 renamed
gaur_PB_polished_mito <- read_tsv(path4,col_names = FALSE)
gaur_PB_polished_mito <- gaur_PB_polished_mito[c(1,3,4,5)]
names(gaur_PB_polished_mito) <- c("name","gap","length","perc_gap")

#sort for length, longest at the top row
gaur_PB_polished_mito <- gaur_PB_polished_mito %>% arrange(desc(length))

sum(gaur_PB_polished_mito$length)
#2717352581
