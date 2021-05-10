#------------------------------------------------------
# Program name: gaur_cafe_filter_for_gaur_branch.R
# Objective: Filter the cafe results for families gain or loss on gaur branch
#           
# Author: Kelly Ren
# Email add: kellydren@gmail.com
#------------------------------------------------------

# Loading packages

library(magrittr)
library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(gtools)
library(tibble)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(GenomicFeatures)
library(biomaRt)


# Find the branch of Bgau
## Load the result from cafe


# there are 7395 families
report_runmodel1_error <- read.table("input/From_CAFE/report_runmodel1_error.cafe", skip = 11)
dim(report_runmodel1_error)



#report_runmodel1_error[,c(1,3,4)]%>%
#  write_csv("report_runmodel1_error.csv")
# can manully edit the csv file and save it to txt or run R script
# report_runmodel1_error <- read_delim("report_runmodel1_error.txt", col_names = F,delim = " ", skip = 1)%>%
#   as.data.frame()%>%
#   set_colnames(c("PTHR_ID","Family_p_value", "Oari", "Chir", "Oari_Chir", "Other_Bbub", "Hbta", "Hbin", "Hbta_Hbin", "Bgau","Other_Bgau", "Bmut", "Other_Bmut","Bbub", "except_pig_human", "Sscr","except_human", "Hsap"))

report_runmodel1_error <- report_runmodel1_error[,c(1,3,4)]

report_runmodel1_error$V4 <- report_runmodel1_error$V4%>%
  gsub("[()]"," ",.)

report_runmodel1_error %<>% mutate_each(funs(str_replace_all(., "\"", "")))

report_runmodel1_error <- report_runmodel1_error[,c(1:2)]%>%cbind(read.table(text = report_runmodel1_error$V4, sep = ",", colClasses = "character"))%>%
  as.data.frame()%>%
  set_colnames(c("PTHR_ID","Family_p_value", "Oari", "Chir", "Oari_Chir", "Other_Bbub", "Hbta", "Hbin", "Hbta_Hbin", "Bgau","Other_Bgau", "Bmut", "Other_Bmut","Bbub", "except_pig_human", "Sscr","except_human", "Hsap"))

# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): (0,2) (1,11) (Hbta,Hbin) (Hbta_Hbin,Bgau) (7,10) (9,12) (3,14) (13,16) 

head(report_runmodel1_error)



Bgau_report_runmodel1_error <- subset(report_runmodel1_error, Bgau < 0.05)

Hbta_report_runmodel1_error <- subset(report_runmodel1_error, Hbta < 0.05)

Hbin_report_runmodel1_error <- subset(report_runmodel1_error, Hbin < 0.05)



# number of sig families in Gaur
Bgau_report_runmodel1_error$PTHR_ID%>%unique()%>%length()

# number of sig families in Hbta
Hbta_report_runmodel1_error$PTHR_ID%>%unique()%>%length()

# number of sig families in Hbin
Hbin_report_runmodel1_error$PTHR_ID%>%unique()%>%length()
 


Select_family <- Bgau_report_runmodel1_error[,c("PTHR_ID","Bgau")]%>%
  set_colnames(c("FAMILY", "Bgau_pvalues"))

dim(Select_family)
hist(as.numeric(Select_family$Bgau_pvalues), main = "p-value distribution")

### Compare to all species 


# read the count table
com_cafe_input_data <- read_delim("Output/CAFE/filtered_cafe_input_all.txt", delim = "\t")%>%as.data.frame()
dim(com_cafe_input_data)

Bgau_com_cafe_input_data <- subset(com_cafe_input_data, FAMILY %in% Select_family$FAMILY)

# In all, 112 no gaur sequence
(Bgau_com_cafe_input_data$Bgau <= 0)%>%table()

Bgau_com_cafe_input_data$cattle_ave <- Bgau_com_cafe_input_data[,c(5,6)]%>%apply(MARGIN = 1,FUN = function(x){sum(x)/2})

# Compare to cattle 
(Bgau_com_cafe_input_data$Bgau < Bgau_com_cafe_input_data$cattle_ave)%>%table()

#FALSE  TRUE 
#  258   129  

# Save the families for 244 gains
Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau > Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%
  write_csv("Output/CAFE/Bgaur_sig_families_gains.csv")

# Save the families for 14 equal
Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau == Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%
  write_csv("Output/CAFE/Bgaur_sig_families_equal.csv")

# Save the families for 129 losses
Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau < Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%
  write_csv("Output/CAFE/Bgaur_sig_families_losses.csv")


### Compare to all cattles 
#### losses

Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau < Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%left_join(Bgau_com_cafe_input_data[,c(2:11)], by= c("FAMILY"= "FAMILY"))%>%write_csv("Output/CAFE/counts_sig_losses_gaur_branch.csv")



branch_losses <- Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau < Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%left_join(Bgau_com_cafe_input_data, by= c("FAMILY"= "FAMILY"))

branch_losses[abs(branch_losses$Bgau - branch_losses$cattle_ave) > 1,]


#### gains

Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau > Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%
  left_join(Bgau_com_cafe_input_data[,c(2:8)], by= c("FAMILY"= "FAMILY"))%>%write_csv("Output/CAFE/counts_sig_gains_gaur_branch.csv")

branch_gains <- Select_family[,c("FAMILY", "Bgau_pvalues")]%>%
  subset(FAMILY %in% c(Bgau_com_cafe_input_data[(Bgau_com_cafe_input_data$Bgau > Bgau_com_cafe_input_data$cattle_ave),]%>%extract2("FAMILY")))%>%
  arrange(Bgau_pvalues)%>%
  left_join(Bgau_com_cafe_input_data, by= c("FAMILY"= "FAMILY"))

branch_gains[abs(branch_gains$Bgau - branch_gains$cattle_ave) > 1,]


