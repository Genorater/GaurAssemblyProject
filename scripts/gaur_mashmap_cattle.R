#------------------------------------------------------
# Program name: gaur_mashmap_cattle.R
# Objective: check mashmap of gaur against cattle 
#           ARS-UCD1.2
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)
library(xlsx)
library(reshape2)
library(ggplot2)

# reading mashmap_gaur_vs_cattle.mashmap
gaur_vs_cattle_mashmap <- 
  read_delim("/Users/lloyd/Documents/lloyd_2019/Research/gaur/ARS_UOA_Gaur_1_align_to_cattle/mashmap_ARSUCD1.2/mashmap_gaur_vs_cattle.mashmap",
             " ", col_names = FALSE, col_types = list(col_character(), col_integer(),col_integer(),
                                                      col_integer(),col_character(),col_character(),
                                                      col_integer(),col_integer(),col_integer(),
                                                      col_double()))

names(gaur_vs_cattle_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                              "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")

#only analyse those in chr
gaur_chr_selc <- gaur_vs_cattle_mashmap$Query_name %in% c("1","2","3","4","5","6","7","8","9","10",
                                             "11","12","13","14","15","16","17",
                                             "18","19","20","21","22","23","24",
                                             "25","26","27","28","29","X")

gaur_vs_cattle_mashmap <- gaur_vs_cattle_mashmap[gaur_chr_selc,]

# cattle_chr_selc <- gaur_vs_cattle_mashmap$Ref_name %in% c("1","2","3","4","5","6","7","8","9","10",
#                                        "11","12","13","14","15","16","17",
#                                        "18","19","20","21","22","23","24",
#                                        "25","26","27","28","29","X")

# create alignL var which has the Ref (in this case cattle ARS-UCD1.2) alignment length
# need to add 1 in alignL calc bcos it is 0-based index
gaur_vs_cattle_mashmap_modi <- gaur_vs_cattle_mashmap %>% 
  arrange(Ref_name,Ref_start) %>% mutate(alignedRefL = Ref_end - Ref_start +1) %>% 
  mutate(alignedQueryL = Query_end - Query_start +1) %>%
  mutate(proportion_alignedRefL = alignedRefL/Ref_length) %>%
  mutate(proportion_alignedQueryL = alignedQueryL/Query_length)

summary(gaur_vs_cattle_mashmap_modi$Percentid)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 80.97   93.03   96.70   94.99   98.06  100.00

gaur_vs_cattle_mashmap_RefQuery <- gaur_vs_cattle_mashmap_modi %>% 
  group_by(Ref_name,Query_name) %>% 
  summarise(sum_alignRefL = sum(alignedRefL), sum_alignQueryL = sum(alignedQueryL),
            sum_proportion_alignedRefL = sum(proportion_alignedRefL),
            sum_proportion_alignedQueryL = sum(proportion_alignedQueryL),
            new_query_length = mean(Query_length),
            new_ref_length = mean(Ref_length), .groups = 'drop')

# sort by where cattle chr matched the most
gaur_vs_cattle_mashmap_sort <- gaur_vs_cattle_mashmap_RefQuery %>% 
  arrange(desc(sum_proportion_alignedQueryL))

#special case with a lot of not matching seqs in gaur chr19
#Query_name is gaur
gaur_vs_cattle_mashmap_sort_chr19 <- gaur_vs_cattle_mashmap_sort %>% filter(Query_name == "19")

#special case with a lot of not matching seqs in gaur chr29
gaur_vs_cattle_mashmap_sort_chr29 <- gaur_vs_cattle_mashmap_sort %>% filter(Query_name == "29")

#which coordinates in chr29 that are not matching
#the non-matching to cattle chr19 in gaur chr29, where else do they match
#the conclusion from this analysis of gaur chr29 is that the genomic regions that do not match 
#cattle chr 19 are matching elsewhere in the genome in small bits for most of them.
#could the non-matching bits be some kinds of repeats?
#Using Bos_gaurus.ARS_UOA_Gaur_1.dna.toplevel.coor, which has gap counts for all gaur chr, it seems
#that chr 29 is quite bad with 102 gaps.
#the dot plot of cattle chr 19 vs gaur chr 29 shows large regions before 4230000 not matching cattle.
#also, regions 4347500 to 5287500, 7343750 to 9047500 are not matching cattle as well.
#these non-matching regions have lots of gaps
gaur_vs_cattle_mashmap_modi_chr29_to_allcattle <- gaur_vs_cattle_mashmap_modi %>% filter(Query_name == "29") %>% arrange(Query_start)
table(gaur_vs_cattle_mashmap_modi_chr29_to_allcattle$Ref_name)

#restricting to just cattle chr 19
gaur_vs_cattle_mashmap_modi_chr29_tocattle_chr19 <- gaur_vs_cattle_mashmap_modi %>% filter(Query_name == "29") %>%
  filter(Ref_name == "19") %>% arrange(Query_start)

# the top30 rows is the corresponding homologous chr
# in which query aligned proportion in the cattle chr is at least >= 0.25
gaur_vs_cattle_mashmap_sort <- gaur_vs_cattle_mashmap_sort[1:31,]
gaur_vs_cattle_mashmap_sort_modi <- gaur_vs_cattle_mashmap_sort %>% arrange(as.numeric(Ref_name))

# create a supp table that shows perc aligned to cattle and vice versa
gaur_vs_cattle_mashmap_tidied1 <- gaur_vs_cattle_mashmap_sort_modi %>% 
  select(Ref_name,Query_name,sum_proportion_alignedRefL,sum_proportion_alignedQueryL,
         new_query_length,new_ref_length) %>%
  mutate(sum_proportion_alignedRefL = round(sum_proportion_alignedRefL * 100,1)) %>%
  mutate(sum_proportion_alignedQueryL = round(sum_proportion_alignedQueryL * 100,1))

write_excel_csv(gaur_vs_cattle_mashmap_tidied1,
                "/Users/lloyd/Documents/lloyd_2019/Research/gaur/ARS_UOA_Gaur_1_align_to_cattle/mashmap_ARSUCD1.2/gaur_vs_cattle_mashmap_tidied1.csv")

#below is for ploting chr length comparison
#####
dir_plot1 <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/ARS_UOA_Gaur_1_align_to_cattle/mashmap_ARSUCD1.2/"
fullpathnameplot1 <- paste0(dir_plot1,"gaur_vs_cattle_mashmap_tidied1.xlsx")

chr_scaff_synteny <- read.xlsx2(fullpathnameplot1,2,stringsAsFactors = FALSE)

#make cattle chr length = 0 when matched to gaur chr29
# chr_scaff_synteny$Cattle_chromosome_length[29] <- "0"

#create DF suitable for plotting
df = melt(data.frame(gaur=as.numeric(chr_scaff_synteny$Gaur_chromosome_length), 
                     cattle=as.numeric(chr_scaff_synteny$Cattle_chromosome_length), 
                     chromosome=c(1:29,"X")),
          variable.name="species")

#plot barplot
orderer <- c("1","2","3","4","5","6","7","8","9","10",
             "11","12","13","14","15","16","17",
             "18","19","20","21","22","23","24",
             "25","26","27","28","29","X")

df$chromosome <- factor(df$chromosome, levels = orderer)

tiff(filename = "Gaur_vs_cattle_length_comparison.tiff",width = 600, height = 600)
ggplot(df, aes(chromosome, value, fill=species)) + 
  geom_bar(position="dodge",stat="identity") +
  labs(x="Gaur chromosome and corrresponding cattle chromosome",
       y="Chromosome length (bp)")
dev.off()
