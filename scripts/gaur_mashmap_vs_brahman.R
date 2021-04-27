#------------------------------------------------------
# Program name: gaur_mashmap_vs_brahman.R
# Objective: Does the gaur scaffolds PB+HIC correspond 
#           well to cattle scaffolds?
# Author: Lloyd Low
# Email add: lloydlow@hotmail.com
#------------------------------------------------------

library(readr)
library(dplyr)

# directory containing mashmap output
dir1 <- "/Users/lloyd/Documents/lloyd_2019/Research/gaur/mashmap_gaur_brahman_genome_correspondence/"

# path to filename
path1 <- paste0(dir1,"gaur_PB_HIC_vs_brahman_frozen_as_ref.split.scaffold_only.mashmap")

# reading gaur_PB_HIC_vs_brahman_frozen_as_ref.split.scaffold_only.mashmap
gaur_vs_brahman_mashmap <- read_delim(path1," ", col_names = FALSE)
names(gaur_vs_brahman_mashmap) <- c("Query_name","Query_length","Query_start","Query_end","Orientation",
                                                     "Ref_name","Ref_length","Ref_start","Ref_end","Percentid")


gaur_vs_brahman_mashmap_modi <- gaur_vs_brahman_mashmap %>% 
  arrange(Ref_name,Ref_start) %>% mutate(alignedL = Ref_end - Ref_start +1) %>%
  mutate(alignProportion = alignedL/Query_length)

#sort by query length
gaur_vs_brahman_mashmap_modi_sortLength <- gaur_vs_brahman_mashmap_modi %>% arrange(desc(Query_length))

query_name <- unique(gaur_vs_brahman_mashmap_modi_sortLength$Query_name)

# loop not working, shut it for now 
for (i in 1:length(query_name)){
  df <- gaur_vs_brahman_mashmap_modi_sortLength %>% filter(Query_name == query_name[i])
  df2 <- df %>% group_by(Ref_name) %>%
    summarise(total_align_prop = sum(alignProportion))
  assign(query_name[i],df2)
}

#Investigate the 2 scaffolds one-by-one manually as a test

query_name[1]
# PGA_scaffold0__82_contigs__length_183742024
gaur_vs_brahman_mashmap_modi_sortLength %>% 
  filter(Query_name == query_name[1]) %>%
  group_by(Ref_name) %>% 
  summarise(total_align_prop = sum(alignProportion))

query_name[2]
# PGA_scaffold1__49_contigs__length_152753188
gaur_vs_brahman_mashmap_modi_sortLength %>% 
  filter(Query_name == query_name[2]) %>%
  group_by(Ref_name) %>% 
  summarise(total_align_prop = sum(alignProportion))

#Investigate the 30 scaffolds one-by-one 
query_name[1]
PGA_scaffold0__82_contigs__length_183742024
#PGA_scaffold0__82_contigs__length_183742024 is 73% chr2 and 25% chr28

query_name[2]
PGA_scaffold1__49_contigs__length_152753188
#PGA_scaffold1__49_contigs__length_152753188 is 99% chr1

query_name[3]
PGA_scaffold2__156_contigs__length_146132534
#PGA_scaffold2__156_contigs__length_146132534 is 95% chrX

query_name[4]
PGA_scaffold4__127_contigs__length_123593481
#PGA_scaffold4__127_contigs__length_123593481 is 96% chr3

query_name[5]
PGA_scaffold3__140_contigs__length_121504598
#PGA_scaffold3__140_contigs__length_121504598 is 96% chr4

query_name[6]
PGA_scaffold5__64_contigs__length_120493815
#PGA_scaffold5__64_contigs__length_120493815 is 98% chr5

query_name[7]
PGA_scaffold6__71_contigs__length_119473097
#PGA_scaffold6__71_contigs__length_119473097 is 97% chr6

query_name[8]
PGA_scaffold7__21_contigs__length_113825074
#PGA_scaffold7__21_contigs__length_113825074 is 99% chr8

query_name[9]
PGA_scaffold8__58_contigs__length_112566739
#PGA_scaffold8__58_contigs__length_112566739 is 96% chr7

query_name[10]
PGA_scaffold9__54_contigs__length_108521854
#PGA_scaffold9__54_contigs__length_108521854 is 97% chr11

query_name[11]
PGA_scaffold10__98_contigs__length_104474243
#PGA_scaffold10__98_contigs__length_104474243 is 96% chr10

query_name[12]
PGA_scaffold11__27_contigs__length_103814667
#PGA_scaffold11__27_contigs__length_103814667 is 99% chr9

query_name[13]
PGA_scaffold12__125_contigs__length_91620959
#PGA_scaffold12__125_contigs__length_91620959 is 91% chr15 and 2.5% chr1

query_name[14]
PGA_scaffold13__96_contigs__length_89800262
#PGA_scaffold13__96_contigs__length_89800262 is 93% chr12

query_name[15]
PGA_scaffold14__28_contigs__length_84794504
#PGA_scaffold14__28_contigs__length_84794504 is 98% chr13

query_name[16]
PGA_scaffold15__32_contigs__length_82824024
#PGA_scaffold15__32_contigs__length_82824024 is 98% chr14

query_name[17]
PGA_scaffold16__27_contigs__length_81199330
#PGA_scaffold16__27_contigs__length_81199330 is 98% chr16

query_name[18]
PGA_scaffold17__34_contigs__length_73199797
#PGA_scaffold17__34_contigs__length_73199797 is 97% chr17

query_name[19]
PGA_scaffold18__201_contigs__length_72908606
#PGA_scaffold18__201_contigs__length_72908606 is 78% chr19, 3.8% chr1, 2.5% chr23 and 2.7% chr29

query_name[20]
PGA_scaffold19__6_contigs__length_71077695
#PGA_scaffold19__6_contigs__length_71077695 is 99% chr20

query_name[21]
PGA_scaffold20__22_contigs__length_67692141
#PGA_scaffold20__22_contigs__length_67692141 is 98% chr21

query_name[22]
PGA_scaffold22__99_contigs__length_67640060
#PGA_scaffold22__99_contigs__length_67640060 is 89% chr18, 1.3% chr27 and 1.3% chrX

query_name[23]
PGA_scaffold21__86_contigs__length_66768343
#PGA_scaffold21__86_contigs__length_66768343 is 93% chr24, 1.3% chr10 and 1.6% chr3

query_name[24]
PGA_scaffold23__15_contigs__length_61205293
#PGA_scaffold23__15_contigs__length_61205293 is 98% chr22

query_name[25]
PGA_scaffold24__189_contigs__length_60321932
#PGA_scaffold24__189_contigs__length_60321932 is 84% chr23, 2.3% chr10, 1.2% chr4 and 1.6% chrX

query_name[26]
PGA_scaffold25__101_contigs__length_54068254
#PGA_scaffold25__101_contigs__length_54068254 is 92% chr26

query_name[27]
PGA_scaffold26__44_contigs__length_49959101
#PGA_scaffold26__44_contigs__length_49959101 is 97% chr29

query_name[28]
PGA_scaffold27__70_contigs__length_46337533
#PGA_scaffold27__70_contigs__length_46337533 is 93% chr27

query_name[29]
PGA_scaffold28__26_contigs__length_42418488
#PGA_scaffold28__26_contigs__length_42418488 is 98% chr25

query_name[30]
PGA_scaffold29__121_contigs__length_12941853
#PGA_scaffold29__121_contigs__length_12941853 is 50% chr19, 2.6% chr10, 3.7% chr20, 3.1% chr21, 4.7% chr26, and 1.4% chr9
