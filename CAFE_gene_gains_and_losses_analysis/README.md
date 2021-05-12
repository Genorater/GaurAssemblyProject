# Analysis of gaur gene gains/losses

## Table of Content
Inside the Buffalo folder, you will find the script in R markdown format and associated input and output files. As some scripts have dependencies on output from other scripts, it may be best to run the script in the order specified below.
### 1. CAFE_Match_gene_family.Rmd
  Group genes in the othrogroups into PANTHER gene families, and make a count table for gene families as CAFE input.
### 2. CAFE_filter_buffalo_branch.Rmd
  Reading CAFE output and select the gene families with significant gene gains/losses for the buffalo branch.
### 3. CAFE_GO_KEGG_Reactome_pathway.Rmd
  Importing CAFE results to perform GO, KEGG and Reactome enrichment analysis based on human and cattle annotations as the references.
