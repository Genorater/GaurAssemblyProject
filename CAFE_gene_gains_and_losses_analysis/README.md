# Analysis of gaur gene gains/losses

CAFE was used to do the gene gains and losses anlaysis for 9 selected species inlcuding gaur.
The analysis was set to focus on the chances on gaur branch.
Three scripts in the scripts folder are related to gene gains/losses anlysis:

### 1. gaur_cafe_match_gene_famlies.R
  Group genes in the othrogroups into PANTHER gene families, and make a count table for gene families as CAFE input.
### 2. gaur_cafe_filter_for_gaur_branch.R
  Reading CAFE output and select the gene families with significant gene gains/losses (p-value < 0.05) for the gaur branch.
### 3. gaur_cafe_go_kegg_pathway_analysis.R
  Importing CAFE results to perform GO and KEGG pathway based on cattle annotations as the references.



