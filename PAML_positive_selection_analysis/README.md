# Analysis results of positively selected sites from PAML

PAML was used to select the positive selected sites in gaur among 9 selected species.
Two scripts in the scripts folder are related to positive selection anlysis:

### 1. gaur_paml_branch_positive_selection.R
Read the log Likelihood vaules from PAML for Null and Alternative test model. 
By calculating the [chi-squared distribution](https://www.researchgate.net/publication/254333395_Continuous_Univariate_Distributions_Volume_1)of the 2(NULL_logLikelihood  - Alternative_logLikelihood), the p-values were calculated for both sites.
The produced p-values were devied for one side.
### 2. gaur_paml_go_kegg_pathway_analysis.R
This script collected the cattle represented genes from the positive selection results to apply GO and KEGG pathway analysis.
### 3. gaur_find_genes_in_GO_terms.R
This script can find the significant genes that related to the GO terms. 
This script has been used for GO results of both positive selection and gene gains/losses.

