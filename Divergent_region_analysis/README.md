# Divergent_region_analysis
Compare genomic region flanking the divergent region in gaur, Brahman and Hereford.
Then determine which genes in this region has gene expression detected in [Cattle Gene Atlas](http://cattlegeneatlas.roslin.ed.ac.uk). Finally plot a heatmap for all genes in this region. I have also included scripts to plot the heat map of gene expression of SLC and lysozyme c families. The script(s) below are used for this part.

gaur_divergent_regions_genes.R
gaur_find_Btau_gene_atlas.R

Check whether the odorant receptors surrounding the divergent region identified in Hereford overlaps with any [QTL](https://www.animalgenome.org/cgi-bin/QTLdb/index). The script(s) below are used for this part.

gaur_divergent_regions_arsucd1_2_qtl.R

Since I reverse complemented gaur and Brahman sequences to compare with Hereford, I need to change the gene annotation coordinates to match the reverse complemented sequences. The script(s) below are used for this part.

gaur_annotation_gtf_flipping_coordinates.R

Since the gene symbol of the gaur is not available in current Ensembl annotation release 101, reciprocal best blast hit was used to find corresponding gene symbols in cattle. The script(s) below are used for this part.

gaur_reciprocal_best_hit_cattle.R
