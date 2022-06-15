# GeneLen/
# Investigating gene length distribution
# First plot CDF of baseline mutation rates => GeneLen/MutRateCDFs.png
Rscript GeneLen.R
# Then, for known ASD/NDD genes, extract their mutation rates (quantiles), CAFs, and effect sizes
# they will be used as illustrating examples for given parameters

# Note on: Mutation rates 
# For mutation rate, we still use baseline LoF rate 
# We did not use re-annotated version because it is only possible to do this annotation
# on nonsense and frameshift variants.


# Selection coefficient
# This can be easily derived from s = mu/CAF
# We compared mu, CAF in reference populations, 
# also compare estimated sel coeff bin with fraction of de novos observed in ASD cases

# We first calculate point estimates for known ASD genes
# they are selected by denovo enrich P<1e-4 and Mega P<0.05 for known genes or novel genes
Rscript SelCoeffs.R

# Draw heatmaps of power and sample size for future studies
Rscript Heatmap.R
