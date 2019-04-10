# tag-based RNA-Seq analysis, version Apr 10, 2019
# Created by Michael Studivan (mstudiva@fau.edu)

#------------------------------
# Create a "raw" directory for all your read counts and associated design spreadsheets

# If your experimental design was complicated, or your library preps included samples from
# multiple projects, you will need to modify your allcounts table to work with
# differential analysis packages

# FIRST, AND MOST IMPORTANTLY
# Duplicate the allcounts.txt file to allcounts_raw.txt, then open allcounts.txt with
# Excel and save as allcounts.xlsx
# ALWAYS keep an unmodified backup of all data

# Using the order of the sample columns, create an design.xlsx table to include the factor
# levels in the SAME order
# Hopefully your sample IDs sort in a way that makes sense to you
# Do note: the allcounts matrix will be organized with genes as rows and samples as
# columns, while the design matrix will be samples as rows and factors as columns
# Use =IF(cell1=cell2, "", "NO") statements to make sure your allcounts and design id
# labels match up perfectly

# If you have a concatenated transcriptome (host and algal symbiont, as in the M. cavernosa and
# Cladocopium sp. transcriptome), create separate datasheets by species-specific isogroups,
# which should be labeled differently if you followed the pipeline in Mcav-Annotated-Transcriptome

# Also double check that each sample has counts for at least one gene by using the =max()
# formula at the bottom of the counts matrix
# If any samples have all zeroes, remove them and make sure to remove the corresponding
# rows in the design file

# Save the counts as allcounts.txt and the design table as design.csv

#------------------------------
# Quantifying Differentially-Expressed Genes (DEGs) with DESeq2

# Create a "DESeq2" directory

# Copy allcounts.txt and design.csv files to the same working directory as your DESeq2.R
# script
# Before running the R script, go through and determine how the experimental factors are
# coded in the template experiment
# Use Find & Replace to change them to your applicable factors
# Single factor experiments will require more work to get the code to run, as the template
# is for two factor
# As a two factor design, the first factor is your batch factor (like site), and serves as
# replicates for your second factor, or treatment
# In the R script example, site is the batch factor and depth is the treatment
# Ideally, your treatment factor should only have TWO levels (a la control, treatment)

# Depending on your experiment size, the rlog transformation may be more suitable than vsd
# Larger experiments will take too long with rlog, and should use vsd instead

# Anytime an .RData package is created, you can safely start a new R session as needed
# The "load(name.RData)" function will bring all the associated dataframes back into R

# The details for each step are in the deseq2.R script, but the overview is as follows:
# 1) reading in data and experiment design
# 2) prefiltering genes for WGCNA
# 3) making a DESeq2 dataset and running data transformations for normalized expression
# 4) running outlier detection, and amending the dataset if needed
# 5) creating heatmaps and PCoAs to visualize sample-to-sample differences
# 6) ANOVA and Discriminant Analysis of Principal Components (DAPC)
# 7) differential expression models with DESeq2
# 8) Venn diagrams of DEGs
# 9) exporting differentially expressed genes

#------------------------------
# Gene Ontologeny (GO) Analysis with Mann-Whitney U-Test

# Create a "GO-MWU" directory

# Download Misha Matz's Github repository from:
git@github.com:z0on/GO_MWU.git
# Follow the instructions in his README
# In Step 1, the GO hierarchy file is go.obo
# the table of GO annotations is mcav_iso2go.tab, found in my Github repository:
https://github.com/mstudiva/Mcav-Annotated-Transcriptome.git
# the -log(p-value) tables were generated at the end of the DESeq2 script
# you should have one _lpv.csv for each of your experimental factors, and one for the
# interaction
# you can also compare fc datasets as well for factors with corresponding dataframes

#------------------------------
# EuKaryotic Orthologous Groups (KOG) Analysis with Mann-Whitney U-Test

# This analysis is meant to compare enrichment of universal KOG annotations
# through gene expression across groups
# For example, subjecting different groups of corals to the same kind of treatment
# Take a look at the example found at ?KOGMWU, describing a heat stress experiment with
# adults and two age groups of larvae

# Create a "KOG-MWU" directory

# Copy mcav_iso2kogClass.tab into the directory
# It is found in my GitHub repository:
https://github.com/mstudiva/Mcav-Annotated-Transcriptome.git

# Use the template included in this repository
# Replace your datasets and factor names with those in your experiment
# The template follows 3 analyses: full model, effect of site, depth, and interaction on
# expression; comparison of shallow vs mesophotic expression across 4 sites; and comparison
# of conserved genes across depths among sites

#------------------------------
# Weighted Correlation Network Analysis (WGCNA)

# This analysis explores correlations between gene expression across treatment groups and
# other measured response variables (physiological monitoring, growth rates, really anything)
# First, genes are clustered into modules based on similar expression patterns
# Then, module membership is correlated to your other response metrics to identify
# interesting patterns
# For example, you could identify if the gene expression of a particular module is
# associated with increased growth rate, or if a particular genotype shows different
# expression patterns

# Create a "WGCNA" directory

# Copy the 3 template R scripts into the working directory
# s.threshold.signed.R, tom_calc_signed.R, wgcna_MEanalysis_signed_heavy.R
# Copy the data4wgcna.RData package from the "DESeq2" directory to "WGCNA"
# Copy mcav_iso2gene.tab included in this repository to "WGCNA"

# Run the R script as directed

#------------------------------
# heatmaps

# Once you have any trends in gene expression across your factors identified using the
# above tests, you can design heatmaps to visually represent the differences
# Use the heatmapEveryWhichWay R script to create custom heatmaps by: GO analysis, KOG
# analysis, gene name, top significant DEGs, WGCNA membership, etc

# Create a "heatmaps" directory
# Copy the mcav_iso2gene.tab and mcav_iso2kogClass.tab included in this repository to
# "heatmaps"
# Copy the data input files as needed for each respective section of code into "heatmaps"
