# What to do in between generating a counts table and running DESeq2
# Created Michael Studivan (mstudiva@fau.edu) 

#------------------------------
# If your experimental design was complicated, or your library preps included samples from multiple projects, you will need
# to modify your allcounts table to work with differential analysis packages

# FIRST, AND MOST IMPORTANTLY
# Rename the allcounts.txt file to allcounts_raw.txt, then open allcounts.txt with Excel and save as allcounts.xlsx
# ALWAYS keep an unmodified backup of all data

# Using the order of the sample columns, create an design.xlsx table to include the factor levels in the SAME order
# Hopefully your sample IDs sort in a way that makes sense to you
# Do note: the allcounts matrix will be organized with genes as rows and samples as columns, while the design matrix will be samples as rows
# and factors as columns
# Use =IF(cell=cell, "", "NO") statements to make sure your allcounts and design id labels match up perfectly

# Also double check that each sample has counts for at least one gene by using the =max() formula at the bottom of the counts matrix
# If any samples have all zeroes, remove them and make sure to remove the corresponding rows in the design file

# Save the counts as allcounts.txt and the design table as design.csv