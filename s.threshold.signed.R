library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()
setwd('.') # replace with your dir name or do command-D (Mac) Ctrl-Shift-H (Rstudio) to browse
lnames=load("data4wgcna.RData") # replace myfile with your file name
datt=t(degs) # replace vmysd with the name of your vsd dataframe

# Choose a set of soft-thresholding powers
powers = c(seq(from = 12, to=32, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")
# Plot the results:
#sizeGrWindow(9, 5)
pdf("soft_threshold_signed.pdf",height=4, width=8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


