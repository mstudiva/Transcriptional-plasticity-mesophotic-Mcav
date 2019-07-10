# installing WGCNA:
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("flashClust")
# install.packages("WGCNA",dependencies=TRUE)
# repos="http://cran.us.r-project.org"
# run these above commands once, then comment out

# always run these before running any of the following script chunks
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
# ?WGCNA
#--------------------------------
getwd()
setwd("~/path/to/local/directory")

lnames=load("data4wgcna.RData")
lnames # "vsd.wg"  "design" # log-transformed variance-stabilized gene expression, and table or experimental conditions
datt=t(vsd.wg)
ncol(datt)
nrow(datt)

head(design)
str(design)

# assembling table of QUANTITATIVE traits
shallow=as.numeric(design$depth=="shallow")
mesophotic=as.numeric(design$depth=="mesophotic")
transplant=as.numeric(design$depth=="transplant")
zero=as.numeric(design$time=="zero")
six=as.numeric(design$time=="six")
twelve=as.numeric(design$time=="twelve")

traits <- cbind(shallow, mesophotic, transplant, zero, six, twelve, design[c(9:14)])
# traits <- design[c(9:14)]
traits

write.csv(traits, file="traits.csv")

save(datt,traits,file="wgcnaData.RData")

####################
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
load("wgcnaData.RData")

# Try different betas ("soft threshold") - power factor for calling connections between genes
powers = c(seq(from = 2, to=26, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")

# Plot the results:
# Run from the line below to dev.off()
sizeGrWindow(9, 5)
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

#####################
# making modules

# take a look at the threshold plots produced above, and the output table from the pickSoftThreshold command
# pick the power that corresponds with a SFT.R.sq value above 0.90

# run from the line below to the save command
s.th=7 # re-specify according to previous section
# the following two lines take a long time, prepare to wait 15-20 min
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules.RData")

#########################
# merging modules:

mm=load('1stPassModules.RData')
mm
lnames=load('wgcnaData.RData')
# traits
# head(datt)

quartz()

MEDissThres = 1 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed.RData")

###################
# plotting correlations with traits:
load(file = "networkdata_signed.RData")
load(file = "wgcnaData.RData");

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# gene-trait correlations - a gene-by-gene heatmap corresponding to the droopy tree
# (to make augmented ridiculogram as in mice-men-embryos paper)
 # quartz()
 # geneTraitCor = cor(datt, traits, use = "p");
 # colnames(geneTraitCor)
 # geneTraitCor=geneTraitCor[geneTree$order,]
 # head(geneTraitCor)
 # labeledHeatmap(Matrix = geneTraitCor,
 # xLabels = colnames(geneTraitCor),
 # xLabelsAngle=90,
 # ySymbols = FALSE,
 # colorLabels = FALSE,
 # colors = blueWhiteRed(50),
 # setStdMargins = FALSE, cex.text = 0.5,
 # zlim = c(-1,1),
 # main = paste("Gene-trait relationships"))

# module-trait correlations
quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(traits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

table(moduleColors) # gives numbers of genes in each module

# if it was first pass with no module merging, this is where you examine your heatmap
# and dendrogram of module eigengenes to see
# where you would like to set cut height (MEDissThres parameter) in the previous section
# to merge modules that are telling the same story for your trait data

# good way to do it is to find a group of similar modules in the heat map and then see
# at which tree height they connect in the dendrogram.

#############
# scatterplots of gene significance (correlation-based) vs kME

load(file = "networkdata_signed.RData")
load(file = "wgcnaData.RData");
traits
table(moduleColors)
whichTrait="chl_ac"

nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(traits[,whichTrait]);
names(selTrait) = whichTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
quartz()
par(mfrow=c(2,3))
counter=0
for(module in modNames[1:length(modNames)]){
counter=counter+1
if (counter>9) {
	quartz()
	par(mfrow=c(3,3))
	counter=1
}
column = match(module, modNames);
moduleGenes = moduleColors==module;
#trr="heat resistance"
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste(module,"module membership"),
ylab = paste("GS for", whichTrait),
col = "grey50",mgp=c(2.3,1,0))
}

################
# eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules

load(file = "networkdata_signed.RData")
load(file = "wgcnaData.RData");

which.module="skyblue"
datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # number of genes in chosen module

#################
# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)

library(WGCNA)
load(file = "networkdata_signed.RData") # moduleColors, MEs
load(file = "wgcnaData.RData") # vsd table
load(file = "data4wgcna.RData") # vsd table

# calculating modul memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs))
names(allkME)=gsub("kME","",names(allkME))

whichModule="blue"
table(moduleColors==whichModule) # how many genes are in it?

# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==whichModule)
combo=data.frame("gene"=row.names(vsd.wg),"Fish_kME"=allkME[,whichModule]*inModuleBinary)
write.csv(combo,file=paste(whichModule,".csv",sep=""),row.names=F,quote=F)

################
# plotting heatmap for named top-kME genes

library(WGCNA)
load(file = "networkdata_signed.RData")
load(file = "data4wgcna.RData")
load(file = "wgcnaData.RData");
allkME =as.data.frame(signedKME(datt, MEs))
gg=read.delim(file="Mcavernosa_Cladocopium_iso2geneName.tab",sep="\t")
library(pheatmap)

whichModule="black"
top=30 # number of named top-kME genes to plot

datME=MEs
datExpr=datt
modcol=paste("kME",whichModule,sep="")
sorted=vsd.wg[order(allkME[,modcol],decreasing=T),]
head(sorted)
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
	if (row.names(sorted)[i] %in% gg$V1) {
		counts=counts+1
		gn=gg[gg$V1==row.names(sorted)[i],2]
		gn=paste(gn,row.names(sorted)[i],sep=".")
		if (gn %in% gnames) {
			gn=paste(gn,counts,sep=".")
		}
		gnames=append(gnames,gn)
		hubs=data.frame(rbind(hubs,sorted[i,]))
		if (counts==top) {break}
	}
}
row.names(hubs)=gnames

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

pheatmap(hubs,scale="row",col=contrasting2,border_color=NA,treeheight_col=0,cex=0.9,cluster_rows=F)
