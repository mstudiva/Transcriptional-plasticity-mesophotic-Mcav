setwd("~/path/to/local/directory")

# run these once, then comment out
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# biocLite("DESeq2",dependencies=T)
# biocLite("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# biocLite("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")

#---------------------
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)

#read in counts
counts = read.table("allcounts_enviro.txt")

# how many genes we have total?
nrow(counts)
# 17901 for host, 13375 for symbiont
ncol(counts)
# 265 samples

# how does the data look?
head(counts)

#---------------------
# the following section prefilters the counts dataframe to remove genes with low counts (counts under 10 across all samples)
# while DESeq2 filters out low-count genes during testing anyways, prefiltering makes data transformations go faster

keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
# 12929 for host, 5014 for symbiont
ncol(countData)
write.csv(countData, file="countData.csv")

# importing a design .csv file
design = read.csv("design_enviro.csv", head=TRUE)
design
str(design)

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+depth+site:depth)

# reorders treatment factor according to "control" vs "treatment" levels
dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and is recommended for datasets with more than 20 or so samples
Vsd=varianceStabilizingTransformation(dds)
# rl=rlog(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))
# e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("site","depth"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers, array 150:
outs=c(150)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
# rl=rl[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+depth+site:depth)
dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,file="initial.RData")
# save(dds,design,countData,rl,file="initial.RData")

#---------------------
# generating normalized variance-stabilized data for PCoA, heatmaps, etc

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# creating normalized dataframe
vsd=assay(Vsd)
# vsd=assay(rl)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2],design[,3],sep=".")
# renames the column names
colnames(vsd)=snames
save(vsd,design,file="vsd.RData")

#-------------------
# exploring similarities among samples

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_enviro_mcav.pdf", width=50, height=50)
pheatmap(cor(vsd))
dev.off()

# Principal coordinate analysis
library(vegan)
library(rgl)
library(ape)

conditions=design

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's
# 3 PCs for host, 4 PCs for symbiont

# plotting PCoA
pdf(file="PCoA_enviro_mcav.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#f6e8c3","#5ab4ac", "#8c510a","#01665e"),pch=c(19,17)[as.numeric((as.factor(conditions$depth)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Site")
# cluster overlay of site
ordiellipse(scores, conditions$site, label=F, draw= "polygon", col=c("#f6e8c3","#01665e", "#5ab4ac","#8c510a"))
legend("topleft", legend=c("BLZ","WFGB","EFGB","PRG-DRT"), fill = c("#f6e8c3","#01665e", "#5ab4ac","#8c510a"), bty="n")
legend("bottomleft", legend=c("mesophotic","shallow"), pch=c(19,17), bty="n")
# cluster overlay of depth
plot(scores[,1], scores[,2],col=c("coral","cyan3"),pch=c(0,10,6,4)[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1", ylab="Coordinate 2", main="Depth")
ordiellipse(scores, conditions$depth, label=F, draw= "polygon", col=c("coral","cyan3"))
legend("bottomleft", legend=c("mesophotic","shallow"), fill = c("coral","cyan3"), bty="n")
legend("topleft", legend=c("BLZ","WFGB","EFGB","PRG-DRT"), pch=c(0,10,6,4), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=15, height=40)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: site and site:depth interaction are significant
ad=adonis(t(vsd)~site*depth,data=conditions,method="manhattan")
ad
# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:4],labels=row.names(ad$aov.tab)[1:4],col=cols,main="site vs depth")
dev.off()

# DAPC analysis: creating discriminant function to tell sites apart, based on shallow from mesophotic
# cannot test depth here since there are only 2 groups (need to remove 1 for the later predictions)
library(adegenet)

# runs simulations on randomly-chosen datasets of 90% of the total dataset to test the number of PCs to retain
set.seed(999)
# by site, excluding mesophotic samples
xvalDapc(t(vsd[,conditions$depth!="mesophotic"]),conditions$site[conditions$depth!="mesophotic"], n.rep=100, parallel="multicore", ncpus= 8)
# 20 PCs for host, 60 PCs for symbiont
# need to test again with a smaller range of possible PCs and more reps
# change the n.pca= statement depending on your target range
xvalDapc(t(vsd[,conditions$depth!="mesophotic"]),conditions$site[conditions$depth!="mesophotic"], n.rep=1000, n.pca=10:30, parallel="multicore", ncpus= 8)
# 23 PCs for host, 60 PCs for symbiont

# now running the dapc without transplants
dp.s=dapc(t(vsd[,conditions$depth!="mesophotic"]),conditions$site[conditions$depth!="mesophotic"],n.pca=23, n.da=3)

# can we predict site for the mesophotic samples based on patterns among the shallow samples?
pred.s=predict.dapc(dp.s,newdata=(t(vsd[,conditions$depth=="mesophotic"])))
pred.s
# look at the posterior section for assignments by probability

# creating a new dapc object to add in mesophotic corals for plotting
dp.meso=dp.s
dp.meso$ind.coord=pred.s$ind.scores
dp.meso$posterior=pred.s$posterior
dp.meso$assign=pred.s$assign
# for host
dp.meso$grp<-as.factor(c("BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","EFGB","EFGB","EFGB","WFGB","EFGB","EFGB","WFGB","WFGB","WFGB","WFGB","EFGB","EFGB","WFGB","EFGB","WFGB","EFGB","WFGB","WFGB","EFGB","EFGB","WFGB","EFGB","EFGB","WFGB","WFGB","EFGB","WFGB","EFGB","EFGB","EFGB","EFGB","WFGB"))
# for symbionts
# dp.meso$grp<-as.factor(c("BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","BLZ","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","PRG-DRT","WFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","WFGB","EFGB","WFGB","EFGB","EFGB","EFGB","WFGB","EFGB","EFGB","WFGB","WFGB","WFGB","WFGB","EFGB","EFGB","EFGB","WFGB","EFGB","WFGB","WFGB","EFGB","EFGB","WFGB","EFGB","EFGB","WFGB","WFGB","EFGB","WFGB","EFGB","EFGB","EFGB","EFGB","WFGB"))

# now exporting side by side figures of shallow vs mesophotic
pdf(file="DAPC_enviro_mcav_site.pdf", width=12, height=6)
par(mfrow=c(1,2))
scatter(dp.s, bg="white",scree.da=FALSE,legend=TRUE,solid=0.6, col= c("#f6e8c3","#01665e", "#5ab4ac","#8c510a"))
scatter(dp.meso, bg="white",scree.da=FALSE,legend=FALSE,solid=0.6, col= c("#f6e8c3","#01665e", "#5ab4ac","#8c510a"))
dev.off()

#----------------------
# FITTING GENE BY GENE MODELS

# with multi-factor, multi-level design - using LRT

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# Running DESeq with LRT for the effect of interaction (using full model dds = same as design argument on line 59, ~ site+depth+site:depth)
# then site+depth are removed, to look just at interaction term
dds.i=DESeq(dds,test="LRT",reduced=~site+depth, parallel=TRUE)

# Creating dataset with design =~site+depth (no interaction), to investigate importance of site and depth
dds1=DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+depth)
dds1$depth <- factor(dds1$depth, levels = c("shallow","mesophotic"))

# model for the effect of depth (2 factor levels => Wald test)
dds.depth=DESeq(dds1, parallel=TRUE)

# model for the effect of site: (>2 factor levels => LRT)
dds.site=DESeq(dds1,test="LRT",reduced=~depth, parallel=TRUE)

# saving all models
save(dds,dds.site,dds.depth,dds.i,file="realModels.RData")

#--------------
# GLM analysis

load("realModels.RData")
library(DESeq2)

# depth treatment - Wald test
depth=results(dds.depth,contrast=c("depth","shallow","mesophotic")) # factor depth, level shallow compared to mesophotic
summary(depth) # look at how many DEGs there are and how many low-expressed genes ("low counts") had to be removed
depth
degs.depth=row.names(depth)[depth$padj<0.1 & !(is.na(depth$padj))]

# site
site=results(dds.site)
summary(site)
site
degs.site=row.names(site)[site$padj<0.1 & !(is.na(site$padj))]

# depth:site interaction
int=results(dds.i)
summary(int)
int
degs.int=row.names(int)[int$padj<0.1 & !(is.na(int$padj))]

# site contrasts
CBC.EFGB=results(dds.site,contrast=c("site","CBC","EFGB"))
summary(CBC.EFGB)
CBC.EFGB
degs.CBC.EFGB=row.names(CBC.EFGB)[CBC.EFGB$padj<0.1 & !(is.na(CBC.EFGB$padj))]

CBC.PRTER=results(dds.site,contrast=c("site","CBC","PRTER"))
summary(CBC.PRTER)
CBC.PRTER
degs.CBC.PRTER =row.names(CBC.PRTER)[CBC.PRTER$padj<0.1 & !(is.na(CBC.PRTER$padj))]

CBC.WFGB=results(dds.site,contrast=c("site","CBC","WFGB"))
summary(CBC.WFGB)
CBC.WFGB
degs.CBC.WFGB=row.names(CBC.WFGB)[CBC.WFGB$padj<0.1 & !(is.na(CBC.WFGB$padj))]

EFGB.PRTER=results(dds.site,contrast=c("site","EFGB","PRTER"))
summary(EFGB.PRTER)
EFGB.PRTER
degs.EFGB.PRTER=row.names(EFGB.PRTER)[EFGB.PRTER$padj<0.1 & !(is.na(EFGB.PRTER$padj))]

EFGB.WFGB=results(dds.site,contrast=c("site","EFGB","WFGB"))
summary(EFGB.WFGB)
EFGB.WFGB
degs.EFGB.WFGB=row.names(EFGB.WFGB)[EFGB.WFGB$padj<0.1 & !(is.na(EFGB.WFGB$padj))]

PRTER.WFGB=results(dds.site,contrast=c("site","PRTER","WFGB"))
summary(PRTER.WFGB)
PRTER.WFGB
degs.PRTER.WFGB=row.names(PRTER.WFGB)[PRTER.WFGB$padj<0.1 & !(is.na(PRTER.WFGB$padj))]

# depth within site contrasts
# For figuring out your contrast statements, the following vignette is helpful
# ?results
# Contrast statements need to be a combination of those listed under resultsNames
resultsNames(dds)

CBC.d=results(dds,contrast=c("depth", "mesophotic", "shallow"))
summary(CBC.d)
CBC.d
degs.CBC=row.names(CBC.d)[CBC.d$padj<0.1 & !(is.na(CBC.d$padj))]

EFGB.d=results(dds,contrast=list("depth_mesophotic_vs_shallow","siteEFGB.depthmesophotic"))
summary(EFGB.d)
EFGB.d
degs.EFGB=row.names(EFGB.d)[EFGB.d$padj<0.1 & !(is.na(EFGB.d$padj))]

PRTER.d=results(dds,contrast=list("depth_mesophotic_vs_shallow","sitePRTER.depthmesophotic"))
summary(PRTER.d)
PRTER.d
degs.PRTER=row.names(PRTER.d)[PRTER.d$padj<0.1 & !(is.na(PRTER.d$padj))]

WFGB.d=results(dds,contrast=list("depth_mesophotic_vs_shallow","siteWFGB.depthmesophotic"))
summary(WFGB.d)
WFGB.d
degs.WFGB=row.names(WFGB.d)[WFGB.d$padj<0.1 & !(is.na(WFGB.d$padj))]

# is the depth effect different across sites?
CBC.EFGB.d=results(dds, name="siteEFGB.depthmesophotic")
summary(CBC.EFGB.d)
CBC.EFGB.d
degs.CBC.EFGB.d=row.names(CBC.EFGB.d)[CBC.EFGB.d$padj<0.1 & !(is.na(CBC.EFGB.d$padj))]

CBC.PRTER.d=results(dds, name="sitePRTER.depthmesophotic")
summary(CBC.PRTER.d)
CBC.PRTER.d
degs.CBC.PRTER.d=row.names(CBC.PRTER.d)[CBC.PRTER.d$padj<0.1 & !(is.na(CBC.PRTER.d$padj))]

CBC.WFGB.d=results(dds, name="siteWFGB.depthmesophotic")
summary(CBC.WFGB.d)
CBC.WFGB.d
degs.CBC.WFGB.d=row.names(CBC.WFGB.d)[CBC.WFGB.d$padj<0.1 & !(is.na(CBC.WFGB.d$padj))]

EFGB.PRTER.d=results(dds, contrast=list("siteEFGB.depthmesophotic","sitePRTER.depthmesophotic"))
summary(EFGB.PRTER.d)
EFGB.PRTER.d
degs.EFGB.PRTER.d=row.names(EFGB.PRTER.d)[EFGB.PRTER.d$padj<0.1 & !(is.na(EFGB.PRTER.d$padj))]

EFGB.WFGB.d=results(dds, contrast=list("siteEFGB.depthmesophotic","siteWFGB.depthmesophotic"))
summary(EFGB.WFGB.d)
EFGB.WFGB.d
degs.EFGB.WFGB.d=row.names(EFGB.WFGB.d)[EFGB.WFGB.d$padj<0.1 & !(is.na(EFGB.WFGB.d$padj))]

PRTER.WFGB.d=results(dds, contrast=list("sitePRTER.depthmesophotic","siteWFGB.depthmesophotic"))
summary(PRTER.WFGB.d)
PRTER.WFGB.d
degs.PRTER.WFGB.d=row.names(PRTER.WFGB.d)[PRTER.WFGB.d$padj<0.1 & !(is.na(PRTER.WFGB.d$padj))]

save(site,depth,int,CBC.EFGB,CBC.PRTER,CBC.WFGB,EFGB.PRTER,EFGB.WFGB,PRTER.WFGB,CBC.d,EFGB.d,PRTER.d,WFGB.d,CBC.EFGB.d,CBC.PRTER.d,CBC.WFGB.d,EFGB.PRTER.d,EFGB.WFGB.d,PRTER.WFGB.d,degs.site,degs.depth,degs.int,degs.CBC.EFGB,degs.CBC.PRTER,degs.CBC.WFGB,degs.EFGB.PRTER,degs.EFGB.WFGB,degs.PRTER.WFGB,degs.CBC,degs.EFGB,degs.PRTER,degs.WFGB,degs.CBC.EFGB.d,degs.CBC.PRTER.d,degs.CBC.WFGB.d,degs.EFGB.PRTER.d,degs.EFGB.WFGB.d,degs.PRTER.WFGB.d,file="pvals.RData")

#-------------------
# density plots: are my DEGs high-abundant or low-abundant?

load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs.site]),col="blue")
lines(density(means[degs.depth]),col="orange")
lines(density(means[degs.int]),col="lightblue")
legend("topright", title = "Factor", legend=c("site","depth","interaction"), fill = c("blue","orange","lightblue"))
dev.off()

#-------------------
# venn diagrams

load("pvals.RData")
library(DESeq2)

# overall factors, full model
candidates=list("depth"=degs.depth,"site"=degs.site, "interaction"=degs.int)

# DEGs across depths within site - are some genes conserved?
sitedepth=list("CBC"=degs.CBC,"PRTER"=degs.PRTER,"WFGB"=degs.WFGB,"EFGB"=degs.EFGB)

# install.packages("VennDiagram")
library(VennDiagram)

# overall factors, full model
fullmodel.venn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "orange", "lightblue"),
  alpha = 0.5,
  label.col = c("darkblue", "white", "darkred", "white", "white", "white", "cornflowerblue"),
  cex = 5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkblue", "darkred", "cornflowerblue"),
  cat.cex = 5,
  cat.fontfamily = "sans",
  cat.dist = c(0.06, 0.06, -0.06),
  cat.pos = 3
)
pdf(file="Venn_enviro_mcav.pdf", height=12, width=12)
grid.draw(fullmodel.venn)
dev.off()

# DEGs across depths within site - are some genes conserved?
sitedepth.venn=venn.diagram(
  x = sitedepth,
  filename=NULL,
  col = "transparent",
  fill = c("#f6e8c3", "#8c510a", "#01665e", "#5ab4ac"),
  alpha = 0.5,
  label.col = c("#01665e","white","#35978f","white","white","black","white", "white","#dfc27d","white","white","white","white","#8c510a","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("#dfc27d", "#8c510a", "#01665e", "#5ab4ac"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_sitedepth_mcav.pdf", height=10, width=12)
grid.draw(sitedepth.venn)
dev.off()

#-------------------
# saving data for GO and KOG analysis
load("realModels.RData")
load("pvals.RData")

# depth response
# log2 fold changes:
head(depth)
source=depth[!is.na(depth$pvalue),]
depth.fc=data.frame("gene"=row.names(source))
depth.fc$lfc=source[,"log2FoldChange"]
head(depth.fc)
write.csv(depth.fc,file="depth_fc.csv",row.names=F,quote=F)
save(depth.fc,file="depth_fc.RData")

# signed log p-values: -log(pvalue)* direction:
depth.p=data.frame("gene"=row.names(source))
depth.p$lpv=-log(source[,"pvalue"],10)
depth.p$lpv[source$stat<0]=depth.p$lpv[source$stat<0]*-1
head(depth.p)
write.csv(depth.p,file="depth_lpv.csv",row.names=F,quote=F)
save(depth.p,file="depth_lpv.RData")

# site response
# since DESeq2 was designed for 2 factor designs (batch vs treatment), only the treatment fold change data are relevant
# signed log p-values: -log(pvalue)* direction:
head(site)
source=site[!is.na(site$pvalue),]
site.p=data.frame("gene"=row.names(source))
site.p$lpv=-log(source[,"pvalue"],10)
site.p$lpv[source$stat<0]=site.p$lpv[source$stat<0]*-1
head(site.p)
write.csv(site.p,file="site_lpv.csv",row.names=F,quote=F)
save(site.p,file="site_lpv.RData")

# site:depth interaction
# signed log p-values: -log(pvalue)* direction:
head(int)
source=int[!is.na(int$pvalue),]
int.p=data.frame("gene"=row.names(source))
int.p$lpv=-log(source[,"pvalue"],10)
int.p$lpv[source$stat<0]=int.p$lpv[source$stat<0]*-1
head(int.p)
write.csv(int.p,file="int_lpv.csv",row.names=F,quote=F)
save(int.p,file="int_lpv.RData")

#-------------------
load("realModels.RData")
load("pvals.RData")

# depth response within CBC
# log2 fold changes:
head(CBC.d)
source=CBC.d[!is.na(CBC.d$pvalue),]
CBC.fc=data.frame("gene"=row.names(source))
CBC.fc$lfc=source[,"log2FoldChange"]
head(CBC.fc)
write.csv(CBC.fc,file="CBC_fc.csv",row.names=F,quote=F)
save(CBC.fc,file="CBC_fc.RData")

# signed log p-values: -log(pvalue)* direction:
CBC.p=data.frame("gene"=row.names(source))
CBC.p$lpv=-log(source[,"pvalue"],10)
CBC.p$lpv[source$stat<0]= CBC.p$lpv[source$stat<0]*-1
head(CBC.p)
write.csv(CBC.p,file="CBC_lpv.csv",row.names=F,quote=F)
save(CBC.p,file="CBC_lpv.RData")

# depth response within EFGB
# log2 fold changes:
head(EFGB.d)
source=EFGB.d[!is.na(EFGB.d$pvalue),]
EFGB.fc=data.frame("gene"=row.names(source))
EFGB.fc$lfc=source[,"log2FoldChange"]
head(EFGB.fc)
write.csv(EFGB.fc,file="EFGB_fc.csv",row.names=F,quote=F)
save(EFGB.fc,file="EFGB_fc.RData")

# signed log p-values: -log(pvalue)* direction:
EFGB.p=data.frame("gene"=row.names(source))
EFGB.p$lpv=-log(source[,"pvalue"],10)
EFGB.p$lpv[source$stat<0]= EFGB.p$lpv[source$stat<0]*-1
head(EFGB.p)
write.csv(EFGB.p,file="EFGB_lpv.csv",row.names=F,quote=F)
save(EFGB.p,file="EFGB_lpv.RData")

# depth response within PRTER
# log2 fold changes:
head(PRTER.d)
source=PRTER.d[!is.na(PRTER.d$pvalue),]
PRTER.fc=data.frame("gene"=row.names(source))
PRTER.fc$lfc=source[,"log2FoldChange"]
head(PRTER.fc)
write.csv(PRTER.fc,file="PRTER_fc.csv",row.names=F,quote=F)
save(PRTER.fc,file="PRTER_fc.RData")

# signed log p-values: -log(pvalue)* direction:
PRTER.p=data.frame("gene"=row.names(source))
PRTER.p$lpv=-log(source[,"pvalue"],10)
PRTER.p$lpv[source$stat<0]= PRTER.p$lpv[source$stat<0]*-1
head(PRTER.p)
write.csv(PRTER.p,file="PRTER_lpv.csv",row.names=F,quote=F)
save(PRTER.p,file="PRTER_lpv.RData")

# depth response within WFGB
# log2 fold changes:
head(WFGB.d)
source=WFGB.d[!is.na(WFGB.d$pvalue),]
WFGB.fc=data.frame("gene"=row.names(source))
WFGB.fc$lfc=source[,"log2FoldChange"]
head(WFGB.fc)
write.csv(WFGB.fc,file="WFGB_fc.csv",row.names=F,quote=F)
save(WFGB.fc,file="WFGB_fc.RData")

# signed log p-values: -log(pvalue)* direction:
WFGB.p=data.frame("gene"=row.names(source))
WFGB.p$lpv=-log(source[,"pvalue"],10)
WFGB.p$lpv[source$stat<0]= WFGB.p$lpv[source$stat<0]*-1
head(WFGB.p)
write.csv(WFGB.p,file="WFGB_lpv.csv",row.names=F,quote=F)
save(WFGB.p,file="WFGB_lpv.RData")

#-------------------
# exporting dataset of conserved genes across all 4 sites
# differs from the previous export steps in that only significant DEGs are exported here, to specifically look at conserved genes
load("realModels.RData")
load("pvals.RData")

library(plyr)

# creates dataframe of significant DEGs for each site, one each for lfc and lpv
# first section is for lpv
source=BLZ.d[BLZ.d$padj<0.1 & !(is.na(BLZ.d$padj)),]
BLZ=data.frame("gene"=row.names(source))
BLZ$lpv=-log(source[,"pvalue"],10)
BLZ$lpv[source$stat<0]= BLZ$lpv[source$stat<0]*-1
head(BLZ)

source=EFGB.d[EFGB.d$padj<0.1 & !(is.na(EFGB.d$padj)),]
EFGB=data.frame("gene"=row.names(source))
EFGB$lpv=-log(source[,"pvalue"],10)
EFGB$lpv[source$stat<0]= EFGB$lpv[source$stat<0]*-1
head(EFGB)

source=PRG_DRT.d[PRG_DRT.d$padj<0.1 & !(is.na(PRG_DRT.d$padj)),]
PRG_DRT=data.frame("gene"=row.names(source))
PRG_DRT$lpv=-log(source[,"pvalue"],10)
PRG_DRT$lpv[source$stat<0]= PRG_DRT$lpv[source$stat<0]*-1
head(PRG_DRT)

source=WFGB.d[WFGB.d$padj<0.1 & !(is.na(WFGB.d$padj)),]
WFGB=data.frame("gene"=row.names(source))
WFGB$lpv=-log(source[,"pvalue"],10)
WFGB$lpv[source$stat<0]= WFGB$lpv[source$stat<0]*-1
head(WFGB)

# finds and outputs common genes and associated DEG reponses across all sites
commongenes<- join_all(list(BLZ,EFGB,PRG_DRT,WFGB), by="gene", type="inner")
str(commongenes)
names(commongenes)<- c("gene","BLZ","EFGB","PRG-DRT","WFGB")
write.csv(commongenes, file="commongenes_lpv.csv")

# takes the corresponding site column for export
BLZ.common=data.frame("gene"=commongenes$gene)
BLZ.common$lpv= commongenes[,"BLZ"]
BLZ.common
write.csv(BLZ.common,file="BLZ_common_lpv.csv",row.names=F,quote=F)
save(BLZ.common,file="BLZ_common_lpv.RData")

EFGB.common=data.frame("gene"=commongenes$gene)
EFGB.common$lpv= commongenes[,"EFGB"]
EFGB.common
write.csv(EFGB.common,file="EFGB_common_lpv.csv",row.names=F,quote=F)
save(EFGB.common,file="EFGB_common_lpv.RData")

PRG_DRT.common=data.frame("gene"=commongenes$gene)
PRG_DRT.common$lpv= commongenes[,"PRG-DRT"]
PRG_DRT.common
write.csv(PRG_DRT.common,file="PRG-DRT_common_lpv.csv",row.names=F,quote=F)
save(PRG_DRT.common,file="PRG-DRT_common_lpv.RData")

WFGB.common=data.frame("gene"=commongenes$gene)
WFGB.common$lpv= commongenes[,"WFGB"]
WFGB.common
write.csv(WFGB.common,file="WFGB_common_lpv.csv",row.names=F,quote=F)
save(WFGB.common,file="WFGB_common_lpv.RData")

#----------------------
# now fc
source=BLZ.d[BLZ.d$padj<0.1 & !(is.na(BLZ.d$padj)),]
BLZ=data.frame("gene"=row.names(source))
BLZ$lfc=source[,"log2FoldChange"]
head(BLZ)

source=EFGB.d[EFGB.d$padj<0.1 & !(is.na(EFGB.d$padj)),]
EFGB=data.frame("gene"=row.names(source))
EFGB$lfc=source[,"log2FoldChange"]
head(EFGB)

source=PRG_DRT.d[PRG_DRT.d$padj<0.1 & !(is.na(PRG_DRT.d$padj)),]
PRG_DRT=data.frame("gene"=row.names(source))
PRG_DRT$lfc=source[,"log2FoldChange"]
head(PRG_DRT)

source=WFGB.d[WFGB.d$padj<0.1 & !(is.na(WFGB.d$padj)),]
WFGB=data.frame("gene"=row.names(source))
WFGB$lfc=source[,"log2FoldChange"]
head(WFGB)

# finds and outputs common genes and associated DEG reponses across all sites
commongenes<- join_all(list(BLZ,EFGB,PRG_DRT,WFGB), by="gene", type="inner")
str(commongenes)
names(commongenes)<- c("gene","BLZ","EFGB","PRG-DRT","WFGB")
write.csv(commongenes, file="commongenes_fc.csv")

# takes the corresponding site column for export
BLZ.common=data.frame("gene"=commongenes$gene)
BLZ.common$lfc= commongenes[,"BLZ"]
BLZ.common
write.csv(BLZ.common,file="BLZ_common_fc.csv",row.names=F,quote=F)
save(BLZ.common,file="BLZ_common_fc.RData")

EFGB.common=data.frame("gene"=commongenes$gene)
EFGB.common$lfc= commongenes[,"EFGB"]
EFGB.common
write.csv(EFGB.common,file="EFGB_common_fc.csv",row.names=F,quote=F)
save(EFGB.common,file="EFGB_common_fc.RData")

PRG_DRT.common=data.frame("gene"=commongenes$gene)
PRG_DRT.common$lfc= commongenes[,"PRG-DRT"]
PRG_DRT.common
write.csv(PRG_DRT.common,file="PRG-DRT_common_fc.csv",row.names=F,quote=F)
save(PRG_DRT.common,file="PRG-DRT_common_fc.RData")

WFGB.common=data.frame("gene"=commongenes$gene)
WFGB.common$lfc= commongenes[,"WFGB"]
WFGB.common
write.csv(WFGB.common,file="WFGB_common_fc.csv",row.names=F,quote=F)
save(WFGB.common,file="WFGB_common_fc.RData")

#------------------------
# exporting dataset of common genes with vsd normalized expression
# first, export the vsd matrix to csv, then reimport to generate a dataframe with gene name as a true column
write.csv(vsd,file="vsd.csv")
vsd.common<- read.csv("vsd.csv", head=TRUE)
colnames(vsd.common)[1]<- "gene"
str(vsd.common)

# next, joins commongenes and vsd.common datasets to produce a dataframe of the conserved genes across all sites
# this includes each sample, so you can look at individual responses across the conserved genes
common<- join_all(list(commongenes,vsd.common), by="gene", type="inner")
str(common)
head(commongenes)
head(common)
# double check the isogroup numbers line up properly
write.csv(common, file="commongenes_vsd.csv")
