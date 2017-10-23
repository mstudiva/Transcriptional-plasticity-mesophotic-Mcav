setwd("~/path/to/local/directory")

# run these once, then comment out
# source("http://bioconductor.org/biocLite.R")
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

#---------------------
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)

#read in counts
counts = read.table("allcounts_enviro.txt")

# how many genes we have total?
nrow(counts) 

# how does the data look? 
head(counts)

#---------------------
# the following section prefilters the counts dataframe to remove genes with low counts
# DESeq2 does this naturally, but prefiltering might make data transformations go faster
# for larger datasets that will take forever anyways, I suggest skipping prefiltering
# and use vsd transformation instead of rlog

# how many genes have mean count >=3?
# means=apply(counts,1,mean)
# table(means>=3)
# plot(density(means),log="x")

# distribution of total counts among samples: are there samples with less the 2SD counts?
# smeans=apply(counts,2,sum)
# lsm=log(smeans,10)
# hist(lsm,breaks=12)
# sdl=sd(lsm)
# which(mean(lsm)-lsm > 2*sdl)

# counts[1:10,1]
# head(counts)
# removing all genes with mean count less than 3
# countData=counts[means>=3,]
# nrow(countData)
# str(countData)
# write.csv(countData, file="countData.csv")

#---------------------
countData=counts
nrow(countData)
ncol(countData)
write.csv(countData, file="countData.csv")

# for WGCNA: removing all genes with mean count less than 10 
# counts4wgcna=counts[means>=10,]
# alternate filtering method below

# for WCGNA: removing all genes with counts in less than 10 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x==0))<ncol(counts)*0.9,]
nrow(counts4wgcna)
ncol(counts4wgcna)

# importing a design .csv file
design = read.csv("design_enviro.csv", head=TRUE)
design
str(design)

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+depth+site:depth)

# pre-filtering at this point will speed up DESeq models without losing critical data
# creates a series of TRUE/FALSE statements for isogroups with total counts above 10 across all samples
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# Note: prefiltering caused my processing time to speed up, but many outliers were detected. Sticking with original data

# reorders treatment factor according to "control" vs "treatment" levels
dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)
# rl=rlog(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))
# e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("site","depth"),force=T)
dev.off()
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples
# samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# the same outliers from Fig 2 are also larger points
# move your mouse over and it will identify the sample name and array number
# use the array number for removal in the following section

## if there were outliers, say, sample number MS_145 (array 138):
outs=c(138, 168, 170)
countData=countData[,-outs]
counts4wgcna=counts4wgcna[,-outs]
Vsd=Vsd[,-outs]
# rl=rl[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+depth+site:depth)
dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,counts4wgcna,Vsd,file="initial.RData")
# save(dds,design,countData,counts4wgcna,rl,file="initial.RData")

#---------------------
# generating normalized variance-stabilized data for PCoA, heatmaps, WGCNA etc

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

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ site+depth+site:depth)
vsd.wg=assay(varianceStabilizingTransformation(wg))
# vsd.wg=assay(rlog(wg))
colnames(vsd.wg)=snames
save(vsd.wg,design,file="data4wgcna.RData")

#-------------------
# EXPLORING SIMILARITIES AMONG SAMPLES

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_enviro.pdf", width=50, height=50)
pheatmap(cor(vsd))
dev.off()

# Principal coordinate analysis
library(vegan)
library(rgl)
library(ape)

conditions=design
# creates a new column in the conditions dataframe with a comination of sample ID and factors
conditions$ibyt=paste(conditions[,1],conditions[,2],conditions[,3],sep=".")

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

# plotting PCoA
pdf(file="PCoA_enviro.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=as.numeric(as.factor(conditions$ind)),pch=19, xlab="Coordinate 1", ylab="Coordinate 2", main="Site")
# cluster overlay of site
ordihull(scores,conditions$site,label=F,draw="polygon",col=c("blue","orange", "red","lightblue"),cex=1.5)
legend("topleft", title = "Site", legend=c("CBC","WFGB","EFGB","PRTER"), fill = c("blue","lightblue","orange","red"))
# cluster overlay of depth
plot(scores[,1], scores[,2],col=as.numeric(as.factor(conditions$ind)),pch=19, xlab="Coordinate 1", ylab="Coordinate 2", main="Depth")
ordihull(scores,conditions$depth,label=F,draw="polygon",col=c("coral","cyan3"),cex=1.5)
legend("topleft", title = "Depth", legend=c("shallow","mesophotic"), fill = c("cyan3","coral"))
dev.off()

# other overlay options
# ordispider(scores,conditions$ibyt,label=T)
# ordiellipse(scores,conditions$trt,label=T,draw="polygon",col="grey90",cex=2)

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=15, height=40)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# interactive 3d plot (package rgl)
plot3d(scores[,1], scores[,2], scores[,3],col=as.numeric(as.factor(conditions$site)),type="s",radius=0.5*as.numeric(as.factor(conditions$depth)))
# this thing is intense

# formal analysis of variance in distance matricies: site and site:depth interaction are significant
ad=adonis(t(vsd)~site*depth,data=conditions,method="manhattan")
ad
# creating pie chart to represent ANOVA results
cols=c("skyblue","green2","coral","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:4],labels=row.names(ad$aov.tab)[1:4],col=cols,main="Site vs Depth")
dev.off()

# DAPC analysis: creating discriminant function to tell shallow from mesophotic, based on factors
library(adegenet)

pdf(file="DAPC_enviro.pdf", width=12, height=6)
par(mfrow=c(1,2))
dp.s=dapc(t(vsd[,conditions$depth!="mesophotic"]),conditions$site[conditions$depth!="mesophotic"],n.pca=3, n.da=1)
scatter(dp.s,bg="white",scree.da=FALSE,legend=TRUE,solid=.4, col = c("blue","lightblue","orange","red")) #discrimination of mesophotic expression by site

dp.d=dapc(t(vsd[,conditions$site!="mesophotic"]),conditions$depth[conditions$site!="mesophotic"],n.pca=3, n.da=1)
scatter(dp.d,bg="white",scree.da=FALSE,legend=TRUE,solid=.4, col= c("cyan3","coral")) #discrimination of mesophotic expression by depth
dev.off()

# can we predict site by mesophotic trends?
pred.s=predict.dapc(dp.s,newdata=(t(vsd[,conditions$depth=="mesophotic"]))) 
pred.s 
# not especially

# can we predict depth by mesophotic trends?
pred.d=predict.dapc(dp.d,newdata=(t(vsd[,conditions$depth=="mesophotic"]))) 
pred.d 
# at certain sites

#----------------------
# FITTING GENE BY GENE MODELS

# with multi-factor, multi-level design - using LRT

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# alternative method to obtain contrast statements
# group all factors into single factor "group"
# dds$group <- factor(paste0(dds$site, dds$depth))
# design(dds) <- ~ group
# then run model same as below

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# Running DESeq with LRT for the effect of interaction (using full model dds = same as design argument on line 63, ~ site+depth+site:depth)
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

plot(density(means))
lines(density(means[degs.site]),col="blue")
lines(density(means[degs.depth]),col="red")
lines(density(means[degs.int]),col="gold")
legend("topright", title = "Factor", legend=c("site","depth","interaction"), fill = c("blue","red","gold"))

#-------------------
# venn diagrams

load("pvals.RData")
library(DESeq2)

# overall factors, full model
#candidates=list("depth"=row.names(depth)[depth$padj<0.1 & !(is.na(depth$padj))],"site"=row.names(site)[site$padj<0.1 & !(is.na(site$padj))],"interaction"=row.names(int)[int$padj<0.1 & !(is.na(int$padj))])

candidates=list("depth"=degs.depth,"site"=degs.site, "interaction"=degs.int)

# DEGs across depths within site - are some genes conserved?
#sitedepth =list("CBC"=row.names(CBC.d)[CBC.d$padj<0.1 & !(is.na(CBC.d$padj))],"WFGB"=row.names(WFGB.d)[WFGB.d$padj<0.1 & !(is.na(WFGB.d$padj))],"EFGB"=row.names(EFGB.d)[EFGB.d$padj<0.1 & !(is.na(EFGB.d$padj))],"PRTER"=row.names(PRTER.d)[PRTER.d$padj<0.1 & !(is.na(PRTER.d$padj))])

sitedepth=list("CBC"=degs.CBC,"PRTER"=degs.PRTER,"WFGB"=degs.WFGB,"EFGB"=degs.EFGB)

# simple venn
library(gplots)
venn(candidates)

#-------------------
# pretty venn diagram

# install.packages("VennDiagram")

library(VennDiagram)

# overall factors, full model
fullmodel.venn=venn.diagram(
	x = candidates,
	filename=NULL,
	col = "transparent",
	fill = c("coral", "cyan3", "khaki3"),
	alpha = 0.5,
	label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
	cex = 5,
	fontfamily = "sans",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col =c("darkred", "darkblue", "darkgreen"),
	cat.cex = 5,
	cat.fontfamily = "sans",
	cat.dist = c(0.06, 0.06, -0.06),
	cat.pos = 3
	)
pdf(file="DEG_enviro_venn.pdf", height=12, width=12)
grid.draw(fullmodel.venn)
dev.off()

# DEGs across depths within site - are some genes conserved?
sitedepth.venn=venn.diagram(
	x = sitedepth,
	filename=NULL,
	col = "transparent",
	fill = c("coral", "khaki3", "cyan3", "lightskyblue"),
	alpha = 0.5,
	label.col = c("darkblue","white","royalblue1","white","white","black","white", "white","darkred","white","white","white","white","darkgreen","white"),
	cex = 3.5,
	fontfamily = "sans",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col =c("coral", "khaki3", "cyan3", "lightskyblue"),
	cat.cex = 3.5,
	cat.fontfamily = "sans",
	cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
	)
pdf(file="DEG_sitedepth_venn.pdf", height=10, width=12)
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
# run each of these twice with the corresponding lines uncommented to export lfc and lpv datasets
source=CBC.d[CBC.d$padj<0.1 & !(is.na(CBC.d$padj)),]
CBC=data.frame("gene"=row.names(source))
#CBC$lfc=source[,"log2FoldChange"]
CBC$lpv=-log(source[,"pvalue"],10)
CBC$lpv[source$stat<0]= CBC$lpv[source$stat<0]*-1
head(CBC)

source=EFGB.d[EFGB.d$padj<0.1 & !(is.na(EFGB.d$padj)),]
EFGB=data.frame("gene"=row.names(source))
#EFGB$lfc=source[,"log2FoldChange"]
EFGB$lpv=-log(source[,"pvalue"],10)
EFGB$lpv[source$stat<0]= EFGB$lpv[source$stat<0]*-1
head(EFGB)

source=PRTER.d[PRTER.d$padj<0.1 & !(is.na(PRTER.d$padj)),]
PRTER=data.frame("gene"=row.names(source))
#PRTER$lfc=source[,"log2FoldChange"]
PRTER$lpv=-log(source[,"pvalue"],10)
PRTER$lpv[source$stat<0]= PRTER$lpv[source$stat<0]*-1
head(PRTER)

source=WFGB.d[WFGB.d$padj<0.1 & !(is.na(WFGB.d$padj)),]
WFGB=data.frame("gene"=row.names(source))
#WFGB$lfc=source[,"log2FoldChange"]
WFGB$lpv=-log(source[,"pvalue"],10)
WFGB$lpv[source$stat<0]= WFGB$lpv[source$stat<0]*-1
head(WFGB)

# finds and outputs common genes and associated DEG reponses (lfc or lpv) across all sites
commongenes<- join_all(list(CBC,EFGB,PRTER,WFGB), by="gene", type="inner")
head(commongenes)
names(commongenes)<- c("gene","CBC","EFGB","PRTER","WFGB")

# takes the corresponding site column of lfc or lpv for export
# run each of these twice with the corresponding lines uncommented to export lfc and lpv datasets
CBC.common=data.frame("gene"=commongenes$gene)
#CBC.common$lfc= commongenes[,"CBC"]
CBC.common$lpv= commongenes[,"CBC"]
CBC.common
#write.csv(CBC.common,file="CBC_common_lfc.csv",row.names=F,quote=F)
#save(CBC.common,file="CBC_common_lfc.RData")
write.csv(CBC.common,file="CBC_common_lpv.csv",row.names=F,quote=F)
save(CBC.common,file="CBC_common_lpv.RData")

EFGB.common=data.frame("gene"=commongenes$gene)
#EFGB.common$lfc= commongenes[,"EFGB"]
EFGB.common$lpv= commongenes[,"EFGB"]
EFGB.common
#write.csv(EFGB.common,file="EFGB_common_lfc.csv",row.names=F,quote=F)
#save(EFGB.common,file="EFGB_common_lfc.RData")
write.csv(EFGB.common,file="EFGB_common_lpv.csv",row.names=F,quote=F)
save(EFGB.common,file="EFGB_common_lpv.RData")

PRTER.common=data.frame("gene"=commongenes$gene)
#PRTER.common$lfc= commongenes[,"PRTER"]
PRTER.common$lpv= commongenes[,"PRTER"]
PRTER.common
#write.csv(PRTER.common,file="PRTER_common_lfc.csv",row.names=F,quote=F)
#save(PRTER.common,file="PRTER_common_lfc.RData")
write.csv(PRTER.common,file="PRTER_common_lpv.csv",row.names=F,quote=F)
save(PRTER.common,file="PRTER_common_lpv.RData")

WFGB.common=data.frame("gene"=commongenes$gene)
#WFGB.common$lfc= commongenes[,"WFGB"]
WFGB.common$lpv= commongenes[,"WFGB"]
WFGB.common
#write.csv(WFGB.common,file="WFGB_common_lfc.csv",row.names=F,quote=F)
#save(WFGB.common,file="WFGB_common_lfc.RData")
write.csv(WFGB.common,file="WFGB_common_lpv.csv",row.names=F,quote=F)
save(WFGB.common,file="WFGB_common_lpv.RData")