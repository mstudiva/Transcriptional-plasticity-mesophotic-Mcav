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
counts = read.table("allcounts_trans_mcav.txt")

# how many genes we have total?
nrow(counts)
# 17901 for host, 13375 for symbiont
ncol(counts)
# 68 samples

# how does the data look?
head(counts)

#---------------------
# the following section prefilters the counts dataframe to remove genes with low counts (counts under 10 across all samples)
# while DESeq2 filters out low-count genes during testing anyways, prefiltering makes data transformations go faster

keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
# 9318 for host, 1628 for symbiont
ncol(countData)
write.csv(countData, file="countData.csv")

# for WCGNA: removing all genes with counts <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna)
# 2189 for host, 105 for symbiont
ncol(counts4wgcna)
write.csv(counts4wgcna, file="counts4wgcna.csv")

# importing a design .csv file
design = read.csv("design_trans.csv", head=TRUE)
design
str(design)

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time+depth+time:depth)

# reorders treatment factor according to "control" vs "treatment" levels
dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic","transplant"))
dds$time <- factor(dds$time, levels = c("zero","six","twelve"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and is recommended for datasets with more than 20 or so samples
Vsd=varianceStabilizingTransformation(dds)
#rl=rlog(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))
#e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("time","depth"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# not removing outliers to preserve experimental design
# if there were outliers, say, arrays 138, 168, and 170:
# outs=c(138, 168, 170)
# countData=countData[,-outs]
# counts4wgcna=counts4wgcna[,-outs]
# Vsd=Vsd[,-outs]
# rl=rl[,-outs]
# design=design[-outs,]

# not necessary if no outliers were found and removed
# remaking model with outliers removed from dataset
# dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time+depth+time:depth)
# dds$depth <- factor(dds$depth, levels = c("shallow","mesophotic","transplant"))

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
# takes the sample IDs, colony IDs, and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,3],design[,5],design[,8],sep=".")
# renames the column names
colnames(vsd)=snames
save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ time+depth+time:depth)
vsd.wg=assay(varianceStabilizingTransformation(wg))
head(vsd.wg)
# vsd.wg=assay(rlog(wg))
colnames(vsd.wg)=snames
save(vsd.wg,design,file="data4wgcna.RData")

#-------------------
# exploring similarities among samples

# heatmap and hierarchical clustering:
load("vsd.RData")

library(pheatmap)
# similarity among samples
pdf(file="heatmap_trans_mcav.pdf", width=20, height=20)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
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
# 1 PC for host, 3 PCs for symbiont

# plotting PCoA
pdf(file="PCoA_trans_mcav.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#225ea8","#41b6c4","#a1dab4"),pch=c(19,17,21)[as.numeric((as.factor(conditions$depth)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
# cluster overlay of time
ordiellipse(scores, conditions$time, label=F, draw= "polygon", col=c("#225ea8","#41b6c4","#a1dab4"))
legend("topleft", legend=c("zero","six","twelve"), fill = c("#225ea8","#41b6c4","#a1dab4"), bty="n")
legend("bottomleft", legend=c("mesophotic","shallow","transplant"), pch=c(19,17,21), bty="n")
# cluster overlay of depth
plot(scores[,1], scores[,2],col=c("coral","cyan3","coral"),pch=c(3,13,12)[as.numeric(as.factor(conditions$time))], xlab="Coordinate 1", ylab="Coordinate 2", main="Depth")
ordiellipse(scores, conditions$depth, label=F, draw= "polygon", col=c("coral","cyan3","coral"))
legend("bottomleft", legend=c("mesophotic","shallow","transplant"), fill = c("coral","cyan3","coral"), bty="n")
legend("topleft", legend=c("zero","six","twelve"), pch=c(3,13,12), bty="n")
dev.off()
# this will require some modification in Adobe Illustrator to make the transplant ellipse more transparent

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=15, height=20)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: time and time:depth interaction are significant
ad=adonis(t(vsd)~time*depth,data=conditions,method="manhattan")
ad
# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:4],labels=row.names(ad$aov.tab)[1:4],col=cols,main="time vs depth")
dev.off()

# DAPC analysis: creating discriminant function to tell transplant from shallow/mesophotic, based on factors
library(adegenet)
conditions$time<-factor(conditions$time, levels=c("zero","six","twelve"))
conditions[order(conditions$time),]

# runs simulations on randomly-chosen datasets of 90% of the total dataset to test the number of PCs to retain
set.seed(999)
# by depth, excluding transplants
xvalDapc(t(vsd[,conditions$depth!="transplant"]),conditions$depth[conditions$depth!="transplant"], n.rep=100, parallel="multicore", ncpus= 8)
# 25 PCs for host, 20 PCs for symbiont
# need to test again with a smaller range of possible PCs and more reps
# change the n.pca= statement depending on your target range
xvalDapc(t(vsd[,conditions$depth!="transplant"]),conditions$depth[conditions$depth!="transplant"], n.rep=1000, n.pca=15:35, parallel="multicore", ncpus= 8)
# 22 PCs for host, 19 PCs for symbiont

# now running the dapc without transplants
dp.d=dapc(t(vsd[,conditions$depth!="transplant"]),conditions$depth[conditions$depth!="transplant"],n.pca=19, n.da=1)

# can we predict depth treatment for the transplants?
pred.d=predict.dapc(dp.d,newdata=(t(vsd[,conditions$depth=="transplant"])))
pred.d
# look at the posterior section for assignments by probability

# creating a new dapc object to add in transplants for plotting
trans.d=dp.d
trans.d$ind.coord=pred.d$ind.scores
trans.d$posterior=pred.d$posterior
trans.d$assign=pred.d$assign
trans.d$grp<-as.factor(c("transplant","transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant"))
# same for symbionts

# now exporting side by side figures of controls vs transplants
# use Adobe Illustrator to overlay transplants curve onto controls
pdf(file="DAPC_trans_mcav_depth.pdf", width=12, height=6)
par(mfrow=c(1,2))
scatter(dp.d, bg="white",scree.da=FALSE,legend=TRUE,solid=0.6, col= c("coral","cyan3"))
scatter(trans.d, bg="white",scree.da=FALSE,legend=FALSE,solid=0.6, col= "coral")
dev.off()

# by time, excluding transplants
xvalDapc(t(vsd[,conditions$depth!="transplant"]),conditions$time[conditions$depth!="transplant"], n.rep=100, parallel="multicore", ncpus= 8)
# 10 PCs for host, 15 PCs for symbiont
# need to test again with a smaller range of possible PCs and more reps
# change the n.pca= statement depending on your target range
xvalDapc(t(vsd[,conditions$depth!="transplant"]),conditions$time[conditions$depth!="transplant"], n.rep=1000, n.pca=1:20, parallel="multicore", ncpus= 8)
# 10 PCs for host, 14 PCs for symbiont

# now running the dapc without transplants
dp.t=dapc(t(vsd[,conditions$depth!="transplant"]),conditions$time[conditions$depth!="transplant"],n.pca=14, n.da=2)

# can we predict depth treatment for the transplants?
pred.t=predict.dapc(dp.t,newdata=(t(vsd[,conditions$depth=="transplant"])))
pred.t
# look at the posterior section for assignments by probability

# exporting for significance testing below
dpc.d=data.frame(rbind(dp.d$ind.coord,pred.d$ind.scores))

write.csv(dpc.d, "DAPC_trans_mcav_DFA_depth.csv", quote=F)
# modify the output CSV to add in the respective columns for time and depth conditions

# then reimport
dpc.d<- read.csv("DAPC_trans_mcav_DFA_depth.csv")

# a little bit of rearranging
dpc.d$depth<-factor(dpc.d$depth, levels=c("transplant","mesophotic","shallow"))
dpc.d[order(dpc.d$depth),]
dpc.d$id<-as.factor(dpc.d$id)
str(dpc.d)

# testing significance of DFA differences with MCMCglmm
# install.packages("MCMCglmm")
library(MCMCglmm)

# sets prior distribution and creates a glm
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
glm.d <-MCMCglmm(LD1~depth,random=~id, family="gaussian", data=dpc.d,prior=prior,nitt=75000,thin=25,burnin=5000)
summary(glm.d)
# check to make sure you don't have autocorrelation with the reps (shown as "walks" in model traces)
plot(glm.d)

# calculating difference in magnitudes of depthmesophotic and depthshallow using sampled sets of parameters
transDelta=abs(glm.d$Sol[,"depthmesophotic"])-abs(glm.d$Sol[,"depthshallow"])
# 95% credible interval
HPDinterval(transDelta)

# MCMC p-value ie are transplants equally different from controls
if (is.na(table(transDelta<0)[2])) {
  cat("p <",signif(1/length(transDelta),1))
} else { cat("p =",signif(table(transDelta<0)[2]/length(transDelta),2)) }

#----------------------
# FITTING GENE BY GENE MODELS

# with multi-factor, multi-level design - using LRT

load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# Running DESeq with LRT for the effect of interaction (using full model dds = same as design argument on line 59, ~ time+depth+time:depth)
# then time+depth are removed, to look just at interaction term
dds.i=DESeq(dds,test="LRT",reduced=~time+depth, parallel=TRUE)

# Creating dataset with design =~time+depth (no interaction), to investigate importance of time and depth
dds1=DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ time+depth)
dds1$depth <- factor(dds1$depth, levels = c("shallow","mesophotic","transplant"))

# model for the effect of depth (>2 factor levels => LRT)
dds.depth=DESeq(dds1,test="LRT",reduced=~time, parallel=TRUE)

# model for the effect of time: (>2 factor levels => LRT)
dds.time=DESeq(dds1,test="LRT",reduced=~depth, parallel=TRUE)

# saving all models
save(dds,dds.time,dds.depth,dds.i,file="realModels.RData")

#--------------
# GLM analysis

load("realModels.RData")
library(DESeq2)

# depth treatment
depth=results(dds.depth)
summary(depth)
depth
degs.depth=row.names(depth)[depth$padj<0.1 & !(is.na(depth$padj))]

# time
time=results(dds.time)
summary(time)
time
degs.time=row.names(time)[time$padj<0.1 & !(is.na(time$padj))]

# depth:time interaction
int=results(dds.i)
summary(int)
int
degs.int=row.names(int)[int$padj<0.1 & !(is.na(int$padj))]

# depth contrasts
depth_m.s=results(dds.depth,contrast=c("depth","mesophotic","shallow")) # factor depth, level shallow compared to mesophotic
summary(depth_m.s) # look at how many DEGs there are and how many low-expressed genes ("low counts") had to be removed
depth_m.s
degs.depth_m.s=row.names(depth_m.s)[depth_m.s$padj<0.1 & !(is.na(depth_m.s$padj))]

depth_t.s=results(dds.depth,contrast=c("depth","transplant","shallow")) # factor depth, level transplant compared to shallow
summary(depth_t.s)
depth_t.s
degs.depth_t.s=row.names(depth_t.s)[depth_t.s$padj<0.1 & !(is.na(depth_t.s$padj))]

depth_t.m=results(dds.depth,contrast=c("depth","transplant","mesophotic")) # factor depth, level transplant compared to mesophotic
summary(depth_t.m)
depth_t.m
degs.depth_t.m=row.names(depth_t.m)[depth_t.m$padj<0.1 & !(is.na(depth_t.m$padj))]

# time contrasts
time_6.0=results(dds.time,contrast=c("time","six","zero"))
summary(time_6.0)
time_6.0
degs.time_6.0=row.names(time_6.0)[time_6.0$padj<0.1 & !(is.na(time_6.0$padj))]

time_12.0=results(dds.time,contrast=c("time","twelve","zero"))
summary(time_12.0)
time_12.0
degs.time_12.0 =row.names(time_12.0)[time_12.0$padj<0.1 & !(is.na(time_12.0$padj))]

time_12.6=results(dds.time,contrast=c("time","twelve","six"))
summary(time_12.6)
time_12.6
degs.time_12.6=row.names(time_12.6)[time_12.6$padj<0.1 & !(is.na(time_12.6$padj))]

# depth within time contrasts
# For figuring out your contrast statements, the following vignette is helpful
# ?results
# Contrast statements need to be a combination of those listed under resultsNames
resultsNames(dds)

zero_m.s=results(dds,contrast=c("depth", "mesophotic", "shallow"))
summary(zero_m.s)
zero_m.s
degs.zero_m.s=row.names(zero_m.s)[zero_m.s$padj<0.1 & !(is.na(zero_m.s$padj))]

zero_t.s=results(dds,contrast=c("depth", "transplant", "shallow"))
summary(zero_t.s)
zero_t.s
degs.zero_t.s=row.names(zero_t.s)[zero_t.s$padj<0.1 & !(is.na(zero_t.s$padj))]

zero_t.m=results(dds,contrast=c("depth", "transplant", "mesophotic"))
summary(zero_t.m)
zero_t.m
degs.zero_t.m=row.names(zero_t.m)[zero_t.m$padj<0.1 & !(is.na(zero_t.m$padj))]

six_m.s=results(dds,contrast=list("depth_mesophotic_vs_shallow","timesix.depthmesophotic"))
summary(six_m.s)
six_m.s
degs.six_m.s=row.names(six_m.s)[six_m.s$padj<0.1 & !(is.na(six_m.s$padj))]

six_t.s=results(dds,contrast=list("depth_transplant_vs_shallow","timesix.depthtransplant"))
summary(six_t.s)
six_t.s
degs.six_t.s=row.names(six_t.s)[six_t.s$padj<0.1 & !(is.na(six_t.s$padj))]

six_t.m=results(dds,contrast=list("timesix.depthtransplant","timesix.depthmesophotic"))
summary(six_t.m)
six_t.m
degs.six_t.m=row.names(six_t.m)[six_t.m$padj<0.1 & !(is.na(six_t.m$padj))]

twelve_m.s=results(dds,contrast=list("depth_mesophotic_vs_shallow","timetwelve.depthmesophotic"))
summary(twelve_m.s)
twelve_m.s
degs.twelve_m.s=row.names(twelve_m.s)[twelve_m.s$padj<0.1 & !(is.na(twelve_m.s$padj))]

twelve_t.s=results(dds,contrast=list("depth_transplant_vs_shallow","timetwelve.depthtransplant"))
summary(twelve_t.s)
twelve_t.s
degs.twelve_t.s=row.names(twelve_t.s)[twelve_t.s$padj<0.1 & !(is.na(twelve_t.s$padj))]

twelve_t.m=results(dds,contrast=list("timetwelve.depthtransplant","timetwelve.depthmesophotic"))
summary(twelve_t.m)
twelve_t.m
degs.twelve_t.m=row.names(twelve_t.m)[twelve_t.m$padj<0.1 & !(is.na(twelve_t.m$padj))]

save(depth,time,int,depth_m.s,depth_t.s,depth_t.m,time_6.0,time_12.0,time_12.6,zero_m.s,zero_t.s,zero_t.m,six_m.s,six_t.s,six_t.m,twelve_m.s,twelve_t.s,twelve_t.m,degs.depth,degs.time,degs.int,degs.depth_m.s,degs.depth_t.s,degs.depth_t.m,degs.time_6.0,degs.time_12.0,degs.time_12.6,degs.zero_m.s,degs.zero_t.s,degs.zero_t.m,degs.six_m.s,degs.six_t.s,degs.six_t.m,degs.twelve_m.s,degs.twelve_t.s,degs.twelve_t.m,file="pvals.RData")

#-------------------
# density plots: are my DEGs high-abundant or low-abundant?

load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs.time]),col="blue")
lines(density(means[degs.depth]),col="orange")
lines(density(means[degs.int]),col="lightblue")
legend("topright", title = "Factor", legend=c("time","depth","interaction"), fill = c("blue","orange","lightblue"))
dev.off()

#-------------------
# venn diagrams

load("pvals.RData")
library(DESeq2)

# overall factors, full model
candidates=list("time"=degs.time, "depth"=degs.depth, "int"=degs.int)

# DEGs across depths - are some genes conserved?
time0=list("m.s"=degs.zero_m.s,"t.s"=degs.zero_t.s,"t.m"=degs.zero_t.m)
time6=list("m.s"=degs.six_m.s,"t.s"=degs.six_t.s,"t.m"=degs.six_t.m)
time12=list("m.s"=degs.twelve_m.s,"t.s"=degs.twelve_t.s,"t.m"=degs.twelve_t.m)

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
	cat.dist = c(0.1, 0.1, 0.05),
	cat.pos = 3
	)
pdf(file="Venn_trans_mcav.pdf", height=12, width=12)
grid.draw(fullmodel.venn)
dev.off()

# DEGs across depths within time - are some genes conserved?
time0.venn=venn.diagram(
	x = time0,
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
pdf(file="Venn_time0_mcav.pdf", height=12, width=12)
grid.draw(time0.venn)
dev.off()

time6.venn=venn.diagram(
	x = time6,
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
pdf(file="Venn_time6_mcav.pdf", height=12, width=12)
grid.draw(time6.venn)
dev.off()

time12.venn=venn.diagram(
	x = time12,
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
pdf(file="Venn_time12_mcav.pdf", height=12, width=12)
grid.draw(time12.venn)
dev.off()

#-------------------
# saving data for GO and KOG analysis
load("realModels.RData")
load("pvals.RData")

# depth response
# signed log p-values: -log(pvalue)* direction:
head(depth)
source=depth[!is.na(depth$pvalue),]
depth.p=data.frame("gene"=row.names(source))
depth.p$lpv=-log(source[,"pvalue"],10)
depth.p$lpv[source$stat<0]=depth.p$lpv[source$stat<0]*-1
head(depth.p)
write.csv(depth.p,file="depth_lpv.csv",row.names=F,quote=F)
save(depth.p,file="depth_lpv.RData")

# time response
# signed log p-values: -log(pvalue)* direction:
head(time)
source=time[!is.na(time$pvalue),]
time.p=data.frame("gene"=row.names(source))
time.p$lpv=-log(source[,"pvalue"],10)
time.p$lpv[source$stat<0]=time.p$lpv[source$stat<0]*-1
head(time.p)
write.csv(time.p,file="time_lpv.csv",row.names=F,quote=F)
save(time.p,file="time_lpv.RData")

# time:depth interaction
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

# differences between depth factor levels
# log2 fold changes:
head(depth_m.s)
source=depth_m.s[!is.na(depth_m.s$pvalue),]
depth_m.s.fc=data.frame("gene"=row.names(source))
depth_m.s.fc$lfc=source[,"log2FoldChange"]
head(depth_m.s.fc)
write.csv(depth_m.s.fc,file="depth_m.s_fc.csv",row.names=F,quote=F)
save(depth_m.s.fc,file="depth_m.s_fc.RData")

head(depth_t.s)
source=depth_t.s[!is.na(depth_t.s$pvalue),]
depth_t.s.fc=data.frame("gene"=row.names(source))
depth_t.s.fc$lfc=source[,"log2FoldChange"]
head(depth_t.s.fc)
write.csv(depth_t.s.fc,file="depth_t.s_fc.csv",row.names=F,quote=F)
save(depth_t.s.fc,file="depth_t.s_fc.RData")

head(depth_t.m)
source= depth_t.m[!is.na(depth_t.m$pvalue),]
depth_t.m.fc=data.frame("gene"=row.names(source))
depth_t.m.fc$lfc=source[,"log2FoldChange"]
head(depth_t.m.fc)
write.csv(depth_t.m.fc,file="depth_t.m_fc.csv",row.names=F,quote=F)
save(depth_t.m.fc,file="depth_t.m_fc.RData")

#-------------------
load("realModels.RData")
load("pvals.RData")

# depth response within time zero
# log2 fold changes:
head(zero_m.s)
source=zero_m.s[!is.na(zero_m.s$pvalue),]
zero_m.s.fc=data.frame("gene"=row.names(source))
zero_m.s.fc$lfc=source[,"log2FoldChange"]
head(zero_m.s.fc)
write.csv(zero_m.s.fc,file="zero_m.s_fc.csv",row.names=F,quote=F)
save(zero_m.s.fc,file="zero_m.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
zero_m.s.p=data.frame("gene"=row.names(source))
zero_m.s.p$lpv=-log(source[,"pvalue"],10)
zero_m.s.p$lpv[source$stat<0]= zero_m.s.p$lpv[source$stat<0]*-1
head(zero_m.s.p)
write.csv(zero_m.s.p,file="zero_m.s_lpv.csv",row.names=F,quote=F)
save(zero_m.s.p,file="zero_m.s_lpv.RData")

head(zero_t.s)
source= zero_t.s[!is.na(zero_t.s$pvalue),]
zero_t.s.fc=data.frame("gene"=row.names(source))
zero_t.s.fc$lfc=source[,"log2FoldChange"]
head(zero_t.s.fc)
write.csv(zero_t.s.fc,file="zero_t.s_fc.csv",row.names=F,quote=F)
save(zero_t.s.fc,file="zero_t.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
zero_t.s.p=data.frame("gene"=row.names(source))
zero_t.s.p$lpv=-log(source[,"pvalue"],10)
zero_t.s.p$lpv[source$stat<0]= zero_t.s.p$lpv[source$stat<0]*-1
head(zero_t.s.p)
write.csv(zero_t.s.p,file="zero_t.s_lpv.csv",row.names=F,quote=F)
save(zero_t.s.p,file="zero_t.s_lpv.RData")

head(zero_t.m)
source= zero_t.m[!is.na(zero_t.m $pvalue),]
zero_t.m.fc=data.frame("gene"=row.names(source))
zero_t.m.fc$lfc=source[,"log2FoldChange"]
head(zero_t.m.fc)
write.csv(zero_t.m.fc,file="zero_t.m_fc.csv",row.names=F,quote=F)
save(zero_t.m.fc,file="zero_t.m_fc.RData")

# signed log p-values: -log(pvalue)* direction:
zero_t.m.p=data.frame("gene"=row.names(source))
zero_t.m.p$lpv=-log(source[,"pvalue"],10)
zero_t.m.p$lpv[source$stat<0]= zero_t.m.p$lpv[source$stat<0]*-1
head(zero_t.m.p)
write.csv(zero_t.m.p,file="zero_t.m_lpv.csv",row.names=F,quote=F)
save(zero_t.m.p,file="zero_t.m_lpv.RData")

# depth response within time six
# log2 fold changes:
head(six_m.s)
source=six_m.s[!is.na(six_m.s$pvalue),]
six_m.s.fc=data.frame("gene"=row.names(source))
six_m.s.fc$lfc=source[,"log2FoldChange"]
head(six_m.s.fc)
write.csv(six_m.s.fc,file="six_m.s_fc.csv",row.names=F,quote=F)
save(six_m.s.fc,file="six_m.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
six_m.s.p=data.frame("gene"=row.names(source))
six_m.s.p$lpv=-log(source[,"pvalue"],10)
six_m.s.p$lpv[source$stat<0]= six_m.s.p$lpv[source$stat<0]*-1
head(six_m.s.p)
write.csv(six_m.s.p,file="six_m.s_lpv.csv",row.names=F,quote=F)
save(six_m.s.p,file="six_m.s_lpv.RData")

head(six_t.s)
source= six_t.s[!is.na(six_t.s$pvalue),]
six_t.s.fc=data.frame("gene"=row.names(source))
six_t.s.fc$lfc=source[,"log2FoldChange"]
head(six_t.s.fc)
write.csv(six_t.s.fc,file="six_t.s_fc.csv",row.names=F,quote=F)
save(six_t.s.fc,file="six_t.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
six_t.s.p=data.frame("gene"=row.names(source))
six_t.s.p$lpv=-log(source[,"pvalue"],10)
six_t.s.p$lpv[source$stat<0]= six_t.s.p$lpv[source$stat<0]*-1
head(six_t.s.p)
write.csv(six_t.s.p,file="six_t.s_lpv.csv",row.names=F,quote=F)
save(six_t.s.p,file="six_t.s_lpv.RData")

head(six_t.m)
source= six_t.m[!is.na(six_t.m $pvalue),]
six_t.m.fc=data.frame("gene"=row.names(source))
six_t.m.fc$lfc=source[,"log2FoldChange"]
head(six_t.m.fc)
write.csv(six_t.m.fc,file="six_t.m_fc.csv",row.names=F,quote=F)
save(six_t.m.fc,file="six_t.m_fc.RData")

# signed log p-values: -log(pvalue)* direction:
six_t.m.p=data.frame("gene"=row.names(source))
six_t.m.p$lpv=-log(source[,"pvalue"],10)
six_t.m.p$lpv[source$stat<0]= six_t.m.p$lpv[source$stat<0]*-1
head(six_t.m.p)
write.csv(six_t.m.p,file="six_t.m_lpv.csv",row.names=F,quote=F)
save(six_t.m.p,file="six_t.m_lpv.RData")

# depth response within time twelve
# log2 fold changes:
head(twelve_m.s)
source=twelve_m.s[!is.na(twelve_m.s$pvalue),]
twelve_m.s.fc=data.frame("gene"=row.names(source))
twelve_m.s.fc$lfc=source[,"log2FoldChange"]
head(twelve_m.s.fc)
write.csv(twelve_m.s.fc,file="twelve_m.s_fc.csv",row.names=F,quote=F)
save(twelve_m.s.fc,file="twelve_m.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
twelve_m.s.p=data.frame("gene"=row.names(source))
twelve_m.s.p$lpv=-log(source[,"pvalue"],10)
twelve_m.s.p$lpv[source$stat<0]= twelve_m.s.p$lpv[source$stat<0]*-1
head(twelve_m.s.p)
write.csv(twelve_m.s.p,file="twelve_m.s_lpv.csv",row.names=F,quote=F)
save(twelve_m.s.p,file="twelve_m.s_lpv.RData")

head(twelve_t.s)
source= twelve_t.s[!is.na(twelve_t.s$pvalue),]
twelve_t.s.fc=data.frame("gene"=row.names(source))
twelve_t.s.fc$lfc=source[,"log2FoldChange"]
head(twelve_t.s.fc)
write.csv(twelve_t.s.fc,file="twelve_t.s_fc.csv",row.names=F,quote=F)
save(twelve_t.s.fc,file="twelve_t.s_fc.RData")

# signed log p-values: -log(pvalue)* direction:
twelve_t.s.p=data.frame("gene"=row.names(source))
twelve_t.s.p$lpv=-log(source[,"pvalue"],10)
twelve_t.s.p$lpv[source$stat<0]= twelve_t.s.p$lpv[source$stat<0]*-1
head(twelve_t.s.p)
write.csv(twelve_t.s.p,file="twelve_t.s_lpv.csv",row.names=F,quote=F)
save(twelve_t.s.p,file="twelve_t.s_lpv.RData")

head(twelve_t.m)
source= twelve_t.m[!is.na(twelve_t.m $pvalue),]
twelve_t.m.fc=data.frame("gene"=row.names(source))
twelve_t.m.fc$lfc=source[,"log2FoldChange"]
head(twelve_t.m.fc)
write.csv(twelve_t.m.fc,file="twelve_t.m_fc.csv",row.names=F,quote=F)
save(twelve_t.m.fc,file="twelve_t.m_fc.RData")

# signed log p-values: -log(pvalue)* direction:
twelve_t.m.p=data.frame("gene"=row.names(source))
twelve_t.m.p$lpv=-log(source[,"pvalue"],10)
twelve_t.m.p$lpv[source$stat<0]= twelve_t.m.p$lpv[source$stat<0]*-1
head(twelve_t.m.p)
write.csv(twelve_t.m.p,file="twelve_t.m_lpv.csv",row.names=F,quote=F)
save(twelve_t.m.p,file="twelve_t.m_lpv.RData")
