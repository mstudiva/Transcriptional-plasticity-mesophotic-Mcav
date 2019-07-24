setwd("/Users/Mike/Documents/HBOI/Data/RNA_Seq/Analysis/KOG-MWU/enviro/mcav")
getwd()
# install.packages("KOGMWU")
library(KOGMWU)

#---------------------------
# full model, assessing changes in KOG expression across site and depth

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

add=load('depth_fc.RData')
add # names of datasets in the package
fc.d=kog.mwu(depth.fc,gene2kog) 
fc.d 

ads=load('site_lpv.RData')
ads # names of datasets in the package
lpv.s=kog.mwu(site.p,gene2kog) 
lpv.s 

adi=load('int_lpv.RData')
adi # names of datasets in the package
lpv.i=kog.mwu(int.p,gene2kog) 
lpv.i

barplot(fc.d$delta.rank,horiz=T)

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Site"=lpv.s,"Depth"=fc.d,"Site:Depth"=lpv.i))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdYlBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="enviro_mcav_heatmap_enviro_fc.pdf", width=6, height=7)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pdf(file="enviro_mcav_corrplot_enviro_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
while (!is.null(dev.list()))  dev.off()
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pdf(file="enviro_mcav_pvalplot_enviro_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
while (!is.null(dev.list()))  dev.off()

# plotting individual delta-rank correlations:
# produces individual plots from the grid above
# corrPlot(x="Site",y="Depth",ktable)
# corrPlot(x="Site",y="Site:Depth",ktable)
# corrPlot(x="Depth",y="Site:Depth",ktable)

#---------------------------
# assessing KOG differences across depth within site

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

blz=load('BLZ_fc.RData')
blz # names of datasets in the package
BLZ=kog.mwu(BLZ.fc,gene2kog) 
BLZ 

wfgb=load('WFGB_fc.RData')
wfgb # names of datasets in the package
WFGB=kog.mwu(WFGB.fc,gene2kog) 
WFGB 

efgb=load('EFGB_fc.RData')
efgb # names of datasets in the package
EFGB=kog.mwu(EFGB.fc,gene2kog) 
EFGB 

prg_drt=load('PRG-DRT_fc.RData')
prg_drt # names of datasets in the package
PRG_DRT=kog.mwu(PRG_DRT.fc,gene2kog) 
PRG_DRT 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("BLZ"=BLZ,"WFGB"=WFGB,"EFGB"=EFGB,"PRG-DRT"=PRG_DRT))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdYlBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="enviro_mcav_heatmap_site_fc.pdf", width=6, height=7)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pdf(file="enviro_mcav_corrplot_site_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
while (!is.null(dev.list()))  dev.off()
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pdf(file="enviro_mcav_pvalplot_site_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
while (!is.null(dev.list()))  dev.off()

#---------------------------
# assessing KOG differences within conserved genes across depth and site

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

blz.common=load('BLZ_common_fc.RData')
blz.common # names of datasets in the package
BLZ.common=kog.mwu(BLZ.common,gene2kog) 
BLZ.common 

wfgb.common=load('WFGB_common_fc.RData')
wfgb.common # names of datasets in the package
WFGB.common=kog.mwu(WFGB.common,gene2kog) 
WFGB.common 

efgb.common=load('EFGB_common_fc.RData')
efgb.common # names of datasets in the package
EFGB.common=kog.mwu(EFGB.common,gene2kog) 
EFGB.common 

prg_drt.common=load('PRG-DRT_common_fc.RData')
prg_drt.common # names of datasets in the package
PRG_DRT.common=kog.mwu(PRG_DRT.common,gene2kog) 
PRG_DRT.common 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("BLZ"=BLZ.common,"WFGB"=WFGB.common,"EFGB"=EFGB.common,"PRG-DRT"=PRG_DRT.common))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdYlBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="enviro_mcav_heatmap_common_fc.pdf", width=6, height=7)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pdf(file="enviro_mcav_corrplot_common_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
while (!is.null(dev.list()))  dev.off()
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pdf(file="enviro_mcav_pvalplot_common_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
while (!is.null(dev.list()))  dev.off()
