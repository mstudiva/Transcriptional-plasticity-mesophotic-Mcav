
# install.fcackages("KOGMWU")
library(KOGMWU)

# assessing KOG differences across depth within times

# loading KOG annotations
gene2kog=read.table("Mcavernosa_Cladocopium_iso2kogClass.tab",sep="\t")
head(gene2kog)

zero_m.s=load('zero_m.s_fc.RData')
zero_m.s # names of datasets in the package
zms=kog.mwu(zero_m.s.fc,gene2kog) 
zms 

zero_t.s=load('zero_t.s_fc.RData')
zero_t.s # names of datasets in the package
zts=kog.mwu(zero_t.s.fc,gene2kog) 
zts 

zero_t.m=load('zero_t.m_fc.RData')
zero_t.m # names of datasets in the package
ztm=kog.mwu(zero_t.m.fc,gene2kog) 
ztm 

six_m.s=load('six_m.s_fc.RData')
six_m.s # names of datasets in the package
sms=kog.mwu(six_m.s.fc,gene2kog) 
sms 

six_t.s=load('six_t.s_fc.RData')
six_t.s # names of datasets in the package
sts=kog.mwu(six_t.s.fc,gene2kog) 
sts 

six_t.m=load('six_t.m_fc.RData')
six_t.m # names of datasets in the package
stm=kog.mwu(six_t.m.fc,gene2kog) 
stm 

twelve_m.s=load('twelve_m.s_fc.RData')
twelve_m.s # names of datasets in the package
tms=kog.mwu(twelve_m.s.fc,gene2kog) 
tms 

twelve_t.s=load('twelve_t.s_fc.RData')
twelve_t.s # names of datasets in the package
tts=kog.mwu(twelve_t.s.fc,gene2kog) 
tts 

twelve_t.m=load('twelve_t.m_fc.RData')
twelve_t.m # names of datasets in the package
ttm=kog.mwu(twelve_t.m.fc,gene2kog) 
ttm 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("0 mo m.s"=zms,"0 mo t.s"=zts,"0 mo t.m"=ztm,"6 mo m.s"=sms,"6 mo t.s"=sts,"6 mo t.m"=stm,"12 mo m.s"=tms,"12 mo t.s"=tts,"12 mo t.m"=ttm))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdYlBu"),"royalblue","darkblue")))(100)
  
# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG-MWU_heatmap_transplant_depth_fc.pdf", width=10, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color) 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pdf(file="KOG-MWU_corrplot_transplant_depth_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
while (!is.null(dev.list()))  dev.off()
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pdf(file="KOG-MWU_pvalplot_transplant_depth_fc.pdf", width=6, height=6)
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
while (!is.null(dev.list()))  dev.off()