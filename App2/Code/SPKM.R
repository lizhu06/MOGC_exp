rm(list=ls())
server <- "/net/wong07/home/liz86/MOGC/"
wd <- "/net/wong07/home/liz86/MOGC/Exp2/App1/Data/TCGA_filterCutoff07_intersectKEGG/"
load(paste0(wd, "d_nodup.RData"))
load(paste0(wd, "featureNames_nodup.RData"))
load(paste0(wd, "type_nodup.RData"))
load(paste0(wd, "pam50.RData"))
load(paste0(wd, "pathway_data_intersected.RData"))

res_wd <- "/net/wong07/home/liz86/MOGC/Exp2/App1/Results/TCGA_filterCutoff07_intersectKEGG/"
system(paste('mkdir -p',res_wd))

# pam50 genes
library("genefu")
library("genefilter")
data(pam50.robust) # pam50 genes
pam50_gene <- rownames(pam50.robust$centroids)
gene_label <- colnames(d_nodup) %in% pam50_gene

# Create U1 and U2
X <- d_nodup
X_s <- scale(X)
feature_types <- sapply(1:ncol(X), function(i) paste( 
  colnames(d_nodup)[i],type_nodup[i],sep="_"))
types <- as.numeric(as.factor(type_nodup))

## create U1 and U2
uni_genes <- unique(colnames(d_nodup))
U1 <- matrix(0, ncol(X_s), length(uni_genes))
colnames(U1) <- uni_genes
for(i in 1:ncol(X_s)){
  U1[i, colnames(X_s)[i]] <- 1
}
U2 <- matrix(0, length(uni_genes), length(pathway_data))
rownames(U2) <- uni_genes
for(i in 1:length(pathway_data)){
  U2[pathway_data[[i]], i]<-1
}

library(MASS)
library(sparcl)
library(mclust)
library(pROC)
library(label.switching)
library(ggplot2)
library(MCMCpack)
library(ISKmeans) 


## SPKM (tune use gap statistics)##
set.seed(20181027)
t0 <- proc.time()
tune <- KMeansSparseCluster.permute(X_s, K=3, nperms=25, 
  wbounds=NULL, silent=TRUE, nvals=10, centers=NULL)
save(tune, file=paste0(res_wd, "res_spkm_tune.RData"))

res_spkm <- KMeansSparseCluster(X_s, K=3, wbounds=tune$bestw, 
  nstart=20, silent =TRUE, maxiter=6, centers=NULL)
save(res_spkm, file=paste0(res_wd, "res_spkm.RData"))

pam50_2 <- as.character(pam50)
pam50_2[which(pam50_2 %in% c("LumA", "LumB", "Normal"))] <- "Luminal"
(ARI_spkm <- adjustedRandIndex(res_spkm[[1]]$Cs, pam50_2))   # 0.39
table(res_spkm[[1]]$Cs, pam50)
(AUC_spkm <- auc(roc(factor(gene_label), res_spkm[[1]]$ws)))   #0.53


