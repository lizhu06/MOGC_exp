rm(list=ls())
server <- "/net/wong07/home/liz86/MOGC/"
wd <- "/net/wong04/home/liz86/MOG_clustering/Exp3/App3_metabric_v1/Data/"
load(paste0(wd, "d_nodup.RData"))
load(paste0(wd, "featureNames_nodup.RData"))
load(paste0(wd, "type_nodup.RData"))
load(paste0(wd, "pam50.RData"))
load(paste0(wd, "pathway_data_intersected.RData"))

res_wd <- "/net/wong04/home/liz86/MOG_clustering/Exp3/App3_metabric_v1/Output/"

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

source(paste0(server, "Code/Function/MOGC_DPT_v6_callC.R"))

## load SPKM results
load(paste0(res_wd, "res_spkm.RData"))

## MOGC_DP
burnInIter <- 1000
keepIter <- 3000
maxIter <- 10000
if(FALSE){
  burnInIter <- 5
  keepIter <- 10
  maxIter <- 30
}
K <- 5
res_mogc_dp <- MOGC_DPT(X_s, U1, U2, K=K, center=FALSE, 
  burnInIter=burnInIter, keepIter=keepIter, maxIter=maxIter,  
  print_int=1000, debug=FALSE, init_mu=NULL, init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
  fix_mu=FALSE, adj_ls=TRUE, 
  BernoulliWeighted=TRUE, MH_ind=1, pi_g_prop_n=10,
  pi_j_prop_n=10, 
  cppfile=paste0(server, "Code/Function/MOGC_DPT_v6.cpp"), 
  recover_overlap=TRUE, types=types)
save(res_mogc_dp, file=paste0(res_wd, "res_mogc_dp.RData"))
#AUC_mog_dp <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))  
#ARI_mog_dp <- adjustedRandIndex(res_mogc_dp$label_map,pam50)








