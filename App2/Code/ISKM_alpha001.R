rm(list=ls())
server <- "/home/liz86/MOGC/"
wd <- paste0(server, "Exp/App2/Data/")
load(paste0(wd, "d_nodup.RData"))
load(paste0(wd, "featureNames_nodup.RData"))
load(paste0(wd, "type_nodup.RData"))
load(paste0(wd, "pam50.RData"))
load(paste0(wd, "pathway_data_intersected.RData"))

res_wd <- paste0(server, "Exp/App2/Results/")
system(paste('mkdir -p',res_wd))

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

## ISKM
groupset <- lapply(1:ncol(U1), function(x) which(U1[,x]==1))

alpha <- 0.01
#gammas <- seq(0.1,0,0.1)
sparseStart <- TRUE
gapStat = gapStatistics(X_s, K=3, B=25, gamma=NULL, 
  alpha=alpha, group=groupset, silence=TRUE)
res_iskm_alpha001 <- ISKmeans(X_s, K=3, gamma=gapStat$bestGamma, 
  alpha=alpha,
  group=groupset, nstart=20, silent=TRUE, maxiter=20,
  sparseStart=sparseStart)
save(res_iskm_alpha001, file=paste0(res_wd, "res_iskm_alpha001.RData"))



