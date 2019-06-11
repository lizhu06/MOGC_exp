rm(list=ls())

library(cluster)
library(foreach)
library(doParallel)

setwd('/net/wong04/home/liz86/MOG_clustering/RawData')

MetaBric_gene <- get(load('MetaBric_gene.rdata'))
MetaBric_CNV <- get(load('MetaBric_CNV.rdata'))

WD <- "/net/wong04/home/liz86/MOG_clustering/Exp3/App3_metabric_v1/Data/"
setwd(WD)
#WD <- "Metabric_filterCutoff07_intersectKEGG"
#system(paste('mkdir -p',WD))
#setwd(WD)

gene_mean <- rowMeans(MetaBric_gene)
geneCutoff <- 0.7

## fileter by mean 
BRCA_Expr_m <- MetaBric_gene[gene_mean>quantile(gene_mean,geneCutoff),]
dim(BRCA_Expr_m) #5788 1981

## filter by sd
gene_sd <- apply(BRCA_Expr_m,1,sd)
BRCA_Expr_ms <- BRCA_Expr_m[gene_sd>quantile(gene_sd,geneCutoff),]
dim(BRCA_Expr_ms) #1737 1981


mRNAnames <- rownames(BRCA_Expr_ms)
CNVindex <- rownames(MetaBric_CNV) %in% mRNAnames
dim(BRCA_Expr_ms) # 1737 1981
BRCA_CNV_ms <- MetaBric_CNV[CNVindex,]
dim(BRCA_CNV_ms) #1587 1981

omicsTyps <- c(rep('gene',nrow(BRCA_Expr_ms)),rep('CNV',nrow(BRCA_CNV_ms)))
featureNames <- c(rownames(BRCA_Expr_ms),rownames(BRCA_CNV_ms))

d <- t(rbind(BRCA_Expr_ms,BRCA_CNV_ms))
dim(d) #1981 3324

## intersect with pathway data
## pathway dataset from Joey
load("/home06/liz86/BayesGL/Data/PathwayDB/KEGG_2016.RData")

length(unique(featureNames)) # 1737
sum(!(unique(featureNames) %in% 
	unique(unlist(KEGG_2016)))) # 960

# intersect with genes in the dataset
uni_genes <- unique(featureNames)
pathway_data_2 <- lapply(1:length(KEGG_2016), 
	function(x)intersect(KEGG_2016[[x]],uni_genes))
names(pathway_data_2) <- names(KEGG_2016)

## Calculate gene numbers in each pathway
pathway_length <- sapply(1:length(pathway_data_2),
	function(x)length(pathway_data_2[[x]]))
#summary(pathway_length)
table(pathway_length)
sum(pathway_length>0)  #285

## select a subset of pathways
pathway_data<-pathway_data_2[which(pathway_length>0)]
length(pathway_data)  #285
uni_genes<-unique(unlist(pathway_data))
length(uni_genes)  #777
sum(featureNames %in% uni_genes) #1490
path_length<-sapply(1:length(pathway_data),
  function(x)length(pathway_data[[x]]))

## data without duplication, only keep pathway genes
d_nodup <- d[,which(featureNames %in% uni_genes)]
featureNames_nodup <- featureNames[which(featureNames %in% uni_genes)]
type_nodup <- omicsTyps [which(featureNames %in% uni_genes)]
dim(d_nodup) #1981 1490
length(featureNames_nodup) #1490
length(type_nodup) #1490

num_dup <- sapply(1:ncol(d_nodup), function(x) 
	sum(unlist(pathway_data)==colnames(d_nodup)[x]))
summary(num_dup)

### load patients clinical data ####
load("/net/wong04/home/liz86/MOG_clustering/RawData/Complete_METABRIC_Clinical_Features_Data.rbin")
dim(Complete_METABRIC_Clinical_Features_Data)
all(rownames(Complete_METABRIC_Clinical_Features_Data) == rownames(d_nodup)) #1981   25
pam50 <- Complete_METABRIC_Clinical_Features_Data$NOT_IN_OSLOVAL_Pam50Subtype
table(pam50)

# remove NC
d_nodup <- d_nodup[which(pam50!="NC"), ]
pam50 <- pam50[which(pam50!="NC")]
dim(d_nodup) #1975 1490
length(pam50) # 1975
table(pam50)

## save data
save(d_nodup,file="d_nodup.RData")
save(featureNames_nodup,file="featureNames_nodup.RData") 
save(type_nodup,file="type_nodup.RData")
save(pam50, file="pam50.RData")
save(pathway_data, file="pathway_data_intersected.RData")



