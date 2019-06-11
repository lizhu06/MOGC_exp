rm(list=ls())

setwd("/net/wong07/home/liz86/MOGC/Exp2/App1/Data/TCGA_fullData_processed")
load("Expr_full_processed.RData")
load("Methyl_full_processed.RData")
load("CNV_full_processed.RData")

#load("/home06/liz86/MetDiffNet/package/MetaDiffNetwork/data/pathway.database.v4.0.Rdata")


WD <- "/net/wong07/home/liz86/MOGC/Exp2/App1/Data/"
setwd(WD)
WD <- "TCGA_filterCutoff07_intersectKEGG"
system(paste('mkdir -p',WD))
setwd(WD)

#### filtering
gene_mean <- rowMeans(Expr_processed)
ncore <- 10
geneCutoff <- 0.7
	
## fileter by mean 
ExprC_m <- Expr_processed[gene_mean>quantile(gene_mean,geneCutoff),]
dim(ExprC_m) ## 6151  770

#quantile(gene_sd,geneCutoff)
gene_sd <- apply(ExprC_m,1,sd)
ExprC_ms <- ExprC_m[gene_sd>quantile(gene_sd,geneCutoff),]
dim(ExprC_ms) ## 1845 770

## filter CNV and methylation by gene expression
mRNAnames <- rownames(ExprC_ms)
CNVindex <- rownames(CNV_processed) %in% mRNAnames
methylindex <- rownames(Methyl_processed) %in% mRNAnames 

CNVC_ms <- CNV_processed[CNVindex,]
methylC_ms <- Methyl_processed[methylindex,]
#Methy_genes_ms <- uni_genes[methylindex]
dim(CNVC_ms) # 1761 770
dim(methylC_ms) # 1816 770

omicsTyps <- c(rep('gene',nrow(ExprC_ms)),
	rep('CNV',nrow(CNVC_ms)),rep('methyl',nrow(methylC_ms)))
featureNames <- c(rownames(ExprC_ms),rownames(CNVC_ms),rownames(methylC_ms))
table(omicsTyps)
#   CNV   gene methyl
#  1761   1845   1816
length(unique(featureNames)) # 1845
d <- t(rbind(ExprC_ms,CNVC_ms,methylC_ms))
dim(d) #770 5422

## intersect with pathway data
## pathway dataset from Joey
load("/home06/liz86/BayesGL/Data/PathwayDB/KEGG_2016.RData")
pathw_names <- names(KEGG_2016)
pathw_names[order(pathw_names)]
length(KEGG_2016[["ErbB signaling pathway_Homo sapiens_hsa04012"]]) # 87
length(KEGG_2016[["Estrogen signaling pathway_Homo sapiens_hsa04915"]]) # 99

length(unique(featureNames)) # 1845
sum(!(unique(featureNames) %in% 
	unique(unlist(KEGG_2016)))) # 999

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
length(uni_genes)  #846
sum(featureNames %in% uni_genes) #2507
path_length<-sapply(1:length(pathway_data),
  function(x)length(pathway_data[[x]]))
summary(path_length)

names(pathway_data)[which(path_length>20)]
pathway_data[["ErbB signaling pathway_Homo sapiens_hsa04012"]]
pathway_data[["Estrogen signaling pathway_Homo sapiens_hsa04915"]]

## data without duplication, only keep pathway genes
d_nodup <- d[,which(featureNames %in% uni_genes)]
featureNames_nodup <- featureNames[which(featureNames %in% uni_genes)]
type_nodup <- omicsTyps [which(featureNames %in% uni_genes)]
dim(d_nodup) #770 2507
length(featureNames_nodup) #2507
length(type_nodup) #2507
table(type_nodup)
#type_nodup
#   CNV   gene methyl
#   825    846    836

num_dup <- sapply(1:ncol(d_nodup), function(x) 
	sum(unlist(pathway_data)==colnames(d_nodup)[x]))
summary(num_dup)

### load patients clinical data ####
patients <- rownames(d_nodup)
pat_names <- substr(patients, start=1, stop=12)
rownames(d_nodup) <- pat_names

load("/home05/liz86/Steffi/Kevin_IDC_ILC_DE/Data/TCGA_PAM50_histology_TumorPurity.rda")
pam50 <- TCGA_PAM50[pat_names, "final_assign"]

## save data
save(d_nodup,file="d_nodup.RData")
save(featureNames_nodup,file="featureNames_nodup.RData")
save(type_nodup,file="type_nodup.RData")
save(pam50, file="pam50.RData")
save(pathway_data, file="pathway_data_intersected.RData")

