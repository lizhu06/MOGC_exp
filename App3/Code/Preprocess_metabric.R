rm(list=ls())
setwd('/net/wong04/home/liz86/MOG_clustering/RawData')

#source("http://bioconductor.org/biocLite.R")
##biocLite("illuminaHumanv3.db")
library("illuminaHumanv3.db")

load("Complete_METABRIC_Expression_Data.rbin")
rawMetaBric = exprs(Complete_METABRIC_Expression_Data)
dim(rawMetaBric)	## 49576  1981

probeNames = rownames(rawMetaBric)
geneNames = unlist(mget(probeNames,illuminaHumanv3SYMBOL))
geneNameSplit = split(1:length(geneNames),geneNames)
length(geneNameSplit) ## 19292


IQR2 <- function(amat){
	if(is.null(dim(amat)))
		return(amat)
	iqr = apply(amat,1,IQR)
	return(amat[which.max(iqr),])
}
MetaBric_gene0 <- sapply(geneNameSplit,function(x) IQR2(rawMetaBric[x,]))
MetaBric_gene = t(MetaBric_gene0)
dim(MetaBric_gene) ## 19292  1981

save(MetaBric_gene,file='MetaBric_gene.rdata')

## load metabric clinical 
load("Complete_METABRIC_Clinical_Features_Data.rbin")
dim(Complete_METABRIC_Clinical_Features_Data)
names(Complete_METABRIC_Clinical_Features_Data)

table(Complete_METABRIC_Clinical_Features_Data$histological_type)

#Complete_METABRIC_Clinical_Survival_Data__DSS.rbin
#Complete_METABRIC_Clinical_Survival_Data_OS.rbin

## load survival
load("Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")
library(survival)
y10 = 365*10
y10index = Complete_METABRIC_Clinical_Survival_Data__DSS[,1]>y10
Complete_METABRIC_Clinical_Survival_Data__DSS10years = Complete_METABRIC_Clinical_Survival_Data__DSS
Complete_METABRIC_Clinical_Survival_Data__DSS10years[y10index,1] = y10
Complete_METABRIC_Clinical_Survival_Data__DSS10years[y10index,2] = 0

save(Complete_METABRIC_Clinical_Survival_Data__DSS10years,file='Complete_METABRIC_Clinical_Survival_Data__DSS10years.rdata')


## load copy number
load('Complete_METABRIC_Copy_Number_Data.rbin')
#source("http://bioconductor.org/biocLite.R")
#biocLite("pd.genomewidesnp.6")
#source("http://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")

library(org.Hs.eg.db)
# wong04: zlibbioc affyio oligoClasses oligo 

CNV.array = exprs(Complete_METABRIC_Copy_Number_Data)
dim(CNV.array)	## 18538  1981
entrezid = unlist(lapply(strsplit(row.names(CNV.array),"_"),function(x) x[1]))
Gene.symbol = unlist(mget(entrezid, org.Hs.egSYMBOL, ifnotfound=NA))
row.names(CNV.array) = unlist(Gene.symbol)

sum(duplicated(row.names(CNV.array))) #44

CNVNameSplit = split(1:length(row.names(CNV.array)),row.names(CNV.array))
length(CNVNameSplit) ## 18493

IQR2 <- function(amat){
	if(is.null(dim(amat)))
		return(amat)
	iqr = apply(amat,1,IQR)
	return(amat[which.max(iqr),])
}

MetaBric_CNV0 <- sapply(CNVNameSplit,function(x) IQR2(CNV.array[x,]))
MetaBric_CNV = t(MetaBric_CNV0)
dim(MetaBric_CNV) ## 18493  1981


save(MetaBric_CNV,file='MetaBric_CNV.rdata')



