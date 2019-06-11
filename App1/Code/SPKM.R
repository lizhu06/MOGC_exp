rm(list=ls())
library(sparcl)
library(ISKmeans)
library(foreach)
library(doParallel)
library(mclust)
library(pROC)

server <- "/home/liz86/MOGC/"
load(paste0(server, "Exp3/App1_v3/Data/input_list.RData" ))

res_wd <- paste0(server, "Exp3/App1_v3/Results/")

#sapply(s_gene,dim)

len <- length(input_list)
cl <- makeCluster(len)
registerDoParallel(cl)

result <- foreach(s = 1:len) %dopar% {
	library(sparcl)

	data_name <- input_list[[s]]$data_name
	data <- t(input_list[[s]]$data) # sample by feature
	data_s <- scale(data)
	label <- input_list[[s]]$label
	pathway_name <- input_list[[s]]$pathway_name
	pathway_data <- input_list[[s]]$pathway_data

	# with standardization
	set.seed(20190314)
	tune <- KMeansSparseCluster.permute(data, K=3, nperms=25, 
	  wbounds=NULL, silent=TRUE, nvals=10, centers=NULL)
	res_spkm <- KMeansSparseCluster(data, K=3, wbounds=tune$bestw, 
	  nstart=20, silent =TRUE, maxiter=6, centers=NULL)
	#adjustedRandIndex(res_spkm[[1]]$Cs, label)
	#table(res_spkm[[1]]$Cs, label)
	#summary(res_spkm[[1]]$ws)
	save(res_spkm, file=paste0(res_wd, "res_spkm_", 
		data_name, "_", pathway_name, ".RData"))
}		
stopCluster(cl)			
