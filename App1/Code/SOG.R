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
len <- length(input_list)
#len <- 3
cl <- makeCluster(3)
registerDoParallel(cl)
source(paste0(server, "Code/Function/SOGC_v16_callC.R"))

result <- foreach(s = 1:len) %dopar% {
	library(sparcl)

	data_name <- input_list[[s]]$data_name
	data <- t(input_list[[s]]$data) # sample by feature
	data_s <- scale(data)
	label <- input_list[[s]]$label
	pathway_name <- input_list[[s]]$pathway_name
	pathway_data <- input_list[[s]]$pathway_data

	# create U1
	U1 <- matrix(0, ncol(data_s), length(pathway_data))
	colnames(U1) <- names(pathway_data)
	for(i in 1:length(pathway_data)){
	  U1[which(colnames(data_s) %in% pathway_data[[i]]), i] <- 1
	}

	# load SPKM data
	load(paste0(res_wd, "res_spkm_", 
		data_name, "_", pathway_name, ".RData"))

	#NSIM <- 30000
	burnInIter <- 1000
	keepIter <- 3000
	maxIter <- 10000
	print_int <- 10000
	K <- 3
	t0 <- proc.time()
	res_sogc <- SOGC(data_s, U1, K=K, 
		center=FALSE, burnInIter=burnInIter, keepIter=keepIter, 
		maxIter=maxIter, 
	  	print_int=print_int, debug=FALSE, init_mu=NULL, 
	  	init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
	  	fix_mu=FALSE, adj_ls=TRUE,
	  	seed=123, BernoulliWeighted=TRUE, 
	  	MH_ind=1, pi_j_prop_n=10, recover_overlap=TRUE,
	   cppfile=paste0(server, "Code/Function/SOGC_v16.cpp"),
	   s2=NULL, types=NULL, report_dup_mean=TRUE)
	t1 <- proc.time()-t0

	save(res_sogc, t1, file=paste0(res_wd, "res_sogc_", 
		data_name, "_", pathway_name, ".RData"))

}		
stopCluster(cl)			
