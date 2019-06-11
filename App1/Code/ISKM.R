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
cl <- makeCluster(len)
registerDoParallel(cl)

result <- foreach(s = 1:len) %dopar% {

	library(ISKmeans) 
	data_name <- input_list[[s]]$data_name
	data <- t(input_list[[s]]$data) # sample by feature
	data_s <- scale(data)
	label <- input_list[[s]]$label
	pathway_name <- input_list[[s]]$pathway_name
	pathway_data <- input_list[[s]]$pathway_data

	# with standardization
	set.seed(20190320)
	groupset <- lapply(1:length(pathway_data), function(x) 
		match(pathway_data[[x]], colnames(data)))
	alpha_vec <- c(0.5, 0.1, 0.01)
	#alpha_vec <- 0.5

	for(av in 1:length(alpha_vec)){
	  alpha <- alpha_vec[av]
	  #gammas <- seq(0.1,0,0.1)

	  sparseStart <- TRUE
	  t0 <- proc.time()
	  gapStat = gapStatistics(data_s, K=3, B=25, gamma=NULL, 
	    alpha=alpha, group=groupset, silence=TRUE)
	  res_iskm <- ISKmeans(data_s, K=3, gamma=gapStat$bestGamma, alpha=alpha,
	    group=groupset, nstart=20, silent=TRUE, maxiter=20,
	    sparseStart=sparseStart)
	  delta_t <- proc.time() - t0
	  save(res_iskm, delta_t, file=paste0(res_wd, "res_iskm_alpha_",alpha, "_", 
	  	data_name, "_", pathway_name, ".RData"))
	}
}		
stopCluster(cl)			
