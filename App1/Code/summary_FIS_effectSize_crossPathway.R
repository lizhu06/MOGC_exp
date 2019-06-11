rm(list=ls())
library(sparcl)
library(ISKmeans)
library(foreach)
library(doParallel)
library(mclust)
library(pROC)
library(GSA)
server <- "/home/liz86/"
res_wd <- paste0(server, "MOGC/Exp3/App1_v3/Results/")
load(paste0(server, "MOGC/Exp3/App1_v3/Data/input_list.RData"))

# load pathways
pathway_wd <- paste0(server, "Dataset/pathway_data_from_Caleb/")
pathwayFileList <- c('c2.cp.biocarta.v5.2.symbols.gmt.txt',
					'c2.cp.kegg.v5.2.symbols.gmt.txt',
					'c2.cp.reactome.v5.2.symbols.gmt.txt')
all_pathways <- list()
for(i in 1:length(pathwayFileList)){
	database = GSA.read.gmt(paste0(pathway_wd, pathwayFileList[i]))
	database2 = database$genesets
	names(database2) = database$geneset.names
	all_pathways[[i]] <- database2
}
names(all_pathways) <- c("biocarta", "kegg", "reactome")
sapply(all_pathways, length)

# pathway analysis function
source(paste0(server, "Functions/Pathway_Fisher_Function.R"))
pathway_enrichment <- function(p, mu=NULL, all_genes, 
	num_genes, pathway_data, 
	pathway_lowerBound=14, pathway_upperBound=101){
	if(FALSE){
		weight <- res_spkm[[1]]$ws
		all_genes <- colnames(d_nodup)
		num_genes <- 1000
	}
	if(is.null(mu)){
		sorted_gene_index <- seq(1,length(all_genes)) [order(
			p, decreasing=TRUE)]
	}else{
		sorted_gene_index <- seq(1,length(all_genes)) [order(
			p, mu, decreasing=TRUE)]
	}
	feature_rank <- match(seq(1, length(all_genes)), sorted_gene_index)
	sig <- all_genes[feature_rank <= num_genes]
	pathwayRes <- pathway_fisher_detail_gene(significant=sig,
	  whole=all_genes, fdr=0.05, database=pathway_data,
	  pathway_lowerBound=pathway_lowerBound, pathway_upperBound=pathway_upperBound)
	return(pathwayRes)
}

ari_mat <- matrix(NA, 3, 6)
num_sig_path <- matrix(NA, 3, 14)
num_genes_select_vec <- rep(NA, 3)
num_sig_path_cross_validation <- list()
for(i in 1:3){
	num_sig_path_cross_validation[[i]] <- matrix(NA, 3, 14)
}
names(num_sig_path_cross_validation) <- names(all_pathways)

for(s in 1:3){
	# load data
	data_name <- input_list[[s]]$data_name
	data <- t(input_list[[s]]$data) # sample by feature
	data_s <- scale(data)
	label <- input_list[[s]]$label
	pathway_name <- input_list[[s]]$pathway_name
	pathway_data <- input_list[[s]]$pathway_data
	num_genes_select <- num_genes_select_vec[s] <- round(ncol(data)/3)

	# SPKM
	load(paste0(res_wd, "res_spkm_", 
			data_name, "_", pathway_name, ".RData"))
	ari_mat[s,1] <- adjustedRandIndex(res_spkm[[1]]$Cs, label)

	spkm_path_res <- pathway_enrichment(p=res_spkm[[1]]$ws, mu=NULL, 
		colnames(data_s), num_genes_select, pathway_data)
	num_sig_path[s,1] <- sum(unlist(spkm_path_res[,2]) < 0.05)
	#write.csv(spkm_path_res, file=paste0(res_wd, "PathwayRes_spkm_", 
	#		data_name, "_", pathway_name, ".csv"))

	for(i in 1:length(all_pathways)){
		pathway_res <- pathway_enrichment(p=res_spkm[[1]]$ws, mu=NULL, 
				colnames(data_s), num_genes_select, all_pathways[[i]],
				14, 101)
		num_sig_path_cross_validation[[i]][s, 1] <- sum(unlist(pathway_res[,2]) < 0.05)
		write.csv(pathway_res, file=paste0(res_wd, "PathwayRes_spkm_", 
				data_name, "_", pathway_name, "_testPath_",
				names(all_pathways)[i], "_crossValidation.csv"))
	}

	# ISKM
	alpha_vec <- c(0.5, 0.1, 0.01)
	for(alpha_index in 1:length(alpha_vec)){
		load(paste0(res_wd, "res_iskm_alpha_",alpha_vec[alpha_index], "_", 
				data_name, "_", pathway_name, ".RData"))
		ari_mat[s,alpha_index+1] <- adjustedRandIndex(res_iskm$Cs, label)
		iskm_path_res <- pathway_enrichment(p=res_iskm$ws, mu=NULL, 
			colnames(data_s), num_genes_select, pathway_data)
		num_sig_path[s,alpha_index+1] <- sum(unlist(iskm_path_res[,2]) < 0.05)
		#write.csv(iskm_path_res, file=paste0(res_wd, "PathwayRes_iskm_alpha_",
		#	alpha_vec[alpha_index], "_",  data_name, "_", pathway_name, ".csv"))

		for(i in 1:length(all_pathways)){
			pathway_res <- pathway_enrichment(p=res_iskm$ws, mu=NULL, 
					colnames(data_s), num_genes_select, all_pathways[[i]],
					14, 101)
			num_sig_path_cross_validation[[i]][s, alpha_index+1] <- sum(unlist(pathway_res[,2]) < 0.05)
			write.csv(pathway_res, file=paste0(res_wd, "PathwayRes_iskm_alpha_",
				alpha_vec[alpha_index], "_", data_name, "_", pathway_name,
				 "_testPath_", names(all_pathways)[i],  "_crossValidation.csv"))
		}
	}

	# SOG
	load(paste0(res_wd, "res_sogc_", 
			data_name, "_", pathway_name, ".RData"))
	ari_mat[s,5] <- adjustedRandIndex(res_sogc$label_map, label)
	ave_p_sog <- sapply(1:3, function(k) 
		apply(res_sogc$MU_store[, k, ]!=0, 2, mean))
	ave_mu_sog <- sapply(1:3, function(k) 
		abs(apply(res_sogc$MU_store[, k, ], 2, mean)))
	for(k in 1:3){
		sog_path_res <- pathway_enrichment(p=ave_p_sog[, k], 
			mu=ave_mu_sog[, k], 
			colnames(data_s), num_genes_select, pathway_data)
		num_sig_path[s,k+4] <- sum(unlist(sog_path_res[, 2]) < 0.05)
	}
	ave_mu_max_sog <- apply(ave_mu_sog, 1, max)
	sog_comb_path_res <- pathway_enrichment(p=res_sogc$feature_select_prob, 
		mu=ave_mu_max_sog, colnames(data_s), num_genes_select, pathway_data)
	#write.csv(sog_comb_path_res, file=paste0(res_wd, "PathwayRes_sogc_",
	#	data_name, "_", pathway_name, ".csv"))
	num_sig_path[s, 8] <- sum(unlist(sog_comb_path_res[,2]) < 0.05)

	for(i in 1:length(all_pathways)){
		pathway_res <- pathway_enrichment(p=res_sogc$feature_select_prob, mu=ave_mu_max_sog,
				colnames(data_s), num_genes_select, all_pathways[[i]],
				14, 101)
		num_sig_path_cross_validation[[i]][s, 8] <- sum(unlist(pathway_res[,2]) < 0.05)
		write.csv(pathway_res, file=paste0(res_wd, "PathwayRes_sogc_", 
				data_name, "_", pathway_name, "_testPath_", names(all_pathways)[i],  "_crossValidation.csv"))
	}

	# SOG_dp
	load(paste0(res_wd, "res_sogc_dp_", 
			data_name, "_", pathway_name, ".RData"))
	ari_mat[s,6] <- adjustedRandIndex(res_sogc_dp$label_map, label)
	ave_p_sog_dp <- sapply(1:5, function(k) 
		apply(res_sogc_dp$MU_store[, k, ]!=0, 2, mean))
	ave_mu_sog_dp <- sapply(1:5, function(k) 
		abs(apply(res_sogc_dp$MU_store[, k, ], 2, mean)))
	for(k in 1:5){
		sog_dp_path_res <- pathway_enrichment(p=ave_p_sog_dp[,k], 
			mu=ave_mu_sog_dp[,k], colnames(data_s),
			num_genes_select, pathway_data)
		num_sig_path[s,k+8] <- sum(unlist(sog_dp_path_res[, 2]) < 0.05)
	}
	ave_mu_max_sog_dp <- apply(ave_mu_sog_dp, 1, max)
	sog_dp_comb_path_res <- pathway_enrichment(p=res_sogc_dp$feature_select_prob, 
		mu=ave_mu_max_sog_dp, colnames(data_s), num_genes_select, pathway_data)
	#write.csv(sog_dp_comb_path_res, file=paste0(res_wd, "PathwayRes_sogc_dp_",
	#	data_name, "_", pathway_name, ".csv"))
	num_sig_path[s, 14] <- sum(unlist(sog_dp_comb_path_res[,2]) < 0.05)

	for(i in 1:length(all_pathways)){
		pathway_res <- pathway_enrichment(p=res_sogc_dp$feature_select_prob, mu=ave_mu_max_sog_dp,
				colnames(data_s), num_genes_select, all_pathways[[i]],
				14, 101)
		num_sig_path_cross_validation[[i]][s, 14] <- sum(unlist(pathway_res[,2]) < 0.05)
		write.csv(pathway_res, file=paste0(res_wd, "PathwayRes_sogc_dp_", 
				data_name, "_", pathway_name,"_testPath_", names(all_pathways)[i],  "_crossValidation.csv"))
	}
}

colnames(ari_mat) <- c("SPKM",
	"ISKM_alpha05", "ISKM_alpha01", "ISKM_alpha001",
	"SOGC", "SOGC_dp")

colnames(num_sig_path) <- c("SPKM",
	"ISKM_alpha05", "ISKM_alpha01", "ISKM_alpha001",
	"SOGC_k1", "SOGC_k2", "SOGC_k3", "SOGC_combined",
	"SOGC_dp_k1", "SOGC_dp_k2", "SOGC_dp_k3", 
	"SOGC_dp_k4", "SOGC_dp_k5", "SOGC_dp_combined")

round(ari_mat, 2)

#num_sig_path
#write.csv(ari_mat, file=paste0(res_wd, "ARI_summary.csv"))
for(i in 1:3){
	colnames(num_sig_path_cross_validation[[i]]) <- c("SPKM",
		"ISKM_alpha05", "ISKM_alpha01", "ISKM_alpha001",
		"SOGC_k1", "SOGC_k2", "SOGC_k3", "SOGC_combined",
		"SOGC_dp_k1", "SOGC_dp_k2", "SOGC_dp_k3", 
		"SOGC_dp_k4", "SOGC_dp_k5", "SOGC_dp_combined")
	write.csv(num_sig_path_cross_validation[[i]], file=paste0(res_wd, 
		"Num_sig_pathways_summary_FIS_effectSize_crossValidation_", names(num_sig_path_cross_validation)[i],".csv"))
}

#sapply(1:length(input_list),function(x) input_list[[x]]$data_name)

#sapply(1:length(input_list),function(x) input_list[[x]]$pathway_name)
