rm(list=ls())
library(mclust)
library(pROC)

# load data
wd <- "/home/liz86/MOGC/Exp/App2/"
load(paste0(wd, "Data/d_nodup.RData"))
load(paste0(wd, "Data/pam50.RData"))
pam50_2 <- as.character(pam50)
pam50_2[which(pam50_2 %in% c("LumA", "LumB", "Normal"))] <- "Luminal"
gene_names <- colnames(d_nodup)
load(paste0(wd, "Data/pathway_data_intersected.RData"))

# load results
res_wd <- paste0(wd, "Results/")
load(paste0(res_wd, "res_spkm.RData"))
load(paste0(res_wd, "res_iskm_alpha05.RData"))
load(paste0(res_wd, "res_iskm_alpha001.RData")) # more group penalty
load(paste0(res_wd, "res_sogc.RData"))
load(paste0(res_wd, "res_sogc_dp.RData"))
load(paste0(res_wd, "res_mogc.RData"))
load(paste0(res_wd, "res_mogc_dp.RData"))
all_res <- list(res_sogc, res_sogc_dp, res_mogc, res_mogc_dp)
names(all_res) <- c("SOGC", "SOGC_dp", "MOGC", "MOGC_dp")

# score results
ARI_table <- matrix(NA, 7, 1)
group_table <- matrix(NA, 23, 5)
rownames(ARI_table) <- 
	c("SPKM", "ISKM_alpha05", "ISKM_alpha001", names(all_res))
rownames(group_table) <- c("SPKM", "ISKM_alpha05", "ISKM_alpha001",
	"SOGC_k1", "SOGC_k2", "SOGC_k3", "SOGC_comb",
	"SOGC_dp_k1", "SOGC_dp_k2", "SOGC_dp_k3", "SOGC_dp_k4", "SOGC_dp_k5", "SOGC_dp_comb",
	"MOGC_k1", "MOGC_k2", "MOGC_k3", "MOGC_comb",
	"MOGC_dp_k1", "MOGC_dp_k2", "MOGC_dp_k3", "MOGC_dp_k4", "MOGC_dp_k5", "MOGC_dp_comb")	
colnames(group_table) <- c("num_feat", "G1", "G2", "G3", "num_sig_path")

# ARI
(ARI_table[1,1] <- adjustedRandIndex(res_spkm[[1]]$Cs, pam50_2))
(ARI_table[2,1] <- adjustedRandIndex(res_iskm_alpha05$Cs, pam50_2))
(ARI_table[3,1] <- adjustedRandIndex(res_iskm_alpha001$Cs, 
	pam50_2)) 
for(i in 1:4){
	ARI_table[i+3,1] <- adjustedRandIndex(all_res[[i]]$label_map, 
		pam50_2)
}

write.csv(ARI_table, file=paste0(res_wd, "ARI_table.csv"))


j <- 0
all_p <- all_mu <- list()
all_p[[1]] <- res_spkm[[1]]$ws
all_p[[2]] <- res_iskm_alpha05$ws
all_p[[3]] <- res_iskm_alpha001$ws
all_mu[[1]] <- NULL
all_mu[[2]] <- NULL
all_mu[[3]] <- NULL
for(i in 1:length(all_res)){
	for(k in 1:dim(all_res[[i]]$MU_store)[2]){
		j <- j + 1
		all_p[[j+3]] <- apply(all_res[[i]]$MU_store[, k, ]!=0, 2, mean)
		all_mu[[j+3]] <- abs(apply(all_res[[i]]$MU_store[, k, ], 2, mean))
	}
	j <- j + 1
	all_p[[j+3]] <- all_res[[i]]$feature_select_prob
	mu_by_k <- sapply(1:dim(all_res[[i]]$MU_store)[2], function(k) 
		abs(apply(all_res[[i]]$MU_store[, k, ], 2, mean)))
	all_mu[[j+3]] <- apply(mu_by_k, 1, max)
}
names(all_p) <- names(all_mu) <- rownames(group_table)

## count gene groups
count_groups <- function(p, mu=NULL){
	if(FALSE){
		p <- all_weights[[1]]
		mu <- all_mu[[1]]
	}
	if(is.null(mu)){
		sorted_gene_index <- seq(1,length(gene_names)) [order(p, decreasing=TRUE)]
	}else{
		sorted_gene_index <- seq(1,length(gene_names)) [order(p, mu, decreasing=TRUE)]
	}
	feature_rank <- match(seq(1, length(gene_names)), sorted_gene_index)
	weight_cut_off <- quantile(seq(1,length(gene_names)), probs=0.33)
	grouping_table <- matrix(NA, length(weight_cut_off), 4)
	colnames(grouping_table) <- c("num_feat", "G1", "G2", "G3")
	for(i in 1:length(weight_cut_off)){
		in_list <- split(feature_rank < weight_cut_off[i], gene_names)
		sum_in_per_gene <- sapply(1:length(in_list), function(x) 
			sum(in_list[[x]]))
		tb <- table(sum_in_per_gene)
		grouping_table[i, 1] <- sum(feature_rank < weight_cut_off[i])
		grouping_table[i, 2] <- tb[match("1", names(tb))]
		grouping_table[i, 3] <- tb[match("2", names(tb))]
		grouping_table[i, 4] <- tb[match("3", names(tb))]
	}	
	return(grouping_table)
}

all_gene_tables <- lapply(1:length(all_p), function(x) 
	count_groups(p=all_p[[x]], mu=all_mu[[x]]))
names(all_gene_tables) <- names(all_p)

for(i in 1:nrow(group_table)){
	group_table[i, 1:4] <- all_gene_tables[[i]][1, ]
}

## pathway enrichment analysis
source("/home/liz86/Functions/Pathway_Fisher_Function.R")
pathway_enrichment <- function(p, mu=NULL, all_genes, 
	num_genes, pathway_data){
	if(FALSE){
		weight <- res_spkm[[1]]$ws
		all_genes <- colnames(d_nodup)
		num_genes <- 1000
	}

	if(is.null(mu)){
		sorted_gene_index <- seq(1,length(gene_names)) [order(
			p, decreasing=TRUE)]
	}else{
		sorted_gene_index <- seq(1,length(gene_names)) [order(
			p, mu, decreasing=TRUE)]
	}
	feature_rank <- match(seq(1, length(gene_names)), sorted_gene_index)

	sig <- gene_names[feature_rank <= num_genes]
	pathwayRes <- pathway_fisher_detail_gene(significant=sig,
	  whole=all_genes, fdr=0.05, database=pathway_data)
	
	return(pathwayRes)
}

all_pathway_res <- lapply(1:length(all_p), function(x) 
	pathway_enrichment(p=all_p[[x]], mu=all_mu[[x]],
		colnames(d_nodup),
		num_genes=round(ncol(d_nodup)/3), pathway_data))

for(i in 1:length(all_pathway_res)){
	group_table[i, 5] <- sum(unlist(all_pathway_res[[i]][,2]) < 0.05)
}
group_table
write.csv(group_table, file=paste0(res_wd, "group_selection_FIS_effectSize.csv"))











