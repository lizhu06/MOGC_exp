rm(list=ls())
wd <- "/net/wong04/home/liz86/MOG_clustering/Exp3/App3_metabric_v1/"
res_wd <- paste0(wd, "Output/")
load(paste0(res_wd, "res_sogc_dp.RData"))
load(paste0(res_wd, "res_mogc_dp.RData"))

all_res <- list(res_sogc_dp, res_mogc_dp)
names(all_res) <- c("SOGC_dp", "MOGC_dp")

K_vec <- matrix(NA, 2, 3000)
prop_table <- matrix(NA, 2, 5)
for(s in 1:2){
	K_vec[s,] <- sapply(1:nrow(all_res[[s]]$Z_store), function(i) 
		length(unique(all_res[[s]]$Z_store[i,])))
	prop_table[s, ] <- prop.table(table(all_res[[s]]$label_map))
}

pdf(file=paste0(res_wd, "prop_per_cluster_DP.pdf"), width=8, height=6)
par(mfrow=c(2, 2))
for(s in 1:2){
	plot(seq(1, 3000), K_vec[s, ], xlab="iteration", ylab="K", 
		main=names(all_res)[s])
}
for(s in 1:2){
	plot(seq(1,5), prop_table[s,]*100, type="b", col=s, 
		xlab="cluster index", ylab="% samples",
		main=names(all_res)[s])
}
dev.off()