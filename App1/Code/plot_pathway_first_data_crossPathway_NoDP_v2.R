rm(list=ls())
library(sparcl)
library(ISKmeans)
library(foreach)
library(doParallel)
library(mclust)
library(pROC)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(latex2exp)
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

p_list <- list()
counter <- 0
for(s in 1:3){
	data_name <- input_list[[s]]$data_name
	data <- t(input_list[[s]]$data) # sample by feature
	data_s <- scale(data)
	label <- input_list[[s]]$label
	pathway_name <- input_list[[s]]$pathway_name
	pathway_data <- input_list[[s]]$pathway_data

	for(i in 1:3){
		counter <- counter + 1
		path_spkm <- read.delim(paste0(res_wd, "PathwayRes_spkm_", 
			data_name, "_", pathway_name,  "_testPath_",
				names(all_pathways)[i], "_crossValidation.csv"))
		pval_spkm <- as.numeric(sapply(1:nrow(path_spkm), function(x) 
			strsplit(as.character(path_spkm[x,1]), split=",") [[1]][3]))

		alpha_vec <- c(0.5, 0.1, 0.01)
		pval_iskm_mat <- list()
		for(alpha_index in 1:length(alpha_vec)){
			path_iskm <- read.delim(paste0(res_wd, "PathwayRes_iskm_alpha_",
				alpha_vec[alpha_index], "_",  data_name, "_", pathway_name,
				 "_testPath_",
					names(all_pathways)[i],  "_crossValidation.csv"))
			pval_iskm_mat[[alpha_index]] <- as.numeric(sapply(1:nrow(path_iskm), 
				function(x) 
				strsplit(as.character(path_iskm[x,1]), split=",") [[1]][3]))
		}

	# SOG
	path_sog <- read.delim(paste0(res_wd, "PathwayRes_sogc_", 
				data_name, "_", pathway_name, "_testPath_", 
				names(all_pathways)[i],  "_crossValidation.csv"))
	pval_sog <- as.numeric(sapply(1:nrow(path_sog), function(x) 
		strsplit(as.character(path_sog[x,1]), split=",") [[1]][3]))
	print(sum(pval_sog<0.05))

	# SOGDP
	path_sog_dp <- read.delim(paste0(res_wd, "PathwayRes_sogc_dp_", 
				data_name, "_", pathway_name,"_testPath_", 
				names(all_pathways)[i],  "_crossValidation.csv"))
	pval_sog_dp <- as.numeric(sapply(1:nrow(path_sog_dp), function(x) 
		strsplit(as.character(path_sog_dp[x,1]), split=",") [[1]][3]))

	nlog10_pval <- (-1) * log10(c(pval_spkm, pval_iskm_mat[[1]], 
		pval_iskm_mat[[3]],pval_sog))
	method <- c(rep("SPKM", length(pval_spkm)), 
		rep("ISKM05", length(pval_spkm)),
		rep("ISKM001", length(pval_spkm)),
		rep("SOGC", length(pval_sog)))

	data_plot <- data.frame(nlog10pval=nlog10_pval, method=method)
	data_plot$method <- factor(data_plot$method, levels=c("SPKM", 
		"ISKM05", "ISKM001", "SOGC"))
	p_list[[counter]] <- ggplot(data=data_plot, aes(method, nlog10pval)) + 
		geom_jitter(aes(colour=method))+
		labs(y="- log10 p-value", x="", 
			title=paste0(pathway_name, "->", names(all_pathways)[i])) + 
		guides(fill=FALSE, color=FALSE) + 
	  #geom_hline(yintercept = -log(cutoff,base=10))+
	  theme(axis.title.x=element_text(size=12))+
	  theme(axis.title.y=element_text(size=12))+
	  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))+
	  theme(axis.text.y = element_text(size=10))+
	  theme(plot.title = element_text(size=14)) +
	  theme(plot.margin = unit(c(0.1,0.2,0,0.2), "cm")) + 
	  geom_hline(yintercept=(-1)*log10(0.05), linetype="dashed") +
	  scale_x_discrete(name="", labels=c("SPKM"="SPKM", 
	    "ISKM05"="ISKM-I", 
	    "ISKM001"="ISKM-II", 
	    "SOGC"="SOGC")) +
	  ylim(-0.2, 6.2)
	#p_list[s] <- p
	}
}

pdf(file=paste0(res_wd, "Plot_pathway_first_data_crossValidation_NoDP.pdf"), 
	width=8, height=10)
ggarrange(plotlist=p_list, nrow=3, ncol=3)
dev.off()







