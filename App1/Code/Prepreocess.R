rm(list=ls())
library(GSA)

server <- "/home/liz86/"
data_wd <- paste0(server, "Dataset/Leukemia_ISKM_copied_from_Caleb/")

load(paste0(data_wd, "s_gene.rdata"))
load(paste0(data_wd, "leukemiaLabel3.Rdata"))

sapply(s_gene, dim)
#     S6891 S17855 S13159
#[1,] 20154  20155  20155
#[2,]    89     74    105
all(rownames(s_gene[[2]]) == rownames(s_gene[[3]]))

sapply(label, table)
#         [,1] [,2] [,3]
#inv(16)    33   27   28
#t(15;17)   21   19   37
#t(8;21)    35   28   40

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

summary(sapply(1:length(all_pathways[[1]]), function(x) length(all_pathways[[1]][[x]])))
summary(sapply(1:length(all_pathways[[2]]), function(x) length(all_pathways[[2]][[x]])))
summary(sapply(1:length(all_pathways[[3]]), function(x) length(all_pathways[[3]][[x]])))

## filtering genes and intersect with pathways
input_list <- list()
num_path <- rep(NA, 9)
for(s in 1:3){
	for(p in 1:3){
		# load gene expression
		adata = s_gene[[s]]
		alabel <- label[[s]]
		adata_name <- names(s_gene)[s]
		cat(paste0("data name: ", adata_name, "\n"))

		# load pathways
		aPathway_database <- all_pathways[[p]]
		aPathway_name <- names(all_pathways)[p]
		cat(paste0("pathway name: ", aPathway_name, "\n"))
		pathway_names <- names(aPathway_database)

		## filter out 30% low expression genes
		geneRowMeans <- rowMeans(adata)
		if(aPathway_name == "biocarta"){
			aquantile <- 0.3
		}else{
			aquantile <- 0.3
		}
		adata <- adata[geneRowMeans >= quantile(geneRowMeans,aquantile),]
		uni_genes <- unique(rownames(adata))
		cat(paste0("data dimension after filtering: ", dim(adata)[1], "\n"))

		## intersect with pathways and filtering based on size
		if(aPathway_name == "biocarta"){
			size_limit <- c(15, 100)
		}else{
			size_limit <- c(15, 100)
		}
		aPathway_database <- lapply(1:length(aPathway_database), 
			function(x) intersect(aPathway_database[[x]], uni_genes))
		path_length <- sapply(aPathway_database, length)
		aPathway_database <- aPathway_database[which(
			path_length >=size_limit[1] & path_length <= size_limit[2])]
		names(aPathway_database) <- pathway_names[which(
			path_length >=size_limit[1] & path_length <= size_limit[2])]
		cat(paste0("pathway length: ", length(aPathway_database), "\n"))

		# intersect gene expression with pathways
		adata <- adata[which(
			rownames(adata) %in% 
			unlist(aPathway_database)), ]
		cat(paste0("data dimension after intersection: ", dim(adata)[1], "\n"))

		input_list[[(s-1)*3+p]] <- list(data_name=adata_name, data=adata, 
			label=alabel, pathway_name=aPathway_name, 
			pathway_data=aPathway_database)


	}
	cat(paste0("===========================", "\n"))
}
save(input_list, file=paste0(server, "MOGC/Exp3/App1_v3/Data/input_list.RData"))
#sapply(1:length(input_list), function(x) dim(input_list[[x]]$data)[1])
#sapply(1:length(input_list), function(x) dim(input_list[[x]]$data)[2])

#sapply(1:length(input_list), function(x) input_list[[x]]$data_name)
#sapply(1:length(input_list), function(x) input_list[[x]]$pathway_name)


