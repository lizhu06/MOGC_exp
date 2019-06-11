ISKM_version <- "new"
code_server <- "/net/wong07/home/liz86/MOGC/"
# detach or remove 
if("ISKmeans" %in% (.packages())){
  detach("package:ISKmeans", unload=TRUE)
}
if("ISKmeans" %in% rownames(installed.packages())){
  remove.packages("ISKmeans")
}
library(devtools)
install_github("Caleb-Huo/IS-Kmeans")

## res in wong04, code(function) in wong07
rm(list=ls())
library(snowfall)

simu <- function(diffmu, weakSignalRatio, ntry){
  options(stringsAsFactors=FALSE)
  if(FALSE){
    diffmu <- 0.8
    weakSignalRatio <- 0.5
    ntry <- 2
  }
  seed <- ntry*1000
  set.seed(seed)  
  code_server <- "/net/wong07/home/liz86/MOGC/Code/"
  res_wd <- "/net/wong04/home/liz86/MOG_clustering/Exp2/Simu3/"

  setwd(res_wd)
  if(!("MASS" %in% rownames(installed.packages())))
  {
    install.packages("MASS", repos='http://cran.us.r-project.org')
    library(MASS)
  }else{
    library(MASS)
  }

  if(!("sparcl" %in% rownames(installed.packages())))
  {
    install.packages("sparcl", repos='http://cran.us.r-project.org')
    library(sparcl)
  }else{
    library(sparcl)
  }

  if(!("mclust" %in% rownames(installed.packages())))
  {
    install.packages("mclust", repos='http://cran.us.r-project.org')
    library(mclust)
  }else{
    library(mclust)
  }

  if(!("pROC" %in% rownames(installed.packages())))
  {
    install.packages("pROC", repos='http://cran.us.r-project.org')
    library(pROC)
  }else{
    library(pROC)
  }

  if(!("label.switching" %in% rownames(installed.packages())))
  {
    install.packages("label.switching", repos='http://cran.us.r-project.org')
    library(label.switching)
  }else{
    library(label.switching)
  }

  if(!("MCMCpack" %in% rownames(installed.packages())))
  {
    install.packages("MCMCpack", repos='http://cran.us.r-project.org')
    library(MCMCpack)
  }else{
    library(MCMCpack)
  }

  if(!("ggplot2" %in% rownames(installed.packages())))
  {
    install.packages("ggplot2", repos='http://cran.us.r-project.org')
    library(ggplot2)
  }else{
    library(ggplot2)
  }

  # install ISKM
  if(!("ISKmeans" %in% rownames(installed.packages()))){
    library(devtools)
    install_github("Caleb-Huo/IS-Kmeans")
  }else{
    library(ISKmeans)
  }

  ##### Simulated data #######
  simu_code_path <- paste0(code_server, "Simu_from_Caleb/MOGC_simu_v3/")
  source(paste0(simu_code_path, "generateMOG.R"))
  source(paste0(code_server, "Function/SOGC_v16_callC.R"))
  source(paste0(code_server, "Function/SOGC_DPT_v11_callC.R"))
  source(paste0(code_server, "Function/MOGC_v8_callC.R"))
  source(paste0(code_server, "Function/MOGC_DPT_v6_callC.R"))

  ## Generate data
  numSamplesPerK <- c(30, 30, 30) 
  m2 <- 4
  m1 <- 40

  U2 <- matrix(0, m1, m2)
  temp <- 0
  for(g in 1:m2){
    U2[temp+seq(1, 10), g] <- 1
    temp <- temp + 10
  }
  U2[1, 2] <- 1
  U2[2, 2] <- 1

  U1 <- matrix(0, m1*3, m1)
  for(i in 1:m1){
    U1[((i-1)*3+1):(i*3), i] <- 1
  }

  numFeatInLevel1Group <- 3
  strongGroupSignal <- diffmu 
  weakGroupSignalRatio <- weakSignalRatio
  strongGroupPerc <- 0.4
  weakGroupPerc <- 0.4  
  noiseMean <- 0 
  sample_sigma <- 1
  sigma_noise <- 1
  percConfounder <- 0
  rho1 <- 0
  rho2 <- 0
  g_index <- c(1, 1, 0, 0)
  k <- 3

  Sdata <- generateMOG(seed=seed, K=k, numSamplesPerK, 
    U2, g_index, numFeatInLevel1Group, 
    strongGroupSignal, weakGroupSignalRatio, 
    strongGroupPerc, weakGroupPerc,   
    noiseMean, sample_sigma, sigma_noise, 
    percConfounder, rho1, rho2)

  ## SPKM ##
  tune <- KMeansSparseCluster.permute(Sdata$data, K=3, nperms=25, 
    wbounds=NULL, silent=TRUE, nvals=10, centers=NULL)
  res_spkm <- KMeansSparseCluster(Sdata$data, K=3, wbounds=tune$bestw, 
    nstart=20, silent =TRUE, maxiter=6, centers=NULL)
  ARI_spkm <- adjustedRandIndex(res_spkm[[1]]$Cs, Sdata$label)  
  AUC_spkm <- auc(roc(factor(Sdata$geneLabel), res_spkm[[1]]$ws))  

  # function to rank features
  sort_feature <- function(p, mu=NULL){
    if(FALSE){
      p <- res$feature_select_prob
      mu <- mu_max
    }
    if(is.null(mu)){
      sorted_gene_index <- seq(1,length(p)) [order(p, 
        decreasing=FALSE)]
    }else{
      sorted_gene_index <- seq(1,length(p)) [order(
        p, mu, decreasing=FALSE)]
    }
    feature_rank <- match(seq(1, length(p)), sorted_gene_index)
    return(feature_rank)
  }
  
  ## SOGC ##
  Y <- Sdata$data 
  K <- 3
  res_sogc <- SOGC(Y, U1, K=K, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=10000, debug=FALSE, init_mu=NULL, 
    init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE,
    seed=123, BernoulliWeighted=TRUE, 
    MH_ind=0, pi_j_prop_n=10, recover_overlap=TRUE, 
    cppfile=paste0(code_server, "Function/SOGC_v16.cpp"))

  ave_mu_sogc <- sapply(1:3, function(k) 
    abs(apply(res_sogc$MU_store[, k, ], 2, mean)))
  mu_max_sogc <- apply(ave_mu_sogc, 1, max)
  sogc_rank <- sort_feature(res_sogc$feature_select_prob, mu_max_sogc)
  AUC_sogc <- auc(roc(factor(Sdata$geneLabel), sogc_rank))
  #AUC_sog <- auc(roc(factor(Sdata$geneLabel), 
  #  res$feature_select_prob))
  ARI_sogc <- adjustedRandIndex(res_sogc$label_map,Sdata$label)

  ## SOGC_DP ##
  Y <- Sdata$data 
  res_sogc_dp <- SOGC_DPT(Y, U1, K=10, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=10000, debug=FALSE, init_mu=NULL, init_z=res_spkm[[1]]$Cs, 
    fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE, seed=123, 
    BernoulliWeighted=TRUE, MH_ind=0,
    pi_j_prop_n=10, recover_overlap=TRUE, 
    cppfile=paste0(code_server, "Function/SOGC_DPT_v11.cpp"))

  #AUC_sog_dp <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))
  #ARI_sog_dp <- adjustedRandIndex(res$label_map, Sdata$label)
  ave_mu_sogc_dp <- sapply(1:10, function(k) 
    abs(apply(res_sogc_dp$MU_store[, k, ], 2, mean)))
  mu_max_sogc_dp <- apply(ave_mu_sogc_dp, 1, max)
  sogc_dp_rank <- sort_feature(res_sogc_dp$feature_select_prob, mu_max_sogc_dp)
  AUC_sogc_dp <- auc(roc(factor(Sdata$geneLabel), sogc_dp_rank))

  #AUC_sogc_dp <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))
  ARI_sogc_dp <- adjustedRandIndex(res_sogc_dp$label_map, Sdata$label)

  ## MOGC ##
  res_mogc <- MOGC(Y, U1, U2, K=3, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=1000, debug=FALSE, init_mu=NULL, 
    init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE, 
    cppfile=paste0(code_server, "Function/MOGC_v8.cpp"), 
    BernoulliWeighted=TRUE, MH_ind=1, pi_g_prop_n=10,
    pi_j_prop_n=10, recover_overlap=TRUE)

  #AUC_mog <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))  
  #ARI_mog <- adjustedRandIndex(res$label_map,Sdata$label)
  ave_mu_mogc <- sapply(1:3, function(k) 
    abs(apply(res_mogc$MU_store[, k, ], 2, mean)))
  mu_max_mogc <- apply(ave_mu_mogc, 1, max)
  mogc_rank <- sort_feature(res_mogc$feature_select_prob, mu_max_mogc)
  AUC_mogc <- auc(roc(factor(Sdata$geneLabel), mogc_rank))
  ARI_mogc <- adjustedRandIndex(res_mogc$label_map,Sdata$label)


  ## MOGC_DP
  K <- 10
  res_mogc_dp <- MOGC_DPT(Y, U1, U2, K=K, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=1000, debug=FALSE, init_mu=NULL, 
    init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE, 
    BernoulliWeighted=TRUE, MH_ind=1, pi_g_prop_n=10,
    pi_j_prop_n=10, 
    cppfile=paste0(code_server, "Function/MOGC_DPT_v6.cpp"), 
    recover_overlap=TRUE)

  #AUC_mog_dp <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))  
  

  ave_mu_mogc_dp <- sapply(1:10, function(k) 
    abs(apply(res_mogc_dp$MU_store[, k, ], 2, mean)))
  mu_max_mogc_dp <- apply(ave_mu_mogc_dp, 1, max)
  mogc_dp_rank <- sort_feature(res_mogc_dp$feature_select_prob, mu_max_mogc_dp)
  AUC_mogc_dp <- auc(roc(factor(Sdata$geneLabel), mogc_dp_rank))
  ARI_mogc_dp <- adjustedRandIndex(res_mogc_dp$label_map, Sdata$label)

  ## save data
  savewd <- paste0(res_wd, "Output_v2/")
  #if (file.exists(savewd)){
  #    setwd(savewd)
  #} else {
  #    dir.create(savewd)
  #    setwd(savewd)
  #}
  save(ARI_sogc, AUC_sogc, ARI_sogc_dp, AUC_sogc_dp, 
    ARI_mogc, AUC_mogc, ARI_mogc_dp, AUC_mogc_dp,
    file=paste0(savewd, "Simu3_diffmu_",
    diffmu, "_weakSignalRatio_", weakSignalRatio, "_ntry_", ntry, 
    "_SPKM_init.RData"))
}

simu_multiple <- function(diffmu_mul, weakSignalRatio_mul, ntry_mul){
  for(i in 1:length(diffmu_mul)){
    foo <- simu(diffmu_mul[i], weakSignalRatio_mul[i], ntry_mul[i])
  }
}


# settings
num_try <- 50
diffmu_uni <- c(0.8, 1)
weak_uni <- c(0.2, 0.5)
diffmu_vec <- rep(diffmu_uni, each=num_try*length(weak_uni))
weak_vec <- rep(rep(weak_uni, each=num_try), times=length(diffmu_uni))
ntry_vec <- rep(seq(1, num_try), times=length(diffmu_uni)*length(weak_uni))
num_cpus <- 20

# create list
diffmu_list <- weak_list <- ntry_list <- list()
num_each <- length(diffmu_vec)/num_cpus
for(i in 1:num_cpus){
  diffmu_list[[i]] <- diffmu_vec[((i-1)*num_each+1) : (i*num_each)]
  weak_list[[i]] <- weak_vec[((i-1)*num_each+1) : (i*num_each)]
  ntry_list[[i]] <- ntry_vec[((i-1)*num_each+1) : (i*num_each)]
}

t0 <- proc.time()

snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=num_cpus) 
snowfall::sfExport("simu", "simu_multiple", 
  "diffmu_list", "weak_list", "ntry_list")
res_return <- snowfall::sfClusterApply(seq(1, num_cpus),function(x) 
  simu_multiple(diffmu_mul=diffmu_list[[x]], 
    weakSignalRatio_mul=weak_list[[x]], ntry_mul=ntry_list[[x]]))
snowfall::sfStop()

delta_t <- proc.time() - t0 
delta_t



