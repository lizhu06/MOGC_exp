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
    ntry <- 1
  }
  seed <- ntry*1000
  set.seed(seed)  
  code_server <- "/net/wong07/home/liz86/MOGC/Code/"
  res_wd <- "/net/wong04/home/liz86/MOG_clustering/Exp2/Simu1/"
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
  simu_code_path <- paste0(code_server, "Simu_from_Caleb/SOGC_simu_v5/")
  source(paste0(simu_code_path, "generateS.R"))
  source(paste0(code_server,"Function/SOGC_v16_callC.R"))
  source(paste0(code_server,"Function/SOGC_DPT_v11_callC.R"))

  ## Generate data
  G <- 8
  P <- G*30
  U1 <- matrix(0, P, G)
  temp <- 0
  for(g in 1:G){
    U1[temp+seq(1, 30), g] <- 1
    temp <- temp + 30
  }
  #U1[1, 2] <- 1
  #U1[2, 2] <- 1
  rho <- 0
  strongSignal <- diffmu
  weakSignalRatio <- weakSignalRatio
  strongSignalPerc <- 0.5
  weakSignalPerc <- 0.3
  g_index <- c(rep(1,2), rep(0, 6))
  Sdata =generateSOG(seed=seed, K=3, 
    numSamplesPerK=c(40,40,40), 
    U1, g_index,  strongSignal, weakSignalRatio, 
    strongSignalPerc, weakSignalPerc,  
    noiseMean=0, sample_sigma=1, 
    sigma_noise=1, percConfounder=0, rho)

if(TRUE){
  ## SPKM ##
  tune <- KMeansSparseCluster.permute(Sdata$data, K=3, nperms=25, 
    wbounds=NULL, silent=TRUE, nvals=10, centers=NULL)
  res_spkm <- KMeansSparseCluster(Sdata$data, K=3, wbounds=tune$bestw, 
    nstart=20, silent =TRUE, maxiter=6, centers=NULL)
  ARI_spkm <- adjustedRandIndex(res_spkm[[1]]$Cs, Sdata$label)  
  AUC_spkm <- auc(roc(factor(Sdata$geneLabel), res_spkm[[1]]$ws))  

  ## ISKM
  groupset <- lapply(1:ncol(U1), function(x) which(U1[,x]==1))
  alpha_vec <- c(0.5, 0.1, 0.01)
  #alpha_vec <- 0.5
  ARI_iskm_vec <- AUC_iskm_vec <- rep(NA, length(alpha_vec))
  names(ARI_iskm_vec) <- names(AUC_iskm_vec) <- paste0("alpha=", alpha_vec)
  for(av in 1:length(alpha_vec)){
    alpha <- alpha_vec[av]
    #gammas <- seq(0.1,0,0.1)
    sparseStart <- TRUE
    gapStat = gapStatistics(Sdata$data, K=3, B=25, gamma=NULL, 
      alpha=alpha, group=groupset, silence=TRUE)
    res <- ISKmeans(Sdata$data, K=3, gamma=gapStat$bestGamma, alpha=alpha,
      group=groupset, nstart=20, silent=TRUE, maxiter=20,
      sparseStart=sparseStart)
    ARI_iskm_vec[av] <- adjustedRandIndex(res$Cs, Sdata$label)  
    AUC_iskm_vec[av] <- auc(roc(factor(Sdata$geneLabel), res$ws)) 
  }

  ## SOGC ##
  Y <- Sdata$data 
  K <- 3
  res <- SOGC(Y, U1, K=K, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=10000, debug=FALSE, init_mu=NULL, 
    init_z=res_spkm[[1]]$Cs, fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE,
    seed=123, BernoulliWeighted=TRUE, 
    MH_ind=0, pi_j_prop_n=10, recover_overlap=TRUE, 
    cppfile=paste0(code_server, "Function/SOGC_v16.cpp"))

  AUC_sog <- auc(roc(factor(Sdata$geneLabel), 
    res$feature_select_prob))
  ARI_sog <- adjustedRandIndex(res$label_map,Sdata$label)

  ## SOGC_DP ##
  Y <- Sdata$data 
  res <- SOGC_DPT(Y, U1, K=10, center=FALSE, 
    burnInIter=1000, keepIter=3000, maxIter=10000, 
    print_int=10000, debug=FALSE, init_mu=NULL, init_z=res_spkm[[1]]$Cs, 
    fix_z=FALSE, 
    fix_mu=FALSE, adj_ls=TRUE, seed=123, 
    BernoulliWeighted=TRUE, MH_ind=0,
    pi_j_prop_n=10, recover_overlap=TRUE, 
    cppfile=paste0(code_server, "Function/SOGC_DPT_v11.cpp"))

  AUC_sog_dp <- auc(roc(factor(Sdata$geneLabel), res$feature_select_prob))
  ARI_sog_dp <- adjustedRandIndex(res$label_map, Sdata$label)
  }

  ## save data
  savewd <- paste0(res_wd, "Output/")
  #if (file.exists(savewd)){
  #    setwd(savewd)
  #} else {
  #    dir.create(savewd)
  #    setwd(savewd)
  #}
  save(ARI_sog, AUC_sog, ARI_sog_dp, AUC_sog_dp, 
    ARI_spkm, AUC_spkm, ARI_iskm_vec, AUC_iskm_vec,
  #save(G,
    file=paste0(savewd, "Simu1_diffmu_",
    diffmu, "_weakSignalRatio_", weakSignalRatio, "_ntry_", 
    ntry, "_SPKM_init.RData"))
}

simu_multiple <- function(diffmu_mul, weakSignalRatio_mul, ntry_mul){
  for(i in 1:length(diffmu_mul)){
    foo <- simu(diffmu_mul[i], weakSignalRatio_mul[i], ntry_mul[i])
  }
}

# setting
num_try <- 50
diffmu_uni <- c(0.6, 0.8)
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
res_return <- snowfall::sfClusterApply(seq(1,num_cpus),
  function(x) 
  simu_multiple(diffmu_mul=diffmu_list[[x]], 
    weakSignalRatio_mul=weak_list[[x]], 
    ntry_mul=ntry_list[[x]]))
snowfall::sfStop()

delta_t <- proc.time() - t0  
delta_t


