rm(list=ls())
library(ggplot2)
library(latex2exp)
library(gridExtra)
server <- "/net/wong04/home/liz86/MOG_clustering/"
wd1 <- paste0(server, "Exp2/Simu1/Output/")
wd2 <- paste0(server, "Exp2/Simu1/Output_v2/")
res_wd <- paste0(server, "Exp2/Simu1/Results/")
#setwd(wd)

num_try <- 50
diffmu_uni <- c(0.6, 0.8)
weak_uni <- c(0.2, 0.5)
diffmu_vec <- rep(diffmu_uni, each=num_try*length(weak_uni))
weak_vec <- rep(rep(weak_uni, each=num_try), times=length(diffmu_uni))
ntry_vec <- rep(seq(1, num_try), times=length(diffmu_uni)*length(weak_uni))
num_cpus <- 10
num_methods <- 6

ARI_m <- AUC_m <- array(NA, dim=c(num_methods, length(diffmu_uni), 
  length(weak_uni), num_try))
ARI_long <- AUC_long <- matrix(NA, num_methods*length(diffmu_vec), 4)
colnames(ARI_long) <- c("method", "strongSignal",
  "weakSignalRatio", "ARI")
colnames(AUC_long) <- c("method", "strongSignal",
  "weakSignalRatio", "AUC")  
index <- 0
for(diffmu_ind in 1:length(diffmu_uni)){
  diffmu <- diffmu_uni[diffmu_ind]

  for(weak_ind in 1:length(weak_uni)){
    weakSignalRatio <- weak_uni[weak_ind]

    for(ntry_ind in 1:num_try){
      ntry <- ntry_vec[ntry_ind]
      load(paste0(wd1, "Simu1_diffmu_",diffmu, 
        "_weakSignalRatio_", weakSignalRatio, "_ntry_", ntry, "_SPKM_init.RData"))
    
      AUC_m[1, diffmu_ind, weak_ind, ntry_ind] <- AUC_spkm
      AUC_m[4:6, diffmu_ind, weak_ind, ntry_ind] <- AUC_iskm_vec
     
      ARI_m[1, diffmu_ind, weak_ind, ntry_ind] <- ARI_spkm
      ARI_m[4:6, diffmu_ind, weak_ind, ntry_ind] <- ARI_iskm_vec

      load(paste0(wd2, "Simu1_diffmu_",diffmu, 
        "_weakSignalRatio_", weakSignalRatio, "_ntry_", ntry, "_SPKM_init.RData"))

      AUC_m[2, diffmu_ind, weak_ind, ntry_ind] <- AUC_sog
      AUC_m[3, diffmu_ind, weak_ind, ntry_ind] <- AUC_sog_dp

      ARI_m[2, diffmu_ind, weak_ind, ntry_ind] <- ARI_sog
      ARI_m[3, diffmu_ind, weak_ind, ntry_ind] <- ARI_sog_dp


      AUC_long[index+seq(1,num_methods), 1] <- c("SPKM", "SOGC", 
        "SOGCDP", 
        "ISKMalpha05", "ISKMalpha01", "ISKMalpha001")
      AUC_long[index+seq(1,num_methods), 2] <- rep(paste0("E=", diffmu), num_methods)
      AUC_long[index+seq(1,num_methods), 3] <- rep(paste0("t=", weakSignalRatio), num_methods)
      AUC_long[index+seq(1,num_methods), 4] <- c(AUC_spkm,
        AUC_sog, AUC_sog_dp, AUC_iskm_vec)

      ARI_long[index+seq(1,num_methods), 1] <- c("SPKM", 
        "SOGC", "SOGCDP",
        "ISKMalpha05", "ISKMalpha01", "ISKMalpha001")
      ARI_long[index+seq(1,num_methods), 2] <- rep(paste0("E=", diffmu), num_methods)
      ARI_long[index+seq(1,num_methods), 3] <- rep(paste0("t=", weakSignalRatio), num_methods)
      ARI_long[index+seq(1,num_methods), 4] <- c(ARI_spkm,
        ARI_sog, ARI_sog_dp, ARI_iskm_vec)

      index <- index + num_methods
    }
  }
}

#save(AUC_m, file="Combined_AUC_m_SPKM_init_ISKM_more_alphas.RData")
#save(ARI_m, file="Combined_ARI_m_SPKM_init_ISKM_more_alphas.RData")

###### plot ARI ######
ARI_long <- as.data.frame(ARI_long)
# only select alpha=0.01 for ISKM
ARI_long <- ARI_long[which(!(ARI_long$method %in% c(  
        "ISKMalpha01"))), ]
ARI_long$method <- factor(ARI_long$method, 
  levels=c("SPKM", "ISKMalpha05", "ISKMalpha001",
         "SOGC", "SOGCDP"))
ARI_long$ARI <- as.numeric(as.character(ARI_long$ARI))
ARI_plot <- ggplot(ARI_long, aes(x=method, y=ARI, group=method)) + 
  geom_boxplot(aes(fill=method)) + 
  facet_grid(weakSignalRatio ~ strongSignal) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") + 
  scale_x_discrete(name="", labels=c("SPKM"="SPKM", 
    "ISKMalpha05"=parse(text=TeX("ISKM $\\alpha$=0.5")), 
    "ISKMalpha001"=parse(text=TeX("ISKM $\\alpha$=0.01")), 
    "SOGC"="SOGC",
    "SOGCDP"=parse(text=TeX("SOGC$_{dp}$"))))+
  labs(title="Single-layer non-overlapping")

pdf(file=paste0(res_wd, "Compare_ARI_SPKM_init_trimed_v2.pdf"), 
  width=4.5, height=4)
par(mfrow=c(2,4), mar=c(5,5,2,1), par(cex.axis=0.8))
plot(ARI_plot)
dev.off()

###### plot AUC ######
AUC_long <- as.data.frame(AUC_long)
AUC_long <- AUC_long[which(!(AUC_long$method %in% c(
        "ISKMalpha01"))), ]
# only select alpha=0.01 for ISKM
AUC_long$method <- factor(AUC_long$method, 
  levels=c("SPKM", "ISKMalpha05", "ISKMalpha001",
   "SOGC", "SOGCDP"))
AUC_long$AUC <- as.numeric(as.character(AUC_long$AUC))
AUC_plot <- ggplot(AUC_long, aes(x=method, y=AUC, group=method)) + 
  geom_boxplot(aes(fill=method)) + 
  facet_grid( weakSignalRatio ~ strongSignal) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") + 
  scale_x_discrete(name="", labels=c("SPKM"="SPKM", 
    "ISKMalpha05"=parse(text=TeX("ISKM $\\alpha$=0.5")), 
    "ISKMalpha001"=parse(text=TeX("ISKM $\\alpha$=0.01")), 
    "SOGC"="SOGC",
    "SOGCDP"=parse(text=TeX("SOGC$_{dp}$"))))+
  labs(title="Single-layer non-overlapping")

pdf(file=paste0(res_wd, "Compare_AUC_SPKM_init_trimed_v2.pdf"), 
  width=4.5, height=4)
par(mfrow=c(2,4), mar=c(5,5,2,1), par(cex.axis=0.8))
plot(AUC_plot)
dev.off()

pdf(file=paste0(res_wd, "Combine_ARI_AUC_trimed_v2.pdf"), width=7, height=3.5)
grid.arrange(ARI_plot, AUC_plot, nrow = 1)
dev.off()

simu1_plot <- list(ARI_plot, AUC_plot)
save(simu1_plot, file=paste0(res_wd, "ARI_AUC_plot_list_v2.RData"))

