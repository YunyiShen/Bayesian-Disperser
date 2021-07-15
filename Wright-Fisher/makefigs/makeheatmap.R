simudata <- read.csv("./res/big_simu_processed.csv")
simudata$mu_ratio <- simudata$mean_genotypenu/simudata$mean_moved

logit <- function(x){
  log(x/(1-x))
}

aggregate_simu <- aggregate(.~risk+pred_temp+sp_auto, data = simudata, FUN = "mean")

#sp0 <- aggregate_simu[aggregate_simu$sp_auto==0,]
#sp1 <- aggregate_simu[aggregate_simu$sp_auto==0.1,]
#sp2 <- aggregate_simu[aggregate_simu$sp_auto==0.2,]
#sp3 <- aggregate_simu[aggregate_simu$sp_auto==0.3,]

aggregate_simu$sp_auto <- as.factor(aggregate_simu$sp_auto)
gene_names <- c("mu0","nu","alpha","beta","thr")



#rm(simudata)
#rm(aggregate_simu)

ggplot(aggregate_simu, aes(risk, pred_temp, fill=log(mu_ratio))) +
  geom_tile() +
  scale_fill_viridis_c() + 
  facet_grid(.~sp_auto, scales = "free")


