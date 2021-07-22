library(ggplot2)
simudata <- read.csv("Wright-Fisher/res/big_simu_processed_mu_mutation1.csv")
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

aggregate_simu <- aggregate_simu %>% 
  mutate(autocorrelation = paste("autocorrelation", sp_auto),
         predictability = factor(pred_temp))

#rm(simudata)
#rm(aggregate_simu)


require(plyr)
require(magrittr)
require(ggthemes)


gglist <- list(geom_tile(), scale_fill_viridis_c(), 
               facet_grid(.~patchiness),  ylab("predictability"))

gglist <- list(geom_line(), geom_point(),
               facet_grid(.~autocorrelation),
               theme_few())

# Genetically encoded knowledge

ggplot(aggregate_simu %>% 
         mutate(genotype = mean_genotypemu0),
       aes(risk, genotype, col=predictability)) + gglist + 
       #aes(risk, pred_temp, fill=genotype)) + gglist + 
  ggtitle("Mean genotype of environmental estimate") 


# Genetically encoded knowledge

ggplot(aggregate_simu %>% 
         mutate(knowledge = log(mu_ratio)), 
       aes(risk, knowledge, col=predictability)) + gglist + 
       #aes(risk, pred_temp, fill=knowledge)) + gglist + 
  ggtitle("Knowledge use")


# Genetically encoded knowledge

ggplot(aggregate_simu %>% 
         mutate(movement = mean_moved), 
       aes(risk, movement, col=predictability)) + 
       # aes(risk, pred_temp, fill=movement)) + 
  gglist + 
  ggtitle("Total movement")

