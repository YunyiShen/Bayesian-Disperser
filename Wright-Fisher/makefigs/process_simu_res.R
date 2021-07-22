# loading simulation results:
load("./res/lesser_big_simu_mu_mutation1.RData")

get_stationary_mean <- function(w){
  if(ncol(w)>1){
    colMeans(w[1001:2000,])
  }
  else{
    mean(w[1001:2000])
  }
}

get_one_repeat <- function(w){
  temp <- lapply(w, get_stationary_mean)
  Reduce(c, temp)
}

get_data_frame <- function(w){
  gene_names <- c("mu0","nu","alpha","beta","thr")
  ress <- data.frame(matrix(NA,5,22))
  colnames(ress) <- c("risk","pred_temp","sp_auto",
                      paste0("mean_genotype", gene_names),
                      paste0("mean_phenotype", gene_names[-5]),
                      paste0("sd_genotype", gene_names),
                      paste0("sd_phenotype", gene_names[-5]),
                      "mean_moved"
                      )
  ress$risk <- w$risk
  ress$pred_temp <- w$pred_temp
  ress$sp_auto <- w$a_sp
  
  temp <- lapply(w$sample, get_one_repeat)
  ress[,4:22] <- Reduce(rbind, temp)
  return(ress)
}


simulation_data <- lapply(res, get_data_frame)
simulation_data <- Reduce(rbind, simulation_data)

write.csv(simulation_data,"./res/big_simu_processed_mu_mutation1.csv", row.names = F)
