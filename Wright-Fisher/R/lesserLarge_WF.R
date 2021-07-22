library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(Matrix)

sourceCpp("./src/Wright-Fisher-simu.cpp")
set.seed(42)
n_iter <- 2000 # number of years
n <- 3000 # population size for WF model
nrep <- 5
risks <- seq(0.005,0.05,0.005)
pred_temp <- seq(0,0.9,0.1)
pred_sp <- seq(0,0.3,0.1)
#nrep <- 1
#risks <- c(0.005,0.05)
#pred_temp <- c(0,0.9)
#pred_sp <- c(0,0.5)
mutation <- c(1,0.1,0.1,0.1,0.1)

n_situations <- length(risks) * length(pred_temp) * length(pred_sp)
i_situ <- 1

res <- list()

for(r in risks){
    for(a_temp in pred_temp){
        for(a_sp in pred_sp){

            res_temp <- list()
            for(i in 1:nrep){
                cat("r =",r, ", a_temp =", a_temp, ", a_sp =", a_sp, ", rep =", i, "\n")
                res_temp[[i]] <- Wright_Fisher_simulation(n=n, max_step = 1500,
#         // max number of steps per year
                                 n_iter = n_iter,
                                 #      // number of years to evalute
                                 mutation = mutation,
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = a_sp,
                                 # // noise level on landscape, env_i+1 = a env_i + (1-a) * epsilon
                                 a_env_temp = a_temp,
                                 # // noise level on time env_t+1 = a env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = r,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
                                  )
                save.image("./res/lesser_big_simu_mu_mutation1.RData",safe = TRUE)
            }
            res[[i_situ]] <- list(risk=r, 
                    pred_temp = a_temp, 
                    a_sp = a_sp, 
                    sample = res_temp)
            i_situ <- i_situ + 1
        }
    }
}



