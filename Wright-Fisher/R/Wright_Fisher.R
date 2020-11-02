library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(Matrix)

sourceCpp("./src/Wright-Fisher-simu.cpp")



test <- Wright_Fisher_simulation(n=1000, max_step = 1500,
#         // max number of steps per year
                                 n_iter = 5000,
                                 #      // number of years to evolute
                                 mutation = matrix(c(1,0.1,0.1,0.1,0.1),5,1),
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = 0,
                                 # // noise level on landscape, env_i+1 = a*env_i + (1-a) epsilon
                                 a_env_temp = 1,
                                 # // noise level on time env_t+1 = a*env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = 0.01,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
)


plot(test$genotype_mean[,5])
points(test$phenotype_mean[,4],col = "blue")