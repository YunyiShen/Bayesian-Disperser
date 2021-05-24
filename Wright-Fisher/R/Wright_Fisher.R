library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(Matrix)

sourceCpp("./src/Wright-Fisher-simu.cpp")

n_iter <- 2000 # number of years
n <- 3000 # population size for WF model

predictable_lowrisk <- Wright_Fisher_simulation(n=n, max_step = 1500,
#         // max number of steps per year
                                 n_iter = n_iter,
                                 #      // number of years to evalute
                                 mutation = c(.1,0.1,0.1,0.1,0.1),
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = 0,
                                 # // noise level on landscape, env_i+1 = a env_i + (1-a) * epsilon
                                 a_env_temp = 0.9,
                                 # // noise level on time env_t+1 = a env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = 0.005,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
)

unpredictable_lowrisk <- Wright_Fisher_simulation(n=n, max_step = 1500,
#         // max number of steps per year
                                 n_iter = n_iter,
                                 #      // number of years to evalute
                                 mutation = c(.1,0.1,0.1,0.1,0.1),
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = 0,
                                 # // noise level on landscape, env_i+1 = a env_i + (1-a) * epsilon
                                 a_env_temp = 0.1,
                                 # // noise level on time env_t+1 = a env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = 0.005,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
)


predictable_highrisk <- Wright_Fisher_simulation(n=n, max_step = 1500,
#         // max number of steps per year
                                 n_iter = n_iter,
                                 #      // number of years to evalute
                                 mutation = c(.1,0.1,0.1,0.1,0.1),
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = 0,
                                 # // noise level on landscape, env_i+1 = a env_i + (1-a) * epsilon
                                 a_env_temp = 0.9,
                                 # // noise level on time env_t+1 = a env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = 0.05,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
)

unpredictable_highrisk <- Wright_Fisher_simulation(n=n, max_step = 1500,
#         // max number of steps per year
                                 n_iter = n_iter,
                                 #      // number of years to evalute
                                 mutation = c(.1,0.1,0.1,0.1,0.1),
                                 #// noise level for mutation, for mu, nv, alpha, beta and thr
                                 a_env_sp = 0,
                                 # // noise level on landscape, env_i+1 = a env_i + (1-a) * epsilon
                                 a_env_temp = 0.1,
                                 # // noise level on time env_t+1 = a env_t + (1-a) epsilon
                                 trend_env_temp = 0,
                                 #// temporal trend of the env
                                 rate = 0.05,
                                 #// cost for dispersion, as death rate each step
                                 progress = T
                                 #// progress bar?
)


#jpeg('prior_mu_time.jpg', width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
plot(predictable_lowrisk$genotype_mean[,1],
    xlab = "time", ylab = "Mean Prior mu",
    ylim = c(-3,3),
    main = "predictable_lowrisk")
plot(unpredictable_lowrisk$genotype_mean[,1],
    xlab = "time", ylab = "Mean Prior mu",
    ylim = c(-3,3),main = "unpredictable_lowrisk")
plot(predictable_highrisk$genotype_mean[,1],
    xlab = "time", ylab = "Mean Prior mu",
    ylim = c(-3,3),main = "predictable high risk")
plot(unpredictable_highrisk$genotype_mean[,1],
    xlab = "time", ylab = "Mean Prior mu",
    ylim = c(-3,3),main = "unpredictable high risk")
#dev.off()


prior_mu <- data.frame(Pred_Low = predictable_lowrisk$genotype_mean[floor(n_iter/2):n_iter,1],
Pred_High = predictable_highrisk$genotype_mean[floor(n_iter/2):n_iter,1],
Unpred_Low = unpredictable_lowrisk$genotype_mean[floor(n_iter/2):n_iter,1],
Unpred_High = unpredictable_highrisk$genotype_mean[floor(n_iter/2):n_iter,1])

jpeg('./res/prior_mu_stat.jpg', width = 6, height = 6, units = "in", res = 300)
boxplot(prior_mu, xlab = "Scheme",ylab = "(Mean) Genotype Prior on mu")
dev.off()


prior_nu <- data.frame(Pred_Low = predictable_lowrisk$genotype_mean[floor(n_iter/2):n_iter,2],
Pred_High = predictable_highrisk$genotype_mean[floor(n_iter/2):n_iter,2],
Unpred_Low = unpredictable_lowrisk$genotype_mean[floor(n_iter/2):n_iter,2],
Unpred_High = unpredictable_highrisk$genotype_mean[floor(n_iter/2):n_iter,2])

jpeg('./res/prior_nu_stat.jpg', width = 6, height = 6, units = "in", res = 300)
boxplot(prior_nu, xlab = "Scheme",ylab = "(Mean) Genotype Prior certainty on mu",outline = T)
dev.off()




par(mfrow = c(2,2))
plot(predictable_lowrisk$mean_move,
    xlab = "time", ylab = "Distance Moved",
    ylim = c(0,300),main = "predictable_lowrisk")
plot(unpredictable_lowrisk$mean_move,
    xlab = "time", ylab = "Distance Moved",
    ylim = c(0,300),main = "unpredictable_lowrisk")
plot(predictable_highrisk$mean_move,
    xlab = "time", ylab = "Distance Moved",
    ylim = c(0,300),main = "predictable high risk")
plot(unpredictable_highrisk$mean_move,
    xlab = "time", ylab = "Distance Moved",
    ylim = c(0,300),main = "unpredictable high risk")


mean_move <- data.frame(Pred_Low = predictable_lowrisk$mean_move[floor(n_iter/2):n_iter],
Pred_High = predictable_highrisk$mean_move[floor(n_iter/2):n_iter],
Upred_Low = unpredictable_lowrisk$mean_move[floor(n_iter/2):n_iter],
Upred_High = unpredictable_highrisk$mean_move[floor(n_iter/2):n_iter])

jpeg('./res/movedist_stat.jpg', width = 6, height = 6, units = "in", res = 300)
boxplot(mean_move, xlab = "Scheme",ylab = "Mean dispersion steps")
dev.off()

