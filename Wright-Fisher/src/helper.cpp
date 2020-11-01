// [[Rcpp::depends(RcppArmadillo)]]
#include "population.h"
#include <RcppArmadillo.h>
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

// we assume a Wright-Fisher type drift and selection, sample gene with multinomial probability
//  with fitness_i = logit(p_i)
void WrightFisher_selection(population &popu, const double &cost)
{
    arma::vec fitness = popu.get_fitness(cost);
    fitness -= mean(fitness); // logit is translation invariant 
    
    arma::vec prob = exp(fitness) / (sum(exp(fitness)));
    arma::ivec res(popu.size);
    
    rmultinom(popu.size, prob.begin(), popu.size, res.begin());
    
    arma::uvec surv(size(res));
    int i_surv = 0;
    for (int i = 0; i < popu.size; ++i)
    {   
        
        if(res(i)>0){
            for(int j = 0 ; j < res(i) ; ++j){
                surv(i_surv) = static_cast<unsigned int> (i); // we need unsigned integer
                i_surv ++ ;
            }
        }
        
    }
    // generate offsprings
    popu.genotype_thr = popu.genotype_thr(surv);
    popu.environment = popu.environment(surv);
    popu.genotype_beta = popu.genotype_beta(surv);
    popu.genotype_alpha = popu.genotype_alpha(surv);
    popu.genotype_mu = popu.genotype_mu(surv);
    popu.genotype_nu = popu.genotype_nu(surv);

    popu.phenotype_alpha = popu.genotype_alpha;
    popu.phenotype_beta = popu.genotype_beta;
    popu.phenotype_mu = popu.genotype_mu;
    popu.phenotype_nu = popu.genotype_nu;

    popu.phenotype_hat_mu = popu.genotype_mu;
    popu.phenotype_hat_tau = popu.genotype_beta % (1 / popu.genotype_alpha);

    // reset movement parameters
    popu.steps_moved *= 0;
    popu.settled *= 0;

    return;
}

// one step move function

void disperse_onestep(population &popu,
                      const double &trend,
                      double a)
{
    a = a<1 ? a : 1;
    arma::uvec unsettled = find(popu.settled == 0); // only move those not yet settled down
    
    popu.steps_moved(unsettled) += 1;               // move one step
    
    
    // to a new environment
    
    popu.environment(unsettled) = a * popu.environment(unsettled) + ((1-a) * randn(size(popu.environment(unsettled))) + trend);
    // get new knowledge, knowledge updating rule was based on normal-gamma conjugate prior of normal, see https://en.wikipedia.org/wiki/Normal-gamma_distribution
    
    
    popu.phenotype_alpha(unsettled) += 0.5;
    popu.phenotype_beta(unsettled) += popu.phenotype_nu(unsettled) %
                                      (1 / (popu.phenotype_nu(unsettled) + 1)) % pow(popu.phenotype_mu(unsettled) - popu.environment(unsettled), 2) / 2;

    popu.phenotype_mu(unsettled) = (popu.phenotype_mu(unsettled) % popu.phenotype_nu(unsettled) + popu.environment(unsettled)) %
                                   (1 / (popu.phenotype_nu(unsettled) + 1));
    popu.phenotype_nu(unsettled) = popu.phenotype_nu(unsettled) + 1;
    
    // make decisions:
    
    popu.phenotype_hat_mu(unsettled) = popu.genotype_mu(unsettled);
    popu.phenotype_hat_tau(unsettled) = popu.genotype_beta(unsettled) % (1 / popu.genotype_alpha(unsettled));
    
    popu.check_settlement();
    
    return;
}

void save_summary_stats(const population &popu,
                        int i, 
                        arma::mat &genotype_mean, 
                        arma::mat &phenotype_mean, 
                        arma::mat &genotype_sd, 
                        arma::mat &phenotype_sd)
{
    genotype_mean(i,0) = mean(popu.genotype_mu);
    genotype_mean(i,1) = mean(popu.genotype_nu);
    genotype_mean(i,2) = mean(popu.genotype_alpha);
    genotype_mean(i,3) = mean(popu.genotype_beta);
    genotype_mean(i,4) = mean(popu.genotype_thr);

    genotype_sd(i,0) = stddev(popu.genotype_mu);
    genotype_sd(i,1) = stddev(popu.genotype_nu);
    genotype_sd(i,2) = stddev(popu.genotype_alpha);
    genotype_sd(i,3) = stddev(popu.genotype_beta);
    genotype_sd(i,4) = stddev(popu.genotype_thr);

    phenotype_mean(i,0) = mean(popu.phenotype_mu);
    phenotype_mean(i,1) = mean(popu.phenotype_nu);
    phenotype_mean(i,2) = mean(popu.phenotype_alpha);
    phenotype_mean(i,3) = mean(popu.phenotype_beta);

    phenotype_sd(i,0) = stddev(popu.phenotype_mu);
    phenotype_sd(i,1) = stddev(popu.phenotype_nu);
    phenotype_sd(i,2) = stddev(popu.phenotype_alpha);
    phenotype_sd(i,3) = stddev(popu.phenotype_beta);
    return;
}
