// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

// stationary AR1 random walk for a normal r.v., will be used for dispersion, next year and mutation
// x_{t+1} = x_t + a epsilon
void normal_AR1(arma::mat & input,
                const double & trend,
                double a){
    input += (a * randn(size(input)) + trend);
    return;
}

void lognormal_AR1(arma::mat & input,
                const double & trend,
                double a){
    input = exp( log(input) + a * randn(size(input)) + trend);
    return;
}

// several population level parameters
class population {
    public:
        int size;
        arma::vec environment;
        arma::vec genotype_mu, genotype_nu, genotype_alpha, genotype_beta;// genotypes, with prior mean mu, number of effective samples k and precision beta
        arma::vec phenotype_mu, phenotype_nu, phenotype_alpha, phenotype_beta;
        arma::vec phenotype_hat_mu, phenotype_hat_tau; // phenotype of mu and precision tau, point estimations
        arma::vec alpha; // how many standard deviations to settle?
        arma::ivec steps_moved; // how may steps moved since this year
        arma::ivec settled;
    public:
        population(int);
        arma::vec get_fitness (double) const;
        void AR1_mutation(double);
        void AR1_nextyear(double,double); 
        void check_settlement(){
            arma::uvec settled = find(environment > phenotype_mu + alpha / (phenotype_hat_tau));
            settled(settled) *= 0;
            settled(settled) += 1;
            return;
        }
};

// random construct
population::population(int n){
    size = n;
    environment = arma::randn<vec>(n);
    genotype_mu = arma::randn<vec>(n);
    genotype_nu = exp(arma::randn<vec>(n));
    genotype_alpha = exp(arma::randn<vec>(n));
    genotype_beta = exp(arma::randn<vec>(n));

    phenotype_alpha = genotype_alpha; 
    phenotype_beta = genotype_beta; 
    phenotype_mu = genotype_mu; 
    phenotype_nu = genotype_nu; 

    phenotype_hat_mu = genotype_mu;
    phenotype_hat_tau = genotype_beta % (1/genotype_alpha);
    alpha = log(arma::randu<vec>(n));
    steps_moved = zeros<ivec>(n);
    settled = zeros<ivec>(n);
}

// calculate fitness for WF update, good environment minus steps moved. 
arma::vec population::get_fitness (double cost) const {
    return(environment-cost*steps_moved);
}

// normal random walk mutation
void population::AR1_mutation(double a){
    lognormal_AR1(alpha,0,a);
    lognormal_AR1(genotype_beta,0,a);
    lognormal_AR1(genotype_alpha,0,a);
    lognormal_AR1(genotype_nu,0,a);
    normal_AR1(genotype_mu,0,a);
    return;
}


// next year's environment
void population::AR1_nextyear(double a,double trend){
    normal_AR1(environment, a, trend);
    return;
}

// we assume a Wright-Fisher type drift and selection, sample gene with multinomial probability 
//  with fitness_i = logit(p_i)
void WrightFisher_selection(population & popu,const double & cost){
    arma::vec fitness = popu.get_fitness(cost);
    arma::vec prob = exp(fitness)/(sum(exp(fitness)));
    arma::ivec res(popu.size);
    rmultinom(popu.size,prob.begin(),popu.size,res.begin());
    arma::uvec surv(size(res));
    for(int i = 0; i<popu.size; ++i){
        surv(i) = static_cast<unsigned int>(res(i)); // we need unsigned integer
    }

    popu.alpha = popu.alpha(surv);
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
    popu.phenotype_hat_tau = popu.genotype_beta % (1/popu.genotype_alpha);
    popu.steps_moved *= 0;
    popu.settled *= 0;


    return;
} 

// one step move function

void disperse_onestep(population & popu,
                const double & trend,
                double a){
    arma::uvec unsettled = find(popu.settled == 0); // only move those not yet settled down
    popu.steps_moved(unsettled) += 1; // move one step
    // to a new environment
    popu.environment(unsettled) += (a * randn(size(popu.environment(unsettled))) + trend);
    // get new knowledge:
    
    popu.phenotype_alpha(unsettled) += 0.5;
    popu.phenotype_beta(unsettled) += popu.phenotype_nu(unsettled) % 
        (1/(popu.phenotype_nu(unsettled) + 1)) % pow(popu.phenotype_mu(unsettled) - popu.environment(unsettled),2)/2;

    popu.phenotype_mu(unsettled) = (popu.phenotype_mu(unsettled) % popu.phenotype_nu(unsettled) + popu.environment(unsettled)) % 
        (1/(popu.phenotype_nu(unsettled) + 1));
    popu.phenotype_nu = popu.phenotype_nu(unsettled) + 1; 

    // make decisions:
    popu.phenotype_hat_mu(unsettled) = popu.genotype_mu(unsettled);
    popu.phenotype_hat_tau(unsettled) = popu.genotype_beta(unsettled) % (1/popu.genotype_alpha(unsettled)); 

    
    popu.check_settlement();

    return;
}
