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


// we assume a Wright-Fisher type drift and selection, sample gene with multinomial probability 
//  with fitness_i = logit(p_i)
void WrightFisher_selection(population & pop,const double & cost){
    arma::vec fitness = pop.get_fitness(cost);
    arma::vec prob = exp(fitness)/(sum(exp(fitness)));
    arma::ivec res(pop.size);
    rmultinom(pop.size,prob.begin(),pop.size,res.begin());
    arma::uvec surv(size(res));
    for(int i = 0; i<pop.size; ++i){
        surv(i) = static_cast<unsigned int>(res(i)); // we need unsigned integer
    }

    pop.alpha = pop.alpha(surv);
    pop.environment = pop.environment(surv);
    pop.genotype_beta = pop.genotype_beta(surv);
    pop.genotype_alpha = pop.genotype_alpha(surv);
    pop.genotype_mu = pop.genotype_mu(surv);
    pop.genotype_nu = pop.genotype_nu(surv);

    pop.phenotype_alpha = pop.genotype_alpha;
    pop.phenotype_beta = pop.genotype_beta;
    pop.phenotype_mu = pop.genotype_mu;
    pop.phenotype_nu = pop.genotype_nu;


    pop.phenotype_hat_mu = pop.genotype_mu;
    pop.phenotype_hat_tau = pop.genotype_beta % (1/pop.genotype_alpha);
    pop.steps_moved *= 0;
    pop.settled *= 0;


    return;
} 

// one step move function

void disperse_onestep(population & pop,
                const double & trend,
                double a){
    arma::uvec unsettled = find(pop.settled == 0); // only move those not yet settled down
    pop.steps_moved(unsettled) += 1; // move one step
    // to a new environment
    pop.environment(unsettled) += (a * randn(size(pop.environment(unsettled))) + trend);
    // get new knowledge:
    
    pop.phenotype_alpha(unsettled) += 0.5;
    pop.phenotype_beta(unsettled) += pop.phenotype_nu(unsettled) % 
        (1/(pop.phenotype_nu(unsettled) + 1)) % pow(pop.phenotype_mu(unsettled) - pop.environment(unsettled),2)/2;

    pop.phenotype_mu(unsettled) = (pop.phenotype_mu(unsettled) % pop.phenotype_nu(unsettled) + pop.environment(unsettled)) % 
        (1/(pop.phenotype_nu(unsettled) + 1));
    pop.phenotype_nu = pop.phenotype_nu(unsettled) + 1; 

    // make decisions:
    pop.phenotype_hat_mu(unsettled) = pop.genotype_mu(unsettled);
    pop.phenotype_hat_tau(unsettled) = pop.genotype_beta(unsettled) % (1/pop.genotype_alpha(unsettled)); 

    
    arma::uvec settled = find(pop.environment > pop.phenotype_mu + pop.alpha / (pop.phenotype_hat_tau));
    pop.settled(settled) *= 0;
    pop.settled(settled) += 1;

    return;
}
