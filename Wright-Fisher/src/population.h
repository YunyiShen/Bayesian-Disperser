#ifndef POPULATION_H
#define POPULATION_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

// several population level parameters

class population
{
public:
    // vars
    int size;
    arma::vec environment;
    arma::vec genotype_mu, genotype_nu, genotype_alpha, genotype_beta; // genotypes, with prior mean mu, number of effective samples k and precision beta
    arma::vec phenotype_mu, phenotype_nu, phenotype_alpha, phenotype_beta;
    arma::vec phenotype_hat_mu, phenotype_hat_tau; // phenotype of mu and precision tau, point estimations
    arma::vec genotype_thr;                        // how many standard deviations to settle?
    arma::ivec steps_moved;                        // how may steps moved since this year
    arma::ivec settled;                            // whether settled
    // methods
    population(int);
    arma::vec get_fitness(double) const;
    void AR1_mutation(arma::vec);
    void AR1_nextyear(double, double);
    void check_settlement();
    arma::uvec survive_migrator(double) const; // this give which died on the way
};

//random constructer
inline population::population(int n)
{
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
    phenotype_hat_tau = genotype_alpha % (1 / genotype_beta);
    genotype_thr = -log(arma::randu<vec>(n));
    steps_moved = zeros<ivec>(n);
    settled = zeros<ivec>(n);
}

// calculate fitness for WF update, good environment minus steps moved.
inline arma::vec population::get_fitness(double cost) const
{
    return (environment - cost * steps_moved);
}

// see which settled based on their standard
inline void population::check_settlement()
{
    //std::cout << genotype_thr << endl;
    //std::cout << phenotype_hat_tau << endl;
    arma::uvec settled_ind = find(environment > (phenotype_hat_mu + genotype_thr /  sqrt(phenotype_hat_tau)));
    //std::cout << settled_ind << endl;
    settled(settled_ind) *= 0;
    settled(settled_ind) += 1;
    return;
}

// stationary AR1 random walk for a normal r.v., will be used for dispersion, next year and mutation
// x_{t+1} = x_t + a epsilon
inline void normal_AR1(arma::mat &input,
                       const double &trend,
                       double a)
{
    input += (a * randn(size(input)) + trend);
    return;
}

inline void normal_AR1_regress(arma::mat &input,
                               const double &trend,
                               double a)
{
    a = a < 1 ? a : 1;
    input = a * input + ((1 - a) * randn(size(input)) + trend);
    return;
}

inline void lognormal_AR1(arma::mat &input,
                          const double &trend,
                          double a)
{
    input = exp(log(input) + a * randn(size(input)) + trend);
    return;
}

// normal random walk mutation
inline void population::AR1_mutation(arma::vec a)
{
    lognormal_AR1(genotype_thr, 0, a(4));
    lognormal_AR1(genotype_beta, 0, a(3));
    lognormal_AR1(genotype_alpha, 0, a(2));
    lognormal_AR1(genotype_nu, 0, a(1));
    normal_AR1(genotype_mu, 0, a(0));
    return;
}

// next year's environment
inline void population::AR1_nextyear(double a, double trend)
{
    normal_AR1_regress(environment, a, trend);
    return;
}

// which migrator survived
inline arma::uvec population::survive_migrator(double rate) const
{
    rate = rate >= 1 ? 1 : rate;
    arma::vec survival_rate(size);
    for (int i = 0; i < steps_moved.n_elem; ++i)
    {
        survival_rate(i) = log(1 - rate) * steps_moved(i);
    }
    return (find(survival_rate >= log(randu(size))));
}

#endif