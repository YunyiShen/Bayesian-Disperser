// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "helper.h"
#include "population.h"


// [[Rcpp::export]]
Rcpp::List Wright_Fisher_simulation(const int &n,                 // population size
                                    const int &max_step,          // max number of steps per year
                                    const int &n_iter,            // number of years to evolute
                                    const arma::vec &mutation,    // noise level for mutation, for mu, nv, alpha, beta and thr
                                    const double &a_env_sp,       // noise level on landscape
                                    const double &a_env_temp,     // noise level on time
                                    const double &trend_env_temp, // temporal trend of the env
                                    const double &rate,           // death rate for dispersion
                                    bool progress                 // progress bar?
)
{
    // store summary statistics
    arma::mat genotype_mean(n_iter, 5); // mu, nu, alpha, beta, thr
    genotype_mean += NA_REAL;
    arma::mat phenotype_mean(n_iter, 4); // mu nu, alpha, beta
    phenotype_mean += NA_REAL;
    arma::mat genotype_sd(n_iter, 5);
    genotype_sd += NA_REAL;
    arma::mat phenotype_sd(n_iter, 4);
    phenotype_sd += NA_REAL;
    arma::mat dist_move(n_iter, 1);
    dist_move += NA_REAL;

    population popu(n);                // inital a population
    Progress prog((n_iter), progress); // progress bar

    for (int i = 0; i < n_iter; ++i)
    {
        if (Progress::check_abort())
        {
            Rcerr << "keyboard abort\n";
            return (
                Rcpp::List::create(
                    Rcpp::Named("genotype_mean") = genotype_mean,
                    Rcpp::Named("phenotype_mean") = phenotype_mean,
                    Rcpp::Named("genotype_sd") = genotype_sd,
                    Rcpp::Named("phenotype_sd") = phenotype_sd,
                    Rcpp::Named("mean_moved") = dist_move));
        }
        // save summary statistics
        save_summary_stats_genotype(popu, i, genotype_mean, genotype_sd);
        // see if settled
        popu.check_settlement();
        // move
        for (int j = 0; j < max_step; ++j)
        {
            
            
            
            if(as_scalar(sum(popu.settled)) > (popu.size-1)){ 
                
                break;
            }
            
            disperse_onestep(popu, 0, a_env_sp);
            
            
        }
        
        dist_move(i) = mean(popu.steps_moved);
        // store phenotype
        save_summary_stats_phenotype(popu, i, phenotype_mean, phenotype_sd);
        // selection
        WrightFisher_selection(popu, rate);
        
        // mutate
        popu.AR1_mutation(mutation);
        
        // new environment
        popu.AR1_nextyear(a_env_temp, trend_env_temp);
        //Rcout << "env:\n" << popu.environment << endl;
        prog.increment();
    }

    return (
        Rcpp::List::create(
            Rcpp::Named("genotype_mean") = genotype_mean,
            Rcpp::Named("phenotype_mean") = phenotype_mean,
            Rcpp::Named("genotype_sd") = genotype_sd,
            Rcpp::Named("phenotype_sd") = phenotype_sd,
            Rcpp::Named("mean_moved") = dist_move));
}