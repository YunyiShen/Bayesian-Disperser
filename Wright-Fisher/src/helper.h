#ifndef HELPER_H
#define HELPER_H

// [[Rcpp::depends(RcppArmadillo)]]
#include "population.h"
#include <RcppArmadillo.h> 
#include <tgmath.h>
using namespace Rcpp;
using namespace arma;

void normal_AR1(arma::mat & input,
                const double & trend,
                double a);

void WrightFisher_selection(population & popu,const double & cost);

void disperse_onestep(population & popu,
                const double & trend,
                double a);

void save_summary_stats(const population &popu,
                        int i, 
                        arma::mat &genotype_mean, 
                        arma::mat &phenotype_mean, 
                        arma::mat &genotype_sd, 
                        arma::mat &phenotype_sd);

#endif