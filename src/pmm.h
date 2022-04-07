#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <typeinfo>
#include <iterator>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PMM_H
#define PMM_H

namespace pmm {

    std::vector<std::string> decay = {"ciao","sono","giuseppe"};

    double inertia(std::string decay, arma::vec pars){
        return 5.0;
    }
}

#endif