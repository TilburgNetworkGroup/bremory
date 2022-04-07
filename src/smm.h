#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <typeinfo>
#include <iterator>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SMM_H
#define SMM_H

namespace smm {
    double inertia(bool intervals, arma::vec K, double max){
        return 3.0;
    }
}

#endif