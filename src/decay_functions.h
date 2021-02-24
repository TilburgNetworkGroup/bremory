#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <iterator>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

//' stepwise 
//'
//' @param x transpired time
//' @param pars is the matrix of (betas,widths)
//'
//' @export
// [[Rcpp::export]]
double stepwiseWeight(double x, arma::mat pars){
    double out = 0.0;
    arma::uword K = pars.n_rows;
    arma::uword k = 0;
    while(k != K){
        if(x <= pars(k,1)){
            out += pars(k,0);
            k = K;
        }
        else{
            k += 1;
        }
    }
    return out;
}

//' linearWeight 
//'
//' @param x transpired time
//' @param pars is the vector with order (alpha,psi) with alpha in (0,Inf) and psi in (0,Inf) 
//'
//' @export
// [[Rcpp::export]]
double linearWeight(double x, arma::mat pars){
    double out = 0.0;
    if(x <= pars(0,0)){
        out += pars(1,0) - (pars(1,0)/pars(0,0))*x;
    }
    return out;
}

//' twoStepsWeight 
//'
//' @param x transpired time
//' @param pars is the vector of parameters OLD=(shape,scale,psi,phi,shape,scale), (shape1,scale1,psi,phi,x0,shape2,scale2)
//'
//' @export
// [[Rcpp::export]]
double twoStepsWeight(double x, arma::mat pars){
    double out = 0.0;
    if(x < pars(4,0)){ // first step
        out += pars(2,0)*exp(-std::pow(x/pars(1,0), pars(0,0)))-pars(2,0)+pars(3,0);
    }
    else{ // second step
        out += (pars(3,0)-pars(2,0))*exp(-std::pow((x-pars(4,0))/pars(6,0),pars(5,0)));
        //(pars(3)-pars(2))*exp(-std::pow((x-pars(1))/(pars(5)+pars(1)), pars(4)));
    }
    return out;
}


//' weibullWeight (it can be either exponential of a smoothed single-step decrease)
//'
//' @param x transpired time
//' @param pars is the vector with order (shape, scale, psi) 
//'
//' @export
// [[Rcpp::export]]
double weibullWeight(double x, arma::mat pars){
    double out = -std::pow(x/pars(1,0), pars(0,0)); //log-density to which the maximum was substracted (shape must be greater than 1, see log())
    return (pars(2,0)*exp(out));
}


//' halflifeWeight
//'
//' @param x transpired time
//' @param pars is a column vector with the halflife
//'
//' @export
// [[Rcpp::export]]
double halflifeWeight(double x, arma::mat pars){
    double pars_loc = log(2)/pars(0,0); 
    return (pars_loc*exp((-1)*x*pars_loc));
}

//' tmaxWeight
//'
//' @param x transpired time
//' @param pars is the matrix of dim 1x1 with the time threshold
//'
//' @export
// [[Rcpp::export]]
double tmaxWeight(double x, arma::mat pars){
    double out = 0.0;
    if(x <= pars(0,0)) out += 1;
    return out;
}


// Mapping decay functions to strings
std::unordered_map<std::string, std::function<double(double, arma::mat)>>  decayMap =
        {
            {"linear", linearWeight},
            {"oneStep", weibullWeight},
            {"twoSteps", twoStepsWeight},
            {"exponential", weibullWeight},
            {"stepwise", stepwiseWeight},
            {"halflife", halflifeWeight},
            {"tmax", tmaxWeight}
        };


//' decay (function that calculates the decay according)
//'
//' @param type string defining the type of decay
//' @param x elapsed time of interest
//' @param pars is the vector with order 
//'
//' @export
// [[Rcpp::export]] 
double decay(std::string type, double x, arma::mat pars)
    {
        double out = decayMap[type](x, pars);
        return out;
    }

#endif

