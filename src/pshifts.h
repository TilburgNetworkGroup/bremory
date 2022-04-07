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

#ifndef PSHIFTS_H
#define PSHIFTS_H


// // Participation Shifts // // 

// Turn receiving 
arma::vec ABBA(arma::uword actor1, arma::uword actor2, arma::umat riskset, arma::vec out)
    {
        arma::uword index_BA = riskset(actor2,actor1);
        out(index_BA) += 1.0;
        return out;
    }

arma::vec ABBY(arma::uword actor1, arma::uword actor2, arma::umat riskset, arma::vec out)
    {
        arma::uvec cols_to_shed  = {actor1,actor2};
        riskset.shed_cols(cols_to_shed);
        arma::uvec indices_BY = riskset.row(actor2).t();
        out(indices_BY) += 1.0;
        return out;
    }

// Turn continuing

arma::vec ABAY(arma::uword actor1, arma::uword actor2, arma::umat riskset, arma::vec out)
    {
        arma::uvec cols_to_shed  = {actor1,actor2};
        riskset.shed_cols(cols_to_shed);
        arma::uvec indices_AY = riskset.row(actor1).t();
        out(indices_AY) += 1.0;
        return out;
    }
  

// Mapping participation shifts functions to strings
std::unordered_map<std::string, std::function<arma::vec(arma::uword,arma::uword,arma::umat,arma::vec)>> pshiftsMap =
         {
            {"ABBA", ABBA},
            {"ABBY", ABBY},
            {"ABAY", ABAY}
         };

#endif

