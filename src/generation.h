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

#ifndef GENERATION_H
#define GENERATION_H

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


// Mapping decay functions to strings
std::unordered_map<std::string, std::function<double(double, arma::mat)>>  decayMap =
        {
            {"linear", linearWeight},
            {"oneStep", weibullWeight},
            {"twoSteps", twoStepsWeight},
            {"exponential", weibullWeight},
            {"stepwise", stepwiseWeight}

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
 
// updateSndSnd
arma::vec updateSndSnd(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword m,d;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::mat dyadicREH = binaryREH["dyadicREH"];
        arma::mat elapsed_time_matrix = dyadicREH(arma::span(0,(M_partial-1)),arma::span::all) % elapsed_time_matrix_loc; // % compute the element-wise multiplication
        for(d = 0; d < N*(N-1); d++){
            for(m = 0; m < M_partial; m++){
                if(dyadicREH(m,d) == 1){
                    out(d) += decay(beta_endo[0],elapsed_time_matrix(m,d),beta_endo[1]);
                }     
            }
        }
        return out;
    }

// updateRecSnd
arma::vec updateRecSnd(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword m,d;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::mat dyadicREH = binaryREH["dyadicREH"];
        arma::mat elapsed_time_matrix = dyadicREH(arma::span(0,(M_partial-1)),arma::span::all) % elapsed_time_matrix_loc; // % compute the element-wise multiplication
        for(d = 0; d < N*(N-1); d++){
            for(m = 0; m < M_partial; m++){
                if(dyadicREH(m,d) == 1){
                    out(riskset_matrix(riskset(d,1),riskset(d,0))) += decay(beta_endo[0],elapsed_time_matrix(m,d),beta_endo[1]);
                }
            }
        }
        return out;
    }

// updateIDSnd
arma::vec updateIDSnd(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword i,j,n;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::vec SndSnd_vec = updateSndSnd(beta_endo, binaryREH, elapsed_time_matrix_loc, edgelist, riskset, riskset_matrix, N, M_partial);
        for(n = 0; n < N; n++){
            double IDSender_n = 0;
            for(i = 0; i < N; i++){
                if(i!=n){
                        IDSender_n += SndSnd_vec(riskset_matrix(i,n));
                }
            }
            arma::urowvec IDsender_row = riskset_matrix.row(n); //
            for(j = 0; j < N; j++){
                if(j!=n){
                    out(IDsender_row(j)) = IDSender_n;
                }
            }
        }
    return out;
    }

// updateIDRec
arma::vec updateIDRec(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword i,j,n;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::vec SndSnd_vec = updateSndSnd(beta_endo, binaryREH, elapsed_time_matrix_loc, edgelist, riskset, riskset_matrix, N, M_partial);
        for(n = 0; n < N; n++){
            double IDReceiver_n = 0;
            for(i = 0; i < N; i++){
                if(i!=n){
                        IDReceiver_n += SndSnd_vec(riskset_matrix(i,n));
                }
            }
            arma::uvec IDReceiver_col = riskset_matrix.col(n); //
            for(j = 0; j < N; j++){
                if(j!=n){
                    out(IDReceiver_col(j)) = IDReceiver_n;
                }
            }
        }
    return out;

    }

// updateODSnd
arma::vec updateODSnd(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword i,j,n;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::vec SndSnd_vec = updateSndSnd(beta_endo, binaryREH, elapsed_time_matrix_loc, edgelist, riskset, riskset_matrix, N, M_partial);
        for(n = 0; n < N; n++){
            double ODSender_n = 0;
            for(i = 0; i < N; i++){
                if(i!=n){
                        ODSender_n += SndSnd_vec(riskset_matrix(n,i));
                }
            }
            arma::urowvec ODSender_row = riskset_matrix.row(n); //
            for(j = 0; j < N; j++){
                if(j!=n){
                    out(ODSender_row(j)) = ODSender_n;
                }
            }
        }
    return out;
    }

// updateODRec
arma::vec updateODRec(Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::uword i,j,n;
        arma::vec out(N*(N-1), arma::fill::zeros);
        arma::vec SndSnd_vec = updateSndSnd(beta_endo, binaryREH, elapsed_time_matrix_loc, edgelist, riskset, riskset_matrix, N, M_partial);
        for(n = 0; n < N; n++){
            double ODReceiver_n = 0;
            for(i = 0; i < N; i++){
                if(i!=n){
                        ODReceiver_n += SndSnd_vec(riskset_matrix(n,i));
                }
            }
            arma::uvec ODreceiver_col = riskset_matrix.col(n); //
            for(j = 0; j < N; j++){
                if(j!=n){
                    out(ODreceiver_col(j)) = ODReceiver_n;
                }
            }
        }
    return out;
    }

// updateCClosure
arma::vec updateCClosure(Rcpp::List beta_endo,
                        Rcpp::Environment binaryREH,
                        arma::mat elapsed_time_matrix_loc,
                        arma::mat edgelist,
                        arma::umat riskset,
                        arma::umat riskset_matrix,
                        arma::uword N,
                        arma::uword M_partial)
    {
        arma::uword m,d;
        arma::vec out(N*(N-1),arma::fill::zeros); 
        arma::vec elapsed_time_vector = elapsed_time_matrix_loc.col(0);
        double epsilon = 0.005;
        arma::vec lb_vec(2);
        arma::vec ub_vec(2);

        for(m = 0; m < (M_partial-1); m++){
            // backward seeking, every event is an (l,i) type of event
            arma::uword li_index = M_partial-1-m;
            double elapsed_li = elapsed_time_vector(li_index);
            double weight_loc = decay(beta_endo[0],elapsed_li,beta_endo[1]);
            double time_li = edgelist(li_index,0);
            double t_star = time_li - elapsed_li;
            if(t_star >= edgelist(0,0)){
                lb_vec(0) = edgelist(0,0);
                lb_vec(1) =  t_star - epsilon;
                double lb = max(lb_vec);
                ub_vec(0) = t_star + epsilon;
                ub_vec(1) = time_li;
                double ub = min(ub_vec);
                arma::mat edgelist_sub = edgelist(arma::span(0,li_index),arma::span::all); //
                arma::urowvec lb_row = arma::find(edgelist_sub.col(0) >= lb).t();
                arma::urowvec ub_row = arma::find(edgelist_sub.col(0) <= ub).t();

                if(lb_row[0] <= ub_row[ub_row.n_elem-1]){
                    arma::mat edgelist_lb_ub = edgelist_sub(arma::span(lb_row[0],ub_row[ub_row.n_elem-1]),arma::span::all);
                    for(d = 0; d < edgelist_lb_ub.n_rows; d++){
                        if((edgelist_lb_ub(d,2) == edgelist(li_index,1)) & (edgelist_lb_ub(d,1) != edgelist(li_index,2))){
                                out(riskset_matrix(edgelist(li_index,2),edgelist_lb_ub(d,1))) += weight_loc; //
                        }
                    }
                }
            }
        }
        return out;
    }


// updateTClosure
arma::vec updateTClosure(Rcpp::List beta_endo,
                        Rcpp::Environment binaryREH,
                        arma::mat elapsed_time_matrix_loc,
                        arma::mat edgelist,
                        arma::umat riskset,
                        arma::umat riskset_matrix,
                        arma::uword N,
                        arma::uword M_partial)
    {
        arma::uword m,d;
        arma::vec out(N*(N-1),arma::fill::zeros); 
        arma::vec elapsed_time_vector = elapsed_time_matrix_loc.col(0);

        for(m = 0; m < M_partial; m++){
            // commented below is the old definition of TClosure (before mid June 2020)
            //// backward seeking, every event is an (l,j) type of event
            //double elapsed_lj = elapsed_time_vector(m);
            //double time_lj = edgelist(m,0);
            //double t_star = time_lj - elapsed_lj;
            //if((t_star > edgelist(0,0)) & (elapsed_lj>0)){
            //    double weight_loc = decay(beta_endo[0],elapsed_lj,beta_endo[1]);
            //    arma::mat edgelist_sub = edgelist(arma::span(0,m),arma::span::all); 
            //    arma::urowvec greater = arma::find(edgelist_sub.col(0) >= t_star).t();
            //    arma::urowvec lower = arma::find(edgelist_sub.col(0) < t_star).t();
            //    arma::uword  upper_index = min(greater);
            //    arma::uword  lower_index = max(lower);
            //    double upper_index_time_diff = (edgelist_sub(upper_index,0) - t_star);
            //    double lower_index_time_diff = (t_star - edgelist_sub(lower_index,0));               
            //    if(upper_index_time_diff < lower_index_time_diff){
            //        if((edgelist_sub(upper_index,2) == edgelist(m,1)) & (edgelist_sub(upper_index,1) != edgelist(m,2))){
            //                    out(riskset_matrix(edgelist_sub(upper_index,1),edgelist(m,2))) += weight_loc; 
            //        }
            //    }
            //    else{
            //        if((edgelist_sub(lower_index,2) == edgelist(m,1)) & (edgelist_sub(lower_index,1) != edgelist(m,2))){
            //                    out(riskset_matrix(edgelist_sub(lower_index,1),edgelist(m,2))) += weight_loc; 
            //        }
            //    }
            //}
            double elapsed_lj = elapsed_time_vector(m);
            double time_lj = edgelist(m,0);
            double t_star = time_lj - elapsed_lj;
            double weight_loc = decay(beta_endo[0],elapsed_lj,beta_endo[1]);
            if((elapsed_lj>0) & (t_star >= edgelist(0,0))){
                    arma::urowvec lower = arma::find(edgelist.col(0) <= t_star).t();
                    arma::uword  lower_index = max(lower);
                    arma::mat edgelist_sub = edgelist(arma::span(lower_index,m),arma::span::all);
                    for(d = 0; d < edgelist_sub.n_rows; d++){
                        if((edgelist_sub(d,2) == edgelist(m,1)) & (edgelist_sub(d,1) != edgelist(m,2))){
                                out(riskset_matrix(edgelist_sub(d,1),edgelist(m,2))) += weight_loc; 
                        }
                    }
            }
            if((t_star >= edgelist(0,0)) & (m>0)){  // when elapsed_lj == 0             
                if((edgelist(m-1,2) == edgelist(m,1)) & (edgelist(m-1,1) != edgelist(m,2))){
                                out(riskset_matrix(edgelist(m-1,1),edgelist(m,2))) += weight_loc; 
                } 
            }

        }
        return out;
    }



// Mapping update endogenous functions to strings (which will be written as input in the .R function)
std::unordered_map<std::string, std::function<arma::vec(Rcpp::List, Rcpp::Environment, arma::mat, arma::mat, arma::umat, arma::umat, arma::uword, arma::uword)>>  updateEndogenousMap =
        {
            {"SndSnd", updateSndSnd},
            {"RecSnd", updateRecSnd},
            {"IDSnd", updateIDSnd},
            {"ODSnd", updateODSnd},
            {"IDRec", updateIDRec},
            {"ODRec", updateODRec},
            {"CClosure", updateCClosure},
            {"TClosure", updateTClosure}
        };


// updateEffect 
arma::vec updateEffect(std::string effect,
                       Rcpp::List beta_endo,
                       Rcpp::Environment binaryREH,
                       arma::mat elapsed_time_matrix_loc,
                       arma::mat edgelist,
                       arma::umat riskset,
                       arma::umat riskset_matrix,
                       arma::uword N,
                       arma::uword M_partial)
    {
        arma::vec out = updateEndogenousMap[effect](beta_endo, binaryREH, elapsed_time_matrix_loc, edgelist, riskset, riskset_matrix, N, M_partial);
        return out;
    }


#endif
