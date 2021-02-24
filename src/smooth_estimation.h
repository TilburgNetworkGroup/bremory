#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <iterator>
#include <string>
#include "decay_functions.h"
#include "pshifts.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SMOOTH_ESTIMATION_H
#define SMOOTH_ESTIMATION_H


#define LOG(x) std::cout << x << "\n"

// SndSnd
arma::mat smoothSndSnd(Rcpp::List pars_endo,
                arma::umat edgelist,
                arma::mat rehBinary, 
                arma::umat risksetMatrix,
                arma::ucube risksetCube, 
                arma::vec t, 
                arma::uword M,
                arma::uword N,
                arma::uword C,
                arma::uword D,
                arma::uword n_cores)
    {   
        arma::uword m,dy,m_loc;
        arma::mat out(M,D,arma::fill::zeros);
        std::string decay_fun = pars_endo[0];
        arma::mat pars_loc = pars_endo[1];

        omp_set_dynamic(0);         // disabling dynamic teams
        omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
        #pragma omp parallel for private(m,dy,m_loc) shared(t,out,rehBinary,M,D)
        for(m = 1; m < M; m++){ // start from m = 1 meaning the second time point (at the first all the statistics are assumed to have value equal to zero)
            // elapsed time matrix
            double current_time = t(m);
            arma::vec current_time_vec(m) ;
            current_time_vec.fill(current_time);
            arma::vec elapsed_time = current_time_vec - t(arma::span(0,m-1)) ;
            arma::mat elapsed_time_matrix_loc(m,D,arma::fill::zeros);
            for(dy = 0; dy < D; dy++){
                elapsed_time_matrix_loc.col(dy) = elapsed_time;
            }
            arma::mat elapsed_time_matrix = rehBinary(arma::span(0,(m-1)),arma::span::all) % elapsed_time_matrix_loc;
            // computing statistics
            arma::vec out_loc(D, arma::fill::zeros); 
            for(dy = 0; dy < D; dy++){
                for(m_loc = 0; m_loc < m; m_loc++){
                    if(rehBinary(m_loc,dy) == 1){
                        out_loc(dy) += decay(decay_fun,elapsed_time_matrix(m_loc,dy),pars_loc); // 'halflife' by default
                    }     
                }
            }
            out.row(m) += out_loc.t();  
        }
        return out;
    }

// RecSnd
arma::mat smoothRecSnd(Rcpp::List pars_endo,
                arma::umat edgelist,
                arma::mat rehBinary, 
                arma::umat risksetMatrix,
                arma::ucube risksetCube,
                arma::vec t, 
                arma::uword M,
                arma::uword N,
                arma::uword C,
                arma::uword D,
                arma::uword n_cores)
    {   
        arma::uword m,d,m_loc;
        arma::mat out(M,D,arma::fill::zeros);
        std::string decay_fun = pars_endo[0];
        arma::mat pars_loc = pars_endo[1];

        omp_set_dynamic(0);         // disabling dynamic teams
        omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
        #pragma omp parallel for private(m,d,m_loc) shared(t,out,rehBinary,risksetMatrix,risksetCube,M,D)
        for(m = 1; m < M; m++){ // start from m = 1 meaning the second time point (at the first all the statistics are assumed to have value equal to zero)
            // elapsed time matrix
            double current_time = t(m);
            arma::vec current_time_vec(m) ;
            current_time_vec.fill(current_time);
            arma::vec elapsed_time = current_time_vec - t(arma::span(0,m-1)) ;
            arma::mat elapsed_time_matrix_loc(m,D,arma::fill::zeros);
            for(d = 0; d < D; d++){
                elapsed_time_matrix_loc.col(d) = elapsed_time;
            }
            arma::mat elapsed_time_matrix = rehBinary(arma::span(0,(m-1)),arma::span::all) % elapsed_time_matrix_loc;
            // computing statistics
            arma::vec out_loc(D, arma::fill::zeros); 
            for(d = 0; d < D; d++){
                arma::uword actor1 = risksetMatrix(d,0);
                arma::uword actor2 = risksetMatrix(d,1);
                arma::uword type_ = risksetMatrix(d,2);
                arma::uword d_ji = risksetCube(actor2,actor1,type_);
                for(m_loc = 0; m_loc < m; m_loc++){
                    if(rehBinary(m_loc,d) == 1){
                        out_loc(d_ji) += decay(decay_fun,elapsed_time_matrix(m_loc,d),pars_loc); // 'halflife' by default
                    }     
                }
            }
            out.row(m) = out_loc.t();  
        }
        return out;
    }


// updateTClosure
arma::mat smoothTClosure(Rcpp::List pars_endo, 
                            arma::umat edgelist,
                            arma::mat rehBinary, 
                            arma::umat risksetMatrix,
                            arma::ucube risksetCube,
                            arma::vec t, 
                            arma::uword M,
                            arma::uword N,
                            arma::uword C,
                            arma::uword D,
                            arma::uword n_cores)
        {
  /*      Rcpp::List pars_endo,
                        Rcpp::Environment binaryREH,
                        arma::mat elapsed_time_matrix_loc,
                        arma::mat edgelist,
                        arma::umat risksetMatrix,
                        arma::ucube risksetCube,
                        arma::uword D,
                        arma::uword M
          */      
        arma::uword m,dy,m_loc;
        arma::mat out(M,D,arma::fill::zeros); 
        std::string decay_fun = pars_endo[0];
        arma::mat pars_loc = pars_endo[1];

        omp_set_dynamic(0);         // disabling dynamic teams
        omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
        #pragma omp parallel for private(m,dy,m_loc) shared(out,t,rehBinary,M,D,risksetCube)
        for(m = 1; m < M; m++){
            // elapsed time matrix
            double current_time = t(m);
            arma::vec current_time_vec(m);
            current_time_vec.fill(current_time);
            arma::vec elapsed_time = current_time_vec - t(arma::span(0,m-1));
            arma::vec out_loc(D, arma::fill::zeros);

            for(m_loc = 0; m_loc < m; m_loc++){
                double elapsed_lj = elapsed_time(m_loc);
                double time_lj = t(m_loc);
                double t_star = time_lj - elapsed_lj;
                if(t_star < t(0)) t_star = t(0);
                double weight_loc = decay(decay_fun,elapsed_lj,pars_loc);
                if(elapsed_lj>0){
                        arma::urowvec lower = arma::find(t <= t_star).t();
                        arma::uword  lower_index = max(lower);
                        arma::umat edgelist_sub = edgelist(arma::span(lower_index,m_loc),arma::span::all);
                        for(dy = 0; dy < edgelist_sub.n_rows; dy++){
                            if((edgelist_sub(dy,2) == edgelist(m_loc,2)) & (edgelist_sub(dy,1) == edgelist(m_loc,0)) & (edgelist_sub(dy,0) != edgelist(m_loc,1))){
                                out_loc(risksetCube(edgelist_sub(dy,0),edgelist(m_loc,1),edgelist_sub(dy,2))) += weight_loc; 
                            }
                        }
                }
                else{
                    if(m_loc>0){  // when elapsed_lj == 0             
                        if((edgelist(m_loc-1,2) == edgelist(m_loc,2)) & (edgelist(m_loc-1,1) == edgelist(m_loc,0)) & (edgelist(m_loc-1,0) != edgelist(m_loc,1))){
                            out_loc(risksetCube(edgelist(m_loc-1,0),edgelist(m_loc,1),edgelist(dy,2))) += weight_loc; 
                        } 
                    }
                }
            }
            out.row(m) = out_loc.t();
        }

        return out;
    }

// smoothPShift
arma::mat smoothPShift(Rcpp::List pars_endo,
                            arma::umat edgelist,
                            arma::mat rehBinary, 
                            arma::umat risksetMatrix,
                            arma::ucube risksetCube,
                            arma::vec t, 
                            arma::uword M,
                            arma::uword N,
                            arma::uword C,
                            arma::uword D,
                            arma::uword n_cores)
    {
        arma::uword c,m;
        arma::mat out(M,D,arma::fill::zeros);
        std::string pshift = pars_endo[0];
        arma::umat riskset_matrix_0 = risksetCube.slice(0);

        
        //omp_set_dynamic(0);         // disabling dynamic teams
        //omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
        //#pragma omp parallel for private(c,m) shared(t,out,rehBinary,M,D,C)

        for(m = 1; m < M; m++){
            arma::uword actor1 = edgelist(m-1,0);
            arma::uword actor2 = edgelist(m-1,1);
            arma::uword type = edgelist(m-1,2);
            
            for(c = 0; c < C; c++){
               if(type == c){ // pshift is computed only between events of the same type
                    arma::uword lb_d = (c*(N*(N-1)));
                    arma::uword ub_d = (c+1)*(N*(N-1))-1;
                    arma::vec out_loc((N*(N-1)),arma::fill::zeros);
                    out_loc = pshiftsMap[pshift](actor1,actor2,riskset_matrix_0,out_loc);
                    out(m,arma::span(lb_d,ub_d)) = out_loc.t();
                }
            }
        }

        return out;
    }

// Mapping participation shifts functions to strings (which will be written as input in the .R function)
std::unordered_map<std::string, std::function<arma::mat(Rcpp::List, arma::umat, arma::mat, arma::umat, arma::ucube, arma::vec, arma::uword, arma::uword, arma::uword, arma::uword, arma::uword)>>  smoothEndogenousMap =
        {
            {"SndSnd", smoothSndSnd},
            {"RecSnd", smoothRecSnd},
            {"TClosure", smoothTClosure},
            {"PShift", smoothPShift}
        };

// computeEffect (internal routine function)
arma::mat computeSmoothEffect(std::string effect,
                        Rcpp::List pars_endo,
                        arma::umat edgelist,
                        arma::mat rehBinary, 
                        arma::umat risksetMatrix,
                        arma::ucube risksetCube, 
                        arma::vec t, 
                        arma::uword M,
                        arma::uword N,
                        arma::uword C,
                        arma::uword D,
                        arma::uword n_cores)
    {
        arma::mat out = smoothEndogenousMap[effect](pars_endo, edgelist, rehBinary, risksetMatrix, risksetCube, t, M, N, C, D, n_cores);
        return out;
    }
#endif