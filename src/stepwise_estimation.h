#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <iterator>
#include <string>
#include <iomanip>


#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef STEPWISE_ESTIMATION_H
#define STEPWISE_ESTIMATION_H


#define LOG(x) std::cout << x << "\n"


// SndSnd
arma::mat SndSnd(arma::mat counts,
                 arma::mat intervals_backward,
                 arma::mat edgelist,
                 arma::umat riskset, 
                 arma::umat riskset_matrix, 
                 arma::vec t, 
                 arma::uword M, 
                 arma::uword k, 
                 arma::uword K_q, 
                 arma::uword n_cores)
    {   
        arma::uword m;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        
        for(m = 1; m < M; m++){ 
          arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
          out.row(m) = counts(lb_ub_row,arma::span(2,counts.n_cols-1)); //((K_q*m)+k)
        }
        return out;
    }

// RecSnd
arma::mat RecSnd(arma::mat counts,
                 arma::mat intervals_backward,
                 arma::mat edgelist,
                 arma::umat riskset,
                 arma::umat riskset_matrix, 
                 arma::vec t, 
                 arma::uword M, 
                 arma::uword k, 
                 arma::uword K_q, 
                 arma::uword n_cores)
    {
        arma::uword i,j,m; 
        arma::uword N = riskset_matrix.n_cols;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        
        for(m = 1; m < M; m++){
            arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
            for(i = 0; i < N; i++){
                for(j = 0; j < N; j++){
                    if(i!=j){
                        arma::uword column_out = riskset_matrix(i,j);
                        arma::uword column_in = riskset_matrix(j,i)+2; // because the first two columns in counts are lower and upper bound of the corresponding interval
                        out(m,column_out) = counts(lb_ub_row,column_in);
                    }
                }
            }
        }
        return out;
    }


// IDSnd
arma::mat IDSnd(arma::mat counts,
                arma::mat intervals_backward,
                arma::mat edgelist, 
                arma::umat riskset,
                arma::umat riskset_matrix, 
                arma::vec t, 
                arma::uword M, 
                arma::uword k, 
                arma::uword K_q, 
                arma::uword n_cores)
    {
        arma::uword i,j,m,n; 
        arma::uword N = riskset_matrix.n_cols;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        for(m = 1; m < M; m++){
            arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
            for(n = 0; n < N; n++){
                // + 2 (!!) ricorda !
                arma::uvec IDsender_in = riskset_matrix.col(n)+2; // indices used in counts (+2 because the first two columns in counts are lower and upper bound of the corresponding interval)
                double IDsender_value = 0;
                for(i = 0; i < N; i++){
                    if(i!=n){
                        IDsender_value += counts(lb_ub_row,IDsender_in(i));
                    }
                }
                arma::urowvec IDsender_ou = riskset_matrix.row(n); // we don't sum +2 here because the output matrix has no (lower,upper) bounds as first two columns
                for(j = 0; j < N; j++){
                    if(j!=n){
                        out(m,IDsender_ou(j)) = IDsender_value;
                    }
                }
            }
        }
        return out;
    }

// IDRec
arma::mat IDRec(arma::mat counts, 
                arma::mat intervals_backward,
                arma::mat edgelist, 
                arma::umat riskset,
                arma::umat riskset_matrix, 
                arma::vec t, 
                arma::uword M, 
                arma::uword k, 
                arma::uword K_q, 
                arma::uword n_cores)
    {
        arma::uword i,j,m,n; 
        arma::uword N = riskset_matrix.n_cols;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        for(m = 1; m < M; m++){
            arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
            for(n = 0; n < N; n++){
                // + 2 (!!) ricorda !
                arma::uvec IDreceiver_in = riskset_matrix.col(n)+2; // indices used in counts (+2 because the first two columns in counts are lower and upper bound of the corresponding interval)
                double IDreceiver_value = 0;
                for(i = 0; i < N; i++){
                    if(i!=n){
                        IDreceiver_value += counts(lb_ub_row,IDreceiver_in(i));
                    }
                }
                arma::uvec IDreceiver_ou = riskset_matrix.col(n); // we don't sum +2 here because the output matrix has no (lower,upper) bounds as first two columns
                for(j = 0; j < N; j++){
                    if(j!=n){
                        out(m,IDreceiver_ou(j)) = IDreceiver_value;
                    }
                }
            }
        }
        return out;
    }

// ODSnd
arma::mat ODSnd(arma::mat counts, 
                arma::mat intervals_backward,
                arma::mat edgelist,
                arma::umat riskset, 
                arma::umat riskset_matrix, 
                arma::vec t, 
                arma::uword M, 
                arma::uword k, 
                arma::uword K_q, 
                arma::uword n_cores)
    {
        arma::uword i,j,m,n; 
        arma::uword N = riskset_matrix.n_cols;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        for(m = 1; m < M; m++){
            arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
            for(n = 0; n < N; n++){
                // + 2 (!!) ricorda !
                arma::urowvec ODsender_in = riskset_matrix.row(n)+2; // indices used in counts (+2 because the first two columns in counts are lower and upper bound of the corresponding interval)
                double ODsender_value = 0;
                for(i = 0; i < N; i++){
                    if(i!=n){
                        ODsender_value += counts(lb_ub_row,ODsender_in(i));
                    }
                }
                arma::urowvec ODsender_ou = riskset_matrix.row(n); // we don't sum +2 here because the output matrix has no (lower,upper) bounds as first two columns
                for(j = 0; j < N; j++){
                    if(j!=n){
                        out(m,ODsender_ou(j)) = ODsender_value;
                    }
                }
            }
        }
        return out;
    }

// ODRec
arma::mat ODRec(arma::mat counts, 
                arma::mat intervals_backward,
                arma::mat edgelist, 
                arma::umat riskset,
                arma::umat riskset_matrix, 
                arma::vec t,  
                arma::uword M, 
                arma::uword k, 
                arma::uword K_q, 
                arma::uword n_cores)
    {
        arma::uword i,j,m,n; 
        arma::uword N = riskset_matrix.n_cols;
        arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
        for(m = 1; m < M; m++){
            arma::uword lb_ub_row = intervals_backward(((K_q*m)+k),2);
            for(n = 0; n < N; n++){
                // + 2 (!!) ricorda !
                arma::urowvec ODreceiver_in = riskset_matrix.row(n)+2; // indices used in counts (+2 because the first two columns in counts are lower and upper bound of the corresponding interval)
                double ODreceiver_value = 0;
                for(i = 0; i < N; i++){
                    if(i!=n){
                        ODreceiver_value += counts(lb_ub_row,ODreceiver_in(i));
                    }
                }
                arma::uvec ODreceiver_ou = riskset_matrix.col(n); // we don't sum +2 here because the output matrix has no (lower,upper) bounds as first two columns
                for(j = 0; j < N; j++){
                    if(j!=n){
                        out(m,ODreceiver_ou(j)) = ODreceiver_value;
                    }
                }
            }
        }
        return out;
    }


// TClosure (Transivity Closure)
arma::mat TClosure(arma::mat counts, 
                   arma::mat intervals_backward,
                   arma::mat edgelist,
                   arma::umat riskset,
                   arma::umat riskset_matrix, 
                   arma::vec t, 
                   arma::uword M, 
                   arma::uword k, 
                   arma::uword K_q, 
                   arma::uword n_cores)
    {
    //#pragma omp parallel for private(m,i,lb_ub,elapsed_lj,time_lj,t_star,upper_index,lower_index,upper_index_time_diff,lower_index_time_diff) shared(M,k,K_q,out,edgelist,riskset_matrix,intervals_backward)
    arma::uword m,i,d,lower_index; //,upper_index;
    double elapsed_lj,time_lj,t_star; //,upper_index_time_diff,lower_index_time_diff;
    arma::rowvec lb_ub(2);
    arma::mat out(M, riskset.n_rows, arma::fill::zeros);

    //omp_set_dynamic(0);         // disabling dynamic teams
    //omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
    //#pragma omp parallel for private(m,i,d,lb_ub,elapsed_lj,time_lj,t_star,lower_index) shared(M,k,K_q,out,edgelist,riskset_matrix,intervals_backward)
    for(m = 1; m < M; m++){
        lb_ub = intervals_backward(((K_q*m)+k),arma::span(0,1));
        if(lb_ub(0)!=(-1)){
            arma::mat edgelist_loc = edgelist(arma::span(lb_ub(0),lb_ub(1)),arma::span::all);
            for(i = 0; i < edgelist_loc.n_rows; i++){
                  // backward seeking, every event is an (l,j) type of event
                elapsed_lj = edgelist(m,0) - edgelist_loc(i,0);
                time_lj = edgelist_loc(i,0);
                t_star = time_lj - elapsed_lj; 
                if(t_star < edgelist(0,0)) t_star = edgelist(0,0);
               if(elapsed_lj>0){
                    arma::urowvec lower = arma::find(edgelist.col(0) <= t_star).t();
                    lower_index = max(lower);    
                    arma::mat edgelist_sub = edgelist(arma::span(lower_index,(lb_ub(0)+i)),arma::span::all);
                    for(d = 0; d < edgelist_sub.n_rows; d++){
                        if((edgelist_sub(d,2) == edgelist_loc(i,1)) & (edgelist_sub(d,1) != edgelist_loc(i,2))){
                                out(m,riskset_matrix(edgelist_sub(d,1),edgelist_loc(i,2))) += 1; 
                        }
                    }                                   
                }
                else{ // when elapsed_lj == 0 
                    if((lb_ub(0)+i)>0){ 
                        if((edgelist((lb_ub(0)+i-1),2) == edgelist((lb_ub(0)+i),1)) & (edgelist((lb_ub(0)+i-1),1) != edgelist((lb_ub(0)+i),2))){
                                out(m,riskset_matrix(edgelist((lb_ub(0)+i-1),1),edgelist((lb_ub(0)+i),2))) += 1; 
                        }
                    }  
                }                 
            }
        } 

    }
        return out;
    }

// CClosure (Cyclic Closure)
arma::mat CClosure (arma::mat counts, 
                   arma::mat intervals_backward,
                   arma::mat edgelist, 
                   arma::umat riskset,
                   arma::umat riskset_matrix, 
                   arma::vec t, 
                   arma::uword M, 
                   arma::uword k, 
                   arma::uword K_q, 
                   arma::uword n_cores)
    {
    arma::uword m,i,l,h,s,r,l_s_loc,r_l_loc,column_out;
    arma::uword N = riskset_matrix.n_rows;
    int lb_ub_row, lb_row, ub_row, m_loc, lb_ub_l_s;
    double CClosure_value_s_r;
    arma::uvec indices(2);
    arma::mat out(M, riskset.n_rows, arma::fill::zeros);


    // create collapsed matrix for the parallelization loop
    arma::umat collapsed_loops((M-1)*N*(N-1), 4, arma::fill::zeros);
    arma::uvec ones_vec(N*(N-1),arma::fill::ones);
        for(m = 0; m < M-1; m++){
            collapsed_loops(arma::span(m*N*(N-1), (m+1)*N*(N-1)-1),arma::span(1,2)) = riskset;
            collapsed_loops(arma::span(m*N*(N-1), (m+1)*N*(N-1)-1),0) = (m+1)*ones_vec;
        }
        for(r = 1; r < collapsed_loops.n_rows; r++){
            collapsed_loops(r,3) = riskset_matrix(collapsed_loops(r,1),collapsed_loops(r,2));
        }
    //

    omp_set_dynamic(0);         // disabling dynamic teams
    omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
    #pragma omp parallel for private(i,m,lb_ub_row,lb_row,ub_row,s,r,column_out,indices,CClosure_value_s_r,l,l_s_loc,r_l_loc,h,m_loc,lb_ub_l_s) shared(M,out,k,K_q,counts,intervals_backward,riskset_matrix,collapsed_loops)
        for(i = 0; i < collapsed_loops.n_rows; i++){
            m = collapsed_loops(i,0);
            lb_ub_row = intervals_backward(((K_q*m)+k),2);
            lb_row = counts(lb_ub_row,0);
            ub_row = counts(lb_ub_row,1);
            s = collapsed_loops(i,1);
            r = collapsed_loops(i,2);
            column_out = collapsed_loops(i,3);
            indices = {s,r};

            // vector of column indices of (l,r)
            arma::umat loc_matrix_l_s = riskset_matrix;
            loc_matrix_l_s.shed_rows(indices);
            arma::uvec l_s = loc_matrix_l_s.col(s);
            // vector of column indices of (s,l)
            arma::umat loc_matrix_r_l = riskset_matrix;
            loc_matrix_r_l.shed_cols(indices);
            arma::uvec r_l = loc_matrix_r_l.row(r).t();

            CClosure_value_s_r = 0;

            for(l = 0; l < r_l.n_elem; l++){
                l_s_loc = l_s(l);
                r_l_loc = r_l(l);
                // backward seeking
                if(counts(lb_ub_row,(l_s_loc+2)) > 0){
                    arma::vec l_s_binary = counts(arma::span(lb_row,ub_row),l_s_loc+2);
                    arma::uvec l_s_ones = arma::find(l_s_binary == 1);
                    for(h = 0; h < l_s_ones.n_elem; h++){
                        m_loc = counts(lb_ub_row,0) + l_s_ones(h);
                        lb_ub_l_s = intervals_backward(((K_q*m_loc)+k),2);
                        if(counts(lb_ub_l_s,(r_l_loc+2)) > 0){
                            CClosure_value_s_r += 1;
                        }
                       }
                }
            }
            out(m,column_out) = CClosure_value_s_r;
        }
        return out;
    }
// Mapping participation shifts functions to strings (which will be written as input in the .R function)
std::unordered_map<std::string, std::function<arma::mat(arma::mat, arma::mat, arma::mat, arma::umat, arma::umat, arma::vec, arma::uword, arma::uword, arma::uword, arma::uword)>>  endogenousMap =
        {
            {"SndSnd", SndSnd},
            {"RecSnd", RecSnd},
            {"IDSnd", IDSnd},
            {"IDRec", IDRec},
            {"ODSnd", ODSnd},
            {"ODRec", ODRec},
            {"TClosure", TClosure},
            {"CClosure", CClosure}
        };


// computeEffect (internal routine function)
arma::mat computeEffect(std::string effect,
                        arma::mat counts, 
                        arma::mat intervals_backward,
                        arma::mat edgelist, 
                        arma::umat riskset,
                        arma::umat riskset_matrix, 
                        arma::vec t, 
                        arma::uword M, 
                        arma::uword k, 
                        arma::uword K_q, 
                        arma::uword n_cores)
    {
        arma::mat out = endogenousMap[effect](counts, intervals_backward, edgelist, riskset, riskset_matrix,  t, M, k, K_q, n_cores);
        return out;
    }
#endif