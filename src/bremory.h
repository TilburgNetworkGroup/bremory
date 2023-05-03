#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <iterator>
#include <string>
#include <remify/utils.h>

// <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef bremory_H
#define bremory_H

// Tie-oriented modeling

int bremory_getDyadIndex(double actor1, double actor2, double type, int N, bool directed) {

    int dyad = -999; // returning impossible index if the dyad is a self-edge (i.e., sender and receiver are the same actor)
    if(actor1 != actor2){
        if(!directed){ // when directed == FALSE we sort actor1 and actor2
            if(actor1 < actor2){
                double dyad_loc = (N*(N-1)/2)*type+(N-1)*actor1+actor2-actor1-1-(actor1*actor1)/2;
                if(actor1>0){
                    dyad_loc += actor1/2;
                }
                dyad = dyad_loc;
            }
            else{
                double dyad_loc = (N*(N-1)/2)*type+(N-1)*actor2+actor1-actor2-1-(actor2*actor2)/2;
                if(actor2>0){
                    dyad_loc += actor2/2;
                }
                dyad = dyad_loc;
            }
        }
        else{ 
            // when directed == TRUE (we do not sort) (actor1 = sender, actor2 = receiver)
            double dyad_loc = N*(N-1)*type+(N-1)*actor1+actor2;
            if(actor2>actor1){
                dyad_loc -= 1;
            }
            dyad = dyad_loc;
        }
    }

    return dyad;
}



// creating a namespace 'stepwiseTie' with all the routines for computing stepwise endogenous statistics (for tie-oriented models)
namespace stepwiseTie
{
    // stepwiseTie::inertia
    arma::mat inertia(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {   
            arma::uword m;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
            
            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m) shared(intervals,M,counts,out)
            for(m = 1; m < M; m++){ 
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                out.row(m) = counts(lb_ub_row,arma::span(2,counts.n_cols-1)); //((K*m)+k)
            }
            return out;
        }

    // stepwiseTie::reciprocity
    arma::mat reciprocity(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M, 
                    int N,
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword m;
            int i,j;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            omp_set_dynamic(0);           // disabling dynamic teams
            omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            #pragma omp parallel for private(m,i,j) shared(intervals,out,M,N,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(i = 0; i < N; i++){
                    for(j = 0; j < N; j++){
                        if(i!=j){
                            int column_out = remify::utils::getDyadIndex(i,j,0,N,true);
                            int column_in = remify::utils::getDyadIndex(j,i,0,N,true)+2; // because the first two columns in counts are lower and upper bound of the corresponding interval
                            out(m,column_out) = counts(lb_ub_row,column_in);
                        }
                    }
                }
            }

            return out;
        }



    // stepwiseTie::indegreeSender
    arma::mat indegreeSender(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            omp_set_dynamic(0);           // disabling dynamic teams
            omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            #pragma omp parallel for private(m,i,j,n) shared(intervals,out,M,N,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating indegree for sender 'n'
                    double IDsender_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            int IDsender_in = remify::utils::getDyadIndex(i,n,0,N,true)+2;
                            //#pragma omp critical
                                IDsender_value += counts(lb_ub_row,IDsender_in);
                        }
                    }
                    // assigning indegree to dyads where 'n' is the sender
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int IDsender_ou = remify::utils::getDyadIndex(n,j,0,N,true);
                            out(m,IDsender_ou) = IDsender_value;
                        }
                    }
                }
            }
            return out;
        }

    // stepwiseTie::indegreeReceiver
    arma::mat indegreeReceiver(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m,i,j,n,K,k) shared(intervals,out,M,N,actor1,actor2,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating indegree for receiver 'n'
                    double IDreceiver_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            int IDreceiver_in = bremory_getDyadIndex(i,n,0,N,true)+2;
                            //#pragma omp critical
                                IDreceiver_value += counts(lb_ub_row,IDreceiver_in);
                        }
                    }
                    // assigning indegree to dyads where 'n' is the receiver
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int IDreceiver_ou = bremory_getDyadIndex(j,n,0,N,true);
                            out(m,IDreceiver_ou) = IDreceiver_value;
                        }
                    }
                }
            }
            return out;
        }

    // stepwiseTie::outdegreeSender
    arma::mat outdegreeSender(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n; 
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m,i,j,n,K,k) shared(intervals,out,M,N,actor1,actor2,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating outdegree of sender 'n'
                    double ODsender_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            int ODsender_in = bremory_getDyadIndex(n,i,0,N,true)+2;
                            ODsender_value += counts(lb_ub_row,ODsender_in);
                        }
                    }
                    // assigning outdegree to dyads where 'n' is the sender
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int ODsender_ou = bremory_getDyadIndex(n,j,0,N,true);
                            out(m,ODsender_ou) = ODsender_value;
                        }
                    }
                }
            }
            return out;
        }
    
    

    // stepwiseTie::outdegreeReceiver
    arma::mat outdegreeReceiver(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m,i,j,n,K,k) shared(intervals,out,M,N,actor1,actor2,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating outdegree of receiver 'n'
                    double ODreceiver_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            int ODreceiver_in = bremory_getDyadIndex(n,i,0,N,true)+2;
                            ODreceiver_value += counts(lb_ub_row,ODreceiver_in);
                        }
                    }
                    // assigning outdegree to dyads where 'n' is the receiver
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int ODreceiver_ou = bremory_getDyadIndex(j,n,0,N,true);
                            out(m,ODreceiver_ou) = ODreceiver_value;
                        }
                    }
                }
            }
            return out;
        }

     // stepwiseTie::totaldegreeSender
    arma::mat totaldegreeSender(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m,i,j,n,K,k) shared(intervals,out,M,N,actor1,actor2,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating indegree for sender 'n'
                    double totaldegree_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            // when 'n' is the receiver
                            int as_receiver = bremory_getDyadIndex(i,n,0,N,true)+2;
                            //#pragma omp critical
                                totaldegree_value += counts(lb_ub_row,as_receiver);

                            // when 'n'is the sender
                            int as_sender = bremory_getDyadIndex(n,i,0,N,true)+2;
                            //#pragma omp critical
                                totaldegree_value += counts(lb_ub_row,as_sender);
                        }
                    }
                    // assigning indegree to dyads where 'n' is the sender
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int sender_ou = bremory_getDyadIndex(n,j,0,N,true);
                            out(m,sender_ou) = totaldegree_value;
                        }
                    }
                }
            }
            return out;
        }

    // stepwiseTie::totaldegreeReceiver
    arma::mat totaldegreeReceiver(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M,
                    int N, 
                    arma::uword k, 
                    arma::uword K, 
                    int ncores)
        {
            arma::uword i,j,m,n;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);

            //omp_set_dynamic(0);           // disabling dynamic teams
            //omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
            //#pragma omp parallel for private(m,i,j,n,K,k) shared(intervals,out,M,N,actor1,actor2,counts)
            for(m = 1; m < M; m++){
                arma::uword lb_ub_row = intervals(((K*m)+k),2);
                for(n = 0; n < N; n++){
                    // calculating indegree for receiver 'n'
                    double totaldegree_value = 0;
                    for(i = 0; i < N; i++){
                        if(i!=n){
                            // when 'n' is the receiver
                            int as_receiver = bremory_getDyadIndex(i,n,0,N,true)+2;
                            //#pragma omp critical
                                totaldegree_value += counts(lb_ub_row,as_receiver);

                            // when 'n'is the sender
                            int as_sender = bremory_getDyadIndex(n,i,0,N,true)+2;
                            //#pragma omp critical
                                totaldegree_value += counts(lb_ub_row,as_sender);
                            
                        }
                    }
                    // assigning indegree to dyads where 'n' is the receiver
                    for(j = 0; j < N; j++){
                        if(j!=n){
                            int receiver_ou = bremory_getDyadIndex(j,n,0,N,true);
                            out(m,receiver_ou) = totaldegree_value;
                        }
                    }
                }
            }
            return out;
        }

    // stepwiseTie::transitivityClosure
    arma::mat transitivityClosure(arma::mat counts, 
                    arma::mat intervals_backward,
                    arma::mat edgelist,
                    arma::umat riskset,
                    arma::umat riskset_matrix, 
                    arma::vec t, 
                    arma::uword M, 
                    arma::uword k, 
                    arma::uword K_q, 
                    int ncores)
        {
        //#pragma omp parallel for private(m,i,lb_ub,elapsed_lj,time_lj,t_star,upper_index,lower_index,upper_index_time_diff,lower_index_time_diff) shared(M,k,K_q,out,edgelist,riskset_matrix,intervals_backward)
        arma::uword m,i,d,lower_index; //,upper_index;
        double elapsed_lj,time_lj,t_star; //,upper_index_time_diff,lower_index_time_diff;
        arma::rowvec lb_ub(2);
        arma::mat out(M, riskset.n_rows, arma::fill::zeros);

        //omp_set_dynamic(0);         // disabling dynamic teams
        //omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
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

    // stepwiseTie::cyclicClosure
    arma::mat cyclicClosure(arma::mat counts, 
                    arma::mat intervals_backward,
                    arma::mat edgelist, 
                    arma::umat riskset,
                    arma::umat riskset_matrix, 
                    arma::vec t, 
                    arma::uword M, 
                    arma::uword k, 
                    arma::uword K_q, 
                    int ncores)
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
        omp_set_num_threads(ncores); // number of threads for all consecutive parallel regions
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
    std::unordered_map<std::string, std::function<arma::mat(arma::mat, arma::mat, arma::vec, arma::vec, arma::vec, arma::uword, int, arma::uword, arma::uword, arma::uword)>>  endogenousMap =
            {
                {"inertia", inertia}
                ,{"reciprocity", reciprocity}
                ,{"indegreeSender", indegreeSender}
                ,{"indegreeReceiver", indegreeReceiver}
                ,{"outdegreeSender", outdegreeSender}
                ,{"outdegreeReceiver", outdegreeReceiver}
                ,{"totaldegreeSender",totaldegreeSender}
                ,{"totaldegreeReceiver", totaldegreeReceiver}
                //,{"transitivityClosure", transitivityClosure}
                //,{"cyclicClosure", cyclicClosure}
    // missing stats: "totaldegreeSender", "totaldegreeReceiver", "sendingClosure", "receivingClosure", "pshift"
            };




    // computeIntervalStatistic (internal routine function)
    arma::mat computeStatistic(std::string statistic,
                            arma::mat counts, 
                            arma::mat intervals, // intervals_backward
                            arma::vec actor1,
                            arma::vec actor2,
                            arma::vec time, // t
                            arma::uword M, 
                            int N,
                            arma::uword k, 
                            arma::uword K, //K_q
                            int ncores)
        {
            arma::mat out = endogenousMap[statistic](counts, intervals, actor1, actor2, time, M, N, k, K, ncores);
            return out;
        }

}


// Actor-oriented modeling

// creating a namespace 'stepwiseActor' with all the routines for computing stepwise endogenous statistics (for actor-oriented models)
namespace stepwiseActor
{
        // stepwiseTie::inertia
    arma::mat inertia(arma::mat counts,
                    arma::mat intervals, 
                    arma::vec actor1,
                    arma::vec actor2,
                    arma::vec time, 
                    arma::uword M, 
                    int N,
                    arma::uword k, 
                    arma::uword K, 
                    bool senderRate,
                    int ncores)
        {   
            arma::uword m;
            arma::mat out(M, counts.n_cols-2, arma::fill::zeros);
            
            for(m = 1; m < M; m++){ 
            arma::uword lb_ub_row = intervals(((K*m)+k),2);
            out.row(m) = counts(lb_ub_row,arma::span(2,counts.n_cols-1)); //((K*m)+k)
            }
            return out;
        }

    // Mapping participation shifts functions to strings (which will be written as input in the .R function)
    std::unordered_map<std::string, std::function<arma::mat(arma::mat, arma::mat, arma::vec, arma::vec, arma::vec, arma::uword, int, arma::uword, arma::uword, bool, arma::uword)>>  endogenousMap =
            {
                {"inertia", inertia}
                // {"reciprocity", reciprocity},
                //  {"indegreeSender", indegreeSender},
                //  {"indegreeReceiver", indegreeReceiver},
                //  {"outdegreeSender", outdegreeSender},
                //  {"outdegreeReceiver", outdegreeReceiver},
                //  {"transitivityClosure", transitivityClosure},
                //   {"cyclicClosure", cyclicClosure}
    // missing stats : "totaldegreeSender", "totaldegreeReceiver", "sendingClosure", "receivingClosure", "pshift"
            };


    // computeIntervalStatistic (internal routine function)
    arma::mat computeStatistic(std::string statistic,
                            arma::mat counts, 
                            arma::mat intervals, // intervals_backward
                            arma::vec actor1,
                            arma::vec actor2,
                            arma::vec time, // t
                            arma::uword M, 
                            int N,
                            arma::uword k, 
                            arma::uword K, //K_q
                            bool senderRate,
                            int ncores)
        {
            arma::mat out = endogenousMap[statistic](counts, intervals, actor1, actor2, time, M, N, k, K, senderRate, ncores);
            return out;
        }
}


#endif