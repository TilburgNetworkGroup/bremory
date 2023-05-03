#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppDist.h>
#include <mvt.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include <iomanip>
#include "bremory.h"


//  BEGIN Preprocessing functions //





//' getExogenousVariablesArray
//'
//' @param exos vector of exogenous variables in the linear predictor (we avoid to process statistics that are not specified in the linear predictor but provided in 'data')
//' @param data list of exogenous variables (dyadic and actor variables), with actor/type names converted to ID's according to reh object
//' @param reh  reh object obtained from remify::reh(). It is one of the input arguments needed by bremory::smm()
//'
//' @return list of arrays of 3 dimensions: reh$M rows (number of events in the sequence), reh$D columns (number of dyads at risk in the network), reh$S slices (number of exogenous statistics). No interactions are accounted for at this stage.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List getExogenousVariablesArray(Rcpp::StringVector exos, 
                                        Rcpp::List data,
                                        Rcpp::List reh) 
{
    arma::uword S = exos.size(); // number of exogenous statistics specified in the linear predictor
    arma::uword M = reh["M"]; // numebr of events
    arma::uword D = reh["D"]; // number of dyads
    arma::uword N = reh["N"]; // number of actors
    arma::uword C = reh["C"]; // number of event types
    arma::uword s,n,c,i,j;
    Rcpp::List out(S);

    // actor-level statistics
    Rcpp::List actor_level = data["actor"];
    Rcpp::List actor_stats = actor_level["variables"]; // data structure with actor-level statistics
    arma::uvec actor_stats_tv = Rcpp::as<arma::uvec>(actor_level["tv"]);
    arma::uvec actors = Rcpp::as<arma::uvec>(actor_level["actorName"]); // data structure with actors' IDs
    Rcpp::StringVector actor_stats_names = Rcpp::as<Rcpp::StringVector>(actor_stats.names()); // names of actor-level statistics


    // dyad-level statistics
    Rcpp::List dyad_level = data["dyadic"];
    Rcpp::List dyad_stats = dyad_level["variables"]; // data structure with dyad-level statistics
    arma::uvec dyad_stats_tv = Rcpp::as<arma::uvec>(dyad_level["tv"]);   
    arma::uvec dyads = Rcpp::as<arma::uvec>(dyad_level["dyad"]); // data structure with dyads' IDs
    Rcpp::StringVector dyad_stats_names = Rcpp::as<Rcpp::StringVector>(dyad_stats.names()); // names of dyad-level statistics
    
    for(s = 0; s < S; s++){
        Rcpp::StringVector exo_s = Rcpp::as<Rcpp::StringVector>(exos[s]); // get variable name
        // if actor-level
        Rcpp::IntegerVector which_s = Rcpp::match(exo_s,actor_stats_names); // check match with actor variables
        arma::uword s_index = which_s[0]-1;
        if(which_s[0] > 0){ // if actor-level variable
           // Rcpp::Rcout << "elaborating exogenous statistic >>" << exo_s << "\n"; // message to console
            arma::cube stats_s = Rcpp::as<arma::cube>(actor_stats[s_index]); // cube of dimensions Nxchange_pointsxlevels
            if(actor_stats_tv(s_index) == 0){ // when the actor-level statistic is constant over time
                arma::cube out_s(M,D,stats_s.n_cols,arma::fill::zeros);
                for(i = 0; i < actors.n_elem; i++){ // for all the actors 
                    for(c = 0; c < C; c++){ // for all the event types
                        for(n = 0; n < N; n++){ // for all the actors conditioned to the actor being sender or receiver
                            if(actors(i) != n){
                                arma::uword s_i = remify::utils::getDyadIndex(actors(i),n,c,N,true); // when actor is a sender
                                arma::uword r_i = remify::utils::getDyadIndex(n,actors(i),c,N,true); // when actor is a receiver
                                for(j = 0; j < stats_s.n_cols; j++){
                                    arma::mat stats_s_mat_j(M,D,arma::fill::zeros);
                                    stats_s_mat_j.col(s_i) += stats_s(i,j,0);
                                    stats_s_mat_j.col(r_i) += stats_s(i,j,0);
                                    out_s.slice(j) += stats_s_mat_j;
                                }
                            }
                        }
                    }
                }
                out[s] = out_s;
            }
            else{
                // when the actor-level statistic is time-varying [not processed now]
                Rcpp::stop("time varying statistics are not supported in this version of the package");
            }
        }
        else{
        // if dyad-level
            Rcpp::IntegerVector which_s = Rcpp::match(exo_s,dyad_stats_names);
            arma::uword s_index = which_s[0]-1;
            if(which_s[0] > 0){
                Rcpp::Rcout << "elaborating exogenous statistic >>" << exo_s << "\n";
                arma::cube stats_s = Rcpp::as<arma::cube>(dyad_stats[s_index]); // cube of dimensions Nxchange_pointsxlevels
            if(dyad_stats_tv(s_index) == 0){ // when the dyad-level statistic is constant over time
                arma::cube out_s(M,D,stats_s.n_cols,arma::fill::zeros);
                for(i = 0; i < dyads.n_rows; i++){
                    // when actor is a sender
                    arma::uword d_i = dyads(i);
                    for(j = 0; j < stats_s.n_cols; j++){
                        arma::mat stats_s_mat_j(M,D,arma::fill::zeros);
                        stats_s_mat_j.col(d_i) += stats_s.at(i,j,0);
                        out_s.slice(j) += stats_s_mat_j; 
                    }          
                }
                out[s] = out_s;
            }
            else{
                // when the dyad-level statistic is time-varying [not processed now]
                Rcpp::stop("time varying statistics are not supported in this version of the package");
            }
            }
            else{
                Rcpp::stop("variable not found. ");
            }
        }
    }
    out.attr("names") = exos;

    return out;
}

//' getBinaryREH
//'
//' @param dyad vector of dyad occurred (reh$edgelist[,2])
//' @param D numbr of possible dyads (reh$D)
//'
//' @return matrix with dimensions given by the number of events (by row) and the dimension of the risk set (by column) .
//'
//' @export
// [[Rcpp::export]]
arma::mat getBinaryREH(arma::uvec dyad, arma::uword D) // output type was arma::sp_mat
{
    arma::uword m;
    arma::uword M = dyad.n_elem;
    arma::mat out(M,D,arma::fill::zeros); // arma::fill::zeros does not work for sp_mat
    for(m = 0; m < M; m++)
        {
            arma::uword dyad_m = dyad(m);
            out(m,dyad_m) = 1;
        }
    return out;
}
//  END Preprocessing functions //

//  BEGIN Statistics functions //

//' getIntervals
//'
//' given vectors of widths and time, the function returns.
//'
//' @param env  is the global environment
//' @param widths time lengths of intervals as to the stepwise model
//' @param time vector of time points
//' @param M number of events
//' @param K_q number of steps for the q-th model
//' @param ncores number of ncores
//'
//' @return matrix.
//'
//' @export
// [[Rcpp::export]]
void getIntervals(Rcpp::Environment env,
                       arma::vec widths, 
                       arma::vec time, 
                       arma::uword M, 
                       arma::uword K_q,
                       int ncores)
{
    Rcpp::Environment statisticsREH = env["statisticsREH"];
    arma::uword i,k;

    // creating an empty matrix of (-1) (allocating memory)
    arma::mat intervals(M*K_q, 3);
    intervals.fill(-1);   
    	
    omp_set_dynamic(0);           // disabling dynamic teams
    omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
    #pragma omp parallel for private(i,k) shared(intervals) 
    // begin for loop
    for (i = 1; i < M; i++) 
    {   
        // Updating intervals to the i-th iteration:
        double t = time(i); //' current time
        arma::vec time_vec = time(arma::span(0,(i-1))); 
        
        for (k = 0; k < K_q; k++)
        {
            // Finding bounds j-th interval
            arma::vec lbs_ubs(2);
            arma::vec bounds_intervals(2,arma::fill::zeros);

            lbs_ubs(0) = t - widths(k+1);
            lbs_ubs(1) = t - widths(k);

            // Finding lower bounds candidate vector
            arma::uvec lb_temp_vec = arma::find(time_vec >= lbs_ubs(0));
            // Finding upper bounds candidate vector
            arma::uvec ub_temp_vec = arma::find(time_vec < lbs_ubs(1));

            if((lb_temp_vec.n_elem != 0) & (ub_temp_vec.n_elem != 0))
            {
                bounds_intervals(0) = min(lb_temp_vec); // lb_temp_vec(0);

                if(k>0){
                    bounds_intervals(1) = max(ub_temp_vec); //ub_temp_vec(ub_temp_vec.n_elem-1);
                }
                else{
                    bounds_intervals(1) = i-1;
                }

                if(bounds_intervals(0) <= bounds_intervals(1))
                {
                intervals(k+i*K_q,arma::span(0,1)) = bounds_intervals.t();
                }
            }

        }
    }
    statisticsREH.assign("intervals",intervals);
}

//' getCountsOMP_old
//'
//' @param binaryREH matrix
//' @param lbs_ubs matrix of lower bounds (lb's, first column) and upper bounds (ub's, second column) of binaryREH
//' @param ncores ncores
//'
//' @return matrix of accumulated counts within ranges of binaryREH
//'
//' @export
// [[Rcpp::export]]
arma::mat getCountsOMP_old(arma::mat binaryREH, arma::mat lbs_ubs, int ncores)  // arma::vec dyad, arma::uword D
{
    arma::uword n_dyads = binaryREH.n_cols;
    arma::uword n_intervals = lbs_ubs.n_rows;
    arma::uword i;
    arma::mat output_counts(n_intervals,n_dyads+2,arma::fill::zeros);
    output_counts.col(0).fill(-1);
    output_counts.col(1).fill(-1);   
    arma::rowvec count_loc(n_dyads,arma::fill::zeros);

    omp_set_dynamic(0);           // disabling dynamic teams
    omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
    #pragma omp parallel for private(i,count_loc) shared(n_dyads,n_intervals,lbs_ubs,output_counts)
        for(i = 0; i < n_intervals; i++){
            if(lbs_ubs(i,0) != (-1)){
                if(lbs_ubs(i,0) == lbs_ubs(i,1)){
                    count_loc = binaryREH.row(lbs_ubs(i,0));
                    
                }
                else{
                    count_loc = arma::sum(binaryREH(arma::span(lbs_ubs(i,0),lbs_ubs(i,1)),arma::span::all),0);
                }
                output_counts(i,0)=lbs_ubs(i,0);
                output_counts(i,1)=lbs_ubs(i,1);
                output_counts(i,arma::span(2,n_dyads+1)) = count_loc;
            }
        }
    return output_counts;
}


//' getCounts
//'
//' @param dyad is the vector of observed dyads 
//' @param D is the number of dyads
//' @param lbs_ubs matrix of lower bounds (lb's, first column) and upper bounds (ub's, second column) of binaryREH
//' @param ncores ncores
//'
//' @return matrix of accumulated counts within ranges of binaryREH
//'
//' @export
// [[Rcpp::export]]
arma::mat getCounts(arma::vec dyad, arma::uword D, arma::mat lbs_ubs, int ncores) 
{
    arma::uword n_intervals = lbs_ubs.n_rows;
    arma::uword i;
    arma::mat counts(n_intervals,D+2,arma::fill::zeros); // first two columns are lower bound and upper bound of the interval, then D columns are for the dyads
    counts.col(0).fill(-1);
    counts.col(1).fill(-1);   
    arma::rowvec count_loc(D,arma::fill::zeros);
    arma::uword which_dyad_loc, d_obs;

    omp_set_dynamic(0);           // disabling dynamic teams
    omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
    #pragma omp parallel for private(i,count_loc,which_dyad_loc,d_obs) shared(D,n_intervals,lbs_ubs,counts)
        for(i = 0; i < n_intervals; i++){
            if(lbs_ubs(i,0) != (-1)){
                if(lbs_ubs(i,0) == lbs_ubs(i,1)){
                    arma::rowvec row_loc(D,arma::fill::zeros);
                    which_dyad_loc = dyad(lbs_ubs(i,0));
                    row_loc(which_dyad_loc) += 1;
                    count_loc = row_loc;      
                }
                else{
                    arma::rowvec row_loc(D,arma::fill::zeros);
                    arma::vec dyad_loc = dyad(arma::span(lbs_ubs(i,0),lbs_ubs(i,1)));
                    for(d_obs = 0; d_obs < dyad_loc.n_elem; d_obs++){
                        row_loc(dyad_loc(d_obs)) += 1.0;
                    }
                    count_loc = row_loc;
                }
                counts(i,0)=lbs_ubs(i,0);
                counts(i,1)=lbs_ubs(i,1);
                counts(i,arma::span(2,D+1)) = count_loc;
            }
        }
    return counts;
}


//' getCountsIndex (function useful for statistics computation)
//'
//' @param intervals matrix of intervals (non-unique intervals, first two columns are lower bound and upper bound respectively)
//' @param counts matrix of intervals (first two columns are lower bound and upper bound respectively)
//' @param ncores number of ncores
//'
//' @return vector of row indices relative to the matrix 'counts'
//'
//' @export
// [[Rcpp::export]]
arma::vec getCountsIndex(arma::mat intervals, arma::mat counts, int ncores)
{
    arma::uword i;
    arma::vec counts_index(intervals.n_rows);

    omp_set_dynamic(0);           // disabling dynamic teams
    omp_set_num_threads(ncores); // number of ncores for all consecutive parallel regions
    #pragma omp parallel for private(i) shared(counts_index) 
    for(i = 0; i < intervals.n_rows; i++){
        arma::uvec loc_find = arma::intersect(arma::find(counts.col(0)==intervals(i,0)),arma::find(counts.col(1)==intervals(i,1)));
        counts_index[i] = counts(loc_find(0),2);
    }
    return counts_index;
}


//' getIntervalStatistics
//'
//' @param stats is a list of statistics (exogenous and endogenous) which will be updated with the endogenous statistics
//' @param counts matrix of dyad counts for unique intervals
//' @param intervals matrix with all intervals
//' @param actor1 time-ordered sequence of actor1 ID's
//' @param actor2 time-ordered sequence of actor2 ID's
//' @param time vector of time
//' @param M number of events
//' @param N number of actors
//' @param D number of dyads (no event type)
//' @param K number of intervals (q-th model)
//' @param statistics string vector of endogenous statistics
//' @param model 'tie' or 'actor'
//' @param senderRate true/false for actor-oriented modeling. if true, then statistics for sender-rate model, otherwise statistics for the receiver-choice model 
//' @param ncores number of ncores to create in case of parallelization
//'
//' @return list of (interval) statistics
//'
// [[Rcpp::export]]
void getIntervalStatistics(Rcpp::Environment env,
                                arma::mat counts,
                                arma::mat intervals,
                                const arma::vec &actor1,
                                const arma::vec &actor2,
                                arma::vec time,
                                arma::uword M,
                                int N,
                                arma::uword D,
                                arma::uword K,
                                std::vector<std::string> statistics,
                                std::string model,
                                bool senderRate = false,
                                int ncores = 1)
{
    Rcpp::Environment stats = env["stats"];
    arma::uword P = statistics.size();
    arma::uword p,k;
    std::vector<std::string> which_model = {"tie","actor"};

    if(model.compare(which_model[0]) == 0){ // tie-oriented modeling
        for(p = 0; p < P; p++){ // p iterator for the endogenous statistics
            arma::cube stats_p(M,D,K);
            for(k = 0; k < K; k++){ // k iterator for the intervals
                arma::mat stats_p_k = stepwiseTie::computeStatistic(statistics[p], counts, intervals, actor1, actor2, time, M, N, k, K, ncores);
                stats_p.slice(k) = stats_p_k;
            }
            stats.assign(statistics[p],stats_p);
        }
    }
    else if(model.compare(which_model[1]) == 0){
        for(p = 0; p < P; p++){ // p iterator for the endogenous statistics
            arma::cube stats_p(M,D,K);
            for(k = 0; k < K; k++){ // k iterator for the intervals
                arma::mat stats_p_k = stepwiseActor::computeStatistic(statistics[p], counts, intervals, actor1, actor2, time, M, N, k, K, senderRate, ncores);
                stats_p.slice(k) = stats_p_k;
            }
            stats.assign(statistics[p],stats_p);
        }
    }  
}

//' arrangeStatistics
//'
//' @param env global environment where the environment of statistics is ('stats')
//' @param effectsMatrix is a matrix of effects from the processed formula object
//' @param names is the vector of row names from the matrix of effects (effectsMatrix)
//' @param intercept is a TRUE/FALSE value whether the intercept is specified or not
//' @param M is the number of events
//' @param D is the number of dyads (without event type)
//' @param S is the final number of predictors in the model
//'
//' @return three dimensional array of statistics and vector of names (slices)
//'
// [[Rcpp::export]]
arma::cube arrangeStatistics(Rcpp::Environment env, arma::mat effectsMatrix, std::vector<std::string> names, bool intercept, int M, int D, int S){

    int L = effectsMatrix.n_cols; // number of columns of effectsMatrix
    int s = 0; //s is the slice counter
    int l; 
    Rcpp::Environment stats = env["stats"];
    arma::cube statistics(M,D,S); // S is the number of statistcs ('statisticsREH$n_pars_q')

    if(intercept){
        statistics.slice(s) = Rcpp::as<arma::cube>(stats["intercept"]);
        s += 1;
    }
    for(l = 0; l < L; l++){ // for each column of the effectsMatrix
        arma::vec statistic_l = effectsMatrix.col(l);
        arma::uvec which_statistics_l = arma::find(statistic_l); // by default it returns the indices of all the non-zero elements
        if(which_statistics_l.n_elem > 1){
            // interaction effect: first variable
            arma::uword p0 = which_statistics_l(0);
            std::string statistic_name_p0 = names[p0];
            arma::cube stats_p0 =  Rcpp::as<arma::cube>(stats[statistic_name_p0]);
            // interaction effect: second variable
            arma::uword p1 = which_statistics_l(1);
            std::string statistic_name_p1 = names[p1];
            arma::cube stats_p1 =  Rcpp::as<arma::cube>(stats[statistic_name_p1]);
            for(int j = 0; j < stats_p0.n_slices; j++){
                for(int z = 0; z < stats_p1.n_slices; z++){
                    statistics.slice(s) = stats_p0.slice(j) % stats_p1.slice(z); // element-wise multiplication
                    s += 1;
                }
            }
            
        }
        else{
            arma::uword p = which_statistics_l(0);
            std::string statistic_name = names[p];
            arma::cube stats_p =  Rcpp::as<arma::cube>(stats[statistic_name]);
            for(int j = 0; j < stats_p.n_slices; j++){
                statistics.slice(s) = stats_p.slice(j);
                s += 1;
            }
        }
    }
    return statistics;
}



//' getStatisticsNames
//'
//' @param stats_names is a (names) list of vectors of names of statistics
//' @param effectsMatrix is a matrix of effects from the processed formula object
//' @param names is the vector of row names from the matrix of effects (effectsMatrix)
//' @param intercept is a TRUE/FALSE value whether the intercept is specified or not
//' @param M is the number of events
//' @param D is the number of dyads (without event type)
//' @param S is the final number of predictors in the model
//'
//' @return three dimensional array of statistics and vector of names (slices)
//'
//' @export
// [[Rcpp::export]]
std::vector<std::string> getStatisticsNames(Rcpp::List stats_names, arma::mat effectsMatrix, std::vector<std::string> names, bool intercept, int M, int D, int S){

    int L = effectsMatrix.n_cols; // number of columns of effectsMatrix
    int s = 0; //s is the slice counter
    int l; 

    std::vector<std::string> output_names(S);
    if(intercept){
        output_names[s] = "intercept";
        s += 1;
    }
    for(l = 0; l < L; l++){ // for each column of the effectsMatrix
        arma::vec statistic_l = effectsMatrix.col(l);
        arma::uvec which_statistics_l = arma::find(statistic_l); // by default it returns the indices of all the non-zero elements
        if(which_statistics_l.n_elem > 1){
            // interaction effect: first variable
            arma::uword p0 = which_statistics_l(0);
            std::string statistic_name_p0 = names[p0];
            std::vector<std::string> stats_names_p0 = Rcpp::as<std::vector<std::string>>(stats_names[statistic_name_p0]);

            // interaction effect: second variable
            arma::uword p1 = which_statistics_l(1);
            std::string statistic_name_p1 = names[p1];
            std::vector<std::string> stats_names_p1 = Rcpp::as<std::vector<std::string>>(stats_names[statistic_name_p1]);
            for(int j = 0; j < stats_names_p0.size(); j++){
                for( int z = 0; z < stats_names_p1.size(); z++){
                    output_names[s] = stats_names_p0[j] + ":" + stats_names_p1[z];
                    s += 1;
                }
            }
            
        }
        else{
            arma::uword p = which_statistics_l(0);
            std::string statistic_name = names[p];
            std::vector<std::string> stats_names_p = Rcpp::as<std::vector<std::string>>(stats_names[statistic_name]);
            for(int j = 0; j < stats_names_p.size(); j++){
                output_names[s] = stats_names_p[j];
                s += 1;
            }
        }
    }
   
    return output_names;

}

//' lpd (Log-Pointwise Density of REM)
//'
//' @param pars is a vector of parameters (note: the order must be aligned witht the column order in 'stats')
//' @param stats is a matrix of dimensions n_dyads*variables with statistics of interest by column and dyads by row.
//' @param event is a vector of 1/0 : 1 indicating the observed dyad and 0 the non observed dyads.
//' @param interevent_time the time difference between the current time point and the previous event time.
//'
//' @return log-pointwise density value of a specific time point
//'
//' @export
// [[Rcpp::export]]
double lpd(arma::vec pars, arma::mat stats, arma::uvec event, double interevent_time){
        arma::uword n_dyads = event.n_elem;
        arma::uword i;
        arma::vec log_lambda = stats * pars;
        double lpd = 0.0;
        for(i = 0; i < n_dyads; i++){
            if(event(i) == 0){
                lpd -= (exp(log_lambda(i))*interevent_time);
            }
            else{
                lpd += (log_lambda(i)-(exp(log_lambda(i))*interevent_time));
            }
        }
        return lpd;
    }


//' getWAIC
//'
//' @param pars parameters values
//' @param stats sub-array of statistics
//' @param events sub-matirx of events
//' @param interevent_time sub-vector of interevent times
//'
//' @return matrix of p_m and lpd_m values for a specific stepwise model
//'
//' @export
// [[Rcpp::export]]
arma::mat getWAIC(arma::mat pars, arma::cube stats, arma::umat events, arma::vec interevent_time){
        arma::uword m,j;
        arma::vec out_vec_loc_log(pars.n_cols,arma::fill::zeros);
        arma::vec out_vec_loc_p(pars.n_cols,arma::fill::zeros);
        arma::mat out(interevent_time.n_elem,2,arma::fill::zeros);
        double lpd_m_j,p_m;
        for(m = 0; m < interevent_time.n_elem; m++)
            {
                for(j = 0; j < pars.n_cols; j++){
                    lpd_m_j = lpd(pars.col(j),stats.slice(m),events.col(m),interevent_time(m));
                    out_vec_loc_log[j] = lpd_m_j;
                    out_vec_loc_p[j] = exp(lpd_m_j);
                }
                p_m = mean(out_vec_loc_p);
                out(m,0) = log(p_m);
                out(m,1) = var(out_vec_loc_log);
            }

        return out;
    }




//' smoothing posterior.
//'
//' This function calculates the posterior distribution give...
//'
//' @param sample_models sample_models
//' @param knots_seq knots_seq
//' @param which_pars which_pars
//' @param post_betas post_betas
//' @param post_gammas post_gammas
//' @param U n_pars
//'
//' @return a cube of posterior draws.
//'
//' @export
// [[Rcpp::export]]
arma::cube smoothing_posterior(arma::uvec sample_models,
                              arma::rowvec knots_seq,    
                              arma::ucube which_pars, 
                              arma::mat post_betas, 
                              arma::mat post_gammas,
                              arma::uword U)
{
    arma::uword i,j,model_j;
    arma::uword n_knots = knots_seq.n_elem;
    arma::uword nsim = post_betas.n_rows;
    arma::uvec beta_indices_loc(U,arma::fill::zeros);
    arma::uvec row_position(1,arma::fill::zeros);
    arma::umat beta_indices(nsim,U,arma::fill::zeros);   
    arma::cube out(nsim, U, n_knots); 
    out.fill(9999);
    for (i = 0; i < n_knots; i++)
    {
        beta_indices = which_pars.slice(i);
        arma::mat post_betas_loc(nsim,U);
        post_betas_loc.fill(9999);
        for (j = 0; j < nsim; j++)
        {
            model_j = sample_models(j);
            beta_indices_loc = beta_indices.row(model_j).t();
            if(beta_indices_loc(0)!=9999){   // condition is only on the first because by default when ==9999 also the other elements will be 
                row_position[0] = j;
                post_betas_loc.row(j) = post_betas(row_position,beta_indices_loc);
            }
        }
        out.slice(i) = post_betas_loc;
    }
    return out;
}

//' getDraws
//'
//' Function to draw from the posterior distribution (here the second step of drawing from Normal distribution is performed;
//' the first consists in drawing from the multinomial distribution for the Q models)
//'
//' @param sample_models models sample with probabilities given by the specific weighting system used (BIC, pseudoBMA+, stacking, WAIC)
//' @param which_pars which parameters
//' @param n_pars number of parameters per model in the set
//' @param n_stats number of statistics
//' @param input list() of K (number of intervals), gammas (time widths), mle_betas (Maximum Likelihood Estimates), vcov_betas (variance and covariance matrix)
//' @param knots_seq sequence of time widths at which compute the posterior distribution
//'
//' @return matrix.
//'
//' @export
// [[Rcpp::export]]
arma::cube getDraws(arma::uvec sample_models,
                    arma::ucube which_pars,
                    arma::vec n_pars,
                    arma::uword n_stats,
                    Rcpp::List input,
                    arma::rowvec knots_seq) 
{
    // defining useful variables :
    arma::uword i, sample_i, n_pars_i; 
    //arma::uword K_i; //K-i and related objects should be removed
    arma::uword nsim = sample_models.n_elem;
    // arranging list elements of 'configurations' in vectors, matrices and cube :
    Rcpp::List intervals  = input["intervals"];
    arma::vec K = intervals["K"];
    arma::mat widths = intervals["widths"];
    arma::mat mle_betas = input["coef"];
    arma::cube vcov_betas = input["vcov"];

    // allocating memory (output) :
    arma::mat out_widths(nsim, widths.n_cols, arma::fill::ones);
    arma::mat out_betas(nsim, mle_betas.n_cols, arma::fill::ones);
    //arma::vec out_K(nsim, arma::fill::zeros);
    for(i = 0; i < nsim; i++)
    {
        // select i-th sampled model
        sample_i = sample_models(i);
        // set number of intervals (K)
        //K_i = K(sample_i);
        n_pars_i = n_pars(sample_i);

        //out_K[i] = K_i;

        // pick gamma's
        out_widths.row(i) = widths(sample_i,arma::span::all);

        // generate beta's
        arma::vec means_i = mle_betas(sample_i,arma::span(0,n_pars_i-1)).t(); 
        arma::mat vcov_all = vcov_betas.slice(sample_i);
        arma::mat vcov_i =  vcov_all(arma::span(0,n_pars_i-1),arma::span(0,n_pars_i-1));
        arma::vec betas_sampled = arma::mvnrnd(means_i,vcov_i);

        out_betas(i,arma::span(0,n_pars_i-1)) = betas_sampled.t();
    }

    arma::cube out = smoothing_posterior(sample_models,knots_seq,which_pars,out_betas,out_widths,n_stats);

    return out;
}


//  END Statistics functions //


//' tryClone (text here)
//'
//' @param input description of input here
//'
//' @return input [dataframe]
//'
//' @export
// [[Rcpp::export]]
double tryClone(double input){
    
double f =3.14159;
  std::cout << std::setprecision(5) << f << '\n';
  std::cout << std::setprecision(9) << f << '\n';
  std::cout << std::fixed;
  std::cout << std::setprecision(5) << f << '\n';
  std::cout << std::setprecision(9) << f << '\n';
  return 0.0;
}