#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppDist.h>
#include <mvt.h>
#include <iostream>
#include <omp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <typeinfo>
#include <map>
#include <iterator>
#include <string>
#include "stepwise_estimation.h"
#include "smooth_estimation.h"
#include "decay_functions.h"



#define LOG(x) std::cout << x << "\n"

//  BEGIN Preprocessing functions //


//' getRiskset_old (obtain permutations of actors' id). (will become getRisksetMatrix) [[ TO REMOVE ]]
//'
//' @param actors vector of actors' id's (either integer or char vector).
//' @param N number of actors in the network.
//' @param selfedges if selfedges are allowed (TRUE) or not (FALSE).
//'
//' @return matrix of possible dyadic events.
//'
//' @export
// [[Rcpp::export]]
arma::mat getRiskset_old(arma::vec actors_id, arma::uword N, bool selfedges){
    arma::uword i,j;
    if(selfedges == 0){
        arma::mat riskset(N*N,2);
        for(i = 0; i < N; i++){
            for(j = 0; j < N ; j++){
                    if(j != i){
                riskset(j+i*N,0) = actors_id(i);
                riskset(j+i*N,1) = actors_id(j);
                }
                else {
                    riskset(j+i*N,0) = NAN;
                    riskset(j+i*N,1) = NAN;
                }
            }
        }
        return riskset; // in R na.omit() function will be applied to this ouptut so as to filter out the rows (NaN,NaN)
    }
    else
    {
    arma::mat riskset((N*N),2);
    for(i = 0; i < N; i++){
        for(j = 0; j < N ; j++){
                riskset(j+i*N,0) = actors_id(i);
                riskset(j+i*N,1) = actors_id(j);
        }
    }
    return riskset;
    }

}


//' getRisksetMatrix_old (will become getRisksetCube) [[ TO REMOVE ]]
//'
//' @param riskset matrix (sender and receiver integer id's by column)
//' @param N number of actors in the network.
//'
//' @return matrix of possible dyadic events arranged in a square matrx (cell value is the row index of the riskset matrix input, which will be the same to the column index in binary REH).
//'
//' @export
// [[Rcpp::export]]
arma::umat getRisksetMatrix_old(arma::mat riskset , arma::uword N) {
    arma::uword i;
    arma::umat riskset_mat(N,N);
    riskset_mat.fill(N*N*N);
    for(i = 0; i < riskset.n_rows; i++){
        // riskset_mat(sender,receiver) = row index in riskset (that is also the column index in binaryREH)
            riskset_mat((riskset(i,0)),(riskset(i,1))) = i;
        }
    return riskset_mat; 
}


//' getBinaryREH()
//'
//'
//' @param M number of relational events observed
//' @param N number of actors
//' @param edgelist matrix of (time,sender,receiver)
//' @param riskset_matrix matrix of possible dyadic events with column_order per each cell (this will be the order in the binaryREH matrices columns)
//'
//' @return Matrix with dimensions given by the number of events (by row) and the dimension of the risk set (by column) .
//'
//' @export
// [[Rcpp::export]]
arma::mat getBinaryREH(arma::uword M,
                        arma::uword N,
                        arma::mat edgelist, 
                        arma::umat riskset_matrix)
{
    arma::uword m;
    arma::mat out(M,N*(N-1),arma::fill::zeros); 
    arma::vec edgelist_loc(3,arma::fill::zeros);

    for(m = 0; m < M; m++)
        {
            edgelist_loc = edgelist.row(m).t();
            arma::rowvec out_loc(N*(N-1),arma::fill::zeros);
            out_loc[riskset_matrix(edgelist_loc(1),edgelist_loc(2))] = 1;
            out.row(m) = out_loc;
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
//'
//' @return matrix.
//'
//' @export
// [[Rcpp::export]]
void getIntervals(Rcpp::Environment env,
                       arma::vec widths)
{
    Rcpp::Environment initializeREH = env["initializeREH"];
    Rcpp::Environment statisticsREH = env["statisticsREH"];
    arma::vec time = initializeREH["t"];
    arma::uword M = initializeREH["M"];
    arma::uword K_q = statisticsREH["K_q"];
    arma::uword i,k;
    // creating an empty matrix of (-1) (allocating memory)
    arma::mat intervals_backward(M*K_q, 3);
    intervals_backward.fill(-1);   

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
                intervals_backward(k+i*K_q,arma::span(0,1)) = bounds_intervals.t();
                }
            }

        }
    }
    statisticsREH.assign("intervals_backward",intervals_backward);
}


//' getCountsOMP
//'
//' @param binaryREH matrix
//' @param lbs_ubs matrix of lower bounds (lb's, first column) and upper bounds (ub's, second column) of binaryREH
//' @param n_cores n_cores
//'
//' @return matrix of accumulated counts within ranges of binaryREH
//'
//' @export
// [[Rcpp::export]]
arma::mat getCountsOMP(arma::mat binaryREH, arma::mat lbs_ubs, int n_cores)
{
    arma::uword n_dyads = binaryREH.n_cols;
    arma::uword n_intervals = lbs_ubs.n_rows;
    arma::uword i;
    arma::mat output_counts(n_intervals,n_dyads+2,arma::fill::zeros);
    output_counts.col(0).fill(-1);
    output_counts.col(1).fill(-1);   
    arma::rowvec count_loc(n_dyads,arma::fill::zeros);

    omp_set_dynamic(0);           // disabling dynamic teams
    omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
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


//' getCountsIndex (function useful for statistics computation)
//'
//' @param intervals matrix of intervals (non-unique intervals, first two columns are lower bound and upper bound respectively)
//' @param counts matrix of intervals (first two columns are lower bound and upper bound respectively)
//'
//' @return vector of row indices relative to the matrix 'counts'
//'
//' @export
// [[Rcpp::export]]
arma::vec getCountsIndex(arma::mat intervals, arma::mat counts)
{
    arma::uword i;
    arma::vec counts_index(intervals.n_rows);
    for(i = 0; i < intervals.n_rows; i++){
        arma::uvec loc_find = arma::intersect(arma::find(counts.col(0)==intervals(i,0)),arma::find(counts.col(1)==intervals(i,1)));
        counts_index[i] = counts(loc_find(0),2);
    }
    return counts_index;
}


//' getEndoEffects
//'
//' @param env environment where the user is currently working (usually the global one which can be accessed via 'globalenv()' function or '.GlobalEnv' object)
//'
//' @return array of Statistics specified in the
//'
//' @export
// [[Rcpp::export]]
arma::cube getEndoEffects(Rcpp::Environment env)
{
    Rcpp::Environment initializeREH = env["initializeREH"];
    Rcpp::Environment statisticsREH = env["statisticsREH"];
    Rcpp::Environment counts = statisticsREH["counts"];
    arma::uword M = initializeREH["M"];
    arma::uword K_q = statisticsREH["K_q"];
    arma::uword P = statisticsREH["P"];
    arma::uword p,k;
    arma::uword n_dyads = initializeREH["n_dyads"];
    arma::uword n_cores = initializeREH["n_cores"];
    arma::vec t = initializeREH["t"];
    std::vector<std::string> endo_effects = statisticsREH["endo_effects"];
    arma::mat edgelist = initializeREH["edgelist"];
    arma::umat riskset = initializeREH["riskset"];
    arma::umat riskset_matrix = initializeREH["riskset_matrix"]; // i should re-arrange the riskset in a matrix (senderxreceiver and the cell has to give the column index in the statistics matrix)
    arma::mat intervals_backward = statisticsREH["intervals_backward"];

    // allocating memory for the output array
    arma::cube stats_array(M,n_dyads,P*K_q); 

    for(p = 0; p < P; p++){
        for(k = 0; k < K_q; k++){
            arma::mat out = computeEffect(endo_effects[p], counts, intervals_backward, edgelist, riskset, riskset_matrix, t, M, k, K_q, n_cores);
            stats_array.slice((p*K_q)+k) = out;
        }
    }
    
    return stats_array;
}

//' getSmoothEndoEffects
//'
//' @param env environment where the user is currently working (usually the global one which can be accessed via 'globalenv()' function or '.GlobalEnv' object)
//' @param endo_effects string vector indicating the endogenous effects
//' @param decay_fun function for the decay of past events influence
//' @param pars matrix with parameters for the decay_fun
//' @param n_cores number of corse to be used in the parallelization
//'
//' @return array of Statistics specified 
//'
//' @export
// [[Rcpp::export]]
arma::cube getSmoothEndoEffects(Rcpp::List reh, std::vector<std::string> endo_effects, Rcpp::List endo_memory_pars, arma::uword n_cores)
{   
    arma::uword M = reh["M"];
    arma::uword N = reh["N"];
    arma::uword C = reh["C"];
    arma::uword D = reh["D"];
    //Rcpp::DataFrame edgelist = reh["edgelist"];
    arma::umat risksetMatrix = reh["risksetMatrix"];
    arma::ucube risksetCube = reh["risksetCube"];
    arma::mat rehBinary = reh["rehBinary"];
    //Rcpp::RObject 
    Rcpp::DataFrame edgelist_loc = Rcpp::as<Rcpp::DataFrame>(reh["edgelist"]);
    arma::uvec actor1 = Rcpp::as<arma::uvec>(edgelist_loc["actor1"]);
    arma::uvec actor2 = Rcpp::as<arma::uvec>(edgelist_loc["actor2"]);
    arma::uvec type = Rcpp::as<arma::uvec>(edgelist_loc["type"]);
    arma::umat edgelist(M,3);
    edgelist.col(0) = actor1;
    edgelist.col(1) = actor2;
    edgelist.col(2) = type;

    arma::uword P = endo_effects.size(); // check if .size() works with std::vector<std::string>
    arma::vec t = edgelist_loc["time"];
    arma::uword p;

    // allocating memory for the output array
    arma::cube stats_array(M,D,P); 

    for(p = 0; p < P; p++){
        arma::mat out = computeSmoothEffect(endo_effects[p], endo_memory_pars[p], edgelist, rehBinary, risksetMatrix, risksetCube, t, M, N, C, D, n_cores);
        stats_array.slice(p) = out;
        //Rcpp::Rcout << "computation of " << endo_effects[p] << " successfully completed!\n"; // printing out some comments (will be removed)
    }
    
    return stats_array;
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


//' nllik (negative log-likelihood rem model)
//'
//' @param pars is a vector of parameters (note: the order must be aligned witht the column order in 'stats')
//' @param stats is a cube of dimensions n_dyads*variables*M with statistics of interest by column and dyads by row.
//' @param event_binary is a matrix of ones and zeros of dimensions M*n_dyads : 1 indicating the observed dyad and 0 the non observed dyads.
//' @param interevent_time the vector of time differences between the current time point and the previous event time.
//' @param n_cores n_cores
//'
//' @return log-pointwise density value of a specific time point
//'
//' @export
// [[Rcpp::export]]  
double nllik(arma::vec pars, arma::cube stats, arma::umat event_binary, arma::vec interevent_time, int n_cores){
        arma::uword n_dyads = event_binary.n_cols;
        arma::uword i,m;
        arma::uword M = event_binary.n_rows;
        arma::vec log_lambda(n_dyads,arma::fill::zeros) ;
        arma::vec llik(M,arma::fill::zeros);

        omp_set_dynamic(0);           // disabling dynamic teams
        omp_set_num_threads(n_cores); // number of threads for all consecutive parallel regions
        #pragma omp parallel for private(m,i,log_lambda) shared(n_dyads,M,stats,event_binary,interevent_time,llik)
        for(m = 0; m < M; m++)
        {
            log_lambda = stats.slice(m) * pars;
            for(i = 0; i < n_dyads; i++){
                if(event_binary(m,i) == 0){
                    llik(m) -= exp(log_lambda(i))*interevent_time(m);
                }
                else{
                    llik(m) += log_lambda(i)-exp(log_lambda(i))*interevent_time(m);
                }
            }
        }
        return -sum(llik);

    }

//' performBSIR
//' 
//' A function that performs the BSIR on the REM model
//'
//' @param nsim number of simulations
//' @param mean vector of model MLEs
//' @param sigma matric of variances and covariances 
//' @param df degrees of freedom for the multivariate Student t (used as proposal distribution)
//' @param stats cube of statistics [D*P*M]
//' @param event_binary matrix of 1/0
//' @param interevent_time vector of interevent times
//' @param n_cores number of cores to be used in the parallelized calculation of the nllik()
//'
//' @return log-pointwise density value of a specific time point
//'
//' @export
// [[Rcpp::export]]  
Rcpp::List performBSIR(arma::uword nsim,
                       arma::vec mean, 
                       arma::mat sigma, 
                       double df,
                       arma::cube stats,
                       arma::umat event_binary, 
                       arma::vec interevent_time,
                       arma::uword n_cores){
    // (0) create output empty list
    Rcpp::List out = Rcpp::List::create();
    arma::uword i;
    arma::vec density_posterior(nsim,arma::fill::ones);

    // (1) generate from a multivariate Student t and save both draw and density value
    arma::mat random_t = rmvt(nsim, mean,sigma, df);
    out["draws"] = random_t;

    arma::vec density_t = dmvt(random_t, mean, sigma, df, false);
    out["densities"] = density_t;

    // (2) evaluate the generated value with nnlik and get the density
    for(i = 0; i < nsim; i++){
        arma::vec draw_loc = random_t.row(i).t();
        density_posterior(i) = nllik(draw_loc,stats,event_binary,interevent_time,n_cores);
    }
    out["densities_post"] = density_posterior;

    return out;
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

//' lpdWAIC
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
arma::mat lpdWAIC(arma::mat pars, arma::cube stats, arma::umat events, arma::vec interevent_time){
        arma::uword m,j;
        arma::mat out(pars.n_cols,interevent_time.n_elem,arma::fill::zeros);
        for(m = 0; m < interevent_time.n_elem; m++)
            {
                for(j = 0; j < pars.n_cols; j++){
                    out(j,m) = lpd(pars.col(j),stats.slice(m),events.col(m),interevent_time(m));
                }
            }

        return out;
    }

//' logpJi 
//'
//' @param pars is a matrix of parameters (note: the order must be aligned witht the column order in 'stats')
//' @param stats is a matrix of dimensions n_dyads*variables with statistics of interest by column and dyads by row.
//' @param event is a vector of 1/0 : 1 indicating the observed dyad and 0 the non observed dyads.
//' @param interevent_time the time difference between the current time point and the previous event time.
//'
//' @return vector of log-pointwise density values of a specific time point at different parameters values
//'
//' @export
// [[Rcpp::export]]
arma::vec logpJi(arma::mat pars, arma::mat stats, arma::uvec event, double interevent_time){
    arma::uword s;
    arma::vec out(pars.n_cols,arma::fill::zeros);
    for(s=0; s<pars.n_cols; s++){
        out[s] = lpd(pars.col(s), stats, event, interevent_time);
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
//' @param n_pars n_pars
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
    arma::vec K = input["K"];
    arma::mat widths = input["widths"];
    arma::mat mle_betas = input["mle_betas"];
    arma::cube vcov_betas = input["vcov_betas"];

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


//' tryClone
//'
//' @param input data.frame object
//'
//' @return input [dataframe]
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame tryClone(Rcpp::DataFrame input){
    return input;
}