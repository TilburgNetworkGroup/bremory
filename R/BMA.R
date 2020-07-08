# Functions to perform Bayesian Model Averaging #

#' getPosterior
#'
#' A function which return a list of weights according to different possible weighting systems
#'
#' @param env environment stepwiseModelsREH where models are stores
#' @param weights which weighting system to use
#' @param method estimation method to use in the bootstrap
#' @param n_knots number of knots of the elapsed time where to evaluate the posterior distribution
#' @param nsim number of posterior draws
#' @param B number of splines
#' @param robustify TRUE/FALSE whether to get more robust estimate by emulating R times a BMA 
#' @param R number of BMAs replications, when the argument robustify is TRUE
#' @param subset_K if user wants specifically to subset the analysis to specific interval dimensions
#' @param sample_size_per_K number of models to sample per K (number of intervals), when robustify is TRUE: stratified sampling
#' @param plot TRUE/FALSE plotting results
#' @param top_K_models find the top K models for each weighting system and return: indices, widths and mles (inside a list)
#'
#' @return  list of if robustify FALSE (normalized weights, posterior distributions across gamma's and per each weight distribution); if robustify is TRUE ...
#'
#' @export
getPosterior <- function(env = NULL,
                         weights = c("BIC","pseudoBMAplus","stacking","WAIC","loglik","elpd_lfo","all"), 
                         method = c("default","spline"),
                         n_knots = 300,
                         nsim = 5e03,
                         B = 500,
                         robustify = FALSE,
                         R = 100,
                         subset_K = NULL,
                         sample_size_per_K = 10,
                         plot = FALSE,
                         top_K_models = 10){

    # subsetting env to a specific subset of stepwise models (either according to specific indices or to the number of intervals)
    if(!is.null(subset_K)){file_name <- paste(method,subset_K,collapse="_")}
    else{file_name <- paste(method,"all_intervals_dimensions",sep="_") }
    subsetting_input_results <- function(env, indices = NULL, K = NULL){
        if(!is.null(K)){
            indices <- env$K%in%K
        }
	    out <- list()
        out$BIC <- env$BIC[indices]
        out$elpd_lfo <- env$elpd_lfo[,indices]
        out$elpd_lfo_old <- env$elpd_lfo_old[,indices]
        out$elpd_waic <- env$elpd_waic[,indices]
        out$K <- env$K[indices]
        out$k_psis <- env$k_psis[,indices]
        out$k_psis_old <- env$k_psis_old[,indices]
        out$L <- env$L
        out$loglik <- env$loglik[indices]
        out$mle_betas <- env$mle_betas[indices,]
        out$n_sims_is <- env$n_sims_is
        out$p_waic <- env$p_waic[,indices]
        out$Q <- length(out$K)
        out$reestimations_psis <- env$reestimations_psis[indices]
        out$reestimations_psis_old <- env$reestimations_psis_old[indices]
        out$stats_names_endo <- env$stats_names_endo
        out$stats_names_exo <- env$stats_names_exo
        out$tau <- env$tau
        out$vcov_betas <- env$vcov_betas[,,indices]
        out$widths <- env$widths[indices,]
        out$true_effects <- env$true_effects
        return(out)
    }

    # subsetting env to a

    if(!is.null(subset_K)) env <- subsetting_input_results(env = env, K = subset_K)  

    # initializing useful vectors
    Q <- env$Q
    P <- length(env$stats_names_endo)
    S <- length(env$stats_names_exo)
    U <- P+S
    # understanding how many different weights are going to be used
    weights_vec <- c("BIC","pseudoBMAplus","stacking","WAIC","loglik","elpd_lfo")
    if(!is.null(weights) && weights[1] == "all"){
        weights <- weights_vec
        W <- length(weights)
    }
    else{
        status <- sapply(1:length(weights_vec), function(x) any(weights_vec[x] == weights))
        weights <- weights_vec[status]
        W <- length(weights)
    }

    # out : object list where to store the output
    out <- list()

    # when the whole ensemble of stepwise models is used (robustify = FALSE)
    if(!robustify){
        
        if(method == "default"){
            out <- getPosteriorDefault(env = env, weights = weights, nsim = nsim, n_knots = n_knots) #output is list(elapsed_time,weights,posterior_draws)
        }

        if(method == "spline"){
            out <- getPosteriorSpline(env = env, weights = weights, nsim = nsim, B = B, n_knots = n_knots)
        }

        out$top_K_models <- list()
        for(w in 1:W){
            indices_w <- order(out$weights[[names(out$weights)[w]]],decreasing = T)[1:top_K_models]
            out$top_K_models[[names(out$weights)[w]]]$indices <- indices_w
            out$top_K_models[[names(out$weights)[w]]]$widths <- env$widths[indices_w,]
            out$top_K_models[[names(out$weights)[w]]]$mle_betas <- env$mle_betas[indices_w,]
        }   
        
    }
    # when sampling from the ensemble in order to perform BMA R times and make results more robust and less dependent by the specific set of models (robustify = TRUE)
    else{
        if(method == "default"){

            out$weights <- list()
            out$posterior_draws <- list()

            for(w in 1:W){
                out$weights[[weights[w]]] <- matrix(NA, nrow = Q, ncol = R)
                out$posterior_draws[[weights[w]]] <- array(NA, dim=c(R*nsim,U,n_knots)) 
                dimnames(out$posterior_draws[[weights[w]]])[[2]] <- c(env$stats_names_endo,env$stats_names_exo)
            }

            for(r in 1:R){
                print(r)
                sample_r <- as.vector(unlist(tapply(X = 1:Q ,INDEX = env$K,FUN =function(x) sample(x = x, size = sample_size_per_K, replace=FALSE))))
                bag_r <- subsetting_results(env = env, indices = sample_r)
                posterior_r <- getPosteriorDefault(env = bag_r, weights = weights, nsim = nsim, n_knots = n_knots)
                if(r==1){ out$elapsed_time<- posterior_r$elapsed_time}
                for(w in 1:W){
                out$weights[[names(posterior_r$weights)[w]]][sample_r,r] <- posterior_r$weights[[names(posterior_r$weights)[w]]]
                out$posterior_draws[[names(posterior_r$weights)[w]]][c(((r-1)*nsim+1):(r*nsim)),,] <- posterior_r$posterior_draws[[names(posterior_r$weights)[w]]] #change nsim with B when with splines
                }
            }
        }

        if(method == "spline"){

            out$weights <- list()
            out$posterior_draws <- list()

            for(w in 1:W){
                out$weights[[weights[w]]] <- matrix(NA, nrow = Q, ncol = R)
                out$posterior_draws[[weights[w]]] <- array(NA, dim=c(R*B,U,n_knots))
                dimnames(out$posterior_draws[[weights[w]]])[[2]] <- c(env$stats_names_endo,env$stats_names_exo)
            }

            for(r in 1:R){
                print(r)
                sample_r <- as.vector(unlist(tapply(X = 1:Q ,INDEX = env$K,FUN = function(x) sample(x = x, size = sample_size_per_K, replace=FALSE))))
                bag_r <- subsetting_results(env = env, indices = sample_r)
                posterior_r <- getPosteriorSpline(env = bag_r, weights = weights, nsim = nsim, B = B, n_knots = n_knots)
                if(r==1){ out$elapsed_time<- posterior_r$elapsed_time}
                for(w in 1:W){
                out$weights[[names(posterior_r$weights)[w]]][sample_r,r] <- posterior_r$weights[[names(posterior_r$weights)[w]]]
                out$posterior_draws[[names(posterior_r$weights)[w]]][c(((r-1)*B+1):(r*B)),,] <- posterior_r$posterior_draws[[names(posterior_r$weights)[w]]] #change nsim with B when with splines
                }
            }
        }

        out$top_K_models <- list()
        for(w in 1:W){
            out$top_K_models[[names(out$weights)[w]]] <- matrix(NA, nrow = sample_size_per_K*length(unique(env$K)), ncol = R)
            for(r in 1:R){
                indices_w_r <- order(na.omit(out$weights[[names(out$weights)[w]]][,r]),decreasing = T) # this vector will have length R
                out$top_K_models[[names(out$weights)[w]]][,r] <- na.omit(out$weights[[names(out$weights)[w]]][,r])[indices_w_r]
            }

        }
    }

    if(plot){
        out$trend <- plotPosterior(posterior_draws = out, type = "trend", true_effects = env$true_effects, file_name = file_name)
        out$credibility_intervals <- plotPosterior(posterior_draws = out, type = "credibility intervals", true_effects = env$true_effects, file_name = file_name)
    }

    return(out)
}


#' getPosteriorSpline
#'
#' A function which return....
#'
#' @param env environment/list stepwiseModelsREH where models are stored (or a subset of it)
#' @param weights type of weighting system to be used
#' @param nsim number of simulations per gamma
#' @param B number of splines
#' @param n_knots number of gamma's
#'
#' @return  list of: normalized weights, posterior distributions across gamma's and per each weight distribution
#'
#' @export
getPosteriorSpline <- function(env = NULL ,
                                weights = c("BIC","pseudoBMAplus","stacking","WAIC","loglik","elpd_lfo","all"),
                                nsim = 1e04,
                                B = 1e03,
                                n_knots = 1e03){
    if(is.null(env)) stop("env is not specified.")

    # initializing useful vectors
    P <- length(env$stats_names_endo)
    S <- length(env$stats_names_exo)
    U <- P+S

    # getPosteriorDefault()
    posterior_default <- getPosteriorDefault(env = env, weights = weights, nsim = nsim, n_knots = n_knots)
    w_list <- posterior_default$weights
    # output 
    out <- list() 
    out$elapsed_time <- elapsed_time <- posterior_default$elapsed_time
    out$weights <- posterior_default$weights
    out$posterior_draws <- list()

    # BEGIN algorithm
    for(w in 1:length(w_list)){
        out$posterior_draws[[names(w_list)[w]]] <- array(NA, dim=c(B,U,n_knots)) # create output array for the specific weighting system
        dimnames(out$posterior_draws[[names(w_list)[w]]])[[2]] <- dimnames(posterior_default$posterior_draws[[names(w_list)[w]]])[[2]]
        for(b in 1:B){
            sample_b <- sample(x = 1:nsim, size = n_knots, replace = FALSE) # n_knots has to be lower than nsim (condition !)
            posterior_b <- t(sapply(1:n_knots, function(x) posterior_default$posterior_draws[[names(w_list)[w]]][sample_b[x],,x]))
            out$posterior_draws[[names(w_list)[w]]][b,c((P+1):dim(posterior_b)[2]),] <- t(posterior_b[,c((P+1):dim(posterior_b)[2])]) # non funziona la posteriori delle esogene (forse selziono male le colonne)
            # spline smoothing is performed only for endogenous statistics which MUST be the first P variables 
            for(p in 1:P){
                spline_p_b <- stats::smooth.spline(x = elapsed_time, y = posterior_b[,p]) 
                out$posterior_draws[[names(w_list)[w]]][b,p,] <- spline_p_b$y
            }
        }
    }
    # END algorithm
    return(out)
    
}


#' getPosteriorDefault
#'
#' A function which return a list of weights according to different possible weighting systems
#'
#' @param env environment/list stepwiseModelsREH where models are stored (or a subset of it)
#' @param weights type of weighting system to be used
#' @param nsim number of simulations per gamma
#' @param n_knots number of gamma's
#'
#' @return  list of: normalized weights, posterior distributions across gamma's and per each weight distribution
#'
#' @export
getPosteriorDefault <- function(env = NULL ,
                        weights = c("BIC","pseudoBMAplus","stacking","WAIC","loglik","all"),
                        nsim = 1e04,
                        n_knots = 1e03){
    if(is.null(env)) stop("env is not specified.")

    max_width <- max(na.omit(as.vector(env$widths)))
    elapsed_time <- seq(from = 1e-05, to = max_width,length = n_knots)

    out <- list()
    out$elapsed_time <- elapsed_time

    # prepare input for the getModelsWeights function
    input_getModelsWeights_loc <- list()
    input_getModelsWeights_loc$Q <- env$Q
    input_getModelsWeights_loc$BIC <- env$BIC
    input_getModelsWeights_loc$elpd_waic <- env$elpd_waic
    input_getModelsWeights_loc$elpd_lfo <- env$elpd_lfo
    input_getModelsWeights_loc$loglik <- env$loglik
    input_getModelsWeights_loc$vcov_betas <- env$vcov_betas
    input_getModelsWeights_loc$K <- env$K

    w_list <- getModelsWeights(input = input_getModelsWeights_loc, weights = weights)
    
    # find posterior for that knot (for each variable in the model)
    input_drawFromPosterior_loc <- list()
    input_drawFromPosterior_loc$Q <- env$Q
    input_drawFromPosterior_loc$stats_names_endo <- env$stats_names_endo
    input_drawFromPosterior_loc$stats_names_exo <- env$stats_names_exo
    input_drawFromPosterior_loc$mle_betas <- env$mle_betas
    input_drawFromPosterior_loc$vcov_betas <- env$vcov_betas
    input_drawFromPosterior_loc$K <- env$K
    input_drawFromPosterior_loc$widths <- env$widths

    # draws from the posterior distribution
    posterior_draws <- drawFromPosterior(input = input_drawFromPosterior_loc, weights = w_list, nsim = nsim, elapsed_time = elapsed_time, max_width = max_width)

    # store weights and posterior draws
    out$weights <- w_list
    out$posterior_draws <- posterior_draws

    return(out)
}



#' getModelsWeights
#'
#' A function which return a list of weights according to different possible weighting systems
#'
#' @param input list of vector, matrices and array useful for weights computation
#' @param weights type of weighting system to be used
#'
#' @return  list of normalized weights
#'
#' @export
getModelsWeights <- function(input = NULL, weights = c("BIC","pseudoBMAplus","stacking","WAIC","loglik","elpd_lfo","all")){

    weights_vec <- c("BIC","pseudoBMAplus","stacking","WAIC","loglik","elpd_lfo")
    if(!is.null(weights) && weights[1] == "all"){
        status <- !logical(length(weights_vec))
    }
    else{
        status <- sapply(1:length(weights_vec), function(x) any(weights_vec[x] == weights))
    }
    if(is.null(input)) stop("input is not specified.")

    out <- list()
    # BIC
    if(status[1]){
        if(is.null(input$BIC)) stop("BIC vector was not found.")
        BIC_weights <- exp(-(input$BIC/2-min(input$BIC/2))) / sum(exp(-(input$BIC/2-min(input$BIC/2))))
        out$BIC <- BIC_weights
    }
    
    # pseudoBMAplus
    if(status[2]){
        if(is.null(input$elpd_lfo)) stop("ELPD matrix was not found.")

        # Performs pseudoBMA+ by using elpd_(psis)_lfo
        elpd_lfo <- input$elpd_lfo

        # algorithm with loo::pseudobma_weights() function
        out$pseudoBMAplus <- as.vector(loo::pseudobma_weights(elpd_lfo))

        ## my algorithm without the use of gtools package
        ## Random Dirichlet values: Dir(1,...,1)
        #alphas <- matrix(rexp(1e03*dim(elpd_lfo)[1],rate=1),nrow=dim(elpd_lfo)[1],ncol=1e03)
        #alphas <- t(apply(alphas,2,function(x) x/sum(x)))
        #z_v_q <- (alphas%*%elpd_lfo)*dim(elpd_lfo)[1] #z*q
        #w_v_q <- apply(z_v_q,1,function(x) exp(x/2-max(x/2))/sum(exp(x/2-max(x/2))))
        ### old correction w_v_q[is.na(w_v_q)] <- 0
        #w_q <- apply(w_v_q,1,mean)
        #w_q <- w_q/sum(w_q)
        #out$pseudoBMAplus <- w_q

        ## algorithm with gtools::rdirichlet() function
        #alphas <- rdirichlet(1e3,rep(1,dim(elpd_lfo)[1]))
        #z_v_q <- (alphas%*%elpd_lfo)*dim(elpd_lfo)[1] #[VxQ]
        #z_v_q <- t(apply(z_v_q,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))
        #w_q <- apply(z_v_q,2,mean)
        #z_q <- apply(z_v_q,2,mean)
        #w_q <- exp(z_q-max(z_q))
        #out$pseudoBMAplus_dirichlet <- w_q/sum(w_q)

    }

    # Stacking of predictive distributions
    if(status[3]){
        if(is.null(input$elpd_lfo)) stop("ELPD matrix was not found.")

        #stacking_opt<- function(pars,matrix_elpd)
        #{
        #   pars <- exp(pars)/sum(exp(pars))
        #   -mean(log(pars%*%matrix_elpd))
        #}

        # Performs SPD by using elpd_(psis)_lfo
        elpd_lfo <- input$elpd_lfo
        
        # algorithm with loo::stacking_weights() function
        if(!is.null(out$pseudoBMAplus)) {start_values <- out$pseudoBMAplus}
        else{
        start_values <- as.vector(loo::pseudobma_weights(elpd_lfo))
        }
        out$stacking <- as.vector(loo::stacking_weights(elpd_lfo,optim_control = list(start_values)))

        # my alrgorithm to perform SPD
        #w_pseudo_bma_plus <- out$pseudoBMAplus # in a future function if out$pseudoBMAplus doesn't exist, estimate it
        #min_value <- optimized_w <- NULL
        #for(i in 1) #:10
        #{
        #   opt0 <- optim(par=log(w_pseudo_bma_plus),fn=stacking_opt,matrix_elpd=t(exp(elpd_lfo)))
        #   min_value <- c(min_value,opt0$value)
        #   optimized_w <- rbind(optimized_w,opt0$par)
        #   cat('Stacking: ',i,'-th finished \n')
        #}
        #w_opt <- optimized_w[which.min(min_value),]
        #stacking_w <- exp(w_opt-max(w_opt))/sum(exp(w_opt-max(w_opt)))
        #out$stacking <- stacking_w
    }

    # WAIC
    if(status[4]){
        if(is.null(input$elpd_waic)) stop("ELPD (WAIC) matrix was not found.")
        #if(dim(na.omit(input$elpd_waic))[2]<input$Q) stop("Not all the stepwise models have an ELPD (WAIC) estimated (vcov matrix was singular).")
        elpd_waic <- apply(rbind(input$elpd_waic),2,sum)
        WAIC_weights <- exp((elpd_waic/2-max(elpd_waic/2))) / sum(exp((elpd_waic/2-max(elpd_waic/2))))
        out$WAIC <- WAIC_weights
    }

    #loglik
    if(status[5]){
        if(is.null(input$loglik)) stop("loglik vector was not found.")
        loglik_weights <- exp((input$loglik/2-max(input$loglik/2))) / sum(exp((input$loglik/2-max(input$loglik/2))))
        out$loglik <- loglik_weights
    }  

    #elpd_lfo
    if(status[6]){
        if(is.null(input$elpd_lfo)) stop("ELPD (LFO) matrix was not found.")
        #if(dim(na.omit(input$elpd_waic))[2]<input$Q) stop("Not all the stepwise models have an ELPD (WAIC) estimated (vcov matrix was singular).")
        elpd_lfo <- apply(rbind(input$elpd_lfo),2,sum)
        elpd_lfo_weights <- exp((elpd_lfo/2-max(elpd_lfo/2))) / sum(exp((elpd_lfo/2-max(elpd_lfo/2))))
        out$elpd_lfo <- elpd_lfo_weights
    }
    return(out)
}


#' drawFromPosterior
#'
#' A function which return a list of draws according to different possible weighting systems
#'
#' @param input stepwiseModelsREH environment
#' @param weights list of vectors of normalized weights (BIC, pseudoBMA+, ...)
#' @param nsim number of simulations from the posterior
#' @param elapsed_time number of knots for the elapsed time axis
#' @param max_width maximum widths of intervals (used in the "stepwise_models" method)
#'
#' @return posterior estimates of memory decay per each type of weights
#'
#' @export
drawFromPosterior <- function(input,weights,nsim,elapsed_time,max_width){

    stats_names <- c(input$stats_names_endo,input$stats_names_exo) 
    P <- length(input$stats_names_endo)
    S <- length(input$stats_names_exo)
    U <- P + S
    Q <- input$Q
    W <- length(weights)
      
    out <- list()
    ## stepwise models ##
    K <- input$K
    n_knots <- length(elapsed_time)
    # function to find which stepwise effect to pick at x (elapsed time value) and give intervals described by the vector 'widths'
    find_step <- function(x,widths){
        k <- length(widths)-1
        if(x <= widths[length(widths)]){
            while(k >= 1){
                if(x >= widths[k]){
                    return(k)
                }
                else{
                    k <- k - 1
                }
            }
        }
        else{return(NaN)}
    }

    # find step (corresponding to which stepwise parameter per each model to consider given the knot)
    which_pars <- array(9999, dim =c(Q,U,n_knots))
    for(i in 1:n_knots){
        for(q in 1:Q) {
                    which_step <- find_step(x = elapsed_time[i], widths = input$widths[q,1:(K[q]+1)])
                    if(!is.na(which_step)){
                        which_pars[q,,i] <- c(which_step+(c(0:(P-1))*K[q]),c((P*K[q]+1):(P*K[q]+S)))-1 # (-1) for Rcpp indexing (!)
                    }
                    rm(which_step)
        }
    }
    

    # BEGIN for loop generation from the posterior for each weighting system in the input object 'weights'
    print("Drawing from posterior...")
    for(w in 1:W){
        #print(names(weights)[w])
        # first stage: sample nsim times from the set of Q models and following probabilities given by weights[[w]]
        sample_models <- sample(x = c(1:Q), size = nsim, replace = TRUE, prob = weights[[names(weights)[w]]])-1

        # second stage : sample from normal distribution with specific parameters given by which_pars and for each model sampled in the previous stage
        local_draws <- getDraws(sample_models = sample_models, 
                                which_pars = which_pars,
                                n_pars = ((K*P)+S),
                                n_stats = U,
                                input = input,
                                knots_seq = elapsed_time)
        local_draws[which(local_draws==9999)] <- NaN
        # create output list object
        out[[names(weights)[w]]] <- local_draws #t(apply(local_draws,c(2,3),function(x) quantile(x,c(0.5),na.rm=TRUE))) # find median was the old command (we want all the draws)
        dimnames(out[[names(weights)[w]]])[[2]] <- stats_names
        rm(local_draws)
    }
    # END for loop
return(out)
}


#' plotPosterior
#'
#' A function which plots posterior trends and eventually returns the ggplot objects
#'
#' @param posterior_draws list of posterior draws from different weighting systems (BIC, pseudoBMA+, ...)
#' @param type type of plot: trend type will just plot the median, credibility intervals type will also plot the credibitlity interval at 0.95
#' @param true_effects effects list of endo effects and exo_effects 
#' @param save TRUE/FALSE whether to plot results right after the creation of the list of ggplot objects
#' @param file_name file name is usually a combination of method and subset_K parameter from the getPosterior function (it is useful just for now in order to save results fastly)
#'
#' @return ggplot objects depending on the specified type
#'
#' @export
plotPosterior <- function(posterior_draws, type = c("trend","credibility intervals"), true_effects = NULL, save = TRUE, file_name=NULL){

    weights <- names(posterior_draws$weights)
    W <- length(weights)
    elapsed_time <- posterior_draws$elapsed_time
    n_knots <- dim(posterior_draws$posterior_draws[[weights[1]]])[3]
    P <- dim(posterior_draws$posterior_draws[[weights[1]]])[2]
    stats_names <- dimnames(posterior_draws$posterior_draws[[weights[1]]])[[2]]
    credibility_intervals <- list()
    for(w in 1:W) credibility_intervals[[weights[w]]] <- apply(posterior_draws$posterior_draws[[weights[w]]],c(2,3),function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))
    
    data <- list()
    for(p in 1:P){
        mat_loc <- matrix(NA, nrow=n_knots*W, ncol=3)
        for(w in 1:W) mat_loc[c((1+(w-1)*n_knots):(n_knots+(w-1)*n_knots)),] <- t(credibility_intervals[[weights[w]]][,p,])

        data[[stats_names[p]]] <- data.frame(timeElapsed = rep(elapsed_time,W), 
                                                weights = rep(weights,each=n_knots),
                                                lb = mat_loc[,1],
                                                effect = mat_loc[,2],
                                                ub = mat_loc[,3])

    }

    out <- list()

    vdecay <- Vectorize(bremory::decay,vectorize.args="x")

    if(type=="trend"){
        out <- list() #list of ggplot objects
        
        for(p in 1:P) #change p with u (u is the index that refers to the statistics)
        {
            ggplot_obj_loc <- ggplot(data[[stats_names[p]]], aes(x=timeElapsed, y=effect, group=weights, color=weights)) + geom_line()

            if(stats_names[p]%in%names(true_effects$exo_effects)){ 
                ggplot_obj_loc <- ggplot_obj_loc + geom_hline(aes_string(yintercept = true_effects$exo_effects[[stats_names[p]]],color = "TRUE"))
                }
            else{
                ggplot_obj_loc <- ggplot_obj_loc + stat_function(fun=vdecay, args=list(type = true_effects$endo_effects[[stats_names[p]]]$decay_type, pars = true_effects$endo_effects[[stats_names[p]]]$pars), geom="line", aes(color="TRUE"),linetype="dashed")
            }

            #ggplot_obj_loc <- ggplot_obj_loc + scale_color_gradientn(colours = rainbow(W)) + ggtitle(paste("Trend of",stats_names[p],sep=" "))
            #scale_colour_manual("Legend", values = c("BIC" ="green","pseudoBMAplus" = "blue","stacking"="purple","WAIC"="yellow","logL"="turquoise","TRUE"="red")) +

            out[[stats_names[p]]] <- ggplot_obj_loc
            rm(ggplot_obj_loc)    

            
        }
        if(save == TRUE) {
            svg(paste(file_name,"_trends.svg",sep=""), width=16, height=16)
            grid.arrange(grobs = out, ncol=ceiling(sqrt(P)), nrow=ceiling(sqrt(P)))
            dev.off()
        }
    }
  
    
    if(type=="credibility intervals"){
        out <- list() #list of ggplot objects
        for(p in 1:P){
            out[[stats_names[p]]] <- list()
            
            for(w in 1:W){
                data_p <- data.frame(timeElapsed = elapsed_time, 
                                    lb = credibility_intervals[[weights[w]]][1,p,],
                                    effect =credibility_intervals[[weights[w]]][2,p,],
                                    ub = credibility_intervals[[weights[w]]][3,p,])
                ggplot_obj_loc <- ggplot(data_p, aes(x=timeElapsed, y=ub)) + 
                geom_line(aes(y = lb), colour="transparent") + 
                geom_line(aes(y = ub), colour="transparent") +
                geom_ribbon(data=subset(data_p, 0 <= timeElapsed & timeElapsed <= max(elapsed_time)), aes(ymin=lb,ymax=ub), fill="purple", alpha=0.2) +
                geom_line(aes(y = effect),size=0.5, linetype = "longdash",colour="purple") +
                ggtitle(paste("Trend of",stats_names[p],"(",weights[w],")",sep=" ")) +
                xlab(expression(gamma)) + 
                ylab(expression(beta)) +
                labs(colour = "Trend") +
                theme(axis.text=element_text(size=12),axis.title = element_text(size=12),
                          title=element_text(size=12,face="italic")) +theme(legend.position = "none") #+
                # command to check scale_fill_manual(values=c(clear,blue)) +
                if(stats_names[p]%in%names(true_effects$exo_effects)){ 
                    ggplot_obj_loc <- ggplot_obj_loc + geom_hline(aes_string(yintercept = true_effects$exo_effects[[stats_names[p]]],color = "TRUE"),size = 0.5)
                    }
                else{
                    ggplot_obj_loc <- ggplot_obj_loc + stat_function(fun=vdecay, args=list(type = true_effects$endo_effects[[stats_names[p]]]$decay_type, pars = true_effects$endo_effects[[stats_names[p]]]$pars), geom="line", size=0.5, aes(color="TRUE"),linetype="dashed")
                }
                out[[stats_names[p]]][[weights[w]]] <- ggplot_obj_loc
                rm(ggplot_obj_loc)
            }
            if(save == TRUE) {
                svg(paste(file_name,"_",stats_names[p],".svg",sep=""), width=16, height=16)
                grid.arrange(grobs = out[[stats_names[p]]], ncol=ceiling(sqrt(W)), nrow=ceiling(sqrt(W)))
                dev.off()
                }
        }
    }
    return(out)

}




#' DoAnalysis
#'
#' A function which plots posterior trends and eventually returns the ggplot objects
#'
#' @param env environment where to work
#' @param do which analysis to perform
#' @param weights which vector of weights to consider
#' @param K which dimension
#' @param n_knots number of knots 
#' @param top best top-x models to select 
#'
#' @return ggplot objects depending on the specified type
#'
#' @export
DoAnalysis <- function(env = NULL ,
                       do = c("mle","sd","sd_increasing_time_lengths","mle_precision_fit","mle_width","mle_elapsed_time","weights_precision","fit_precision_interval","mle_transparent_lines","beta_loglik","widths_weights"),
                       weights = c("BIC","pseudoBMAplus","stacking","WAIC","loglik","all"),
                       K = 2,
                       n_knots = 300,
                       top = 10){

    if(is.null(env)) stop("env is not specified.")                 
    subsetting_intervals <- function(env,K){
        out <- list()
        out$BIC <- env$BIC[env$K%in%K]
        out$elpd_lfo <- env$elpd_lfo[,env$K%in%K]
        out$elpd_waic <- env$elpd_waic[,env$K%in%K]
        out$k_psis <- env$k_psis[,env$K%in%K]
        out$loglik <- env$loglik[env$K%in%K]
        out$mle_betas <- env$mle_betas[env$K%in%K,]
        out$p_waic <- env$p_waic[,env$K%in%K]
        out$reestimations_psis <- env$reestimations_psis[env$K%in%K]
        out$vcov_betas <- env$vcov_betas[,,env$K%in%K]
        out$widths <- env$widths[env$K%in%K,]
        out$K <- env$K[env$K%in%K]
        out$Q <- length(out$K)
        out$tau <- env$tau
        out$n_sims_is <- env$n_sims_is
        out$L <- env$L 
        out$stats_names_endo <- env$stats_names_endo
        out$stats_names_exo <- env$stats_names_exo
        return(out)
    }
    getStepBeta <- function(x,widths,pars){
        k <- length(widths)
        while(k >= 1){
            if(x >= widths[k]){
                return(pars[k])
            }
            else{
                k <- k - 1
            }
        }
    }
    vdecay <- Vectorize(decay,vectorize.args="x")


    if(do == "mle"){
        ###
        subset_models <- subsetting_intervals(env = env, K = K)
        ###
        max_width <- max(na.omit(subset_models$widths[1,]))
        elapsed_time <- seq(1e-04,max_width-1e-04,length=n_knots)
        n_endo_effects <- length(subset_models$stats_names_endo)
        # prepare input for the getModelsWeights function
        input_getModelsWeights_loc <- list()
        input_getModelsWeights_loc$Q <- subset_models$Q
        input_getModelsWeights_loc$BIC <- subset_models$BIC
        input_getModelsWeights_loc$elpd_waic <- subset_models$elpd_waic
        input_getModelsWeights_loc$elpd_lfo <- subset_models$elpd_lfo
        input_getModelsWeights_loc$loglik <- subset_models$loglik
        w_list <- getModelsWeights(input = input_getModelsWeights_loc, weights = weights)
        W <- length(w_list)
        weights <- names(w_list)
        for(w in 1:W){
            top_w <- order(w_list[[weights[w]]],decreasing=T)[1:top]
            x11()
            par(mfrow=c(ceiling(sqrt(n_endo_effects)),ceiling(sqrt(n_endo_effects))))
            for(p in 1:n_endo_effects){
                output_endo_p <- array(NA, dim=c(length(top_w),length(elapsed_time),3))                                                    
                for(i in 1:length(top_w)){
                    mles <- subset_models$mle_betas[top_w[i],(1+(p-1)*K):(p*K)]
                    sds <- diag(subset_models$vcov_betas[(1+(p-1)*K):(p*K),(1+(p-1)*K):(p*K),top_w[i]])**0.5
                    lbs <- mles - sds*qnorm(0.975)
                    ubs <- mles + sds*qnorm(0.975)
                    output_endo_p[i,,2] <- sapply(elapsed_time,function(x) getStepBeta(x=x,widths=na.omit(subset_models$widths[top_w[i],]),pars=mles)) # mle
                    output_endo_p[i,,1] <- sapply(elapsed_time,function(x) getStepBeta(x=x,widths=na.omit(subset_models$widths[top_w[i],]),pars=lbs)) # lower bound
                    output_endo_p[i,,3] <- sapply(elapsed_time,function(x) getStepBeta(x=x,widths=na.omit(subset_models$widths[top_w[i],]),pars=ubs)) # upper bound                                              
                }
                plot(1:10,1:10,type="n",ylim=c(min(output_endo_p),max(output_endo_p)),xlim=c(min(elapsed_time),max(elapsed_time)),
                            xlab=expression(gamma),
                            ylab=expression(beta),
                            main=paste("effect",subset_models$stats_names_endo[p],"with",weights[w],sep=" ")) 
                for(i in 1:length(top_w)){                                                       
                    for(value in 1:3){
                        col_value <- ifelse(value == 2,2,1)
                        points(elapsed_time,output_endo_p[i,,value],pch="-",col=col_value)
                    }
                } 
            }                                             
        }
    }
    if(do == "sd"){
        ###
        subset_models <- subsetting_intervals(env = env, K = K)
        ###
        max_width <- max(na.omit(subset_models$widths[1,]))
        elapsed_time <- seq(1e-04,max_width-1e-04,length=n_knots)
        n_endo_effects <- length(subset_models$stats_names_endo)
        time_lengths <- t(apply(subset_models$widths[,1:(K+1)],1,function(x) diff(x)))
        x11()
        par(mfrow=c(ceiling(sqrt(n_endo_effects)),ceiling(sqrt(n_endo_effects))))
        for(p in 1:n_endo_effects){
            sds <- t(sapply(1:dim(subset_models$vcov_betas)[3], function(x) diag(subset_models$vcov_betas[(1+(p-1)*K):(p*K),(1+(p-1)*K):(p*K),x])**0.5))
            matplot(time_lengths,sds,type="p",main=subset_models$stats_names_endo[p])
        }
    }

    if(do == "sd_increasing_time_lengths"){
        ###
        subset_models <- subsetting_intervals(env = env, K = K)
        ###
        max_width <- max(na.omit(subset_models$widths[1,]))
        elapsed_time <- seq(1e-04,max_width-1e-04,length=n_knots)
        n_endo_effects <- length(subset_models$stats_names_endo)
        time_lengths <- t(apply(subset_models$widths[,1:(K+1)],1,function(x) diff(x)))
        selection <- apply(time_lengths,1,function(x) !is.unsorted(x))
        if(length(selection)!=0){
            time_lengths <- time_lengths[selection,]
            vcov_selection <- subset_models$vcov_betas[,,selection]
            x11()
            par(mfrow=c(ceiling(sqrt(n_endo_effects)),ceiling(sqrt(n_endo_effects))))
            for(p in 1:n_endo_effects){
                sds <- t(sapply(1:dim(vcov_selection)[3], function(x) diag(vcov_selection[(1+(p-1)*K):(p*K),(1+(p-1)*K):(p*K),x])**0.5))
                matplot(time_lengths,sds,type="p",main=subset_models$stats_names_endo[p])
            }

        input_getModelsWeights_loc <- list()
        input_getModelsWeights_loc$Q <- subset_models$Q
        input_getModelsWeights_loc$BIC <- subset_models$BIC
        input_getModelsWeights_loc$elpd_waic <- subset_models$elpd_waic
        input_getModelsWeights_loc$elpd_lfo <- subset_models$elpd_lfo
        w_list <- getModelsWeights(input = input_getModelsWeights_loc, weights = weights)
        W <- length(w_list)
        weights <- names(w_list)
        x11()
        par(mfrow=c(ceiling(sqrt(W)),ceiling(sqrt(W))))
        for(w in 1:W){
            plot(1:subset_models$Q,w_list[[weights[w]]],xlab="q",ylab="w",main=weights[w])
            points(which(selection==TRUE),w_list[[weights[w]]][selection],pch=19,col=2)
        }

        }
        else{stop("No models found !")}
    }

    if(do == "mle_precision_fit"){
         P <- length(env$stats_names_endo)
        for(K in min(env$K):max(env$K)){ 
            mle_K <- env$mle_betas[env$K==K,1:(K*P)]
            vcov_K <- env$vcov_betas[1:(K*P),1:(K*P),env$K==K]
            sds <- t(apply(vcov_K,3,function(x) diag(x)**0.5))
            trace_sds <- apply(sds,1,function(x) sum(x**2))
            widths_K <- env$widths[env$K==K,]
            loglik_K <- env$loglik[env$K==K]
            pseudo_fic_llik <- exp(((1/trace_sds+loglik_K)-max(1/trace_sds+loglik_K)))/sum(exp(((1/trace_sds+loglik_K)-max(1/trace_sds+loglik_K))))


            # first: creating ggplot objects
            plot_K <- list()
            for(p in 1:P){
                for(k in 1:K){
                    if(k==K){
                        widths_loc <- widths_K[,k] # this just because the last interval would have always the same gamma
                    }
                    else{
                        widths_loc <- widths_K[,k+1]
                    }
                    precision <- (1/sds[,k+K*(p-1)]) 
                    fit <- exp(loglik_K-max(loglik_K)) 
                    data <- data.frame(gamma=widths_loc, beta=mle_K[,k+K*(p-1)] , precision = precision, fit = fit, pseudo_fic = pseudo_fic_llik)
                    plot_loc <- ggplot(data, aes(gamma, beta)) + geom_point(colour = "orange", aes(alpha = precision, size = pseudo_fic_llik))
                    plot_loc <- plot_loc  + ggtitle(paste(env$stats_names_endo[p],"(",k,")",sep="")) + xlab(expression(gamma)) + ylab(expression(beta)) + theme(legend.position = "none") 
                    plot_K[[k+K*(p-1)]] <- plot_loc
                }
            }
            # second: arrange and save the plot
            svg(paste(K,"intervals_mle_fit_precision.svg",sep=""), width=20, height=15)
            grid.arrange(grobs = plot_K, ncol=K)
            dev.off()
            }
        }

    if(do == "mle_width"){
        P <- length(env$stats_names_endo)
        widths <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) diff(env$widths[x,1:(env$K[x]+1)])))
        interval <- as.factor(unlist(sapply(1:dim(env$mle_betas)[1], function(x) c(1:env$K[x]))))
        number_of_intervals <- as.factor(unlist(sapply(1:dim(env$mle_betas)[1], function(x) rep(env$K[x],env$K[x]))))
        loglik <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) rep(env$loglik[x],env$K[x]) ))
        loglik_diff <- (loglik - max(loglik)) /10 # /10 scaling factor

        # first: creating ggplot objects
        for(p in 1:P){
            plot_gg <- list()
            betas <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) env$mle_betas[x,(1+env$K[x]*(p-1)):(env$K[x]*p)]))
            sds <- unlist(sapply(1:dim(env$mle_betas)[1],function(x) diag(env$vcov_betas[(1+env$K[x]*(p-1)):(env$K[x]*p),(1+env$K[x]*(p-1)):(env$K[x]*p),x])))
            data <- data.frame(widths = widths, betas = betas, interval = interval, K = number_of_intervals, sds = sds, loglik_diff = loglik_diff)

            plot_loc <- ggplot(data, aes(widths, betas)) + geom_point(aes(size = sds, alpha = loglik_diff, colour = K)) + facet_wrap(~interval) 
            plot_loc <- plot_loc + ggtitle(paste(env$stats_names_endo[p],sep="")) + xlab("interval width") + ylab(expression(beta)) 
            plot_gg[[1]] <- plot_loc

            # second: arrange and save the plot
            svg(paste(env$stats_names_endo[p],"_mle_width.svg",sep=""), width=15, height=15)
            grid.arrange(grobs = plot_gg, ncol = 1)
            dev.off()
        }

    }

    if(do == "mle_elapsed_time"){
        P <- length(env$stats_names_endo)
        elapsed_time <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) c(env$widths[x,2:(env$K[x])],env$widths[x,env$K[x]]))) # (we always take the second last gamma for the last interval)
        interval <- as.factor(unlist(sapply(1:dim(env$mle_betas)[1], function(x) c(1:env$K[x]))))
        number_of_intervals <- as.factor(unlist(sapply(1:dim(env$mle_betas)[1], function(x) rep(env$K[x],env$K[x]))))
        loglik <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) rep(env$loglik[x],env$K[x]) ))
        # first: creating ggplot objects
        
        for(p in 1:P){
            plot_gg <- list()
            betas <- unlist(sapply(1:dim(env$mle_betas)[1], function(x) env$mle_betas[x,(1+env$K[x]*(p-1)):(env$K[x]*p)]))
            sds <- unlist(sapply(1:dim(env$mle_betas)[1],function(x) diag(env$vcov_betas[(1+env$K[x]*(p-1)):(env$K[x]*p),(1+env$K[x]*(p-1)):(env$K[x]*p),x]))) * 10
            data <- data.frame(time = elapsed_time, betas = betas, interval = interval, K = number_of_intervals, sds = sds, loglik = loglik)

            plot_loc <- ggplot(data, aes(time, betas)) + geom_point(aes(size = sds, alpha = loglik, colour = K)) + stat_function(fun=vdecay, args=list(type = env$true_effects$endo_effects[[env$stats_names_endo[p]]]$decay_type, pars = env$true_effects$endo_effects[[env$stats_names_endo[p]]]$pars), geom="line", size=0.2, aes(color="TRUE")) + facet_wrap(~interval)
            plot_loc <- plot_loc + ggtitle(paste(env$stats_names_endo[p],sep="")) + xlab("interval width") + ylab(expression(beta)) 
            plot_gg[[1]] <- plot_loc

            # second: arrange and save the plot
            svg(paste(env$stats_names_endo[p],"_mle_elapsed_time.svg",sep=""), width=15, height=15)
            grid.arrange(grobs = plot_gg, ncol = 1)
            dev.off()
        }

    }
    

    if(do == "weights_precision"){
        input_getModelsWeights_loc <- list()
        input_getModelsWeights_loc$Q <- env$Q
        input_getModelsWeights_loc$BIC <- env$BIC
        input_getModelsWeights_loc$elpd_waic <- env$elpd_waic
        input_getModelsWeights_loc$elpd_lfo <- env$elpd_lfo
        input_getModelsWeights_loc$K <- env$K 
        input_getModelsWeights_loc$loglik <- env$loglik
        input_getModelsWeights_loc$vcov_betas <- env$vcov_betas
        w_list <- getModelsWeights(input = input_getModelsWeights_loc, weights = weights)
        W <- length(w_list)
        names_w <- names(w_list)
        P <- length(env$stats_names_endo)
        interval <- as.factor(env$K)
        precision <- unlist(sapply(1:dim(env$mle_betas)[1],function(x) 1/sum(diag(env$vcov_betas[1:(env$K[x]*P),1:(env$K[x]*P),x]))))
        # first: creating ggplot objects
        plot_gg <- list()
        for(w in 1:W){
            weights <- w_list[[w]]
            data <- data.frame(precision = precision, weights = weights, interval = interval)

            plot_loc <- ggplot(data, aes(precision,weights)) + geom_point(aes(colour = interval)) 
            plot_loc <- plot_loc + ggtitle(paste(names_w[w],sep="")) + theme(legend.position = "none") 
            plot_gg[[w]] <- plot_loc
        }
        # second: arrange and save the plot
        svg("weights_precision.svg", width=15, height=15)
        grid.arrange(grobs = plot_gg, ncol=ceiling(sqrt(W)),nrow=ceiling(sqrt(W)))
        dev.off()
    }

    if(do == "fit_precision_interval"){
        loglik <- env$loglik
        interval <- as.factor(env$K)
        P <- length(env$stats_names_endo)
        S <- length(env$stats_names_exo)
        precision <- unlist(sapply(1:dim(env$mle_betas)[1],function(x) (env$K[x]*P+S)/sum(diag(env$vcov_betas[1:(env$K[x]*P+S),1:(env$K[x]*P+S),x]))))
        data <- data.frame(precision = log(precision), loglik = loglik, interval = interval)

        plot_gg <- list()
        plot_loc <- ggplot(data, aes(precision,loglik)) + geom_point(aes(colour = interval)) 
        plot_loc <- plot_loc + ggtitle("log(precision) vs. log(L) (color = number of intervals)") + xlab("log(precision)")
        plot_gg[[1]] <- plot_loc
        
        # second: arrange and save the plot
        svg("loglik_precision.svg", width=9, height=9)
        grid.arrange(grobs = plot_gg, ncol=1)
        dev.off()
    }

    if(do == "mle_transparent_lines"){
        subset_models <- subsetting_intervals(env = env, K = K)
        P <- length(env$stats_names_endo)

        for(p in 1:P){
            step_gg <- ggplot(data.frame(x = c(0,max(subset_models$widths))), aes(x))
            for(q in 1:subset_models$Q){ 
                data_q <- data.frame(elapsed_time = subset_models$widths[q,1:(subset_models$K[q]+1)], beta = c(subset_models$mle_betas[q,(1+(p-1)*subset_models$K[q])],subset_models$mle_betas[q,c((1+(p-1)*subset_models$K[q]):(p*subset_models$K[q]))]) )
                step_gg <- step_gg + geom_step(data=data_q, mapping=aes(x=elapsed_time, y=beta), alpha = 0.05, size = 5, direction="vh")      
            }
            step_gg <- step_gg + ggtitle(paste("Step functions ",env$stats_names_endo[p],sep="")) + xlab(expression(gamma)) + ylab(expression(beta)) + stat_function(fun=vdecay, args=list(type = env$true_effects$endo_effects[[env$stats_names_endo[p]]]$decay_type, pars = env$true_effects$endo_effects[[env$stats_names_endo[p]]]$pars), colour= "red")
            plot_gg <- list()
            plot_gg[[1]] <- step_gg
            # second: arrange and save the plot
            svg(paste(env$stats_names_endo[p],"_",K,"_step_functions.svg",sep=""), width=9, height=9)
            grid.arrange(grobs = plot_gg, ncol=1)
            dev.off()
        }
    }

    if(do == "beta_loglik"){ #add precision = size and colour = number of intervals
        
        P <- length(env$stats_names_endo)
        quantili <- quantile(stepwiseModelsREH$loglik,c(0)) 
        for(p in 1:P){
            plot_gg <- list()
            for(k in 1:max(env$K)){
                subset_models <- subsetting_intervals(env = env, K = (2*I(k==1) + I(k>1)*k):max(env$K) )
                select_k <- which((subset_models$widths[,2]-subset_models$widths[,1]) <= (subset_models$widths[,3]-subset_models$widths[,2]))
                #which(subset_models$loglik >= quantili) # which(apply(subset_models$widths,1,function(x) max(diff(na.omit(x))))<=2)
                # #which((subset_models$mle_betas[,1]>=0.09) & (subset_models$mle_betas[,1]<=0.11))
  
                K_loc <- subset_models$K[select_k]
                beta_k_p <- subset_models$mle_betas[select_k,]
                beta_k_p <- as.vector(sapply(1:dim(beta_k_p)[1], function(x) beta_k_p[x,k+(p-1)*K_loc[x]]))
                vcov_beta_k_p <- subset_models$vcov_betas[,,select_k]
                precision_beta_k_p <- as.vector(sapply(1:dim(vcov_beta_k_p)[3], function(x) 1/sqrt(vcov_beta_k_p[k+(p-1)*K_loc[x],k+(p-1)*K_loc[x],x])))
                width_k_p <-subset_models$widths[select_k,k+1]-subset_models$widths[select_k,k]
                alpha_feature <- (precision_beta_k_p) #*(1/width_k_p)
                data_gg <- data.frame(loglik = subset_models$loglik[select_k] ,
                beta = beta_k_p,
                K = as.factor(K_loc),
                alpha =  alpha_feature,
                width = width_k_p)
                plot_gg[[k]] <-ggplot(data_gg, aes(loglik,beta)) + geom_point(aes(colour = K)) + ggtitle(paste("interval",k,sep=" ")) #, alpha = alpha, size = width
                # second: arrange and save the plot
                svg(paste(env$stats_names_endo[p],"beta_loglik.svg",sep=" "), width=9, height=9)
                grid.arrange(grobs = plot_gg, ncol=ceiling(sqrt(max(env$K))), nrow=ceiling(sqrt(max(env$K))))
                dev.off()
            }
        }
        
    }

    if(do=="widths_weights"){
        subset_models <- subsetting_intervals(env = env, K = K)
        w_list <- getModelsWeights(input = subset_models, weights = weights)
        W <- length(w_list)
        widths_vec <- NULL
        for(q in 1:subset_models$Q){
            widths_vec <- c(widths_vec, subset_models$widths[q,2:subset_models$K[q]])
        # step_gg <- step_gg + geom_vline(xintercept=subset_models$widths[q,2:subset_models$K[q]]) + scale_colour_gradientn(colours = weights_w[q])
        }
        for(w in 1:W){
            weights_w <- as.vector(unlist(sapply(1:length(w_list[[names(w_list)[w]]]), function(x) rep(w_list[[names(w_list)[w]]][x],(subset_models$K[x]-1)))))
            data_w <- data.frame(widths = widths_vec, weights = weights_w)
            step_gg <- ggplot(data_w, aes(x = widths, y = weights)) + geom_point() + geom_vline(xintercept = env$true_effects$endo_effects[[1]]$pars[,2], color = "red")

            plot_gg <- list()
            plot_gg[[1]] <- step_gg
            svg(paste(names(w_list)[w],"widths_weights.svg",sep=" "), width=9, height=9)
            grid.arrange(grobs = plot_gg, ncol=1)
            dev.off()
        }
        return(w_list)
    }



}

