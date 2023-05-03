# Functions to perform Bayesian Model Averaging #


#' bma
#'
#' A function which return an object of class 'bma' with the estimate of the posterior distribution of effects according to specific weighting systems
#'
#' @param smm can be either an object of class 'smm' or a list of objects of class 'smm'. The second option allows to combine multiple smm objects provided that they refer to the same model specification but to different sets of stepwise memory models.
#' @param weights vector of weighting systems to be used in the averaging. "BIC" performs the method were smm's are considered to be generative models for the posterior distribution (which might be an unrealistic assumption since the posterior trend can be smooth); The "pseudoBMAplus" and "stacking" are based on the measure of ELPD. "WAIC" is an approximation of the ELPD and will be . If "all", then all the available method will be used. Also a few methods can be specified
#' @param n_knots xx
#' @param nsim xx
#'
#' @return  object of class 'bma' 
#'
#' @export
bma <- function(smm = NULL,
                weights = c("BIC","WAIC","pseudoBMAplus","stacking","all"),
                n_knots = 300,
                nsim = 5e03){  
    # 'smm' input processing
    if(!is.null(smm)){
        if(class(smm) == "list"){
            list_class <- unlist(lapply(smm,class))
            if(all(list_class) == "smm"){
                smm <- merge.smm(smm)
                # ... combine here the list of smm's objects and assign it to 'smm' again
                # don't forget also to add smm$stats_names_endo, smm$stats_names_exo
            }
            else{
                stop("argument smm must be either an object of class 'smm' or a list of objects of class 'smm'")
            }
        }
        else{
            if(class(smm) != "smm"){
                stop("argument smm must be either an object of class 'smm' or a list of objects of class 'smm'")
            }
        }
    }
    else{
        stop("missing argument 'smm'")
    }


    # 'weights' input processing
    if(!is.null(weights)){
        if(is.vector(weights)){
            weights_all <- c("BIC","WAIC","pBMAplus","SPD")
            if(any(weights %in% c("all"))){
                status <- !logical(length(weights_all))
                weights <- weights_all
                W <- length(weights)
            }
            else{
                status <- sapply(1:length(weights_all), function(x) any(weights_all[x] == weights))
                if(all(!status)) stop("Not available or mistyped argument 'weights'")
                weights <- weights_all[status]
                W <- length(weights)
            }
        }
        else{
            stop("argument weights must be a vector")
        }
    }
    else{
        stop("missing argument 'weights'")
    }

    # computing weights

    weights_list <- list() # initializing list of weights

    # BIC weights
    if(status[1]){
        if(is.null(smm$BIC)) stop("BIC vector was not found.")
        bic_w <- (-smm$BIC/2)
        weights_list$BIC <- exp(bic_w-max(bic_w)) / sum(exp(bic_w-max(bic_w))) #-(input$BIC/2-min(input$BIC/2))
    }

    # WAIC
    if(status[2]){
        if(is.null(smm$WAIC)) stop("WAIC object was not found.")
        #if(dim(na.omit(input$elpd_waic))[2]<input$Q) stop("Not all the stepwise models have an ELPD (WAIC) estimated (vcov matrix was singular).")
        elpd_waic <- apply(rbind(smm$WAIC$elpdWAIC),2,sum) 
        weights_list$WAIC <- exp((elpd_waic-max(elpd_waic))) / sum(exp((elpd_waic-max(elpd_waic)))) 
    }

    # pseudoBMAplus #temporarily changed elpd_lfo
    if(status[3]){
        if(is.null(smm$ELPD)) stop("ELPD LFO matrix was not found.")

        # Performs pseudoBMA+ by using elpd_(psis)_lfo
        elpd_lfo <- smm$ELPD$elpdLFO

        # algorithm with loo::pseudobma_weights() function
        weights_list$pseudoBMAplus <- as.vector(loo::pseudobma_weights(elpd_lfo))

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
    if(status[4]){
        if(is.null(smm$ELPD)) stop("ELPD LFO matrix was not found.")

        #stacking_opt<- function(pars,matrix_elpd)
        #{
        #   pars <- exp(pars)/sum(exp(pars))
        #   -mean(log(pars%*%matrix_elpd))
        #}

        # Performs SPD by using elpd_(psis)_lfo
        elpd_lfo <- smm$ELPD$elpdLFO
        
        # algorithm with loo::stacking_weights() function
        if(!is.null(weights_list$pseudoBMAplus)) {start_values <- weights_list$pseudoBMAplus}
        else{
        start_values <- as.vector(loo::pseudobma_weights(elpd_lfo))
        }
        weights_list$stacking <- as.vector(loo::stacking_weights(elpd_lfo,optim_control = list(start_values)))

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


    # generating from the posterior
    max_width <- max(stats::na.omit(as.vector(smm$intervals$widths)))
    elapsedTime <- seq(from = 1e-04, to = max_width-1e-04 ,length = n_knots)
    posterior_draws <- drawFromPosterior(input = smm,
                                        weights = weights_list, 
                                        nsim = nsim, 
                                        elapsed_time = elapsedTime, 
                                        max_width = max_width)

    # returning bma object
    out <- structure(list(weights = weights_list, elapsedTime = elapsedTime, draws = posterior_draws,memory_stats = smm$stats_names_endo),class="bma")
    return(out)
}

#######################################################################################
#######################################################################################
##########(START)           Methods for `bma` object           (START)#################
#######################################################################################
#######################################################################################

#' @title plot.bma
#' @rdname plot.bma
#' @description A function that returns a summary of the temporal network.
#' @param x is a \code{bma} object .. draws
#' @param ci ...
#' @param widthMax ...gamma_max
#' @param xlim ...
#' @param ylim ...
#' @param width ...
#' @param height ...
#' @param ncol ...
#' @param nrow ...
#' @param save ...
#' @param output_folder_name ...
#' @param true_pars ... list of ...
#' @method plot bma
#'
#' @export
plot.bma <- function(x, ci = FALSE, widthMax = NULL, xlim, ylim, width, height, ncol, nrow, save = TRUE, output_folder_name = "output_plots", true_pars){
    # out object list with ggplot's inside
    #p_out <- list()

    var_names <- dimnames(x$draws[[1]])[[2]]
    P <- length(x$memory_stats) 
    S <- length(var_names) - P
    weights <- names(x$draws)
    method_plots <- list()
    if(is.null(widthMax)) widthMax <- max(x$elapsedTime)
    else{
        if(widthMax > max(x$elapsedTime))
        stop(cat('widthMax must be lower or equal than',max(x$elapsedTime)))
    }

    p_out <- list()
    for(u in 1:length(var_names)){
        for(w in 1:length(weights)){
            p_out <- list()
            # calculating (2.5,50,97.5) percentiles
            #df_loc <- apply(draws$posterior_draws[[weights[w]]],c(2,3),function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE))

            # first we estimate the density objects for the intercept

            # intercept/exogenous/endogenous without memory histograms 
            if(u > P){
                data_loc <- data.frame(effect =x$draws[[weights[w]]][,u,1]) # 4th is the column of the intercept
                density_intercept <- stats::density(data_loc$effect)
                lower_upper_intercept <- as.vector(HDInterval::hdi(density_intercept, credMass=0.95))
                
                map_intercept <- as.numeric(bayestestR::map_estimate(x=data_loc$effect))
                
                p_out[[1]] <- ggplot(data=data_loc, aes(x=effect, y= ..density..)) + 
                            geom_histogram(bins=101) + 
                            geom_vline(xintercept=map_intercept,size=0.8) +  # posterior mode
                            theme_classic() +
                            labs(title = paste(weights[w],"\n\n",var_names[u],sep=""), x = expression(beta), y = "density")
                if(ci) p_out[[1]] <- p_out[[1]] + geom_vline(xintercept=c(lower_upper_intercept[1],lower_upper_intercept[2]),size=0.8,lty=2) +
                                scale_x_continuous(breaks=c(lower_upper_intercept[1],map_intercept,lower_upper_intercept[2]),
                                labels=as.character(round(c(lower_upper_intercept[1],map_intercept,lower_upper_intercept[2]),3))) # ci quantiles
            }
            # for loop over endogenous statistics
            #for(u in 1:(dim(draws$posterior_draws[[weights[w]]])[2]-1)){

            else{
                df_loc <- x$draws[[weights[w]]][,u,]
                df_data_loc <- apply(df_loc,2,function(x) as.vector(HDInterval::hdi(x, credMass=0.95)))
                data_median <- apply(df_loc,2,function(x) stats::quantile(x,c(0.5),na.rm=TRUE))
                data_mean <- apply(df_loc,2,function(x) mean(x,na.rm=TRUE))
                map_loc <- apply(df_loc,2,function(x) as.numeric(bayestestR::map_estimate(x))) 
                data_loc <- data.frame(time = stats::na.omit(x$elapsedTime), 
                                        lb = df_data_loc[1,],
                                        mode=map_loc, 
                                        ub = df_data_loc[2,],
                                        median = data_median,
                                        mean= data_mean)
                rm(df_data_loc)
                p_loc <- ggplot(data=data_loc, aes(x=time)) +  #effect
                        geom_line(aes(y=mode,color="red"),size = 0.8) + # trend is transparent
                        geom_line(aes(y=median,color="black"),size = 0.8) +
                        geom_line(aes(y=mean,color="green"),size = 0.8) +
                        geom_hline(yintercept=0,size=0.8,col="gray",lty=2,alpha=0.6) +
                        theme_classic() +
                      #  scale_x_continuous(breaks=round(seq(0,widthMax,length=5)),labels=as.character(round(seq(0,widthMax,length=5))))# +
                        labs(title = var_names[u], x = "transpired time", y = bquote(paste(beta[.(var_names[u])]))) #+
       #                 coord_cartesian(xlim=xlim)

                if(ci) p_loc <- p_loc + geom_ribbon(data=data_loc,aes(ymin=lb,ymax=ub),alpha=0.3) + scale_x_continuous(breaks = pretty)
                p_out[[1]] <- p_loc
            }
       if(u <= P){method_plots <- do.call(c, list(method_plots, p_out)) }
        }
    }


    if(save){
    name <- "all_methods_overview"
    if(ci) name <- paste(name, "ci", sep = "_")
    svg(paste(getwd(),"./",output_folder_name,"/",name,".svg", sep = ""), width=width, height=height) # or width=6, height=10
    grid.arrange(grobs = method_plots, ncol=ncol, nrow=nrow, as.table=FALSE)
    dev.off()
    }

    return(method_plots)

}


####
##summary.bma  <- function(){

##}

#######################################################################################
#######################################################################################
##########(END)             Methods for `bma` object             (END)#################
#######################################################################################
#######################################################################################




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

    stats_names <- c(input$stats_names_endo,input$stats_names_exo) #temporarily added by the user to the input object
    P <- length(input$stats_names_endo)
    S <- length(input$stats_names_exo)
    U <- P + S
    Q <- input$Q
    W <- length(weights)
      
    out <- list()
    ## stepwise models ##
    K <- input$intervals$K # vector of number of intervals
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
                    which_step <- find_step(x = elapsed_time[i], widths = input$intervals$widths[q,1:(K[q]+1)])
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
