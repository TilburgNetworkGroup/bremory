#' getStepwiseModels
#'
#' A function which run the estimation of Q models and updates the environment 'stepwiseModelsREH' 
#'
#' @param env environment where the user is currently working
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param K_range range as to the number of intervals that are going to be generated.
#' @param max_width maximum for width
#' @param nsim_per_K how many simulations per number of intervals
#' @param L number of events to use as history to estimate the predictive pointwise density of future events (set to M-10 by default)
#' @param A_SAP A-Steps-Ahead-Predictions is the number of steps ahead we want consider in the prediction (set to 1 by default)
#' @param file_name a string indicating the name of the file .RData that will save results each 20 models.
#' @param WAIC default is TRUE and is needed where WAIC weights have to be estimated
#' @param n_cores number of threads to create
#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModelsREH'
#'
#' @export
getStepwiseModels <- function(env = globalenv(),
                              effects = NULL, 
                              K_range = c(2:6), 
                              max_width = NULL, 
                              nsim_per_K = 2e02, 
                              L = NULL, 
                              A_SAP = 1, 
                              file_name = "simulation_0", 
                              WAIC = TRUE, 
                              n_cores = 10,
                              intervals = c("default","increasing","equal"))
{
    # check for initializeREH environment
    if(is.null(env$initializeREH)){stop("initializeREH() function must be run first.")}

    # create a new environment 'statisticsREH' (the largest in the estimation process and continuously updated according to the characteristic of the stepswise model
    # at a specific iteration)
    statisticsREH(env = env, effects = effects, K_range = K_range)

    # create environment 'stepwiseModelsREH' where, widths are stored as well as the output of each model (information useful to the estimation stage)
    env$stepwiseModelsREH <- new.env()

    # storing dimensions and variables name
    env$stepwiseModelsREH$Q <- length(K_range)*nsim_per_K
    env$stepwiseModelsREH$K <- rep(K_range, each = nsim_per_K)
    env$stepwiseModelsREH$L <- max(which(env$initializeREH$t<(max(env$initializeREH$t)-max_width))) # M-10 was the FIXED VALUE (previously)
    env$stepwiseModelsREH$stats_names_endo <- env$statisticsREH$endo_effects # or effects$endogenous
    env$stepwiseModelsREH$stats_names_exo <- dimnames(env$statisticsREH$exogenous_stats)[[3]]

    ## BEGIN rem supplist 
    supplist <- list() 
    supplist[[1]]<-matrix(TRUE, nrow = env$initializeREH$M, ncol = env$initializeREH$n_dyads)
    ## END rem supplist

    # allocating space for each model result (Q models in total)
    env$stepwiseModelsREH$mle_betas <- matrix(NA, nrow = env$stepwiseModelsREH$Q, ncol = (( max(K_range) * env$statisticsREH$P) + env$statisticsREH$S)) 
    env$stepwiseModelsREH$vcov_betas <- array(NA , dim = c((( max(K_range) * env$statisticsREH$P) + env$statisticsREH$S),
    (( max(K_range) * env$statisticsREH$P) + env$statisticsREH$S), env$stepwiseModelsREH$Q))
    env$stepwiseModelsREH$loglik <- env$stepwiseModelsREH$BIC <- rep(NA, env$stepwiseModelsREH$Q)
    env$stepwiseModelsREH$reestimations_psis <- rep(NA, env$stepwiseModelsREH$Q)
    env$stepwiseModelsREH$elpd_lfo <- matrix(NA, nrow = (env$initializeREH$M - env$stepwiseModelsREH$L), ncol = env$stepwiseModelsREH$Q)
    env$stepwiseModelsREH$k_psis <- matrix(NA,nrow=((env$initializeREH$M - env$stepwiseModelsREH$L)-1),ncol=env$stepwiseModelsREH$Q)

    if(WAIC){
      # for multi-steps models  :env$initializeREH$M
      env$stepwiseModelsREH$elpd_waic <- matrix(NA, nrow = (env$initializeREH$M-env$stepwiseModelsREH$L), ncol= env$stepwiseModelsREH$Q) #env$stepwiseModelsREH$L
      env$stepwiseModelsREH$p_waic <- matrix(NA, nrow = (env$initializeREH$M-env$stepwiseModelsREH$L), ncol= env$stepwiseModelsREH$Q) #env$stepwiseModelsREH$L
    }
    ###### START FIXED VALUES (MAYBE TO CHANGE) #####
    env$stepwiseModelsREH$n_sims_is <- 1e03
    env$stepwiseModelsREH$tau <- 0.7
    ###### END FIXED VALUES (MAYBE TO CHANGE) #####

    # generating random widths according to K_range and nsim_per_K parameters
    generateWidths(env = env, K_range = K_range, nsim_per_K = nsim_per_K, max_width = max_width)
    # shuffle intervals
    shuffled_indices <- sample(x = c(1:env$stepwiseModelsREH$Q) ,size = env$stepwiseModelsREH$Q ,replace = FALSE)
    env$stepwiseModelsREH$K <- env$K_reh #env$stepwiseModelsREH$K[shuffled_indices]
    env$stepwiseModelsREH$widths <- env$widths_reh #env$stepwiseModelsREH$widths[shuffled_indices,]

    # fixing interval to the whole event sequence where there is no dyadic/triadic statistics
    if(length(names(env$statisticsREH$binaryREH)) == 0) {
      env$stepwiseModelsREH$widths <- rbind(c(0,Inf))
      env$stepwiseModelsREH$K <- 1
    }

    #### START PARALLELIZATION SETTINGS
    # set Cores for parallelization purpose
    # n_cores <- floor(detectCores()) - 2
    # assign object 'n_cores' to environment initializeREH
    env$initializeREH$n_cores <- n_cores
    #### END PARALLELIZATION SETTINGS
    env$statisticsREH$edgelist <- data.matrix(convertEdgelist(edgelist = env$initializeREH$edgelist, 
                                                              riskset = env$initializeREH$riskset,  
                                                              actors = env$initializeREH$actors,
                                                              rem = TRUE))
    env$statisticsREH$counts <- new.env()
    

    ##
    #env$stepwiseModelsREH$widths[1,] <- c(0,0.2,1.3,3.5,7,NA) #only for the stepwise case (change interval widths with the true ones)
    #env$stepwiseModelsREH$K[1] <- 4
    for(q in 1:dim(env$stepwiseModelsREH$widths)[1]) 
    {  
        start <- Sys.time()     
        env$statisticsREH$widths_q <-  na.omit(env$stepwiseModelsREH$widths[q,]) 
        env$statisticsREH$K_q <- env$stepwiseModelsREH$K[q] 

        #### START PARALLELIZATION WITH 'PARALLEL' PKG (by using parLapply() function)
        ## list of range indices in order to go parallel
        #lbs_ubs_list <- list()
        #lbs_ubs <- intervals_loc
        #seq_rows<-trunc(seq(1,dim(lbs_ubs)[1],length=(n_cores+1)))
        #for(i in 2:(n_cores+1)){lbs_ubs_list[[i-1]] <- lbs_ubs[c((seq_rows[i-1]+1*I((i-1)>1)):seq_rows[i]),]}
        ## defining cluster and exporting objects from global environment
        #cl <- makeCluster(n_cores) 
        #loc_env <- new.env()
        #loc_env$binaryREH <- as.matrix(env$initializeREH$binaryREH)
        #loc_env$lbs_ubs_list <- lbs_ubs_list
        #clusterExport(cl = cl, envir = loc_env, c("binaryREH","lbs_ubs_list")) 
        #clusterExport(cl = cl, c("getCounts"))
        #clusterEvalQ(cl = cl, library(Matrix)) # exporting dgCMatrix class for sparse matrices
        ## starting parallelization
        #counts_new <- parLapply(cl=cl, as.list(1:(n_cores)), function(x){
        #getCounts(binaryREH = binaryREH, lbs_ubs = lbs_ubs_list[[x]])
        #})
        #stopCluster(cl)
        ## end parallelization
        ## concatenating matrices of counts coming from the threads (got with the parallelization)
        #counts <- concatCount(env = globalenv(), matrix_list = counts_new)
        ## assigning object counts to environment object insid statisticsREH
        #env$statisticsREH$counts <- counts 
        #### END PARALLELIZATION WITH 'PARALLEL' PKG (by using parLapply() function)

        #getCountsOMP (OpenMP parallelization in C++)
        if(length(env$statisticsREH$endo_effects)>0){
          # getIntervals() will find per each time point the lower and the upper bound for each interval and will automatically assign it to the object 'intervals_backwards'
          getIntervals(env = globalenv(), widths = env$statisticsREH$widths_q) 
          env$statisticsREH$intervals_backward[(env$statisticsREH$K_q+1),1:2] <- c(0,0) #only for the first interval at the second time point
          intervals_loc <- unique(env$statisticsREH$intervals_backward[,1:2])

          for(l in 1:length(names(env$statisticsREH$binaryREH))){
            env$statisticsREH$counts[[names(env$statisticsREH$binaryREH)[l]]] <- getCountsOMP(binaryREH = env$statisticsREH$binaryREH[[names(env$statisticsREH$binaryREH)[l]]], 
                                                                                              lbs_ubs = intervals_loc, 
                                                                                              n_cores = env$initializeREH$n_cores)
          }

          indices <- getCountsIndex(intervals = env$statisticsREH$intervals_backward, 
                                    counts = cbind(env$statisticsREH$counts[[names(env$statisticsREH$binaryREH)[1]]][,1:2],
                                    c(0:(dim(env$statisticsREH$counts[[names(env$statisticsREH$binaryREH)[1]]])[1]-1))))
          env$statisticsREH$intervals_backward[,3] <-  indices

          rm(intervals_loc,indices)

          ### calculating effects ###
          endogenous_stats <- getEndoEffects(env = env)

          # create names for endongenous and exogenous variables
          rep_names <- rep(c(env$statisticsREH$endo_effects), each = env$statisticsREH$K_q)
          rep_int_index <- rep(1:env$statisticsREH$K_q,times = env$statisticsREH$P)
          env$statisticsREH$names_endogenous <- sapply((1:(env$statisticsREH$P*env$statisticsREH$K_q)),
              function(x) { paste(rep_names[x],rep_int_index[x],sep="_")
          })
          dimnames(endogenous_stats)[[3]] <- env$statisticsREH$names_endogenous


          rm(list=ls(envir = env$statisticsREH$counts), envir = env$statisticsREH$counts)
          env$statisticsREH$intervals_backward <- NULL
        }
        else{endogenous_stats <- NULL}
        
        stats <- abind(endogenous_stats, env$statisticsREH$exogenous_stats, along = 3)

        ##### relevent::rem function #####   
        model_loc <- relevent::rem(eventlist = env$statisticsREH$edgelist, 
                                  statslist = stats,
                                  supplist = supplist,
                                  timing = "interval",
                                  estimator = "MLE") # "BSIR"    
        # WAIC routine 
        if(WAIC){
          # ROUTINE FOR THE COMPLETE MODEL 
          start_waic <- Sys.time()  
          # generate from posterior multivariate normal approximation
          Sigma <- tryCatch(as.matrix(solve(model_loc$hessian)),
                                    error = function(error_message) {matrix(-1,nrow=dim(stats)[2],ncol=dim(stats)[2])})
          time_points <- (env$stepwiseModelsREH$L+1):env$initializeREH$M           


          interevent_time <- diff(c(0,env$initializeREH$t))[time_points]  
          stats_waic <- aperm(stats, perm = c(2,3,1))  
          events <- env$statisticsREH$binaryREH$dyadicREH[time_points,]                     
          if(Sigma[1,1] != (-1)){
            post_pars <- t(mvtnorm::rmvnorm(n = env$stepwiseModelsREH$n_sims_is, mean = model_loc$coef, sigma = Sigma)) # dim = [pars x n_sim_is]

            #lpd <- matrix(NA,nrow=env$stepwiseModelsREH$n_sims_is,ncol=length(time_points))
            #for(m in 1:length(time_points))
            #{
            #event <- numeric(env$initializeREH$n_dyads)
            #event[env$statisticsREH$edgelist[time_points[m],1]] <- 1
            #lpd[,m] <- apply(post_pars,2,function(y) lpd(pars = y,
            #                                                stats = stats[time_points[m],,],
            #                                                event = event,
            #                                                interevent_time = env$initializeREH$t[time_points[m]]-env$initializeREH$t[time_points[m]-1]))
          	#}  
            lpd <- bremory::lpdWAIC(pars = post_pars, 
                                    stats = stats_waic[,,time_points], 
                                    events = t(events), 
                                    interevent_time = interevent_time) 
            elpd_waic_obj <- loo::waic(lpd)
            env$stepwiseModelsREH$p_waic[,q] <- elpd_waic_obj$pointwise[,2] #lpd_and_p_waic[,2]
            env$stepwiseModelsREH$elpd_waic[,q] <- elpd_waic_obj$pointwise[,1] #lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
            
            #lpd_and_p_waic <- getWAIC(pars = post_pars, stats = stats_waic[,,time_points], events = t(events), interevent_time = interevent_time)
            #for(m in 1:length(time_points))
            #{
            #  # (lpd_m - p_waic_m)
            #  event <- numeric(env$initializeREH$n_dyads)
            #  event[env$statisticsREH$edgelist[time_points[m],1]] <- 1
            #  p_m <- apply(post_pars,2,function(y) exp(lpd(pars = y,
            #                                                stats = stats[time_points[m],,],
            #                                                event = event,
            #                                                interevent_time = env$initializeREH$t[time_points[m]]-env$initializeREH$t[time_points[m]-1])))
            #  lpd_m <- log(mean(p_m))
            #  p_waic_m <- var(log(p_m))
            #  env$stepwiseModelsREH$p_waic[m,q] <- p_waic_m 
            #  env$stepwiseModelsREH$elpd_waic[m,q] <- lpd_m - p_waic_m
            #}
          }
          end_waic<- Sys.time()  
          cat("\n WAIC routine executed in.... " ,end_waic-start_waic, attr(end_waic-start_waic,'units'))
        }
                  
        ### run psis-lfo ###
        ### PSIS-LFO-CV MSAP (1-step-ahead-predictions) ###
        psis_1_sap <- PSIS_LFO_CV_1_SAP(n_dyads = env$initializeREH$n_dyads,
                                        eventlist = env$statisticsREH$edgelist, 
                                        statslist = stats,
                                        t = env$initializeREH$t,
                                        M = env$initializeREH$M, 
                                        L = env$stepwiseModelsREH$L, 
                                        n_sim_is = env$stepwiseModelsREH$n_sims_is,
                                        tau = env$stepwiseModelsREH$tau,
                                        cores = env$initializeREH$n_cores)                   
        end <- Sys.time()
        # print message 
        cat("\n ",file_name,": Iteration ",q," with ",env$statisticsREH$K_q," intervals and ", psis_1_sap$counts, "  no. re-estimations completed in --> ", end-start, attr(end-start,'units'),
        ".\n")  
        # print message
        

        #### store results ####

        n_pars_loc <- dim(stats)[3]
        env$stepwiseModelsREH$mle_betas[q,c(1:n_pars_loc)] <- as.vector(model_loc$coef)
        env$stepwiseModelsREH$vcov_betas[c(1:n_pars_loc),c(1:n_pars_loc),q] <- tryCatch(as.matrix(solve(model_loc$hessian)),
                                    error = function(error_message) {matrix(-1,nrow=n_pars_loc,ncol=n_pars_loc)})
        env$stepwiseModelsREH$loglik[q] <- model_loc$loglik
        env$stepwiseModelsREH$BIC[q] <- model_loc$BIC 
        env$stepwiseModelsREH$reestimations_psis[q] <- psis_1_sap$counts               
        env$stepwiseModelsREH$elpd_lfo[,q] <- psis_1_sap$elpd_lfo
        env$stepwiseModelsREH$k_psis[,q] <- psis_1_sap$k

        rm(endogenous_stats,n_pars_loc,model_loc,stats,psis_1_sap) 

        if(q%%5==0)
        {
            cat('\n ... Saving partial results ... \n')
            save(initializeREH,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

            # save intervals table
            partial_result_1 <- data.frame(table(env$stepwiseModelsREH$K[1:q]),iteration = q)
            write.csv(partial_result_1, file = paste(file_name,".csv",sep=""))
        }
    }

    cat('\n ... Saving final results ... \n')
    save(initializeREH,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

    rm(statisticsREH, envir = env)
}


#' PSIS_LFO_CV_1_SAP
#'
#' A function which run th PSIS with LFO CV by using 1-Step-Ahead-Predictions 
#'
#' @param n_dyads ...
#' @param eventlist ...
#' @param statslist ...
#' @param t ...
#' @param M ...
#' @param L ...
#' @param n_sim_is ...
#' @param tau ...
#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModelsREH'
#' 
#' @export
PSIS_LFO_CV_1_SAP <- function(n_dyads,eventlist,statslist,t,M,L,n_sim_is,tau,cores)
{
  J_i <- L+1 # set of time points contributing in the calculation of the importance ratios. 
  time_points <- ((L+1):M) 
  
  log_p_J_i <- matrix(NA, nrow = length(time_points), ncol = n_sim_is)
  elpd_lfo <- vector(mode = "numeric", length = length(time_points))
  count_no_estimation <- 0 # know how many times the model is estimated considering a larger subset than the previous step
  
  supplist_L <- list() 
  supplist_L[[1]]<-matrix(TRUE, L, n_dyads)

  # helper functions taken from the tutorial: https://mc-stan.org/loo/articles/loo2-lfo.html
  # more stable than log(sum(exp(x))) 
  log_sum_exp <- function(x) {
    max_x <- max(x)  
    max_x + log(sum(exp(x - max_x)))
  }

  # more stable than log(mean(exp(x)))
  log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
  }

  # estimate the model subsetting in 1:L (first L events)
  #### relevent::rem ####  
  model_subset <- relevent::rem(eventlist = eventlist[1:L,],
                                statslist = statslist[1:L,,],
                                supplist = supplist_L,
                                timing = "interval",
                                estimator = "MLE")
                               
  #### relevent::rem ####

  # store MLEs and standard errors
  mle_subset <- model_subset$coef
  vcov_subset <- as.matrix(solve(model_subset$hessian))
  #draw a sample from a multivariate normal distribution (approximation of the posterior distribution)
  pars_draws <- t(mvtnorm::rmvnorm(n = n_sim_is, mean = mle_subset, sigma = vcov_subset)) # dim = [pars x n_sim_is]

  model_subset_old <- model_subset 
  mle_subset_old <- mle_subset
  vcov_subset_old <- vcov_subset
  pars_draws_old <- pars_draws

  # first step ahead (L+1) is predicted with the exact LFO
  event <- numeric(n_dyads)
  event[eventlist[time_points[1],1]] <- 1
  log_p_J_i[1,] <- logpJi(pars = pars_draws,
                      stats = statslist[time_points[1],,],
                      event = event,
                      interevent_time = t[time_points[1]]-t[time_points[1]-1]) 
  # exact LFO prediction                                   
  elpd_lfo[1] <-  log_mean_exp(log_p_J_i[1,]) #log(mean(exp(log_p_J_i[1,]))) 

  i_star <- 1 # everything starts from i = L+1
  i <- 2
  k_vec <- NULL
  # progress bar
  c('\n Progress of PSIS-LFO-1SAP (latest algorithm')
  pb <- txtProgressBar(min = 1, max = length(time_points), style = 3)
  while(i <= length(time_points)) 
  {
    J_i <- c(J_i,time_points[i])
    event <- numeric(n_dyads)
    event[eventlist[time_points[i],1]] <- 1
    log_p_J_i[i,] <- logpJi(pars = pars_draws,
                      stats = statslist[time_points[i],,],
                      event = event,
                      interevent_time = t[time_points[i]]-t[time_points[i]-1])    

    log_ratios <- as.vector(apply(rbind(log_p_J_i[i_star:i,]),2,sum))

    psis_smooth <- suppressWarnings(loo::psis(log_ratios = log_ratios, cores = cores))
    k <- psis_smooth$diag$par   
    k_vec <- c(k_vec,k)
    if(k < tau)
    {   
      psis_weights <-  weights(psis_smooth, normalize = TRUE)[, 1]

      elpd_lfo[i] <- log_sum_exp(psis_weights + log_p_J_i[i,]) #log(sum(exp(log_p_J_i[i,])*psis_weights))         
    } else
    {
      cat('... reestimating \n')
      i_star <- i # useful to filter out from log_p_J_i those events that become part of the fitted model
      J_i <- time_points[i] 

      supplist_L <- list() 
      supplist_L[[1]] <- matrix(TRUE, (time_points[i]-1), n_dyads)

      #re-estimating the model until (time_points[i]-1)
      #### relevent::rem ####                        
      model_subset <- relevent::rem(eventlist = eventlist[1:(time_points[i]-1),],
                                    statslist = statslist[1:(time_points[i]-1),,],
                                    supplist = supplist_L,
                                    timing = "interval",
                                    estimator = "MLE")
      #### relevent::rem ####

      mle_subset <- model_subset$coef
      vcov_subset <- as.matrix(solve(model_subset$hessian))
      pars_draws <- t(mvtnorm::rmvnorm(n = n_sim_is, mean = mle_subset, sigma = vcov_subset)) # dim = [pars x n_sim_is]
      
      # estimating with exact LFO (using the new "pars_draws" matrix)
      event <- numeric(n_dyads)
      event[eventlist[time_points[i],1]] <- 1
      log_p_J_i[i,] <- logpJi(pars = pars_draws,
                      stats = statslist[time_points[i],,],
                      event = event,
                      interevent_time = t[time_points[i]]-t[time_points[i]-1]) 
      elpd_lfo[i] <- log_mean_exp(log_p_J_i[i,]) #log(mean(exp(log_p_J_i[i,]))) 
      count_no_estimation <- count_no_estimation + 1                                      
    }
    i <- i + 1
    setTxtProgressBar(pb, i)
  }

  return(list(elpd_lfo = elpd_lfo, counts = count_no_estimation, k = k_vec))
}

