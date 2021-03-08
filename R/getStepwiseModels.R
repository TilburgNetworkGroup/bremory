#' getStepwiseModels
#'
#' A function which run the estimation of Q models and updates the environment 'stepwiseModelsREH' 
#'
#' @param reh either a 'reh' object (see remify::reh() function) or the edgelist
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param K range as to the number of intervals that are going to be generated.
#' @param max_width maximum for width
#' @param min_diff minimum width (diff will be changed with width)
#' @param nsim_per_K how many simulations per number of intervals
#' @param L number of events to use as history to estimate the predictive pointwise density of future events (set to M-10 by default)
#' @param A_SAP A-Steps-Ahead-Predictions is the number of steps ahead we want consider in the prediction (set to 1 by default)
#' @param file_name a string indicating the name of the file .RData that will save results each 20 models.
#' @param ELPD default is FALSE and is TRUE if ELPD-LFO has to be calculated (PSIS approximation is used [ref. to literature] for one-step-ahead predictions). This procedure might slow down the algorithm.
#' @param WAIC default is TRUE and is needed where WAIC weights have to be estimated
#' @param n_threads number of threads to create
#' @param first_event integer indicating the event at which start the event sequence and statistic in the estimation stage
#' @param last_event integer indicating the event at which stop the event sequence and statistic in the estimation stage
#' @param interval_input object "stepwise intervals" with widths, K and types
#' @param env environment where the user is currently working
#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModelsREH'
#'
#' @export
getStepwiseModels <- function(reh = NULL,
                              effects = NULL, 
                              K = c(2:6), 
                              max_width = NULL, 
                              min_diff = NULL,
                              nsim_per_K = 101,
                              L = NULL, 
                              A_SAP = 1, 
                              file_name = "simulation_0", 
                              ELPD = FALSE,
                              WAIC = TRUE, 
                              n_threads = parallel::detectCores()-2,
                              first_event = NULL,
                              last_event = NULL,
                              interval_input = NULL,
                              env = globalenv(),
                              ...)
{
    # check for input 'reh'
    if(is.null(reh)) stop("Input 'reh' cannot be left NULL.")
    else{
        if(class(reh)!="reh"){
            if(is.data.frame(reh) | is.matrix(reh)){
                reh <- remify::reh(edgelist = reh, ...)
            }
            else{
                stop("Input 'reh' must be either a 'reh' object (from 'remify' package) or the event sequence (as matrix or data.frame).")
            }
        }
    }

    # if events have weights/intensities (remify should process this stage and return a rehBinary with weights instead of ones)
    #if(!is.null(env$initializeREH$weights))
    #{
    #    weights_loc <- matrix(rep(env$initializeREH$weights,env$initializeREH$n_dyads), nrow = env$initializeREH$M,ncol = env$initializeREH$n_dyads)
    #    env$statisticsREH$binaryREH$dyadicREH <- weights_loc*env$statisticsREH$binaryREH$dyadicREH
    #}


    # create a new environment 'statisticsREH' (the largest in the estimation process and continuously updated according to the characteristic of the stepswise model
    # at a specific iteration)
    statisticsREH(reh = reh, effects = effects, env = env)

    # create environment 'stepwiseModelsREH' where, widths are stored as well as the output of each model (information useful to the estimation stage)
    env$stepwiseModelsREH <- new.env()

    # generating random widths according to K and nsim_per_K parameters
    if(is.null(interval_input)){
      env$stepwiseModelsREH$widths <- matrix(NA, nrow = length(K)*nsim_per_K, ncol = (max(K)+1)) # length(K)*(nsim_per_K-1)*2+length(K) is specific of the case where  for example we want to generate 150 increasing 150 decreasing and 1 equal size , thus for K=c(3,4,5) 3*(151-1)+3 where 151 is nsim_per_K
      env$stepwiseModelsREH$K <- rep(K,each = nsim_per_K)
      env$stepwiseModelsREH$interval_type <- rep(c("equal",rep("increasing",(nsim_per_K-1)/2),rep("decreasing",(nsim_per_K-1)/2)),length(K))
      # shuffling order [now] : this is done not to increasingly overload the cpu and/or the ram (in this way sometimes a complex model is estimated, other times a much simpler one is optimized)
      shuffled_indices <- sample(x=1:dim(env$stepwiseModelsREH$widths)[1],size=dim(env$stepwiseModelsREH$widths)[1])
      env$stepwiseModelsREH$K <- env$stepwiseModelsREH$K[shuffled_indices]
      env$stepwiseModelsREH$interval_type <- env$stepwiseModelsREH$interval_type[shuffled_indices]
    }
    else{
      #if(class(interval_input)!="stepwise intervals") stop("input argument 'interval_input' must be of class 'stepwise intervals'. ")
      env$stepwiseModelsREH$widths <- interval_input$widths
      env$stepwiseModelsREH$K <- interval_input$K
      env$stepwiseModelsREH$interval_type <- interval_input$type
      max_width <- max(na.omit(as.vector(env$stepwiseModelsREH$widths)))
    }
    env$stepwiseModelsREH$widths_rejected <- env$stepwiseModelsREH$K_rejected <- env$stepwiseModelsREH$type_rejected <- NULL

    env$stepwiseModelsREH$first_event <- first_event
    env$stepwiseModelsREH$last_event <- last_event

    if(!is.null(first_event) & !is.null(last_event)){
      if(first_event > last_event){
        stop("first event integer value must be lower than the last event iteger value.")
      }
    }
    else{
      env$stepwiseModelsREH$first_event <- first_event <- 1
      env$stepwiseModelsREH$last_event <- last_event <- reh$M
    }

    ### TO REMOVE ###
    # print we could also remove
    print(c(env$stepwiseModelsREH$first_event,env$stepwiseModelsREH$last_event,max_width))
    ### TO REMOVE ###
    
    # storing dimensions and variables name
    env$stepwiseModelsREH$Q <- dim(env$stepwiseModelsREH$widths)[1]    
    # M-L are the number of events to predict in WAIC and/or ELPD and their time span correspond to the maximum time width observed in the intervals
    if(is.null(L)){
      env$stepwiseModelsREH$L <- c(first_event:last_event)[max(which(reh$edgelist$time[first_event:last_event]<(max(reh$edgelist$time[first_event:last_event])-max_width)))] 
    }
    else{
      if(L>=reh$M) stop("'L' must be smaller than the length of the event sequence.")
      else{
        env$stepwiseModelsREH$L <- L
      }
    }
    reh$edgelist$time <- cumsum(as.numeric(reh$intereventTime)) # with remstimate::remstimate this column won't be part of the estimation stage, instead the intereventTime will be used. Therefore this line WILL BE REMOVED

    # names of endogenous and exogenous statistics
    env$stepwiseModelsREH$stats_names_endo <- env$statisticsREH$endo_effects 
    env$stepwiseModelsREH$stats_names_exo <- dimnames(env$statisticsREH$exogenous_stats)[[3]]

    ## BEGIN rem supplist 
    supplist <- list() 
    supplist[[1]]<-matrix(TRUE, nrow = length(first_event:last_event), ncol = reh$D)
    ## END rem supplist

    # allocating space for each model result (Q models in total)
    env$stepwiseModelsREH$mle_betas <- matrix(NA, nrow = env$stepwiseModelsREH$Q, ncol = (( max(K) * env$statisticsREH$P) + env$statisticsREH$S)) 
    env$stepwiseModelsREH$vcov_betas <- array(NA , dim = c((( max(K) * env$statisticsREH$P) + env$statisticsREH$S),
    (( max(K) * env$statisticsREH$P) + env$statisticsREH$S), env$stepwiseModelsREH$Q))
    env$stepwiseModelsREH$loglik <- env$stepwiseModelsREH$BIC <- rep(NA, env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$loglik_no_endo

    ## ELPD LFO memory allocation
    #if(ELPD){
    #env$stepwiseModelsREH$reestimations_psis <- rep(NA, env$stepwiseModelsREH$Q)
    #env$stepwiseModelsREH$elpd_lfo <- matrix(NA, nrow = (length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L)) , ncol = env$stepwiseModelsREH$Q) 
    #env$stepwiseModelsREH$k_psis <- matrix(NA, nrow = nrow = ((length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L))-1), ncol = env$stepwiseModelsREH$Q) 
    #}

    if(WAIC){
      env$stepwiseModelsREH$elpd_waic <- matrix(NA, nrow = (length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L)), ncol= env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$elpd_waic_no_endo
      env$stepwiseModelsREH$p_waic <- matrix(NA, nrow = (length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L)), ncol= env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$p_waic_no_endo
    }

    ###### START FIXED VALUES (MAYBE TO CHANGE) #####
    env$stepwiseModelsREH$n_sims_is <- 1e03
    env$stepwiseModelsREH$tau <- 0.7
    ###### END FIXED VALUES (MAYBE TO CHANGE) #####

    #### START PARALLELIZATION SETTINGS
    # set Cores for parallelization purpose
    if(is.null(n_threads)) n_threads <- 1
    else{
      if(n_threads > (floor(parallel::detectCores()) - 2)) stop(cat("'n_threads' is recommended to be set at most to ", (floor(parallel::detectCores()) - 2),"."))
    }
    #### END PARALLELIZATION SETTINGS

    # converting the edgelist to a suitable format for relevent::rem()
    env$statisticsREH$edgelist <- data.matrix(convertToReleventEdgelist(reh = reh, risksetCube = reh$risksetCube))

    q <- 1
    q_rej <- 0
    cat("\n Starting estimating the bag of models ! \n")
    while(q <= env$stepwiseModelsREH$Q)
    { 
        model_status <- FALSE
        while(!model_status){

          if(is.null(interval_input)){
            env$statisticsREH$widths_q <-  generateWidths(K= env$stepwiseModelsREH$K[q],
                                                        min_diff = min_diff,
                                                        max_width= max_width,
                                                        intervals= env$stepwiseModelsREH$interval_type[q]) #na.omit(env$stepwiseModelsREH$widths[q,])
            env$stepwiseModelsREH$widths[q,1:(env$stepwiseModelsREH$K[q]+1)] <- env$statisticsREH$widths_q
            env$statisticsREH$K_q <- env$stepwiseModelsREH$K[q] 
          }
          else{
            env$statisticsREH$widths_q <-  env$stepwiseModelsREH$widths[q,1:(env$stepwiseModelsREH$K[q]+1)]
            env$statisticsREH$K_q <- env$stepwiseModelsREH$K[q] 
          }

          
          start <- Sys.time()

          #getCountsOMP (OpenMP parallelization in C++)
          if(length(env$statisticsREH$endo_effects)>0){
            # getIntervals() will find per each time point the lower and the upper bound for each interval and will automatically assign it to the object 'intervals_backwards'
            getIntervals(env = globalenv(), widths = env$statisticsREH$widths_q, time = reh$edgelist$time, M = reh$M, K_q = env$statisticsREH$K_q) # reh$edgelist$time

            ### TO REMOVE ###
            #env$statisticsREH$intervals[(env$statisticsREH$K_q+1),1:2] <- c(0,0) #only for the first interval at the second time point
            ### TO REMOVE ###

            intervals_loc <- unique(env$statisticsREH$intervals[,1:2])
            
            env$statisticsREH$counts <- getCountsOMP(binaryREH = reh$rehBinary, 
                                                    lbs_ubs = intervals_loc, 
                                                    n_threads = n_threads)
            

            indices <- getCountsIndex(intervals = env$statisticsREH$intervals, 
                                      counts = cbind(env$statisticsREH$counts[,1:2],
                                      c(0:(dim(env$statisticsREH$counts)[1]-1)))
                                      )
            env$statisticsREH$intervals[,3] <-  indices

            rm(intervals_loc,indices)

            ### calculating effects ###
            endogenous_stats <- getEndoEffects(env = env,
                                              M = reh$M, 
                                              D = reh$D,
                                              time = reh$edgelist$time,
                                              edgelist = data.matrix(reh$edgelist), 
                                              risksetMatrix = reh$risksetMatrix, 
                                              risksetCube0 = reh$risksetCube[,,1], # still one event type is used (need to make a change)
                                              n_threads = n_threads) 
            # create names for endongenous and exogenous variables
            rep_names <- rep(c(env$statisticsREH$endo_effects), each = env$statisticsREH$K_q)
            rep_int_index <- rep(1:env$statisticsREH$K_q,times = env$statisticsREH$P)
            env$statisticsREH$names_endogenous <- sapply((1:(env$statisticsREH$P*env$statisticsREH$K_q)),
                function(x) { paste(rep_names[x],rep_int_index[x],sep="_")
            })
            dimnames(endogenous_stats)[[3]] <- env$statisticsREH$names_endogenous


            env$statisticsREH$counts  <- NULL
            env$statisticsREH$intervals <- NULL
          }
          else{endogenous_stats <- NULL}
          stats <- abind::abind(endogenous_stats, env$statisticsREH$exogenous_stats, along = 3)
          
          ##### relevent::rem function #####  
          # estimating model with relevent package #
          stats_rem <- stats[first_event:last_event,,]
          
          model_loc <- tryCatch(relevent::rem(eventlist = env$statisticsREH$edgelist[first_event:last_event,], 
                                    statslist = stats_rem,
                                    supplist = supplist,
                                    timing = "interval",
                                    estimator = "MLE"), error = function(error_message) {NULL})  
                            
          # update model_status
          if(!is.null(model_loc)){
              Sigma <- tryCatch(as.matrix(solve(model_loc$hessian)), 
                                        error = function(error_message) {matrix(-1,nrow=dim(stats_rem)[3],ncol=dim(stats_rem)[3])})
            if((isSymmetric(model_loc$hessian) & (Sigma[1,1] != (-1)) )){
              ###-### if(isSymmetric(as.matrix(solve(model_loc$hessian))))
              ###-###   {
                    model_status <- TRUE
              ###-###   }
                            
            }
          }
          if(model_status == FALSE) {print(q); return(model_loc)}
          
          if(model_status){
                  # old logical condition
                  #if(!is.null(model_loc)){
                  # env$stepwiseModelsREH$model_status[q] <- tryCatch(isSymmetric(as.matrix(solve(model_loc$hessian))), 
                  #
                  #                            error = function(error_message) {FALSE})
            # WAIC routine 
            if(WAIC){
              # ROUTINE FOR THE COMPLETE MODEL 
              start_waic <- Sys.time()  
              # generate from posterior multivariate normal approximation
              Sigma <- as.matrix(solve(model_loc$hessian))
              time_points <- (env$stepwiseModelsREH$L+1):last_event           
              interevent_time <- reh$intereventTime[time_points]  
              stats_waic <- aperm(stats, perm = c(2,3,1))  
              events <- reh$rehBinary[time_points,]                     
              post_pars <- t(mvtnorm::rmvnorm(n = env$stepwiseModelsREH$n_sims_is, mean = model_loc$coef, sigma = Sigma)) # dim = [pars x n_sim_is]

              #lpd <- matrix(NA,nrow=env$stepwiseModelsREH$n_sims_is,ncol=length(time_points))
              #for(m in 1:length(time_points))
              #{
              #event <- numeric(reh$D)
              #event[env$statisticsREH$edgelist[time_points[m],1]] <- 1
              #lpd[,m] <- apply(post_pars,2,function(y) lpd(pars = y,
              #                                                stats = stats[time_points[m],,],
              #                                                event = event,
              #                                                interevent_time = reh$edgelist$time[time_points[m]]-reh$edgelist$time[time_points[m]-1]))
              #}  
              lpd <- bremory::lpdWAIC(pars = post_pars, 
                                      stats = stats_waic[,,time_points], 
                                      events = t(events), 
                                      interevent_time = interevent_time) 
              elpd_waic_obj <- loo::waic(lpd)
              env$stepwiseModelsREH$p_waic[,q] <- elpd_waic_obj$pointwise[,2] #lpd_and_p_waic[,2]
              env$stepwiseModelsREH$elpd_waic[,q] <- elpd_waic_obj$pointwise[,1] #lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
              rm(elpd_waic_obj)
              
              #lpd_and_p_waic <- getWAIC(pars = post_pars, stats = stats_waic[,,time_points], events = t(events), interevent_time = interevent_time)
              #for(m in 1:length(time_points))
              #{
              #  # (lpd_m - p_waic_m)
              #  event <- numeric(reh$D)
              #  event[env$statisticsREH$edgelist[time_points[m],1]] <- 1
              #  p_m <- apply(post_pars,2,function(y) exp(lpd(pars = y,
              #                                                stats = stats[time_points[m],,],
              #                                                event = event,
              #                                                interevent_time = reh$edgelist$time[time_points[m]]-reh$edgelist$time[time_points[m]-1])))
              #  lpd_m <- log(mean(p_m))
              #  p_waic_m <- var(log(p_m))
              #  env$stepwiseModelsREH$p_waic[m,q] <- p_waic_m 
              #  env$stepwiseModelsREH$elpd_waic[m,q] <- lpd_m - p_waic_m
              #}
              end_waic<- Sys.time()  
              cat("\n WAIC routine executed in.... " ,end_waic-start_waic, attr(end_waic-start_waic,'units'))
            }
                      
            ### run psis-lfo ###
            ### PSIS-LFO-CV MSAP (1-step-ahead-predictions) ###
            #if(ELPD){
            #psis_1_sap <- PSIS_LFO_CV_1_SAP(n_dyads = reh$D,
            #                                eventlist = env$statisticsREH$edgelist, 
            #                                statslist = stats,
            #                                t = reh$edgelist$time,
            #                                first_event = first_event,
            #                                last_event = last_event, #reh$M, 
            #                                L = last_event - 10, #reh$M-10, #env$stepwiseModels$L #$L is in 
            #                                n_sim_is = env$stepwiseModelsREH$n_sims_is,
            #                                tau = env$stepwiseModelsREH$tau,
            #                                cores = n_threads)  

            #env$stepwiseModelsREH$reestimations_psis[q] <- psis_1_sap$counts               
            #env$stepwiseModelsREH$elpd_lfo[,q] <- psis_1_sap$elpd_lfo
            #env$stepwiseModelsREH$k_psis[,q] <- ifelse(is.null(psis_1_sap$k),rep(NA,9),psis_1_sap$k) #(reh$M - env$stepwiseModelsREH$L)
            #rm(psis_1_sap)
            #}                
            end <- Sys.time()
            # print message 
            cat("\n ",file_name,": Iteration ",q,"(rejected up to now ",q_rej," ) with ",env$statisticsREH$K_q," intervals completed in --> ", end-start, attr(end-start,'units'),
            ".\n")  ## and ", psis_1_sap$counts, "  no. re-estimations ##
            # print message
            

            #### storing results ####
            n_pars_loc <- dim(stats)[3]
            env$stepwiseModelsREH$mle_betas[q,c(1:n_pars_loc)] <- as.vector(model_loc$coef)
            env$stepwiseModelsREH$vcov_betas[c(1:n_pars_loc),c(1:n_pars_loc),q] <- as.matrix(solve(model_loc$hessian))
            env$stepwiseModelsREH$loglik[q] <- model_loc$loglik
            env$stepwiseModelsREH$BIC[q] <- model_loc$BIC 



            rm(endogenous_stats,n_pars_loc,model_loc,stats) 
                    # saving partial results every 5 models
            if(q%%5==0)
            {
                cat('\n ... Saving partial results ... \n')
                env$stepwiseModelsREH$q_rej <- q_rej
                save(reh,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

                # save intervals table
                partial_result_1 <- data.frame(table(env$stepwiseModelsREH$K[1:q]),iteration = q)
                write.csv(partial_result_1, file = paste(file_name,".csv",sep=""))
            }
          }
          else{
            widths_q_rej_loc <- rep(NA,dim(env$stepwiseModelsREH$widths)[2])
            widths_q_rej_loc[1:(env$statisticsREH$K_q+1)] <- env$statisticsREH$widths_q
            env$stepwiseModelsREH$widths_rejected <-  rbind(env$stepwiseModelsREH$widths_rejected, widths_q_rej_loc) # widths_rejected
            env$stepwiseModelsREH$K_rejected <-  c(env$stepwiseModelsREH$K_rejected, env$statisticsREH$K_q) # K_rejected
            env$stepwiseModelsREH$type_rejected <-  c(env$stepwiseModelsREH$type_rejected, env$stepwiseModelsREH$interval_type[q]) # type_rejected
            q_rej <- q_rej + 1
            print("rejected")
          }

        }

        # class and structure definition
        print(q)
        # update q index
        q <- q +1
    }

    cat('\n ... Saving final results ... \n')
    env$stepwiseModelsREH$q_rej <- q_rej
    save(reh,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

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
#' @param first_event ...
#' @param last_event ...
#' @param L ...
#' @param n_sim_is ...
#' @param tau ...
#' @param cores ...

#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModelsREH'
#' 
#' @export
PSIS_LFO_CV_1_SAP <- function(n_dyads,eventlist,statslist,t,first_event,last_event,L,n_sim_is,tau,cores)
{
  J_i <- L+1 # set of time points contributing in the calculation of the importance ratios. 
  time_points <- ((L+1):last_event) 
  
  log_p_J_i <- matrix(NA, nrow = length(time_points), ncol = n_sim_is)

  # output list objects
  elpd_lfo <- rep(NA,length(time_points))
  count_no_estimation <- 0 # know how many times the model is estimated considering a larger subset than the previous step
  k_vec <- NULL
  
  supplist_L <- list() 
  supplist_L[[1]]<-matrix(TRUE, length(first_event:L), n_dyads)

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
  model_subset <- relevent::rem(eventlist = eventlist[first_event:L,],
                                statslist = statslist[first_event:L,,],
                                supplist = supplist_L,
                                timing = "interval",
                                estimator = "MLE")
                               
  #### relevent::rem ####

  # store MLEs and standard errors
  mle_subset <- model_subset$coef
  vcov_subset <- tryCatch(as.matrix(solve(model_subset$hessian)),
                            error = function(error_message) {matrix(-1,nrow=dim(statslist)[2],ncol=dim(statslist)[2])})                  
  if(vcov_subset[1,1] != (-1)){
    #draw a sample from a multivariate normal distribution (approximation of the posterior distribution)
    pars_draws <- t(mvtnorm::rmvnorm(n = n_sim_is, mean = mle_subset, sigma = vcov_subset)) # dim = [pars x n_sim_is]

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
        supplist_L[[1]] <- matrix(TRUE, length(first_event:(time_points[i]-1)) , n_dyads)

        #re-estimating the model until (time_points[i]-1)
        #### relevent::rem ####                        
        model_subset <- relevent::rem(eventlist = eventlist[first_event:(time_points[i]-1),],
                                      statslist = statslist[first_event:(time_points[i]-1),,],
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
  }
  else{cat("PSIS-LFO did not run because of a singular hessian matrix")}

  return(list(elpd_lfo = elpd_lfo, counts = count_no_estimation, k = k_vec))
}

