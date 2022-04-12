#' smm (Stepwise Memory Models)
#'
#' A function which run the estimation of Q models and updates the environment 'stepwiseModels' 
#'
#' @param formula [not working]formula specifing the linear predictor (this argument is still being coded, plase use the argument 'effects')
#' @param reh either a 'reh' object (see remify::reh() function) or the edgelist
#' @param effects [working]list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param data [not working]list of exogenous variables (time varying or not) (this argument is still being coded, please use the argument 'effects')
#' @param intervals it can be specified in two different ways: (1) as a list of intervals with "widths", "K"; (2) as list specifying: a vector of integers indicating the number of steps per stepwise model ("K"), the maximum time width for intervals ("maxWidth"), the minimum width for intervals ("minDiff"), the number of simulations per each K ("nsimK"). With the second specification, the function bremory::intervals() will be used.
#' @param WAIC default is TRUE and is needed where WAIC weights have to be estimated
#' @param ELPD default is FALSE and is TRUE if ELPD-LFO has to be calculated (PSIS approximation is used [ref. to literature] for one-step-ahead predictions). This procedure might slow down the algorithm.
#' @param firstEvent integer indicating the event at which start the event sequence and statistic in the estimation stage
#' @param lastEvent integer indicating the event at which stop the event sequence and statistic in the estimation stage
#' @param nthreads number of threads to create
#' @param filename a string indicating the name of the file .RData that will save results each 20 models.
#' @param env environment where the user is currently working (default is set to .GlobalEnv, globalenv())
#' @param ... other arguments referring to those in bremory::elpd()
#'
#' @return  object of class 'ssm'
#'
#' @export

smm <- function(formula = NULL,
                reh = NULL,
                data = NULL,
                effects = NULL, 
                intervals = NULL,
                WAIC = TRUE,
                ELPD = FALSE,
                firstEvent = NULL,
                lastEvent = NULL,
                nthreads = parallel::detectCores()-2,
                filename = "simulation_0", # probably this argument will be removed
                env = globalenv(),
                ...)
{
    more.args <- list(...)
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

    #### START PARALLELIZATION SETTINGS
    # set Cores for parallelization purpose
    if(is.null(nthreads) | (nthreads<=0)) nthreads <- 1
    else{
      if(nthreads > (floor(parallel::detectCores()) - 2)) stop(cat("'nthreads' is recommended to be set at most to ", (floor(parallel::detectCores()) - 2),"."))
    }
    #### END PARALLELIZATION SETTINGS

    # create environment 'stepwiseModels' where, widths are stored as well as the output of each model (information useful to the estimation stage)
    env$stepwiseModels <- new.env()

    # ... firstEvent and lastEvent ...
    env$stepwiseModels$firstEvent <- firstEvent
    env$stepwiseModels$lastEvent <- lastEvent
    if(!is.null(firstEvent) & !is.null(lastEvent)){
      if(firstEvent > lastEvent){
        stop("first event integer value must be lower than the last event iteger value.")
      }
    }
    else{
      env$stepwiseModels$firstEvent <- firstEvent <- 1
      env$stepwiseModels$lastEvent <- lastEvent <- reh$M
    }

    # WEIGHTS/INTENSITIES of events should be handled when statistics are to be calculated with them, if not specified, statistics won't be calculated by using event weights
    # if events have weights/intensities (remify should process this stage and return a rehBinary with weights instead of ones)
    #if(!is.null(env$initializeREH$weights))
    #{
    #    weights_loc <- matrix(rep(env$initializeREH$weights,env$initializeREH$n_dyads), nrow = env$initializeREH$M,ncol = env$initializeREH$n_dyads)
    #    env$statisticsREH$binaryREH$dyadicREH <- weights_loc*env$statisticsREH$binaryREH$dyadicREH
    #}


    # create a new environment 'statisticsREH' (the largest in the estimation process and continuously updated according to the characteristic of the stepswise model
    # at a specific iteration)
    bremory::statisticsREH(formula = formula, data = data, reh = reh, effects = effects, env = env)

    # Handling argument 'intervals' ...
    generateIntervalsFlag <- FALSE # creating a disabled flag. This will be activated if in the estimation loop the generation o a sequence of widths will be required
    if(is.null(intervals)) stop("argument 'intervals' must be specified in order to run the estimation")
    else{
      if(!is.list(intervals)) stop("argument 'intervals' must be a list according to the documentation")
      else{
        if(class(intervals)=="intervals"){
          if(all(names(intervals) %in% c("widths","type","K"))){
            env$stepwiseModels$intervals <- intervals 
            intervals$maxWidth <- max(na.omit(as.vector(env$stepwiseModels$intervals$widths)))
          } 
          else{stop("wrong structure for object of class 'intervals'")} 
        }
        else{
          if(all(names(intervals) %in% c("widths","K"))){
            # add several checks on widths and K:
            # 1. widths must be increasing and different from one another (throw ERROR if not)
            # 2. the minDiff should be checked with the reh object such that there is no width smaller than the observed min(intereventTime) (throw an ERROR here). Statistics for that  interval will be always zero 
            # 3. consistency between widths and K 
            # 4. check 1-3 also in the case of an object of class 'intervals'
            env$stepwiseModels$intervals <- intervals 
            intervals$maxWidth <- max(na.omit(as.vector(env$stepwiseModels$intervals$widths)))
          } 
          else{
            if(all(names(intervals) %in% c("K","maxWidth","minDiff","nsimK"))){
              generateIntervalsFlag <- TRUE # activating flag for 
              if((intervals$nsimK %% 2) != 0) intervals$nsimK <- (intervals$nsimK - 1)
              env$stepwiseModels$intervals <- list()
              env$stepwiseModels$intervals$widths <- matrix(NA, nrow = length(intervals$K)*(intervals$nsimK+1), ncol = (max(intervals$K)+1)) 
              env$stepwiseModels$intervals$K <- rep(intervals$K,each = (intervals$nsimK+1))
              env$stepwiseModels$intervals$type <- rep(c("equal",rep("increasing",intervals$nsimK/2),rep("decreasing",intervals$nsimK/2)),length(intervals$K))
              # shuffling order of models: this is done not to increasingly overload the cpu (in this way sometimes a complex model is estimated, other times a much simpler one)
              shuffled_indices <- sample(x=1:dim(env$stepwiseModels$intervals$widths)[1],size=dim(env$stepwiseModels$intervals$widths)[1])
              env$stepwiseModels$intervals$K <- env$stepwiseModels$intervals$K[shuffled_indices]
              env$stepwiseModels$intervals$type <- env$stepwiseModels$intervals$type[shuffled_indices]
            }
            else{stop("structure of list 'intervals' is not defined according to the documentation")}
          } 

        }
      }
    }

    ## Handling Rejected models list ...
    env$stepwiseModels$rejected <- list()
    env$stepwiseModels$rejected$widths <- env$stepwiseModels$rejected$K <- env$stepwiseModels$rejected$type <- NULL
    
    # Storing dimensions and variables name
    env$stepwiseModels$Q <- dim(env$stepwiseModels$intervals$widths)[1]    

    # creating empty vector where storing timing information for estimating each model in the bag
    env$stepwiseModels$time_sim <- list()
    env$stepwiseModels$time_sim$value <- rep(NA,env$stepwiseModels$Q)
    env$stepwiseModels$time_sim$unit <- rep(NA,env$stepwiseModels$Q)
    
    ### --RELEVENT-- ###
    # this line WILL BE REMOVED 
    #reh$edgelist$time <- cumsum(as.numeric(reh$intereventTime)) # with remstimate::remstimate this column won't be part of the estimation stage, instead the intereventTime will be used. Therefore this line WILL BE REMOVED
    ### --RELEVENT-- ###

    # names of endogenous and exogenous statistics (will change with the new data and formula input)
    env$stepwiseModels$stats_names_endo <- env$statisticsREH$endo_effects 
    env$stepwiseModels$stats_names_exo <- dimnames(env$statisticsREH$exogenous_stats)[[3]]

    ### --RELEVENT-- ###
    ## BEGIN rem supplist (will be removed once we integrate 'remstimate')
    #supplist <- list() 
    #supplist[[1]]<-matrix(TRUE, nrow = length(firstEvent:lastEvent), ncol = reh$D)
    ## END rem supplist
    ### --RELEVENT-- ###

    # allocating space for each model result (Q models in total)
    env$stepwiseModels$coef  <- matrix(NA, nrow = env$stepwiseModels$Q, ncol = (( max(env$stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S)) #env$stepwiseModels$coef_remstimate
    env$stepwiseModels$vcov <- array(NA , dim = c((( max(env$stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S),
    (( max(env$stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S), env$stepwiseModels$Q))   #env$stepwiseModels$vcov_remstimate
    env$stepwiseModels$loglik <- env$stepwiseModels$BIC <- rep(NA, env$stepwiseModels$Q) # env$stepwiseModels$loglik_no_endo

   
    if(ELPD | WAIC){
      env$stepwiseModels$elpdSettings <- list()
      # ... L parameter ...
      # L is the number of events to use as history to estimate the predictive pointwise density of future events (set to M-10 by default)
      # M-L are the number of events to predict in WAIC and/or ELPD and their time span correspond to the maximum time width observed in the intervals
      if(!("L" %in% names(more.args))){
        which_L <- tryCatch(max(which(cumsum(reh$intereventTime[firstEvent:lastEvent])<(max(cumsum(reh$intereventTime[firstEvent:lastEvent]))-intervals$maxWidth))),warning = function(x) return(NA), erorr = function(x) return(NA))
        env$stepwiseModels$elpdSettings$L <- ifelse(!is.na(which_L),c(firstEvent:lastEvent)[which_L], lastEvent - 10) 
        if(is.na(which_L)) warning("'L' is set by default to 'lastEvent - 10'.")
      }
      else{
        if(more.args$L>=reh$M) stop("'L' must be smaller than the length of the event sequence.")
        else{
          env$stepwiseModels$elpdSettings$L <- more.args$L
        }
      }

      # ... A_SAP parameter ...
      # A_SAP (A-Steps-Ahead-Predictions) is the number of steps ahead we want consider in the prediction (set to 1 by default) 
      if(!("A_SAP" %in% names(more.args))){
        env$stepwiseModels$elpdSettings$A_SAP <- 1
      }
      else{
        env$stepwiseModels$elpdSettings$A_SAP <- more.args$A_SAP
      }

      # ... n_sims_is ... (change name to nsim because they refer both to WAIC and ELPD-PSIS procedure)
      if(!("n_sims_is" %in% names(more.args))){
        env$stepwiseModels$elpdSettings$n_sims_is <- 1e03
      }
      else{
        env$stepwiseModels$elpdSettings$n_sims_is <- more.args$n_sims_is
      }

      if(ELPD){ ## ELPD LFO PSIS memory allocation
        if(!("tau" %in% names(more.args))){
          env$stepwiseModels$elpdSettings$tau <- 0.7 # default value according to literature (cite which paper here)
        }
        else{
          env$stepwiseModels$elpdSettings$tau <- more.args$tau
        }
        #env$stepwiseModels$ELPD <- list()
        #env$stepwiseModels$ELPD$reestimations_psis <- rep(NA, env$stepwiseModels$Q)
        #env$stepwiseModels$ELPD$elpd_lfo <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L)) , ncol = env$stepwiseModels$Q) 
        #env$stepwiseModels$ELPD$k_psis <- matrix(NA, nrow = nrow = ((length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L))-1), ncol = env$stepwiseModels$Q) 
      }
      else{
        env$stepwiseModels$ELPD <- NULL
      }

      if(WAIC){ # WAIC memory allocation
        #env$stepwiseModels$WAIC <- list()
        #env$stepwiseModels$WAIC$elpdWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L)), ncol= env$stepwiseModels$Q)
        #env$stepwiseModels$WAIC$pWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L)), ncol= env$stepwiseModels$Q)

        # start WAIC_mine
        env$stepwiseModels$WAIC_mine <- list()
        env$stepwiseModels$WAIC_mine$pWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L)), ncol= env$stepwiseModels$Q)
        env$stepwiseModels$WAIC_mine$elpdWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:env$stepwiseModels$elpdSettings$L)), ncol= env$stepwiseModels$Q)
      }
      else{
        env$stepwiseModels$WAIC <- NULL
      }
    }

    # adapted temporarily to the current version of the package "bremory"
    # start #
    reh$rehBinary <- getBinaryREH(dyad=reh$edgelist[,2],D=reh$D)
    reh$risksetMatrix <- as.matrix(expand.grid(0:(reh$N-1),0:(reh$N-1)))
    reh$risksetMatrix <- reh$risksetMatrix[-which(reh$risksetMatrix[,1]==reh$risksetMatrix[,2]),c(2,1)]
    reh$risksetCube <- array(NA,dim=c(reh$N,reh$N,1))
    idx <- 0
    for(n_row in 1:reh$N){
      for(n_col in 1:reh$N){
        if(n_row != n_col){
          reh$risksetCube[n_row,n_col,1] <- idx
          idx <- idx + 1
        }
      }
    }
    edgelist_conv <- matrix(NA,nrow=reh$M,ncol=2)
    for(m in 1:reh$M){
      edgelist_conv[m,] <- reh$risksetMatrix[reh$edgelist[m,2]+1,]
    }
    reh_estimate <- reh
    reh$edgelist <- data.frame(time=reh$edgelist[,1],actor1=edgelist_conv[,1],actor2=edgelist_conv[,2],type=rep(0,reh$M),weight=reh$edgelist[,3])
    rm(n_row,n_col,idx,m,edgelist_conv)
    # end #
    # adapted temporarily to the current version of the package "bremory"
  
    ### --RELEVENT-- ###
    # converting the edgelist to a suitable format for relevent::rem()
    #env$statisticsREH$edgelist <- data.matrix(convertToReleventEdgelist(reh = reh))
    ### --RELEVENT-- ###

    q <- 1
    q_rej <- 0
    # Progress bar settings (based on the number of generated events)
    pb <- txtProgressBar(min = 1, max = env$stepwiseModels$Q, style = 3)
    cat("\n Starting estimating the bag of models ! \n")
    while(q <= env$stepwiseModels$Q)
    {   
        start_time_q <- Sys.time()
        modelStatus <- FALSE
        while(!modelStatus){
          if(generateIntervalsFlag){
            env$statisticsREH$widths_q <-  bremory::intervals(K= env$stepwiseModels$intervals$K[q],
                                                        minDiff = intervals$minDiff,
                                                        maxWidth= intervals$maxWidth,
                                                        intervals= env$stepwiseModels$intervals$type[q]) #na.omit(env$stepwiseModels$intervals$widths[q,])
            env$stepwiseModels$intervals$widths[q,1:(env$stepwiseModels$intervals$K[q]+1)] <- env$statisticsREH$widths_q
            env$statisticsREH$K_q <- env$stepwiseModels$intervals$K[q] 
          }
          else{
            env$statisticsREH$widths_q <-  env$stepwiseModels$intervals$widths[q,1:(env$stepwiseModels$intervals$K[q]+1)]
            env$statisticsREH$K_q <- env$stepwiseModels$intervals$K[q] 
          }

          #getCountsOMP (OpenMP parallelization in C++)
          if(length(env$statisticsREH$endo_effects)>0){

            # getIntervals() will find per each time point the lower and the upper bound for each interval and will automatically assign it to the object 'intervals_backwards'
            start_time <-  Sys.time()
            getIntervals(env = globalenv(), 
                          widths = env$statisticsREH$widths_q, 
                          time = cumsum(reh$intereventTime), 
                          M = reh$M, 
                          K_q = env$statisticsREH$K_q,
                          nthreads = nthreads) # reh$edgelist$time

            ### TO CHECK FIRST AND THEN REMOVE ###
            #env$statisticsREH$intervals[(env$statisticsREH$K_q+1),1:2] <- c(0,0) #only for the first interval at the second time point
            ### TO CHECK FIRST AND THEN REMOVE ###

            intervals_loc <- unique(env$statisticsREH$intervals[,1:2])
            
            env$statisticsREH$counts <- getCountsOMP(binaryREH = reh$rehBinary, # make sure that this matrix is 1/0 only 
                                                    lbs_ubs = intervals_loc, 
                                                    nthreads = nthreads)

            indices <- getCountsIndex(intervals = env$statisticsREH$intervals, 
                                      counts = cbind(env$statisticsREH$counts[,1:2],
                                      c(0:(dim(env$statisticsREH$counts)[1]-1))),
                                      nthreads = nthreads)
            env$statisticsREH$intervals[,3] <-  indices
            rm(intervals_loc,indices)
         
            ### calculating effects ###
            endogenous_stats <- getEndoEffects(env = env,
                                              M = reh$M, 
                                              D = reh$D,
                                              time = reh$edgelist[,1],
                                              edgelist = data.matrix(reh$edgelist), 
                                              risksetMatrix = reh$risksetMatrix, 
                                              risksetCube0 = reh$risksetCube[,,1], # still one event type is used (need to make a change)
                                              nthreads = nthreads) 

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
          stats <- abind::abind(endogenous_stats,env$statisticsREH$exogenous_stats, along = 3) #env$Pshift_loc
          
          
          # estimating model with relevent package #
          stats_rem <- stats[firstEvent:lastEvent,,]

          # adapted temporarily for the current version of the package "bremory"
          # start #
          if(q==1) {
            env$stats_nomemory <- stats_rem
            stats_rem[,,-dim(stats_rem)[3]] <- stats_rem[,,-dim(stats_rem)[3]]/reh$M}
          # end #
          # adapted temporarily for the current version of the package "bremory"

          ### --RELEVENT-- ###
          ##### relevent::rem function #####  
          #model_q <- tryCatch(relevent::rem(eventlist = env$statisticsREH$edgelist[firstEvent:lastEvent,], 
          #                          statslist = stats_rem,
          #                          supplist = supplist,
          #                          timing = "interval",
          #                          estimator = "MLE"), error = function(error_message) {NULL}) 
          ### --RELEVENT-- ###
          model_q <- tryCatch(remstimate::remstimate(reh = reh_estimate,
                                          stats = stats_rem,
                                          method =  "MLE",
                                          model = "tie",
                                          ncores = nthreads), error = function(error_message) {NULL}) 
                                          
          # update modelStatus 
          if(!is.null(model_q)){
              Sigma <- tryCatch(as.matrix(solve(model_q$hessian)), 
                                        error = function(error_message) {matrix(-1,nrow=dim(stats_rem)[3],ncol=dim(stats_rem)[3])})
            if(isSymmetric(model_q$hessian) & (Sigma[1,1] != (-1)) ){
                    modelStatus <- TRUE               
            }
          }                                                          
          
          if(modelStatus){
            # WAIC routine 
            if(WAIC){
              # generate from posterior multivariate normal approximation
              Sigma <- as.matrix(solve(model_q$hessian))
              time_points <- (env$stepwiseModels$elpdSettings$L+1):lastEvent           
              interevent_time <- reh$intereventTime[time_points]  
              stats_waic <- aperm(stats, perm = c(2,3,1))  
              events <- reh$rehBinary[time_points,]                     
              post_pars <- t(mvtnorm::rmvnorm(n = env$stepwiseModels$elpdSettings$n_sims_is, mean = model_q$coef, sigma = Sigma)) # dim = [pars x n_sim_is]

              #lpd <- matrix(NA,nrow=env$stepwiseModels$elpdSettings$n_sims_is,ncol=length(time_points))
              #for(m in 1:length(time_points))
              #{
              #event <- numeric(reh$D)
              #event[env$statisticsREH$edgelist[time_points[m],1]] <- 1
              #lpd[,m] <- apply(post_pars,2,function(y) lpd(pars = y,
              #                                                stats = stats[time_points[m],,],
              #                                                event = event,
              #                                                interevent_time = reh$edgelist$time[time_points[m]]-reh$edgelist$time[time_points[m]-1]))
              #}

              ### temporarily commented  
              #lpd <- bremory::lpdWAIC(pars = post_pars, 
              #                        stats = stats_waic[,,time_points], 
              #                        events = t(events), 
              #                        interevent_time = interevent_time) 
              #elpd_waic_obj <- loo::waic(lpd)
              #env$stepwiseModels$WAIC$pWAIC[,q] <- elpd_waic_obj$pointwise[,2] #lpd_and_p_waic[,2]
              #env$stepwiseModels$WAIC$elpdWAIC[,q] <- elpd_waic_obj$pointwise[,1] #lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
              #env$stepwiseModels$WAIC$WAIC <- elpd_waic_obj$pointwise[,3]
              #rm(elpd_waic_obj)
              ### temporarily commented

              ## START my routine for WAIC
              lpd_and_p_waic <- bremory::getWAIC(pars = post_pars, stats = stats_waic[,,time_points], events = t(events), interevent_time = interevent_time)
              env$stepwiseModels$WAIC_mine$pWAIC[,q] <- lpd_and_p_waic[,2]
              env$stepwiseModels$WAIC_mine$elpdWAIC[,q] <- lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
              ## END my routine for WAIC
            }
                      
            ### run psis-lfo ###
            ### PSIS-LFO-CV MSAP (1-step-ahead-predictions) ###
            #if(ELPD){
            #psis_1_sap <- PSIS_LFO_CV_1_SAP(n_dyads = reh$D,
            #                                eventlist = env$statisticsREH$edgelist, 
            #                                statslist = stats,
            #                                t = reh$edgelist$time,
            #                                firstEvent = firstEvent,
            #                                lastEvent = lastEvent, #reh$M, 
            #                                L = lastEvent - 10, #reh$M-10, #env$stepwiseModels$elpdSettings$L #$L is in 
            #                                n_sim_is = env$stepwiseModels$elpdSettings$n_sims_is,
            #                                tau = env$stepwiseModels$elpdSettings$tau,
            #                                cores = nthreads
            #)  

            #env$stepwiseModels$ELPD$reestimations_psis[q] <- psis_1_sap$counts               
            #env$stepwiseModels$ELPD$elpd_lfo[,q] <- psis_1_sap$elpd_lfo
            #env$stepwiseModels$ELPD$k_psis[,q] <- ifelse(is.null(psis_1_sap$k),rep(NA,9),psis_1_sap$k) #(reh$M - env$stepwiseModels$elpdSettings$L)
            #rm(psis_1_sap)
            #}
            

            #### storing results ####
            n_pars_loc <- dim(stats)[3]
            env$stepwiseModels$coef[q,c(1:n_pars_loc)] <- as.vector(model_q$coef)
            env$stepwiseModels$vcov[c(1:n_pars_loc),c(1:n_pars_loc),q] <- as.matrix(solve(model_q$hessian))
            env$stepwiseModels$loglik[q] <- model_q$loglik
            env$stepwiseModels$BIC[q] <- model_q$BIC 


            ## remstimate saving coef and vcov ##
            #env$stepwiseModels$coef_remstimate[q,c(1:n_pars_loc)] <- as.vector(model_q_remstimate$coef)
            #env$stepwiseModels$vcov_remstimate[c(1:n_pars_loc),c(1:n_pars_loc),q] <- as.matrix(solve(model_q_remstimate$hessian))
            ## remstimate saving coef and vcov ##



            rm(endogenous_stats,n_pars_loc,model_q,stats) 
            # saving partial results every 5 models
            if(q%%250==0)
            {
                cat('\n ... Saving partial results ... \n')
                env$stepwiseModels$q_rej <- q_rej
                save(stepwiseModels,envir = env, file = paste(filename,".RData",sep="")) # also save 'reh'

                # save intervals table
                # partial_result_1 <- data.frame(table(env$stepwiseModels$intervals$K[1:q]),iteration = q)
                # write.csv(partial_result_1, file = paste(filename,".csv",sep=""))
            }
          }
          else{
            if(!generateIntervalsFlag) modelStatus <- TRUE # in this way if intervals are not supposed to be generated we force to go on to the next q
            q_rej <- q_rej +1
            widths_q_rej_loc <- rep(NA,dim(env$stepwiseModels$intervals$widths)[2])
            widths_q_rej_loc[1:(env$statisticsREH$K_q+1)] <- env$statisticsREH$widths_q
            env$stepwiseModels$rejected$widths <-  rbind(env$stepwiseModels$rejected$widths, widths_q_rej_loc) # width rejected
            env$stepwiseModels$rejected$K <-  c(env$stepwiseModels$rejected$K, env$statisticsREH$K_q) # K rejected
            env$stepwiseModels$rejected$type <-  c(env$stepwiseModels$rejected$type, env$stepwiseModels$intervals$type[q]) # type rejected
          }

        }
        end_time_q <- Sys.time()
        env$stepwiseModels$time_sim$value[q] <- as.numeric(end_time_q - start_time_q)
        env$stepwiseModels$time_sim$unit[q] <- attr(end_time_q - start_time_q,"units")
        # update q index
        q <- q+1
        # Progress bar         
        setTxtProgressBar(pb, (q-1))
    }
    env$stepwiseModels$q_rej <- q_rej

    cat('\n ... Saving final results ... \n')

    
    # Check for rejected model when 'intervals' argument consist of an external set of intervals
    if(!generateIntervalsFlag){ 
      modelNA <- na.omit(env$stepwiseModels$BIC)
      modelNA <- as.vector(attr(modelNA,"na.action"))
      if(!is.null(modelNA)){
        env$stepwiseModels$Q <- (env$stepwiseModels$Q - length(modelNA))
        env$stepwiseModels$intervals$widths <- env$stepwiseModels$intervals$widths[-modelNA,]
        env$stepwiseModels$intervals$K <- env$stepwiseModels$intervals$K[-modelNA]
        if(!is.null(env$stepwiseModels$intervals$type)){
          env$stepwiseModels$intervals$type <- env$stepwiseModels$intervals$type[-modelNA]
        }
        env$stepwiseModels$coef <- env$stepwiseModels$coef[-modelNA,]
        env$stepwiseModels$vcov <- env$stepwiseModels$vcov[,,-modelNA]
        env$stepwiseModels$loglik <- env$stepwiseModels$loglik[-modelNA]
        env$stepwiseModels$BIC <- env$stepwiseModels$BIC[-modelNA]


        ## remstimate finalizing results ##
        #env$stepwiseModels$coef_remstimate <- env$stepwiseModels$coef_remstimate[-modelNA,]
        #env$stepwiseModels$vcov_remstimate <- env$stepwiseModels$vcov_remstimate[,,-modelNA]
        ## remstimate finalizing results ##

        # WAIC
        if(WAIC){
          #env$stepwiseModels$WAIC$pWAIC <- env$stepwiseModels$WAIC$pWAIC[,-modelNA]
          #env$stepwiseModels$WAIC$elpdWAIC <- env$stepwiseModels$WAIC$elpdWAIC[,-modelNA]
        }
        # ELPD
        if(ELPD){
          #env$stepwiseModels$ELPD$reestimations_psis <- env$stepwiseModels$ELPD$reestimations_psis[-modelNA]             
          #env$stepwiseModels$ELPD$elpd_lfo <- env$stepwiseModels$ELPD$elpd_lfo[,-modelNA]
          #env$stepwiseModels$ELPD$k_psis <- env$stepwiseModels$ELPD$k_psis[,-modelNA]
        }
        warning("x models were rejected. Therefore, the number of stepwise models per number of intervals might differ")
      }
    }

    # class and structure definition 
    str_out <- structure(list(
                            Q = env$stepwiseModels$Q,
                            intervals = env$stepwiseModels$intervals,
                            coef = env$stepwiseModels$coef,
                            vcov = env$stepwiseModels$vcov,
                            #coef_remstimate = env$stepwiseModels$coef_remstimate, ## remstimate
                            #vcov_remstimate = env$stepwiseModels$vcov_remstimate, ## remstimate
                            loglik = env$stepwiseModels$loglik,
                            BIC = env$stepwiseModels$BIC,
                            #WAIC = env$stepwiseModels$WAIC,
                            WAIC_mine = env$stepwiseModels$WAIC_mine,
                            ELPD = env$stepwiseModels$ELPD,
                            stats_names_endo = env$stepwiseModels$stats_names_endo, ## to change
                            stats_names_exo = env$stepwiseModels$stats_names_exo, ## to change
                            time_sim = env$stepwiseModels$time_sim
                            ), class="smm")

    attr(str_out, "rejected") <- env$stepwiseModels$rejected
    attr(str_out, "window") <- c(env$stepwiseModels$firstEvent, env$stepwiseModels$lastEvent)
    attr(str_out, "elpdSettings") <- env$stepwiseModels$elpdSettings
    attr(str_out, "formula") <- formula # this will substitute 'stats_names_endo' and 'stats_names_exo'

    rm(statisticsREH,stepwiseModels, envir = env)
    return(str_out)
}

#######################################################################################
#######################################################################################
##########(START)             Methods for `smm` object             (START)#############
#######################################################################################
#######################################################################################


#' @title merge.smm
#' @description A function that returns a summary of the temporal network.
#' @param smm a list of \code{smm} objects
#'
#' @export
merge.smm <- function(smm){
    if(!is.null(smm)){
        if(class(smm) == "list"){
            list_class <- unlist(lapply(smm,class))
            if(all(list_class) == "smm"){
              S <- length(smm)
              # create empty data structures
              for(s in 1:S){
              # check.smm() this is a function that checks wether an object is of class smm. In case it is not, the iteration will be skipped and finally a warning will be thrown
              # merge here
              }
              
              print("sono qui !")
            }
            else{
                stop("argument smm must be a list of objects of class 'smm'")
            }
        }
        else{
            if(class(smm) != "smm"){
                stop("argument smm must be a list of objects of class 'smm'")
            }
            else{
              return(smm)
            }
        }
    }
    else{
        stop("missing argument 'smm'")
    }
}

#######################################################################################
#######################################################################################
##########(END)             Methods for `smm` object             (END)#################
#######################################################################################
#######################################################################################



#' PSIS_LFO_CV_1_SAP
#'
#' A function which run th PSIS with LFO CV by using 1-Step-Ahead-Predictions 
#'
#' @param n_dyads ...
#' @param eventlist ...
#' @param statslist ...
#' @param t ...
#' @param firstEvent ...
#' @param lastEvent ...
#' @param L ...
#' @param n_sim_is ...
#' @param tau ...
#' @param cores ...

#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModels'
#' 
#' @export
PSIS_LFO_CV_1_SAP <- function(n_dyads,eventlist,statslist,t,firstEvent,lastEvent,L,n_sim_is,tau,cores)
{
  J_i <- L+1 # set of time points contributing in the calculation of the importance ratios. 
  time_points <- ((L+1):lastEvent) 
  
  log_p_J_i <- matrix(NA, nrow = length(time_points), ncol = n_sim_is)

  # output list objects
  elpd_lfo <- rep(NA,length(time_points))
  count_no_estimation <- 0 # know how many times the model is estimated considering a larger subset than the previous step
  k_vec <- NULL
  
  supplist_L <- list() 
  supplist_L[[1]]<-matrix(TRUE, length(firstEvent:L), n_dyads)

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
  model_subset <- relevent::rem(eventlist = eventlist[firstEvent:L,],
                                statslist = statslist[firstEvent:L,,],
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
        supplist_L[[1]] <- matrix(TRUE, length(firstEvent:(time_points[i]-1)) , n_dyads)

        #re-estimating the model until (time_points[i]-1)
        #### relevent::rem ####                        
        model_subset <- relevent::rem(eventlist = eventlist[firstEvent:(time_points[i]-1),],
                                      statslist = statslist[firstEvent:(time_points[i]-1),,],
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

#' formula_lab
#'
#' A function where do build up a formula object for smm and pmm
#'
#' @param input input of the function is a
#'
#' @return  the function returns a NULL
#' 
#' @export
formula_lab <- function(input){

  # parseEffectChoice()

  # Inertia
  inertia <- function(param = NULL, scaling = c("raw", "std", "prop","log"), start = 0, end = 0)
  {
    scaling <- match.arg(scaling)
    prepEndoVar("inertia", param, scaling, start, end) # need to add a way to swap between smm and pmm 
  }
 # prepEndoVar() // prepExoVar() // prepInteractVar() // should address both Actor Oriented and 
  # Reciprocity

  # InDegreeSender 
  
  # InDegreeReceiver
  
  # OutDegreeSender
  
  # OutDegreeReceiver
  
  # TransitivityClosure
  
  # CyclicClosure          
}


#' smmTemp (Stepwise Memory Models [Temporary Function])
#'
#' A function which run the estimation of Q models and updates the environment 'stepwiseModels' 
#'
#' @param formula formula specifing the linear predictor (this argument is still being coded, plase use the argument 'effects')
#' @param reh either a 'reh' object (see remify::reh() function) or the edgelist
#' @param data list of exogenous variables (time varying or not) (this argument is still being coded, please use the argument 'effects')
#' @param subset old {firstEvent, lastEvent}
#' @param intervals it can be specified in two different ways: (1) as a list of intervals with "widths", "K" and "type"; (2) as list specifying: a vector of integers indicating the number of steps per stepwise model ("K"), the maximum time width for intervals ("maxWidth"), the minimum width for intervals ("minDiff"), the number of simulations per each K ("nsimK"). With the second specification, the function bremory::intervals() will be used.
#' @param WAIC default is TRUE and is needed where WAIC weights have to be estimated
#' @param ELPD default is FALSE and is TRUE if ELPD-LFO has to be calculated (PSIS approximation is used [ref. to literature] for one-step-ahead predictions). This procedure might slow down the algorithm.
#' @param threads number of threads to create
#' @param env environment where the user is currently working (default is set to .GlobalEnv, globalenv())
#' @param ... other arguments referring to those in bremory::elpd()
#'
#' @return  object of class 'ssm'
#'
#' @export
smmTemp <- function(formula,
                reh,
                data,
                subset, #firstEvent LastEvent
                intervals = NULL,
                WAIC = TRUE,
                ELPD = FALSE,
                threads = parallel::detectCores()-2,
                env = globalenv(),
                ...)

  # formula1 = 1 + inertia + reciprocity + ABBA + Exo1 + Exo2 + Exo1:inertia + Exo1:Exo2
  # formula2  = 1 + inertia + reciprocity + ABBA + Exo1*inertia + Exo1*Exo2
  # where Exo1*inertia = inertia + Exo1 + Exo1:inertia
  # and Exo1*Exo2 = Exo1 + Exo2 + Exo1:Exo2
  # Therefore : the processing of 'formula2' should bring to the same processing with 'formula1'
{

}