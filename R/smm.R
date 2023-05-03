#' Step-wise Memory Model
#'
#' A function for the estimation of one or more step-wise models (at the moment working without remstats package)
#'
#' @param formula formula specifing the linear predictor (list of formulas if actor-oriented model)
#' @param reh either a 'reh' object (see remify::remify() function) or the edgelist to process
#' @param intervals interval(s) defining the step-wise function(s). It can be specified in two different ways: (1) as a list of intervals with "widths", "K"; (2) as list specifying: a vector of integers indicating the number of steps per stepwise model ("K"), the maximum time width for intervals ("maxTime"), the minimum width for intervals ("minWidth"), the number of simulations per each K ("nsimK"). With the second specification, the function bremory::intervals() will be used.
#' @param data [\emph{optional}] list of exogenous variables (time varying or not)
#' @param WAIC [\emph{optional}] default is TRUE and is needed where WAIC weights have to be estimated
#' @param ELPD [\emph{optional}] default is FALSE and is TRUE if ELPD-LFO has to be calculated (PSIS approximation is used [ref. to literature here] for one-step-ahead predictions). This procedure might slow down the algorithm. [please do not run. Not yet implemented]
#' @param ncores [\emph{optional}] number of ncores to use for the parallelization
#' @param subset [\emph{optional}] first and last event index defining the sub-history on which to estimate the stepwise memory model(s)
#' @param ... [\emph{optional}] other arguments referring to those in bremory::elpd()
#'
#' @return  object of class 'smm' with attributes and methods explained in the vignette(topic = "smm", package = "bremory")
#'
#' @export
#'
#' @examples 
#' 
#' # examples here
smm <- function(formula,
                reh,
                intervals,
                data = NULL,
                WAIC = TRUE,
                ELPD = FALSE,
                ncores = 1L,
                subset = NULL, #firstEvent LastEvent [[ TO REMOVE ]]
                ...)
{
    more.args <- list(...)

    # set 'env' (environment) to .GlobalEnv 
    env <- globalenv()

    # ... remify object ('reh' input argument)
    if(!inherits(reh,"remify")){
        stop("'reh' must be a 'remify' object (see ?remify::remify).")
    }

    # ... check for model == "tie" [[TEMPORARY]]
    if(attr(reh,"model") != "tie"){
      stop("actor-oriented modeling framework not (yet) available.")
    }

    # ... 
    if(is.null(reh$C)){ # local change
        reh$C <- 1 # we set it to one 
    }
    reh$D <- reh$D/reh$C

    # ... ncores
    if(is.null(ncores) | (ncores<=0)) ncores <- 1L
    else{
        if((parallel::detectCores() == 2L) & (ncores > 1L))
            stop("'ncores' is recommended to be set at most to 1.")
        else if((parallel::detectCores() > 2L) & (ncores > floor(parallel::detectCores()-2L)))
            stop("'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)")
    }

    # ... Creating list 'stepwiseModels' (list where all the useful information of all the step-wise models will be stored)
    stepwiseModels <- list() 

    # ... Checking 'subset' argument and defining 'firstEvent' and 'lastEvent' for estimation on sub-history [[TO REMOVE]]
    firstEvent <- lastEvent <- NULL
    if(!is.null(subset)){
      if(is.numeric(subset)){
        if(length(subset==2)){
          if(is.unsorted(subset)){
            subset <- sort(subset)
            warning("Sorting the input 'subset'. Check whether the input provided is correct.")
          }
          firstEvent <- subset[1]
          lastEvent <- subset[2]
        }
        else {
          stop("The input 'subset' must be a numeric vector of length 2: first and last event index of the sub-history.")
        }
      }
      else {
        stop("The input 'subset' must be a numeric vector of length 2: first and last event index of the sub-history.")
      }
    }
    else{
      firstEvent <- 1
      lastEvent <- reh$M
    }
    stepwiseModels$firstEvent <- firstEvent
    stepwiseModels$lastEvent <- lastEvent

    # ... weights of events
    stepwiseModels$weighted <- FALSE
    if(attr(reh,"weighted")){
      stepwiseModels$weighted <- TRUE
    }
    # WEIGHTS/INTENSITIES of events should be handled when statistics are to be calculated with them, if not specified, statistics won't be calculated by using event weights
    # if events have weights/intensities (remify should process this stage and return a rehBinary (or some other structure) with weights instead of ones)

    
    # ... types of events
    stepwiseModels$with_type <- FALSE
    if(attr(reh,"with_type")){
      stepwiseModels$with_type <- TRUE
    }

    # ... Defining vector of endogenous statistics available for the step-wise modeling approach (the names are similar - if not the same - to the ones in the 'remstats' package)
    first_order_statistics <- c("inertia",
                                "reciprocity",
                                "indegreeSender",
                                "indegreeReceiver",
                                "outdegreeSender",
                                "outdegreeReceiver",
                                "totaldegreeSender",
                                "totaldegreeReceiver")
    endo_vec <- c(first_order_statistics,
                  "transitivityClosure",
                  "cyclicClosure", 
                  "sendingClosure",
                  "receivingClosure",
                  "pshift")

    # ... Parsing 'formula'
    smmf <- NULL 
    if(!inherits(formula,"formula")) stop("argument 'formula' must be an object of class 'formula'")
    # ... Saving formula
    stepwiseModels$formula <- formula
    # ... Parsing formula
    smmf <- stats::terms(x = formula, keep.order = TRUE)
    if(sum(attr(smmf,"factors")[1,])==0) stop("The argument 'formula' doesn't need a response variable definition")
    vars <- rownames(attr(smmf, "factors")) #  variables names, alternatively also: as.character(attr(smmf,"variables"))[-1]
    attr(smmf,"endos") <- vars[which(!is.na(match(vars,endo_vec)))]
    attr(smmf,"exos") <- vars[-which(!is.na(match(vars,endo_vec)))]
    if(any(apply(attr(smmf,"factors"),2,sum)>2)) stop('three-terms interactions not supported.')
    # ... Understanding transformation of variables [not coded yet]
    #[ ... only if protected by the function I() user can define transformation of variables, such as I(inertia**0.5), I(Exo2^2), I(inertia - mean(inertia)) and so forth]

    # ... Creating the environment 'statisticsREH' (this environment will temporary store the statistics and processed intervals of each step-wise model and it will occupy the largest memory space)
    env$statisticsREH <- new.env()
    env$statisticsREH$stats_names <- list()
    env$statisticsREH$intervals <- NULL # creating a NULL object inside the environment 'statisticsREH' where the intervals of each step-wise function will be processed

    # ... Processing 'data'
    process_data(data = data, parsed_formula = smmf, reh = reh) # in future updates the processing will be operated by the 'remstats' package
    data  <- NULL # free-ing space
  

    # ... Saving vector of endogenous statistics names and length from 'formula'
    env$statisticsREH$endogenous_statistics <- attr(smmf, "endos")

    # ... Saving vector of exogenous statistics names and length from 'formula'
    env$statisticsREH$exogenous_statistics <- attr(smmf, "exos")

    # ... argument 'intervals' 
    generateIntervalsFlag <- FALSE # creating a disabled flag. This will be activated if in the estimation loop the generation o a sequence of widths will be required
    if(!is.list(intervals)) stop("argument 'intervals' must be a list according to the documentation")
    else{
      if(inherits(intervals,"intervals")){
        if(all(names(intervals) %in% c("widths","type","K"))){
          # add checks about appropriate values of maxTime and minWidth given the data at hand (reh)
          intervals$maxTime <- max(na.omit(as.vector(intervals$widths)))
          stepwiseModels$intervals <- intervals 
        } 
        else{stop("wrong structure for object of class 'intervals'")} 
      }
      else{
        if(all(names(intervals) %in% c("widths","K"))){
          # add several checks on widths and K:
          # 1. widths must be increasing and different from one another (throw ERROR if not)
          # 2. the minWidth should be checked with the reh object such that there is no width smaller than the observed min(intereventTime) (throw an ERROR here). Statistics for that  interval will be always zero 
          # 3. consistency between widths and K 
          # 4. check 1-3 also in the case of an object of class 'intervals'
          intervals$maxTime <- max(na.omit(as.vector(intervals$widths)))
          stepwiseModels$intervals <- intervals 
        } 
        else{
          if(all(names(intervals) %in% c("K","maxTime","minWidth","nsimK"))){
            # add checks about appropriate values of maxTime and minWidth given the data at hand (reh)
            generateIntervalsFlag <- TRUE # activating flag for 
            if((intervals$nsimK %% 2) != 0) intervals$nsimK <- (intervals$nsimK - 1)
            stepwiseModels$intervals$widths <- matrix(NA, nrow = length(intervals$K)*(intervals$nsimK+1), ncol = (max(intervals$K)+1)) 
            stepwiseModels$intervals$K <- rep(intervals$K,each = (intervals$nsimK+1))
            stepwiseModels$intervals$type <- rep(c("equal",rep("increasing",intervals$nsimK/2),rep("decreasing",intervals$nsimK/2)),length(intervals$K))
            # shuffling order of models: this is done not to increasingly overload the cpu (in this way sometimes a complex model is estimated, other times a much simpler one)
            shuffled_indices <- sample(x=1:dim(stepwiseModels$intervals$widths)[1],size=dim(stepwiseModels$intervals$widths)[1])
            stepwiseModels$intervals$K <- stepwiseModels$intervals$K[shuffled_indices]
            stepwiseModels$intervals$type <- stepwiseModels$intervals$type[shuffled_indices]
          }
          else{stop("structure of list 'intervals' is not defined according to the documentation")}
        } 

      }
    }
    
    ## ... Handling Rejected models list 
    stepwiseModels$rejected <- list()
    stepwiseModels$rejected$widths <- stepwiseModels$rejected$K <- stepwiseModels$rejected$type <- NULL
    
    # ... Storing dimensions and variables name
    stepwiseModels$Q <- dim(stepwiseModels$intervals$widths)[1]    

    # ... Creating empty vector where storing timing information for estimating each model in the bag
    stepwiseModels$time_sim <- list()
    stepwiseModels$time_sim$value <- rep(NA,stepwiseModels$Q)
    stepwiseModels$time_sim$unit <- rep(NA,stepwiseModels$Q)

    stepwiseModels$coefficients  <- list() # matrix(NA, nrow = stepwiseModels$Q, ncol = (( max(stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S))  # to correct number of columns ##[[OLD COMMENT]] stepwiseModels$coefficients_remstimate
    stepwiseModels$vcov <- list() # array(NA , dim = c((( max(stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S), (( max(stepwiseModels$intervals$K) * env$statisticsREH$P) + env$statisticsREH$S), stepwiseModels$Q))   ##[[OLD COMMENT]] stepwiseModels$vcov_remstimate
    stepwiseModels$loglik <- stepwiseModels$BIC <- rep(NA, stepwiseModels$Q) ##[[OLD COMMENT]] stepwiseModels$loglik_no_endo

    # ... Handling WAIC and ELPD settings
    if(ELPD | WAIC){
      stepwiseModels$elpdSettings <- list()
      # ... L parameter ...
      # L is the number of events to use as history to estimate the predictive pointwise density of future events (set to M-10 by default)
      # M-L are the number of events to predict in WAIC and/or ELPD and their time span correspond to the maximum time width observed in the intervals
      if(!("L" %in% names(more.args))){
        which_L <- tryCatch(max(which(cumsum(reh$intereventTime[firstEvent:lastEvent])<(max(cumsum(reh$intereventTime[firstEvent:lastEvent]))-intervals$maxTime))),warning = function(x) return(NA), erorr = function(x) return(NA))
        stepwiseModels$elpdSettings$L <- ifelse(!is.na(which_L),c(firstEvent:lastEvent)[which_L], lastEvent - 10) 
        if(is.na(which_L)) warning("'L' is set by default to 'lastEvent - 10'.")
      }
      else{
        if(more.args$L>=reh$M) stop("'L' must be smaller than the length of the event sequence.")
        else{
          stepwiseModels$elpdSettings$L <- more.args$L
        }
      }

      # ... A_SAP parameter ...
      # A_SAP (A-Steps-Ahead-Predictions) is the number of steps ahead we want consider in the prediction (set to 1 by default) 
      if(!("A_SAP" %in% names(more.args))){
        stepwiseModels$elpdSettings$A_SAP <- 1
      }
      else{
        stepwiseModels$elpdSettings$A_SAP <- more.args$A_SAP
      }

      # ... n_sims_is ... (change name to nsim because they refer both to WAIC and ELPD-PSIS procedure)
      if(!("n_sims_is" %in% names(more.args))){
        stepwiseModels$elpdSettings$n_sims_is <- 1e03
      }
      else{
        stepwiseModels$elpdSettings$n_sims_is <- more.args$n_sims_is
      }

      if(ELPD){ ## ELPD LFO PSIS memory allocation
        if(!("tau" %in% names(more.args))){
          stepwiseModels$elpdSettings$tau <- 0.7 # default value according to literature (cite which paper here)
        }
        else{
          stepwiseModels$elpdSettings$tau <- more.args$tau
        }
        #stepwiseModels$ELPD <- list()
        #stepwiseModels$ELPD$reestimations_psis <- rep(NA, stepwiseModels$Q)
        #stepwiseModels$ELPD$elpd_lfo <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L)) , ncol = stepwiseModels$Q) 
        #stepwiseModels$ELPD$k_psis <- matrix(NA, nrow = nrow = ((length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L))-1), ncol = stepwiseModels$Q) 
      }
      else{
        stepwiseModels$ELPD <- NULL
      }

      if(WAIC){ # WAIC memory allocation
        #stepwiseModels$WAIC <- list()
        #stepwiseModels$WAIC$elpdWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L)), ncol= stepwiseModels$Q)
        #stepwiseModels$WAIC$pWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L)), ncol= stepwiseModels$Q)

        # start WAIC_mine
        stepwiseModels$WAIC_mine <- list()
        stepwiseModels$WAIC_mine$pWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L)), ncol= stepwiseModels$Q)
        stepwiseModels$WAIC_mine$elpdWAIC <- matrix(NA, nrow = (length(firstEvent:lastEvent)-length(firstEvent:stepwiseModels$elpdSettings$L)), ncol= stepwiseModels$Q)
      }
      else{
        stepwiseModels$WAIC <- NULL
      }
    }

        # processing attr(reh,"dyad")
    if(reh$C>1){
      for(m in 1:reh$M){
        attr(reh,"dyad")[m] <- remify:::getDyadIndex(actor1 = reh$edgelist$actor1_ID[m]-1,actor2 = reh$edgelist$actor2_ID[m]-1,type = 0, N = reh$N, directed = attr(reh,"directed"))+1
      }
    }

    # adapted temporarily to the current version of the package "bremory"
    # start #
    reh$rehBinary <- getBinaryREH(dyad=attr(reh,"dyad")-1,D=reh$D)

    #reh$risksetMatrix <- as.matrix(expand.grid(0:(reh$N-1),0:(reh$N-1)))
    #reh$risksetMatrix <- reh$risksetMatrix[-which(reh$risksetMatrix[,1]==reh$risksetMatrix[,2]),c(2,1)]
    #reh$risksetCube <- array(NA,dim=c(reh$N,reh$N,1))
    #idx <- 0
    #for(n_row in 1:reh$N){
    #  for(n_col in 1:reh$N){
    #    if(n_row != n_col){
    #      reh$risksetCube[n_row,n_col,1] <- idx
    #      idx <- idx + 1
    #    }
    #  }
    #}
    #edgelist_conv <- matrix(NA,nrow=reh$M,ncol=2)
    #for(m in 1:reh$M){
    #  edgelist_conv[m,] <- reh$risksetMatrix[reh$edgelist[m,2]+1,]
    #}
    #reh_estimate <- reh
    #reh$edgelist <- data.frame(time=reh$edgelist[,1],actor1=edgelist_conv[,1],actor2=edgelist_conv[,2],type=rep(0,reh$M),weight=reh$edgelist[,3])
    #rm(n_row,n_col,idx,m,edgelist_conv)
    # end #
    # adapted temporarily to the current version of the package "bremory"

    # ... Starting step-wise model counter
    q <- 1
    # ... Starting rejected models counter
    q_rej <- 0
            

    # ... Setting up progress bar settings (based on the number of generated step-wise functions)
    pb <- txtProgressBar(min = 1, max = stepwiseModels$Q, style = 3)

    cat("\n Starting estimating the bag of models ! \n")

    while(q <= stepwiseModels$Q){  
        start_time_q <- Sys.time()
        modelStatus <- FALSE
        while(!modelStatus){
          if(generateIntervalsFlag){
            env$statisticsREH$widths_q <-  bremory::intervals(K= stepwiseModels$intervals$K[q],
                                                        minWidth = intervals$minWidth,
                                                        maxTime= intervals$maxTime,
                                                        intervals= stepwiseModels$intervals$type[q]) #na.omit(stepwiseModels$intervals$widths[q,])
            stepwiseModels$intervals$widths[q,1:(stepwiseModels$intervals$K[q]+1)] <- env$statisticsREH$widths_q
            env$statisticsREH$K_q <- stepwiseModels$intervals$K[q] 
          }
          else{
            env$statisticsREH$widths_q <-  stepwiseModels$intervals$widths[q,1:(stepwiseModels$intervals$K[q]+1)]
            env$statisticsREH$K_q <- stepwiseModels$intervals$K[q] 
          }

          # ... Calculating number of statistics (S) in the q-th model
          pars_q <- rep(NA,length(vars))
          for(j in 1:length(vars)){
            if(vars[j] %in% attr(smmf,"endos")){
              pars_q[j] <- env$statisticsREH$K_q
            }
            else{
              pars_q[j] <- dim(env$stats[[vars[j]]])[3]
            }
          }
          S <- ifelse(attr(smmf,"intercept"),1,0)
          S <- S + sum(apply(matrix(rep(pars_q,dim(attr(smmf,"factors"))[2]),ncol = dim(attr(smmf,"factors"))[2]) * attr(smmf,"factors"),2,function(x) prod(x[which(x>0)])))
     
          # if at least one endogenous statistic is specified in the linear predictor
          if(length(attr(smmf,"endos"))>0){
            # (0) arrange names of endogenous statistics 
            for(j in 1:length(attr(smmf,"endos"))){
              env$statisticsREH$stats_names[[attr(smmf,"endos")[j]]] <- paste(attr(smmf,"endos")[j],".",1:env$statisticsREH$K_q,sep="")
            }
            # (1) calculate endogenous statistics
            
            # getIntervals() will find per each time point the lower and the upper bound for each interval
            start_time <-  Sys.time()
            getIntervals(env = globalenv(), 
                          widths = env$statisticsREH$widths_q, 
                          time = cumsum(reh$intereventTime), 
                          M = reh$M, 
                          K_q = env$statisticsREH$K_q,
                          ncores = ncores)


            ### TO CHECK FIRST AND THEN REMOVE ###
            #env$statisticsREH$intervals[(env$statisticsREH$K_q+1),1:2] <- c(0,0) #only for the first interval at the second time point
            ### TO CHECK FIRST AND THEN REMOVE ###

            intervals_loc <- unique(env$statisticsREH$intervals[,1:2])

            

            #getCounts (with OpenMP parallelization in C++)
            if(any(attr(smmf,"endos") %in% first_order_statistics)){ # at least one first order endogenous statistics is specified, 
                                                                    # then run the computation of 'counts' to speed up the calculation of statistics in intervals                      
            env$statisticsREH$counts <- getCounts(dyad = attr(reh,"dyad")-1,
                                                              D = reh$D,
                                                              lbs_ubs = intervals_loc, 
                                                              ncores = ncores)



            indices <- getCountsIndex(intervals = env$statisticsREH$intervals, 
                                      counts = cbind(env$statisticsREH$counts[,1:2],
                                      c(0:(dim(env$statisticsREH$counts)[1]-1))),
                                      ncores = ncores)
            env$statisticsREH$intervals[,3] <-  indices
            rm(intervals_loc,indices)
            }
         
            ### updating endogenous statistics for q-th model ###
            #endogenous_stats <- 
            getIntervalStatistics(env = globalenv(),
                                  counts = env$statisticsREH$counts,
                                  intervals = env$statisticsREH$intervals,
                                  actor1 = reh$edgelist$actor1_ID-1,
                                  actor2 = reh$edgelist$actor2_ID-1,
                                  time = reh$edgelist$time,
                                  M = reh$M,
                                  N = reh$N,
                                  D = reh$D,
                                  K = env$statisticsREH$K_q,
                                  statistics = attr(smmf,"endos"),
                                  model = attr(reh,"model"), # at the moment only working for tie-oriented modeling
                                  senderRate = FALSE,
                                  ncores = ncores)
            env$statisticsREH$counts  <- NULL
            env$statisticsREH$intervals <- NULL
          }
          else{endogenous_stats <- list()}

          # (2) ARRANGE ARRAY "STATS"
          statistics <- arrangeStatistics(env = globalenv(), effectsMatrix = attr(smmf,"factors"), names = rownames(attr(smmf,"factors")), intercept = attr(smmf,"intercept"), M = reh$M, D = reh$D, S = S)
          # statistics <- statistics[firstEvent:lastEvent,,]
          dimnames(statistics)[[3]] <- getStatisticsNames(stats_names =  env$statisticsREH$stats_names, effectsMatrix = attr(smmf,"factors"), names = rownames(attr(smmf,"factors")), intercept = attr(smmf,"intercept"), M = reh$M, D = reh$D, S = S)
          class(statistics) <- c("remstats","tomstats")
       

          # (4) ESTIMATE MODEL PARAMETERS
          model_q <- tryCatch(remstimate::remstimate(reh = reh,
                                          stats = statistics,
                                          method =  "MLE",
                                          model = "tie",
                                          ncores = ncores), error = function(error_message) {NULL}) 
          statistics <- NULL      
          # update modelStatus 
          if(!is.null(model_q)){
              #Sigma <- tryCatch(as.matrix(solve(model_q$hessian)), 
              #                          error = function(error_message) {matrix(-1,nrow=dim(stats_rem)[3],ncol=dim(stats_rem)[3])})
            #if(isSymmetric(model_q$hessian) & (Sigma[1,1] != (-1)) ){
            #        modelStatus <- TRUE               
            #}
            modelStatus <- TRUE
          }                                                          
          
          if(modelStatus){
            # WAIC routine 
            if(WAIC){
              # generate from posterior multivariate normal approximation
              Sigma <- as.matrix(solve(model_q$hessian))
              time_points <- (stepwiseModels$elpdSettings$L+1):lastEvent           
              interevent_time <- reh$intereventTime[time_points]  
              stats_waic <- aperm(stats, perm = c(2,3,1))  
              events <- reh$rehBinary[time_points,]                     
              post_pars <- t(mvtnorm::rmvnorm(n = stepwiseModels$elpdSettings$n_sims_is, mean = model_q$coefficients, sigma = Sigma)) # dim = [pars x n_sim_is]

              #lpd <- matrix(NA,nrow=stepwiseModels$elpdSettings$n_sims_is,ncol=length(time_points))
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
              #stepwiseModels$WAIC$pWAIC[,q] <- elpd_waic_obj$pointwise[,2] #lpd_and_p_waic[,2]
              #stepwiseModels$WAIC$elpdWAIC[,q] <- elpd_waic_obj$pointwise[,1] #lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
              #stepwiseModels$WAIC$WAIC <- elpd_waic_obj$pointwise[,3]
              #rm(elpd_waic_obj)
              ### temporarily commented

              ## START my routine for WAIC
              lpd_and_p_waic <- bremory::getWAIC(pars = post_pars, stats = stats_waic[,,time_points], events = t(events), interevent_time = interevent_time)
              stepwiseModels$WAIC_mine$pWAIC[,q] <- lpd_and_p_waic[,2]
              stepwiseModels$WAIC_mine$elpdWAIC[,q] <- lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
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
            #                                L = lastEvent - 10, #reh$M-10, #stepwiseModels$elpdSettings$L #$L is in 
            #                                n_sim_is = stepwiseModels$elpdSettings$n_sims_is,
            #                                tau = stepwiseModels$elpdSettings$tau,
            #                                cores = ncores
            #)  

            #stepwiseModels$ELPD$reestimations_psis[q] <- psis_1_sap$counts               
            #stepwiseModels$ELPD$elpd_lfo[,q] <- psis_1_sap$elpd_lfo
            #stepwiseModels$ELPD$k_psis[,q] <- ifelse(is.null(psis_1_sap$k),rep(NA,9),psis_1_sap$k) #(reh$M - stepwiseModels$elpdSettings$L)
            #rm(psis_1_sap)
            #}
            

            #### storing results ####
            n_pars_loc <- dim(stats)[3]
            stepwiseModels$coefficients[[q]] <- as.vector(model_q$coefficients)
            stepwiseModels$vcov[[q]] <- as.matrix(solve(model_q$hessian))
            stepwiseModels$loglik[q] <- model_q$loglik
            stepwiseModels$BIC[q] <- model_q$BIC 

            rm(model_q) 
            # saving partial results every 5 models
            #if(q%%250==0)
            #{
                # cat('\n ... Saving partial results ... \n')
            #    stepwiseModels$q_rej <- q_rej
            #    save(stepwiseModels,envir = env, file = paste(filename,".RData",sep="")) # also save 'reh'

                # save intervals table
                # partial_result_1 <- data.frame(table(stepwiseModels$intervals$K[1:q]),iteration = q)
                # write.csv(partial_result_1, file = paste(filename,".csv",sep=""))
            #}
          }
          else{
            if(!generateIntervalsFlag) modelStatus <- TRUE # in this way if intervals are not supposed to be generated we force to go on to the next q
            q_rej <- q_rej +1
            widths_q_rej_loc <- rep(NA,dim(stepwiseModels$intervals$widths)[2])
            widths_q_rej_loc[1:(env$statisticsREH$K_q+1)] <- env$statisticsREH$widths_q
            stepwiseModels$rejected$widths <-  rbind(stepwiseModels$rejected$widths, widths_q_rej_loc) # width rejected
            stepwiseModels$rejected$K <-  c(stepwiseModels$rejected$K, env$statisticsREH$K_q) # K rejected
            stepwiseModels$rejected$type <-  c(stepwiseModels$rejected$type, stepwiseModels$intervals$type[q]) # type rejected
          }

        }
        end_time_q <- Sys.time()
        stepwiseModels$time_sim$value[q] <- as.numeric(end_time_q - start_time_q)
        stepwiseModels$time_sim$unit[q] <- attr(end_time_q - start_time_q,"units")
        # update q index
        q <- q+1
        # Progress bar         
        setTxtProgressBar(pb, (q-1))
    }
    stepwiseModels$q_rej <- q_rej

    cat('\n ... Saving final results ... \n')

    
    # Check for rejected model when 'intervals' argument consist of an external set of intervals
    if(!generateIntervalsFlag){ 
      modelNA <- na.omit(stepwiseModels$BIC)
      modelNA <- as.vector(attr(modelNA,"na.action"))
      if(!is.null(modelNA)){
        stepwiseModels$Q <- (stepwiseModels$Q - length(modelNA))
        stepwiseModels$intervals$widths <- stepwiseModels$intervals$widths[-modelNA,]
        stepwiseModels$intervals$K <- stepwiseModels$intervals$K[-modelNA]
        if(!is.null(stepwiseModels$intervals$type)){
          stepwiseModels$intervals$type <- stepwiseModels$intervals$type[-modelNA]
        }
        stepwiseModels$coefficients <- stepwiseModels$coefficients[-modelNA,]
        stepwiseModels$vcov <- stepwiseModels$vcov[,,-modelNA]
        stepwiseModels$loglik <- stepwiseModels$loglik[-modelNA]
        stepwiseModels$BIC <- stepwiseModels$BIC[-modelNA]


        ## remstimate finalizing results ##
        #stepwiseModels$coefficients_remstimate <- stepwiseModels$coefficients_remstimate[-modelNA,]
        #stepwiseModels$vcov_remstimate <- stepwiseModels$vcov_remstimate[,,-modelNA]
        ## remstimate finalizing results ##

        # WAIC
        if(WAIC){
          #stepwiseModels$WAIC$pWAIC <- stepwiseModels$WAIC$pWAIC[,-modelNA]
          #stepwiseModels$WAIC$elpdWAIC <- stepwiseModels$WAIC$elpdWAIC[,-modelNA]
        }
        # ELPD
        if(ELPD){
          #stepwiseModels$ELPD$reestimations_psis <- stepwiseModels$ELPD$reestimations_psis[-modelNA]             
          #stepwiseModels$ELPD$elpd_lfo <- stepwiseModels$ELPD$elpd_lfo[,-modelNA]
          #stepwiseModels$ELPD$k_psis <- stepwiseModels$ELPD$k_psis[,-modelNA]
        }
        warning("x models were rejected. Therefore, the number of stepwise models per number of intervals might differ")
      }
    }

    # class and structure definition 
    str_out <- structure(list(
                            Q = stepwiseModels$Q,
                            intervals = stepwiseModels$intervals,
                            coef = stepwiseModels$coefficients,
                            vcov = stepwiseModels$vcov,
                            #coef_remstimate = stepwiseModels$coefficients_remstimate, ## remstimate
                            #vcov_remstimate = stepwiseModels$vcov_remstimate, ## remstimate
                            loglik = stepwiseModels$loglik,
                            BIC = stepwiseModels$BIC,
                            #WAIC = stepwiseModels$WAIC,
                            WAIC_mine = stepwiseModels$WAIC_mine,
                            ELPD = stepwiseModels$ELPD,
                            stats_names_endo = stepwiseModels$stats_names_endo, ## to change
                            stats_names_exo = stepwiseModels$stats_names_exo, ## to change
                            time_sim = stepwiseModels$time_sim
                            ), class="smm")

    attr(str_out, "rejected") <- stepwiseModels$rejected
    attr(str_out, "window") <- c(stepwiseModels$firstEvent, stepwiseModels$lastEvent)
    attr(str_out, "elpdSettings") <- stepwiseModels$elpdSettings
    attr(str_out, "formula") <- formula # this will substitute 'stats_names_endo' and 'stats_names_exo'

    rm(stepwiseModels)
    rm(statisticsREH, stats, envir = env)
    return(str_out)
}

#######################################################################################
#######################################################################################
##########(START)             Methods for `smm` object             (START)#############
#######################################################################################
#######################################################################################
 
#' @title merge.smm
#' @rdname merge.smm
#' @description A function that given a list of smm object, it merges into a single 'smm' object
#' @param x a list of \code{smm} objects
#' @method merge smm
#' 
#' @export
merge.smm <- function(x){
  if(inherits(x,"list")){
      list_class <- unlist(lapply(x,class))
      if(all(list_class) == "smm"){
        S <- length(x)
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
      if(!inherits(x,"smm")){
          stop("argument smm must be a list of objects of class 'smm'")
      }
      else{
        return(x) 
      }
  }
}

#######################################################################################
#######################################################################################
##########(END)             Methods for `smm` object             (END)#################
#######################################################################################
#######################################################################################