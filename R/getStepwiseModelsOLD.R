#' getStepwiseModelsOLD
#'
#' A function which run the estimation of Q models and updates the environment 'stepwiseModelsREH' 
#'
#' @param env environment where the user is currently working
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param K range as to the number of intervals that are going to be generated.
#' @param max_width maximum for width
#' @param min_diff minimum width (diff will be changed with width)
#' @param nsim_per_K how many simulations per number of intervals
#' @param L number of events to use as history to estimate the predictive pointwise density of future events (set to M-10 by default)
#' @param A_SAP A-Steps-Ahead-Predictions is the number of steps ahead we want consider in the prediction (set to 1 by default)
#' @param file_name a string indicating the name of the file .RData that will save results each 20 models.
#' @param WAIC default is TRUE and is needed where WAIC weights have to be estimated
#' @param n_cores number of threads to create
#' @param first_event integer indicating the event at which start the event sequence and statistic in the estimation stage
#' @param last_event integer indicating the event at which stop the event sequence and statistic in the estimation stage
#' @param widths_input matirx of sequences of widths from the user
#' @param K_input vector of number of steps per each sequence of widths from the user
#' @param interval_type_input vector of interval types per each sequence of widths from the user
#'
#' @return  the function doesn't return any output but updates objects inside the environment 'stepwiseModelsREH'
#'
#' @export
getStepwiseModelsOLD <- function(env = globalenv(),
                              effects = NULL, 
                              K = c(2:6), 
                              max_width = NULL, 
                              min_diff = NULL,
                              nsim_per_K = 101,
                              L = NULL, 
                              A_SAP = 1, 
                              file_name = "simulation_0", 
                              WAIC = TRUE, 
                              n_cores = 10,
                              first_event = NULL,
                              last_event = NULL,
                              widths_input = NULL,
                              K_input = NULL,
                              interval_type_input = NULL)
{
    # check for initializeREH environment
    if(is.null(env$initializeREH)){stop("initializeREH() function must be run first.")}

    # create a new environment 'statisticsREHOLD' (the largest in the estimation process and continuously updated according to the characteristic of the stepswise model
    # at a specific iteration)
    statisticsREHOLD(env = env, effects = effects)

    # create environment 'stepwiseModelsREH' where, widths are stored as well as the output of each model (information useful to the estimation stage)
    env$stepwiseModelsREH <- new.env()

    # generating random widths according to K and nsim_per_K parameters
    if(is.null(widths_input)){
      env$stepwiseModelsREH$widths <- matrix(NA, nrow = length(K)*nsim_per_K, ncol = (max(K)+1)) # length(K)*(nsim_per_K-1)*2+length(K) is specific of the case where  for example we want to generate 150 increasing 150 decreasing and 1 equal size , thus for K=c(3,4,5) 3*(151-1)+3 where 151 is nsim_per_K
      env$stepwiseModelsREH$K <- rep(K,each = nsim_per_K)
      env$stepwiseModelsREH$interval_type <- rep(c("equal",rep("increasing",(nsim_per_K-1)/2),rep("decreasing",(nsim_per_K-1)/2)),length(K))
      # shuffling order [now] : this is done not to increasingly overload the cpu and/or the ram (in this way sometimes a complex model is estimated, other times a much simpler one is optimized)
      shuffled_indices <- sample(x=1:dim(env$stepwiseModelsREH$widths)[1],size=dim(env$stepwiseModelsREH$widths)[1])
      env$stepwiseModelsREH$K <- env$stepwiseModelsREH$K[shuffled_indices]
      env$stepwiseModelsREH$interval_type <- env$stepwiseModelsREH$interval_type[shuffled_indices]
    }
    else{
      # a check over the dimensions of the three object should be coded but this feature is temporary: in future the user will switch between methods of interval generation
      env$stepwiseModelsREH$widths <- widths_input
      env$stepwiseModelsREH$K <- K_input
      env$stepwiseModelsREH$interval_type <- interval_type_input
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
      env$stepwiseModelsREH$last_event <- last_event <- env$initializeREH$M
    }
    # print we could also remove
    print(c(env$stepwiseModelsREH$first_event,env$stepwiseModelsREH$last_event,max_width))
    # storing dimensions and variables name
 
    env$stepwiseModelsREH$Q <- dim(env$stepwiseModelsREH$widths)[1]    
    env$stepwiseModelsREH$L <- c(first_event:last_event)[max(which(env$initializeREH$t[first_event:last_event]<(max(env$initializeREH$t[first_event:last_event])-max_width)))] 

    # M-10 was the FIXED VALUE (previously)
    # max(env$initializeREH$t[first_event:last_event]) could be avoided and changed with env$initializeREH$t[last_event] since the time variable is supposed to be ordered

    env$stepwiseModelsREH$stats_names_endo <- env$statisticsREHOLD$endo_effects # or effects$endogenous
    env$stepwiseModelsREH$stats_names_exo <- dimnames(env$statisticsREHOLD$exogenous_stats)[[3]]

    ## BEGIN rem supplist 
    supplist <- list() 
    supplist[[1]]<-matrix(TRUE, nrow = length(first_event:last_event), ncol = env$initializeREH$n_dyads)
    ## END rem supplist

    # allocating space for each model result (Q models in total)
    env$stepwiseModelsREH$mle_betas <- matrix(NA, nrow = env$stepwiseModelsREH$Q, ncol = (( max(K) * env$statisticsREHOLD$P) + env$statisticsREHOLD$S)) 
    env$stepwiseModelsREH$vcov_betas <- array(NA , dim = c((( max(K) * env$statisticsREHOLD$P) + env$statisticsREHOLD$S),
    (( max(K) * env$statisticsREHOLD$P) + env$statisticsREHOLD$S), env$stepwiseModelsREH$Q))
    env$stepwiseModelsREH$loglik <- env$stepwiseModelsREH$BIC <- rep(NA, env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$loglik_no_endo

    ## ELPD LFO memory allocation
    #env$stepwiseModelsREH$reestimations_psis <- rep(NA, env$stepwiseModelsREH$Q)
    #env$stepwiseModelsREH$elpd_lfo <- matrix(NA, nrow = 10, ncol = env$stepwiseModelsREH$Q) # nrow = (env$initializeREH$M - env$stepwiseModelsREH$L)
    #env$stepwiseModelsREH$k_psis <- matrix(NA, nrow = 9, ncol = env$stepwiseModelsREH$Q) # nrow = ((env$initializeREH$M - env$stepwiseModelsREH$L)-1)

    if(WAIC){
      # for multi-steps models  :env$initializeREH$M
      env$stepwiseModelsREH$elpd_waic <- matrix(NA, nrow = (length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L)), ncol= env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$elpd_waic_no_endo
      env$stepwiseModelsREH$p_waic <- matrix(NA, nrow = (length(first_event:last_event)-length(first_event:env$stepwiseModelsREH$L)), ncol= env$stepwiseModelsREH$Q) # env$stepwiseModelsREH$p_waic_no_endo
    }
    ###### START FIXED VALUES (MAYBE TO CHANGE) #####
    env$stepwiseModelsREH$n_sims_is <- 1e03
    env$stepwiseModelsREH$tau <- 0.7
    ###### END FIXED VALUES (MAYBE TO CHANGE) #####

    #### START PARALLELIZATION SETTINGS
    # set Cores for parallelization purpose
    # n_cores <- floor(detectCores()) - 2
    # assign object 'n_cores' to environment initializeREH
    env$initializeREH$n_cores <- n_cores
    #### END PARALLELIZATION SETTINGS
    env$statisticsREHOLD$edgelist <- data.matrix(convertEdgelist(edgelist = env$initializeREH$edgelist, 
                                                              riskset = env$initializeREH$riskset,  
                                                              actors = env$initializeREH$actors,
                                                              rem = TRUE))
    env$statisticsREHOLD$counts <- new.env()
    

    ## when widths are provided from the .GlobalEnv
    ## env$stepwiseModelsREH$widths <- globalenv()$widths_loc # just for now
    ## env$stepwiseModelsREH$K <- globalenv()$K_loc # just for now
    ## env$stepwiseModelsREH$Q <- length(env$stepwiseModelsREH$K)
    #env$stepwiseModelsREH$widths[1,] <- c(0,0.2,1.3,3.5,7,NA,NA,NA) #only for the stepwise case (change interval widths with the true ones)
    #env$stepwiseModelsREH$K[1] <- 4

    ###### START ONLY EXOGENOUS PART
    #model_only_exo <- relevent::rem(eventlist = env$statisticsREHOLD$edgelist, 
    #                      statslist = env$statisticsREHOLD$exogenous_stats,
    #                      supplist = supplist,
    #                      timing = "interval",
    #                      estimator = "MLE")
    #env$stepwiseModelsREH$loglik_only_exo <- model_only_exo$loglik 
    #env$stepwiseModelsREH$model_only_exo <- model_only_exo

    #if(WAIC){
    #      # generate from posterior multivariate normal approximation
    #      Sigma <- tryCatch(as.matrix(solve(model_only_exo$hessian)),
    #                                error = function(error_message) {matrix(-1,nrow=dim(env$statisticsREHOLD$exogenous_stats)[2],ncol=dim(env$statisticsREHOLD$exogenous_stats)[2])})
    #      time_points <- (env$stepwiseModelsREH$L+1):env$initializeREH$M           
    #      interevent_time <- diff(c(0,env$initializeREH$t))[time_points]  
    #      stats_waic_exo <- aperm(env$statisticsREHOLD$exogenous_stats, perm = c(2,3,1))  
    #      events <- env$statisticsREHOLD$binaryREH$dyadicREH[time_points,]                     
    #      if(Sigma[1,1] != (-1)){
    #        post_pars_exo <- t(mvtnorm::rmvnorm(n = env$stepwiseModelsREH$n_sims_is, mean = model_only_exo$coef, sigma = Sigma)) # dim = [pars x n_sim_is]
    #
    #        lpd <- bremory::lpdWAIC(pars = post_pars_exo, 
    #                                stats = stats_waic_exo[,,time_points], 
    #                                events = t(events), 
    #                                interevent_time = interevent_time) 
    #        elpd_waic_obj <- loo::waic(lpd)
    #        env$stepwiseModelsREH$p_waic_exo <- elpd_waic_obj$pointwise[,2] #lpd_and_p_waic[,2]
    #        env$stepwiseModelsREH$elpd_waic_exo <- elpd_waic_obj$pointwise[,1] #lpd_and_p_waic[,1] - lpd_and_p_waic[,2]
    #        rm(lpd)
    #      
    #    }}
    ###### END ONLY EXOGENOUS PART
    q <- 1
    q_rej <- 0
    cat("\n Starting estimating the bag of models ! \n")
    while(q <= env$stepwiseModelsREH$Q)
    { 
    

        #if(length(env$statisticsREHOLD$widths_q)==0) {
        #  env$statisticsREHOLD$widths_q <- c(0,max_width) #Inf
        #  env$statisticsREHOLD$K_q <- 1
        #}
        model_status <- FALSE
        while(!model_status){

        if(is.null(widths_input)){
          env$statisticsREHOLD$widths_q <-  generateWidths(K= env$stepwiseModelsREH$K[q],
                                                      min_diff = min_diff,
                                                      max_width= max_width,
                                                      intervals= env$stepwiseModelsREH$interval_type[q]) #na.omit(env$stepwiseModelsREH$widths[q,])
          env$stepwiseModelsREH$widths[q,1:(env$stepwiseModelsREH$K[q]+1)] <- env$statisticsREHOLD$widths_q
          env$statisticsREHOLD$K_q <- env$stepwiseModelsREH$K[q] 
        }
        else{
          env$statisticsREHOLD$widths_q <-  env$stepwiseModelsREH$widths[q,1:(env$stepwiseModelsREH$K[q]+1)]
          env$statisticsREHOLD$K_q <- env$stepwiseModelsREH$K[q] 
        }

          
          start <- Sys.time() 
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
          ## assigning object counts to environment object insid statisticsREHOLD
          #env$statisticsREHOLD$counts <- counts 
          #### END PARALLELIZATION WITH 'PARALLEL' PKG (by using parLapply() function)

          #getCountsOMP (OpenMP parallelization in C++)
          if(length(env$statisticsREHOLD$endo_effects)>0){
            # getIntervals() will find per each time point the lower and the upper bound for each interval and will automatically assign it to the object 'intervals_backwards'
            getIntervals_old(env = globalenv()) 

            #env$statisticsREHOLD$intervals_backward[(env$statisticsREHOLD$K_q+1),1:2] <- c(0,0) #only for the first interval at the second time point
            intervals_loc <- unique(env$statisticsREHOLD$intervals_backward[,1:2])

            for(l in 1:length(names(env$statisticsREHOLD$binaryREH))){
              env$statisticsREHOLD$counts[[names(env$statisticsREHOLD$binaryREH)[l]]] <- getCountsOMP_old(binaryREH = env$statisticsREHOLD$binaryREH[[names(env$statisticsREHOLD$binaryREH)[l]]], 
                                                                                                lbs_ubs = intervals_loc, 
                                                                                                n_cores = env$initializeREH$n_cores)
            }

            indices <- getCountsIndex_old(intervals = env$statisticsREHOLD$intervals_backward, 
                                      counts = cbind(env$statisticsREHOLD$counts[[names(env$statisticsREHOLD$binaryREH)[1]]][,1:2],
                                      c(0:(dim(env$statisticsREHOLD$counts[[names(env$statisticsREHOLD$binaryREH)[1]]])[1]-1))))
            env$statisticsREHOLD$intervals_backward[,3] <-  indices

            rm(intervals_loc,indices)

            ### calculating effects ###
            endogenous_stats <- getEndoEffects_old(env = env)
            # create names for endongenous and exogenous variables
            rep_names <- rep(c(env$statisticsREHOLD$endo_effects), each = env$statisticsREHOLD$K_q)
            rep_int_index <- rep(1:env$statisticsREHOLD$K_q,times = env$statisticsREHOLD$P)
            env$statisticsREHOLD$names_endogenous <- sapply((1:(env$statisticsREHOLD$P*env$statisticsREHOLD$K_q)),
                function(x) { paste(rep_names[x],rep_int_index[x],sep="_")
            })
            dimnames(endogenous_stats)[[3]] <- env$statisticsREHOLD$names_endogenous


            rm(list=ls(envir = env$statisticsREHOLD$counts), envir = env$statisticsREHOLD$counts)
            env$statisticsREHOLD$intervals_backward <- NULL
          }
          else{endogenous_stats <- NULL}
          stats <- abind(endogenous_stats, env$statisticsREHOLD$exogenous_stats, along = 3)
          
          ##### relevent::rem function #####  
          # estimating model with relevent package #
          stats_rem <- stats[first_event:last_event,,] 
          model_loc <- tryCatch(relevent::rem(eventlist = env$statisticsREHOLD$edgelist[first_event:last_event,], 
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
               ###-### {
                    model_status <- TRUE
               ###-### }
                            
            }
          }
          
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
              Sigma <- tryCatch(as.matrix(solve(model_loc$hessian)), 
                                        error = function(error_message) {matrix(-1,nrow=dim(stats)[3],ncol=dim(stats)[3])})
              time_points <- (env$stepwiseModelsREH$L+1):last_event           


              interevent_time <- diff(c(0,env$initializeREH$t))[time_points]  
              stats_waic <- aperm(stats, perm = c(2,3,1))  
              events <- env$statisticsREHOLD$binaryREH$dyadicREH[time_points,]                     
              if((Sigma[1,1] != (-1)) ){###-### & isSymmetric(Sigma)
                post_pars <- t(mvtnorm::rmvnorm(n = env$stepwiseModelsREH$n_sims_is, mean = model_loc$coef, sigma = Sigma)) # dim = [pars x n_sim_is]

                #lpd <- matrix(NA,nrow=env$stepwiseModelsREH$n_sims_is,ncol=length(time_points))
                #for(m in 1:length(time_points))
                #{
                #event <- numeric(env$initializeREH$n_dyads)
                #event[env$statisticsREHOLD$edgelist[time_points[m],1]] <- 1
                #lpd[,m] <- apply(post_pars,2,function(y) lpd(pars = y,
                #                                                stats = stats[time_points[m],,],
                #                                                event = event,
                #                                                interevent_time = env$initializeREH$t[time_points[m]]-env$initializeREH$t[time_points[m]-1]))
                #}  
                lpd <- bremory::lpdWAIC_old(pars = post_pars, 
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
                #  event[env$statisticsREHOLD$edgelist[time_points[m],1]] <- 1
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
            #psis_1_sap <- PSIS_LFO_CV_1_SAP(n_dyads = env$initializeREH$n_dyads,
            #                                eventlist = env$statisticsREHOLD$edgelist, 
            #                                statslist = stats,
            #                                t = env$initializeREH$t,
            #                                first_event = first_event,
            #                                last_event = last_event, #env$initializeREH$M, 
            #                                L = last_event - 10, #env$initializeREH$M-10, #env$initializeREH$L 
            #                                n_sim_is = env$stepwiseModelsREH$n_sims_is,
            #                                tau = env$stepwiseModelsREH$tau,
            #                                cores = env$initializeREH$n_cores)                   
            end <- Sys.time()
            # print message 
            cat("\n ",file_name,": Iteration ",q,"(rejected up to now ",q_rej," ) with ",env$statisticsREHOLD$K_q," intervals completed in --> ", end-start, attr(end-start,'units'),
            ".\n")  ## and ", psis_1_sap$counts, "  no. re-estimations ##
            # print message
            

            #### storing results ####
            n_pars_loc <- dim(stats)[3]
            env$stepwiseModelsREH$mle_betas[q,c(1:n_pars_loc)] <- as.vector(model_loc$coef)
            env$stepwiseModelsREH$vcov_betas[c(1:n_pars_loc),c(1:n_pars_loc),q] <- as.matrix(solve(model_loc$hessian))
            env$stepwiseModelsREH$loglik[q] <- model_loc$loglik
            env$stepwiseModelsREH$BIC[q] <- model_loc$BIC 
            #env$stepwiseModelsREH$reestimations_psis[q] <- psis_1_sap$counts               
            #env$stepwiseModelsREH$elpd_lfo[,q] <- psis_1_sap$elpd_lfo
            #env$stepwiseModelsREH$k_psis[,q] <- ifelse(is.null(psis_1_sap$k),rep(NA,9),psis_1_sap$k) #(env$initializeREH$M - env$stepwiseModelsREH$L)

            rm(endogenous_stats,n_pars_loc,model_loc,stats) #,psis_1_sap
                    # saving partial results every 5 models
            if(q%%5==0)
            {
                cat('\n ... Saving partial results ... \n')
                env$stepwiseModelsREH$q_rej <- q_rej
                save(initializeREH,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

                # save intervals table
                partial_result_1 <- data.frame(table(env$stepwiseModelsREH$K[1:q]),iteration = q)
                write.csv(partial_result_1, file = paste(file_name,".csv",sep=""))
            }
          }
          else{
            env$stepwiseModelsREH$widths_rejected <-  rbind(env$stepwiseModelsREH$widths_rejected, env$statisticsREHOLD$widths_q) # widths_rejected
            env$stepwiseModelsREH$K_rejected <-  c(env$stepwiseModelsREH$K_rejected, env$statisticsREHOLD$K_q) # K_rejected
            env$stepwiseModelsREH$type_rejected <-  c(env$stepwiseModelsREH$type_rejected, env$stepwiseModelsREH$interval_type[q]) # type_rejected
            print("rejected")
            q_rej <- q_rej + 1
          }

        }



        # update q index
        q <- q +1
    }

    cat('\n ... Saving final results ... \n')
    env$stepwiseModelsREH$q_rej <- q_rej
    save(initializeREH,stepwiseModelsREH,envir = env, file = paste(file_name,".RData",sep=""))

    rm(statisticsREHOLD, envir = env)
}
