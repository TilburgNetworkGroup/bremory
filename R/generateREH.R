#' generateREH
#'
#' A function which generate a relational event history according to network characteristics and statistics that are specified by input arguments
#'
#' @param N number of actors in the network.
#' @param max_time this parameter is not used in this function anymore. By setting the maximum value for the variable time we cannot define the number of final events in advance.
#' @param M number of events we want to generate.
#' @param exo_stats matrix of N*(N-1) rows with as many columns as the number of exogenous statistics.
#' @param beta_endo named list of endogenous statistics with decay type (decay_type) and parameters (pars). Names must match with the ones used within the function (see vector of names).
#' @param beta_exo vector of effects (its length must be equal to the number of columns of exo_stats).
#' @param seed used internally. If not specified a random seed is used and finally returned.
#'
#' @return  list of final dgelist, seed used to generate it and vector of offset values.
#'
#' @export
generateREH <- function(N, max_time = NULL, M, exo_stats = NULL, beta_endo, beta_exo = NULL, seed = NULL)
{
        # Check input arguments
        # if(any()) stop("user must define at least N, M and beta_endo")

        # defining riskset (no selfedges) and getting riskset arranged in a matrix with binaryREH column 
        # index corresponding to the specific dyad
        riskset <- matrix(na.omit(getRiskset(actors_id = c(0:(N-1)), N = N, selfedges = FALSE)),ncol=2)
        riskset_matrix <- getRisksetMatrix(riskset = riskset, N = N)  

        # Defining vector of endogenous effects names (useful for the Rcpp functions to update the endogenous statistics specified by the user)
        dyadic_endogenous_effects_names <- c("SndSnd","RecSnd","IDSnd","IDRec","ODSnd","ODRec","TClosure","CClosure")

        # Handling the endogenous stats input (beta_endo) 
        endogenous_stats <- NULL
        if(!is.null(beta_endo)){

                endogenous_stats <- names(beta_endo)
                names_dyadic <- endogenous_stats[sapply(1:length(endogenous_stats), function(x){any(dyadic_endogenous_effects_names==endogenous_stats[x])})]

                # binaryREH (= "binary Relational Event History") is an environment where to store dyadicREH matrix, which is useful for
                # the computation of endogenous statistics
                .GlobalEnv$binaryREH <- new.env() # using the global environment

                if(length(names_dyadic) > 0)
                {.GlobalEnv$binaryREH$dyadicREH <- matrix(0, nrow = M, ncol = N*(N-1))} # the riskset considers N*(N-1) possible dyads (no event type here)

                cat("Generating with endogenous statistics (!!)")
        }

        # Calculating the constant contribute of exogenous variables (non-time-varying) to the log(lambda) that is the logarithm of the event rate
        if(is.null(exo_stats)){
                log_lambda_exo <- cbind(rep(0,N*(N-1)))
        }
        else{
                log_lambda_exo <- exo_stats%*%beta_exo
        }

        # Initializing the logarithm of the event rate as regards endogenous effects
        log_lambda_endo <- cbind(rep(0,N*(N-1)))
        if(is.null(endogenous_stats)){
        endo_stats <- cbind(rep(0,N*(N-1)))
        }
        
        # Pre-allocating output of generated relational events, by row [time,sender,receiver]
        out <- matrix(0, nrow = M, ncol = 3)

        # Initializing counter variables
        m <- 1 # current m-th relational event

        # Initializing time variable
        t <- 0   

        # Progress bar settings (based on the number of generated events)
        pb <- txtProgressBar(min = 0, max = M, style = 3)

        # Initializing vector where to store the offset variable over time
        offset0 <- numeric(M)

        # Setting the seed and running the simulation
        if(is.null(seed)){seed  <- sample(x = c(1:1e04), size = 1)}
        set.seed(seed)
        while((m-1) < M){ 
                
                # t<(max_time+1e-03): this is just a reminder that we can also set the max_time as input and generate until the algorithm reaches the fixed time value. In this case we won't be able to control the final number of generated events anymore.

                # Calculating the log(lambda) per each dyad
                log_lambda <- log_lambda_endo + log_lambda_exo # sum of log-event-rates for endo_stats and exo_stats
                lambda <- exp(log_lambda) # event rates

                # Calculating the offset
                offset0[m] <- 1/max(lambda) 

                # Generating waiting time
                delta <- rexp(n = 1, rate = sum(lambda)) # *offset0[m] # the offset multiplies the sum(lambda)

                # Updating occurrence probabilities
                p <- (lambda/sum(lambda))       # dyad-specific probability 
                p[is.na(p)] <- 1e-30            # correcting for too low probabilities

             
                # Sampling the m-th dyadic event
                ij <- sample(c(1:(N*(N-1))), size = 1, prob  = p)

                # Updating the time variable 
                t <- t + delta 
                
                
                # Updating output matrix: [time,sender,receiver]
                out[m,] <- c(t,riskset[ij,])

                # this is a cat() that serves to me in order to monitor the trend of some quantities during the generation
                cat('\n',c(round(m),                            # m-th event
                round(as.numeric(out[m,1]),5),                  # time variable value at the m-th event
                round(as.numeric(as.vector(out[m,2:3]))),       # dyadic event generated at the m-th iteration
                max(p),                                         # max probability of occurrence at t_m
                round(median(lambda),3),                        # median of event-rates at t_m
                round(sum(lambda),3),                           # actual parameter for generating the waiting time
                round(offset0[m],3),                            # offset value at t_m 
                round(sum(lambda)*offset0[m],3)                 # modified parameter for generating the waiting time (offset)
                ))
                
                ### BEGIN ENDOGENOUS STATISTICS UPDATES 
                if(length(endogenous_stats) > 0){
                        # Updating dyadicREH matrix: 1/0 matrix (1 assigned to the dyad that occurred, 0 otherwise)
                        env$binaryREH[["dyadicREH"]][m,ij] <- 1
                        # update endo_stats (updateEndoEffects returns a vector of length N*(N-1) )
                        endo_stats <- updateEndoEffects(effects = endogenous_stats,
                                                        beta_endo = beta_endo,
                                                        binaryREH = .GlobalEnv$binaryREH, #rbind(binaryREH[1:m,]), 
                                                        elapsed_time = cbind(out[m,1]-out[1:m,1]),
                                                        edgelist = rbind(out[1:m,]),
                                                        riskset = riskset,
                                                        riskset_matrix = riskset_matrix, 
                                                        N = N,
                                                        M_partial = m)  #I could change name from M_partial to partial_sequence_length
                }
                ### END ENDOGENOUS STATISTICS UPDATES ###
                       
                # Updating contribute of endogenous statistics for the next event (the (m+1)-th) to be generated
                log_lambda_endo <- cbind(endo_stats)

                # Updating event counter
                m <- m + 1  # we go for generating the (m+1)-th event

                # Progress bar         
                # setTxtProgressBar(pb, (m-1)) #uncomment if you want to see the progress bar (but comment the cat() first !)

                # Should the generated waiting time (delta) be too short, the algorithm must stop
                if(delta<1e-10) {print("too short waiting time !");break}
        }
        
                
        # processing the output :
        rm(binaryREH, envir = .GlobalEnv) # removing the envirnoment 'binaryREH'

        out<- data.frame(out[1:(m-1),])
        names(out) <-c("time","sender","receiver")
        return(list(edgelist = out, seed = seed, offset = offset0))
}