#' statisticsREH
#'
#' A function which creates an environment to...
#'
#' @param env environment where the user is currently working (it will become and  reh object and env will be set as last input because it is still important since statisticsREH is an environment)
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param K_range range as to the number of intervals that are going to be generated (argument NO MORE USED)
#'
#' @return  size of the specified environment.
#'
statisticsREH <- function(env = globalenv(), effects = NULL, K_range = c(2:6)){
    
    # Some initials check that can stop the function (to add further one as to 'K_range' argument)
    if(is.null(effects)){stop("the argument 'effects' must be specified \n (see documentation for further details).")}
        else{
            if(!is.list(effects)) {stop("user must provide a list of two objects: \n --'endogenous' (vector of strings with effects names); \n -- 'exogenous' (array with exogenous statistics).\n")}
                else{
                        if(!setequal(names(effects),c("endogenous","exogenous"))){stop(" 'effects' must be a list of two objects \n --'endogenous' (vector of strings with endogenous effects names); \n -- 'exogenous' (array with exogenous statistics).\n")}
                    }
            }
    # Do not change these vector of names (!!)
    dyadic_endogenous_effects_names <- c("SndSnd","RecSnd","IDSnd","IDRec","ODSnd","ODRec")
    triadic_endogenous_effects_names <- c("CClosure","TClosure")
    effects_names <- c(dyadic_endogenous_effects_names, triadic_endogenous_effects_names)


    # Check effects names (they must be the same as the names given by default (see ' effects_names'  vector))
    status <- as.double(all(sapply(1:length(effects$endogenous), function(x){any(effects_names==effects$endogenous[x])})))  
    if(status == 0){stop("At least one of the effects is not typed according to the names in the documentation \n (see documentation for details).")}
    #names_dyadic <- effects$endogenous[sapply(1:length(effects$endogenous), function(x){any(dyadic_endogenous_effects_names==effects$endogenous[x])})]

    # Creating environment 'statisticsREH' where temporary infos about the q-th model in the set will be stored
    env$statisticsREH <- new.env()
    env$statisticsREH$endo_effects <- effects$endogenous 
    env$statisticsREH$P <- length(effects$endogenous) # number of endogenous variables

    # Creating Environment 'binaryREH' with binary (1/0) matrices
    env$statisticsREH$binaryREH <- new.env()
    
    # dyadicREH is an M*[N*(N-1)] matrix with 1 where the event happened, 0 otherwise 
    # (useful for the remCpp function and for the computation of dyadic endogenous statistics)
    env$statisticsREH$binaryREH$dyadicREH <- getBinaryREH(M = env$initializeREH$M,
                                                            N = env$initializeREH$N,
                                                            edgelist = env$initializeREH$edgelist,
                                                            riskset_matrix = env$initializeREH$riskset_matrix)
    if(!is.null(env$initializeREH$weights))
    {
        weights_loc <- matrix(rep(env$initializeREH$weights,env$initializeREH$n_dyads), nrow = env$initializeREH$M,ncol = env$initializeREH$n_dyads)
        env$statisticsREH$binaryREH$dyadicREH <- weights_loc*env$statisticsREH$binaryREH$dyadicREH
    }
    # create names for endongenous and exogenous variables (maybe useful for a new arrangement of mles and vcov betas output)
    #rep_names <- rep(c(effects$endogenous), each =  max(K_range))
    #rep_int_index <- rep(1:max(K_range),times = env$statisticsREH$P)
    #env$statisticsREH$names_endogenous <- sapply((1:(env$statisticsREH$P*max(K_range))),
    #    function(x) { paste(rep_names[x],rep_int_index[x],sep="_")
    #})

    if(!is.null(effects$exogenous)){
        if(!is.array(effects$exogenous)) {stop("the object 'exogenous' must be an array with the exogenous statistics of interest.")}
            else{
                if(is.null(dimnames(effects$exogenous)[[3]])){dimnames(effects$exogenous)[[3]] <- sapply(1:dim(effects$exogenous)[3],function(x) paste("V",x,sep=""))}
                env$statisticsREH$exogenous_stats <- effects$exogenous[,env$initializeREH$dyad_position_array,] #dyad_position_array is important to rearrage columns according to the new dyads order
                env$statisticsREH$S <- dim(effects$exogenous)[3]
            }
    }
    else{env$statisticsREH$exogenous_stats <- NULL
         env$statisticsREH$S <- 0}


    # Creating useful NULL objects (for other routines)
    env$statisticsREH$intervals_backward <- NULL
}



#' statisticsREHOLD
#'
#' A function which creates an environment to...
#'
#' @param env environment where the user is currently working
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is an array [M*dyads*exogenous_effects]
#' @param K_range range as to the number of intervals that are going to be generated (argument NO MORE USED)
#'
#' @return  size of the specified environment.
#'
statisticsREHOLD <- function(env = globalenv(), effects = NULL, K_range = c(2:6)){
    
    # Some initials check that can stop the function (to add further one as to 'K_range' argument)
    if(is.null(effects)){stop("the argument 'effects' must be specified \n (see documentation for further details).")}
        else{
            if(!is.list(effects)) {stop("user must provide a list of two objects: \n --'endogenous' (vector of strings with effects names); \n -- 'exogenous' (array with exogenous statistics).\n")}
                else{
                        if(!setequal(names(effects),c("endogenous","exogenous"))){stop(" 'effects' must be a list of two objects \n --'endogenous' (vector of strings with endogenous effects names); \n -- 'exogenous' (array with exogenous statistics).\n")}
                    }
            }
    # Do not change these vector of names (!!)
    dyadic_endogenous_effects_names <- c("SndSnd","RecSnd","IDSnd","IDRec","ODSnd","ODRec")
    triadic_endogenous_effects_names <- c("CClosure","TClosure")
    effects_names <- c(dyadic_endogenous_effects_names, triadic_endogenous_effects_names)


    # Check effects names (they must be the same as the names given by default (see ' effects_names'  vector))
    status <- as.double(all(sapply(1:length(effects$endogenous), function(x){any(effects_names==effects$endogenous[x])})))  
    if(status == 0){stop("At least one of the effects is not typed according to the names in the documentation \n (see documentation for details).")}
    #names_dyadic <- effects$endogenous[sapply(1:length(effects$endogenous), function(x){any(dyadic_endogenous_effects_names==effects$endogenous[x])})]

    # Creating environment 'statisticsREHOLD' where temporary infos about the q-th model in the set will be stored
    env$statisticsREHOLD <- new.env()
    env$statisticsREHOLD$endo_effects <- effects$endogenous 
    env$statisticsREHOLD$P <- length(effects$endogenous) # number of endogenous variables

    # Creating Environment 'binaryREH' with binary (1/0) matrices
    env$statisticsREHOLD$binaryREH <- new.env()
    
    # dyadicREH is an M*[N*(N-1)] matrix with 1 where the event happened, 0 otherwise 
    # (useful for the remCpp function and for the computation of dyadic endogenous statistics)
    env$statisticsREHOLD$binaryREH$dyadicREH <- getBinaryREH_old(M = env$initializeREH$M,
                                                            N = env$initializeREH$N,
                                                            edgelist = env$initializeREH$edgelist,
                                                            riskset_matrix = env$initializeREH$riskset_matrix)
    if(!is.null(env$initializeREH$weights))
    {
        weights_loc <- matrix(rep(env$initializeREH$weights,env$initializeREH$n_dyads), nrow = env$initializeREH$M,ncol = env$initializeREH$n_dyads)
        env$statisticsREHOLD$binaryREH$dyadicREH <- weights_loc*env$statisticsREHOLD$binaryREH$dyadicREH
    }
    # create names for endongenous and exogenous variables (maybe useful for a new arrangement of mles and vcov betas output)
    #rep_names <- rep(c(effects$endogenous), each =  max(K_range))
    #rep_int_index <- rep(1:max(K_range),times = env$statisticsREHOLD$P)
    #env$statisticsREHOLD$names_endogenous <- sapply((1:(env$statisticsREHOLD$P*max(K_range))),
    #    function(x) { paste(rep_names[x],rep_int_index[x],sep="_")
    #})

    if(!is.null(effects$exogenous)){
        if(!is.array(effects$exogenous)) {stop("the object 'exogenous' must be an array with the exogenous statistics of interest.")}
            else{
                if(is.null(dimnames(effects$exogenous)[[3]])){dimnames(effects$exogenous)[[3]] <- sapply(1:dim(effects$exogenous)[3],function(x) paste("V",x,sep=""))}
                env$statisticsREHOLD$exogenous_stats <- effects$exogenous[,env$initializeREH$dyad_position_array,] #dyad_position_array is important to rearrage columns according to the new dyads order
                env$statisticsREHOLD$S <- dim(effects$exogenous)[3]
            }
    }
    else{env$statisticsREHOLD$exogenous_stats <- NULL
         env$statisticsREHOLD$S <- 0}


    # Creating useful NULL objects (for other routines)
    env$statisticsREHOLD$intervals_backward <- NULL
}