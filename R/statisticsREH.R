#' statisticsREH
#'
#' A function which creates an environment to...
#'
#' @param formula formula from getStepwiseModels(...)
#' @param data data from getStepwiseModels
#' @param reh 'reh' object
#' @param effects list of two objects: effects$endogenous is a vector of strings with specified the effects of interest; effects$exogenous is a matrix [dyads*exogenous_effects]
#' @param env environment where the user is currently working
#'
#' @return  size of the specified environment.
#'
#' @export
statisticsREH <- function(formula, data, reh, effects, env){
    
    # Creating environment 'statisticsREH' where temporary infos about the linear predictor will be stored
    env$statisticsREH <- new.env()

    # Some initials check that can stop the function (to add further one as to 'K_range' argument)
    if(is.null(effects)){stop("the argument 'effects' must be specified \n (see documentation for further details)")}
        else{
            if(!is.list(effects)) {stop("user must provide a list of two objects: \n --'endogenous' (vector of strings with effects names) \n -- 'exogenous' (array with exogenous statistics)\n")}
                else{
                        if(!setequal(names(effects),c("endogenous","exogenous"))){stop(" 'effects' must be a list of two objects \n --'endogenous' (vector of strings with endogenous effects names) \n -- 'exogenous' (array with exogenous statistics)\n")}
                    }
            }
    # Do not change these vector of names (!!)
    effects_names <- c("SndSnd","RecSnd","IDSnd","IDRec","ODSnd","ODRec","CClosure","TClosure") # change vector of names to c("inertia","reciprocity","indegreeSender","indegreeReceiver","outdegreeSender","outdegreeReceiver","totaldegreeSender","totaldegreeReceiver","pshift","transitivityClosure","cyclicClosure","sendingClosure","receivingClosure"
    # "pshift" doesn't have memory 


    # Check effects names (they must be the same as the names given by default (see ' effects_names'  vector))
    status <- as.double(all(sapply(1:length(effects$endogenous), function(x){any(effects_names==effects$endogenous[x])})))  
    if(status == 0){stop("At least one of the effects is not typed according to the names in the documentation \n (see documentation for details)")}

    #############################################
    #### START EXOGENOUS VARIABLES PROCESSING ###
    #############################################
    exogenous_internal_names <- NULL
    if(!is.null(data)){
        if(!is.list(data)) stop("'data' must be a list")
        exogenous_dyadic <- exogenous_actor <- NULL
        # checking 'data$dyadic' ...
        if(!is.null(data$dyadic)){
            if(!is.list(data$dyadic)) stop("'data$dyadic' must be a list")
            else{
                if(!is.null(data$dyadic$dyad) & !is.null(data$dyadic$variables)){
                    if(!is.matrix(data$dyadic$dyad) | !is.data.frame(data$dyadic$dyad)) stop("'data$dyadic$dyad' must be either a matrix or a data.frame with columns indicating [actor1,actor2,type]")
                    match_names <- match(c("actor1","actor2","type"),dimnames(data$dyadic$dyad)[[2]])
                    if(is.unsorted(match_names) | anyNA(match_names)) stop("column names of 'data$dyadic$dyad' must be the following (in order): c('actor1','actor2','type')")
                    if(!is.array(data$dyadic$variables) | !is.data.frame(data$dyadic$variables)) stop("'data$dyadic$variables' must be either an array (with two [matrix] or more dimensions [3-d array]) or a data.frame in the simple case of constant exogenous dyadic variables")   
                    if(!is.null(dimnames(data$dyadic$variables)[[2]])) exogenous_dyadic <- dimnames(data$dyadic$variables)[[2]]  # columns in a matrix/data.frame/3-d array
                    else{stop("column names of 'dyadic$variables' are missing")}       
                }

            }
        }   
        # checking 'data$actor' ...
        if(!is.null(data$actor)){
            if(!is.list(data$actor)) stop("'data$actor' must be a list.")
            else{
                if(!is.null(data$actor$actorName) & !is.null(data$actor$variables)){
                    if(!is.matrix(data$actor$actorName) | !is.data.frame(data$actor$actorName)) stop("'data$actor$actorName' must be a column vector or a data.frame with one column containing the name of actors")
                    if(!is.array(data$actor$variables) | !is.data.frame(data$actor$variables)) stop("'data$actor$variables' must be either an array (with two [matrix] or more dimensions [3-d array]) or a data.frame in the simple case of constant exogenous dyadic variables")
                    if(!is.null(dimnames(data$actor$variables)[[2]])) exogenous_actor <- dimnames(data$actor$variables)[[2]]  # columns in a matrix/data.frame/3-d array
                    else{stop("column names of 'actor$variables' are missing")}                  
                }

            }
        }
    if(!all(unique(c(data$dyadic$dyad$actor1,data$dyadic$dyad$actor2,data$actor$actorName)) %in% attr(reh,"dictionary")$actors$actorName)) stop("actor names in 'data' must match with actor names in 'attr(reh,'dictionary')$actors$actorName' ")
    if(!all(unique(data$dyadic$dyad$type) %in% attr(reh,"dictionary")$types$typeName)) stop("event types in 'data' must match with event types in 'attr(reh,'dictionary')$types$typeName' ")
    exogenous_internal_names <- c(exogenous_dyadic,exogenous_actor)

    }
    effects_names <- c(effects_names,exogenous_internal_names)

    # add check on actor names and type names


    # getExogenousVariablesArray(data = data, reh = reh) Rcpp::List data .. Rcpp::List dyadic = data["dyadic"]; Rcpp::List actor = data["actor"]; Rcpp::List dictionary = reh.attr("dictionary");

    ##############################################
    ####  END EXOGENOUS VARIABLES PROCESSING  ####
    ##############################################

    #################################
    #### START FORMULA PROCESSING ###
    #################################
    if(!is.null(formula)){ # this if() is used just temporarily before integrating the use of formulas in the package
        if(class(formula)!="formula") stop("Argument 'formula' must be an object of class 'formula'")
        terms_formula <- terms(formula)
        if(attr(terms_formula,"response")==1) stop("The left element (response) of 'formula' must be left undefiend, e.g. '~ inertia + reciprocity + ...' ")
        all_variables <- all.vars(formula)
        status <- as.double(all(sapply(1:length(all_variables), function(x){any(effects_names==all_variables[x])})))  
        if(status == 0){stop("At least one of either the endogenous or the exogenous statistics is not typed according to the documentation (endogenous variables) or the 'data' input (exogenous variables)")}

        # get endogenous names from formula
        # substitute $endo_effects with the line below
        env$statisticsREH$endogenous_statistics <- endogenous_internal_names[na.omit(match(all_variables,endogenous_internal_names))] # endogenous_internal_names = effects_names

        # get exogenous names from formula
        # substitute $exogenous_stats with the line below
        env$statisticsREH$exogenous_statistics <- exogenous_internal_names[na.omit(match(all_variables,exogenous_internal_names))]

        # get intercept
        intercept <- NULL
        if(attr(terms_formula,"intercept")) intercept <- array(1,dim=c(reh$M,reh$D,1))
        # processing conditioning intercept to some factor variable or event type HERE ...

        # define interactions and/or transformation of variables
        model_frame <- attr(terms_formula,"factors") # regressors are by column, variables implied are by row
        # _//.. yet to be coded ..\\_
    }
    ##################################
    ####  END FORMULA PROCESSING  ####
    ##################################

    # processFormula() # handling interactions between regressors and mathematical transformations of variables (MOVE THIS FUNCTION TO getStepwiseModelsREH)

    env$statisticsREH$endo_effects <- effects$endogenous 
    env$statisticsREH$P <- length(effects$endogenous) # number of endogenous variables

    if(!is.null(effects$exogenous)){
        if(!is.matrix(effects$exogenous) & !is.data.frame(effects$exogenous)) {stop("The object 'exogenous' inside 'effects' must be either a matrix or a data.frame with the exogenous statistics of interest")}
            else{
                if(is.null(colnames(effects$exogenous))){colnames(effects$exogenous) <- sapply(1:dim(effects$exogenous)[2],function(x) paste("V",x,sep=""))}
                
                rearranged_exo_stats <- function(reh,exo_stats){  
                    position_rearranged <- NULL
                        for(i in 1:reh$D){
                            sender_old <- exo_stats$actor1[i] #reh$risksetMatrix[i,1]
                            receiver_old <- exo_stats$actor2[i] #reh$risksetMatrix[i,2]
                            type_old <- exo_stats$type[i] #reh$risksetMatrix[i,3]
                            
                            dict_loc <- attr(reh,"dictionary")
                            sender_new <- as.numeric(dict_loc$actors$actorName[which(dict_loc$actors$actorID == sender_old)])+1
                            receiver_new <- as.numeric(dict_loc$actors$actorName[which(dict_loc$actors$actorID== receiver_old)])+1
                            type_new <- as.numeric(dict_loc$types$typeName[which(dict_loc$types$typeID == type_old)])+1

                            position_new <- remify::getDyadIndex()+1 #reh$risksetCube[sender_new,receiver_new,type_new]+1
                            position_rearranged <- c(position_rearranged,position_new)
                        }
                    exo_stats_r <- cbind(exo_stats[position_rearranged,-c(1:3)])

                    array_exo <- array(NA,dim=c(reh$M,reh$D,dim(exo_stats_r)[2])) # empty array for exogenous statistics
                    for(s in 1:dim(exo_stats_r)[2]){
                        array_exo[,,s] <- matrix(rep(exo_stats_r[,s],reh$M),reh$M,reh$D,byrow=TRUE)
                    }
                    return(array_exo)
                }

                exogenous_rearranged <- rearranged_exo_stats(reh = reh, exo_stats = cbind(effects$exogenous))
                env$statisticsREH$exogenous_stats <- exogenous_rearranged
                env$statisticsREH$S <- dim(exogenous_rearranged)[3] #+ 1 # this +1 is temporarily herebecause of the PShift variable
            }
    }
    else{env$statisticsREH$exogenous_stats <- array(1,dim=c(reh$M,reh$D,1))
        dimnames(env$statisticsREH$exogenous_stats)[[3]] <- c("Intercept")
         env$statisticsREH$S <- 1}


    # Creating NULL objects for other routines
    env$statisticsREH$intervals <- NULL
}