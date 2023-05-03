# 'process_data' is a function for processing exogenous statistics in light of the formula provided as input
# (1) it processes the input data according to the reh object of the event sequence
# (2) it creates a list object inside the statisticsREH environment with all the exogenous statistics needed to define the linear predictor (only for the marginal effects)
# (3) it frees space when exogenous statistics are not defined in the linear predictor
# (4-TO UPGRADE) in the current version, the function can only process actor variables that are not sender or receiver specific. They will then be assigned to the dyads where the specific actor is either sender or receiver of the interaction
process_data <- function(data, parsed_formula, reh){
    env <- globalenv()

    if(!is.null(data)){
        if(!is.list(data)) stop("'data' must be a list of two objects: 'dyadic' for (dyad-level exogenous statistics), 'actor' for actor-level exogenous statistics.")
        # ... Checking 'data$dyadic' 
        if(!is.null(data$dyadic)){
            if(!is.list(data$dyadic)) stop("'data$dyadic' must be a list")
            else{
                if(!is.null(data$dyadic$dyad) & !is.null(data$dyadic$variables)){
                    if(!is.data.frame(data$dyadic$dyad)){
                        stop("'data$dyadic$dyad' must be a 'data.frame' with columns indicating [actor1,actor2]")
                    }
                    else{   
                        match_names <- match(c("actor1","actor2"),dimnames(data$dyadic$dyad)[[2]])
                        if(anyNA(match_names)){
                            stop("column names of 'data$dyadic$dyad' must be the following : 'actor1', 'actor2'")
                        }
                        else{
                            # ... Checking actors' names
                            if(!all(unique(c(data$dyadic$dyad$actor1,data$dyadic$dyad$actor2)) %in% attr(reh,"dictionary")$actors$actorName)) stop("actor names in 'data$dyadic$dyad' must match with actor names in 'attr(reh,'dictionary')$actors$actorName' ")
                        }
                        
                        # ... Checking list 'data$dyadic$variables'
                        if(!is.list(data$dyadic$variables)){
                            stop("'data$dyad$variables' must be a list of matrices (for time varying dyadic variables) or one-column matrices (for non-time-varying dyadic variables). ")
                        }
                        else if(length(data$dyadic$variables) > 0){
                            if(is.null(names(data$dyadic$variables))){
                                stop("The exogenous dyadic variables supplied must be named.")
                            }
                            dataframe_check_dyadic_variables <- unlist(lapply(data$dyadic$variables, function(x) !is.data.frame(x)))
                            if(any(dataframe_check_dyadic_variables)){
                                stop(paste("Some objects of the list 'data$dyadic$variables' must be either a data.frame with reh$M columns (for time-varying dyadic variables) or a one-column data.frame (for constant dyadic variables). Check variables:\n\n,",names(data$dyadic$variables)[dataframe_check_dyadic_variables],"\n",paste=" "))
                            }
                            dim_dyadic_variables <- unlist(lapply(data$dyadic$variables, function(x) dim(x)[1]))
                            if(any(dim_dyadic_variables != dim(data$dyadic$dyad)[1])){
                                stop("one or more variables in 'data$dyadic$variables' have a different number of dyads than those provided in the object 'data$dyadic$dyad'.")
                            }         
                        }
                    }        
                }
                else{
                    stop("'data$dyadic' must be a list of two objects: 'dyad' and 'variables'")
                } 
                ## (??) check for NA's and throw warnings/errors (this shouldn't throw errors when dyads are not in the riskset and the statistic value is set to NA)       
            }    
        }   
        # ... Checking 'data$actor'
        if(!is.null(data$actor)){
            if(!is.list(data$actor)) stop("'data$actor' must be a list.")
            else{
                if(!is.null(data$actor$actorName) & !is.null(data$actor$variables)){
                    if(!is.data.frame(data$actor$actorName)){
                        stop("'data$actor$actorName' must be a 'data.frame' with one column containing the name of actors")
                    }
                    # ... Checking actor names
                    else if(!all(unique(c(data$actor$actorName[,1])) %in% attr(reh,"dictionary")$actors$actorName)){
                        stop("actor names in 'data$actor$actorName' must match with actor names in 'attr(reh,'dictionary')$actors$actorName' ")
                    }
                    # ... Checking list 'data$actor$variables'
                    if(!is.list(data$actor$variables)){
                        stop("'data$actor$variables' must be a list of matrices (for time varying actor variables) or vectors (for non-time-varying actor variables). ")
                    }
                    else if(length(data$actor$variables) > 0){
                        if(is.null(names(data$actor$variables))){
                            stop("The exogenous actor variables supplied must be named")
                        }  
                        dataframe_check_actor_variables <- unlist(lapply(data$actor$variables, function(x) !is.data.frame(x)))
                        if(any(dataframe_check_actor_variables)){
                            stop(paste("Some objects of the list 'data$actor$variables' must be either a data.frame with reh$M columns (for time-varying actor variables) or a one-column data.frame (for constant actor variables). Check variables:\n\n,",names(data$actor$variables)[dataframe_check_actor_variables],"\n",paste=" "))
                            }
                        dim_actor_variables <- unlist(lapply(data$actor$variables, function(x) dim(x)[1]))
                        if(any(dim_actor_variables != dim(data$actor$actorName)[1])){
                            stop("one or more variables in 'data$actor$variables' have a different number of dyads than those provided in the object 'data$actor$actorName'.")
                        } 
                        if(any(unlist(lapply(data$actor$variables,function(x) is.null(attr(x,"which")))))){
                            stop("one or more variables in 'data$actor$variables' have no attribute 'which' specified. Please, provide the attribute for each actor variable. The attribute can assume value 'sender' or 'receiver'")
                        }  
                        else if(!all(unlist(lapply(data$actor$variables,function(x) attr(x,"which"))) %in% c("sender","receiver"))){
                            stop("attribute 'which' for the actor variables must be set either to 'sender' or to 'receiver'")
                        }                                             
                    }             
                }
                else{
                    stop("'data$actor' must be a list of two objects: 'actorName' and 'variables'")
                }

                ## (??) check for NA's and throw warnings/errors (this shouldn't throw errors when actors are not in the riskset and the statistic value is set to NA)    
            }
        }
        # ... Saving names of variables from 'data'
        exogenous_stats_dyadic <- names(data$dyadic$variables)
        exogenous_stats_actor <- names(data$actor$variables)
        if(anyNA(match(attr(parsed_formula,"exos"),c(exogenous_stats_actor,exogenous_stats_dyadic)))){
            stop("one or more exogenous variables specified in 'formula' are not provided in the object 'data'")
        }

        # ... Computing vectors of rearranged id's for ...
        # ... dyadic variables
        rearranged_ids_dyadic <- rep(NA,reh$D) # for dyadic variables, is a vector of new dyad positions 
        dict_loc <- attr(reh,"dictionary")
        if(!is.null(data$dyadic$dyad)){ # [note: this function can be faster if written in Rcpp]
            for(d in 1:dim(data$dyadic$dyad)[1]){
                actor1_id <- dict_loc$actors$actorID[which(dict_loc$actors$actorName == data$dyadic$dyad$actor1[d])]-1
                actor2_id <- dict_loc$actors$actorID[which(dict_loc$actors$actorName == data$dyadic$dyad$actor2[d])]-1
                # ... Finding dyad index and saving it
                rearranged_ids_dyadic[d] <- remify:::getDyadIndex(actor1 = actor1_id, actor2 = actor2_id, type = 0, N = reh$N, directed = attr(reh,"directed")) + 1 # we sum one because we are working at R level and vector positions start at 1
            }
            rearranged_ids_dyadic <- order(rearranged_ids_dyadic)
        }

        # ... actor variables
        rearranged_ids_actor <- list(sender=0,receiver=0)
        if(!is.null(data$actor$actorName)){      
            if(any(unlist(lapply(data$actor$variables,function(x) attr(x,"which")=="sender")))){
                rearranged_ids_actor[["sender"]] <- matrix(NA, nrow = reh$N,ncol = (reh$N-1))
                for(n in 1:dim(data$actor$actorName)[1]){
                    actor_id <- dict_loc$actors$actorID[which(dict_loc$actors$actorName == data$actor$actorName[n,1])]-1
                    counter_idx <- 1
                    for(r in 1:reh$N){
                        if((r-1)!=actor_id){ # avoiding self-loops
                            # dyad when actor_id is the sender
                            rearranged_ids_actor[["sender"]][n,counter_idx] <- remify:::getDyadIndex(actor1 = actor_id, actor2 = (r-1), type = 0, N = reh$N, directed = attr(reh,"directed")) + 1
                            counter_idx <- counter_idx + 1
                        }
                    }
                    rm(counter_idx)
                }
                rearranged_ids_actor[["sender"]] <- order(as.vector(t(rearranged_ids_actor[["sender"]])))  
            }
            if(any(unlist(lapply(data$actor$variables,function(x) attr(x,"which")=="receiver")))){
                rearranged_ids_actor[["receiver"]] <- matrix(NA, nrow = reh$N,ncol = (reh$N-1))
                for(n in 1:dim(data$actor$actorName)[1]){
                    actor_id <- dict_loc$actors$actorID[which(dict_loc$actors$actorName == data$actor$actorName[n,1])]-1
                    counter_idx <- 1
                    for(s in 1:reh$N){
                        if((s-1)!=actor_id){ # avoiding self-loops
                            # dyad when actor_id is the receiver
                            rearranged_ids_actor[["receiver"]][n,counter_idx] <- remify:::getDyadIndex(actor1 = (s-1), actor2 = actor_id, type = 0, N = reh$N, directed = attr(reh,"directed")) + 1
                            counter_idx <- counter_idx + 1
                        }
                    }
                    rm(counter_idx)
                }
                rearranged_ids_actor[["receiver"]] <- order(as.vector(t(rearranged_ids_actor[["receiver"]])))
            }                         
        }
        
        # ... Processing exogenous variables in the formula, attr(parse_formula,"exos")
        exos_names <- attr(parsed_formula, "exos")
        env$stats <- new.env() # empty environment where to save statistics
        # an alternative would be to save each processed endogenous statistic in a separate .RData [next UPGRADE]
        if(attr(parsed_formula,"intercept")){
            env$stats[["intercept"]] <- array(1,dim=c(reh$M,reh$D,1))
            env$statisticsREH$stats_names[["intercept"]] <- "intercept"
        }
        for(s in 1:length(exos_names)){
            find_var_actor <- which(exogenous_stats_actor == exos_names[s])
            find_var_dyad <- which(exogenous_stats_dyadic == exos_names[s])
            # ... Procesing dyadic variable
            if(length(find_var_dyad)>0){ 
                if(is.logical(data$dyadic$variables[[exos_names[s]]][,1])){ # logical variable (TRUE/FALSE coded as 1/0)
                    var_s <- NULL
                    # time-varying variable with M time points as columns
                    var_s <- t(data$dyadic$variables[[exos_names[s]]][rearranged_ids_dyadic,]) # dimensions [M x D]
                    if(dim(data$dyadic$variables[[exos_names[s]]])[2]==1){
                        var_s <- array(rep(as.vector(var_s),each=reh$M),dim=c(reh$M,reh$D,1))# creating a matrix of dimensions [M x D]
                    }
                    env$stats[[exos_names[s]]] <- (var_s*1)
                    env$statisticsREH$stats_names[[exos_names[s]]] <- exos_names[s]
                    rm(var_s)
                }
                else if(is.factor(data$dyadic$variables[[exos_names[s]]][,1]) | is.character(data$dyadic$variables[[exos_names[s]]][,1])){ # factor variable
                    factor_s <- factor(unique(as.vector(as.matrix(data$dyadic$variables[[exos_names[s]]])))) # because there can be levels observed at some time points and not in others 
                    ctr_s <- stats::contrasts(factor_s)
                    names_s <- paste(exos_names[s],".",colnames(ctr_s),sep="")
                    var_s <- data.frame(data$dyadic$variables[[exos_names[s]]][rearranged_ids_dyadic,])
                    if(dim(data$dyadic$variables[[exos_names[s]]])[2]>1){ # time-varying: this step might require lots of memory
                        env$stats[[exos_names[s]]] <- array(0,dim=c(reh$M,reh$D,(length(levels(factor_s))-1))) # dimensions [M x D x (levels-1)]
                        dimnames(env$stats[[exos_names[s]]])[[3]] <- names_s 
                        for(i in 1:(reh$D)){
                            for(j in 1:reh$M){
                                k <- which(rownames(ctr_s)==var_s[i,j])
                                env$stats[[exos_names[s]]][j,i,] <- as.numeric(ctr_s[k,])
                            }
                        }
                    }
                    else{
                        env$stats[[exos_names[s]]] <- array(0,dim=c(reh$M,reh$D,(length(levels(factor_s))-1))) # dimensions [M x D x (levels-1)]
                        dimnames(env$stats[[exos_names[s]]])[[3]] <- names_s
                        for(i in 1:dim(var_s)[1]){
                            k <- which(rownames(ctr_s)==var_s[i,1])
                            env$stats[[exos_names[s]]][,i,] <- matrix(rep(as.numeric(ctr_s[k,]),reh$M),nrow=reh$M,ncol=length(levels(factor_s))-1,byrow=TRUE)
                        }
                    }
                    env$statisticsREH$stats_names[[exos_names[s]]] <- names_s
                    rm(factor_s,var_s,ctr_s,names_s)
                }
                else if(is.numeric(data$dyadic$variables[[exos_names[s]]][,1])){ # numeric variable
                    var_s <- t(data$dyadic$variables[[exos_names[s]]][rearranged_ids_dyadic,]) # [M x D] or a [1 x D]
                    if(dim(data$dyadic$variables[[exos_names[s]]])[2]==1){
                        env$stats[[exos_names[s]]] <- array(rep(as.vector(var_s),each=reh$M),dim=c(reh$M,reh$D,1))
                    }
                    else{
                        env$stats[[exos_names[s]]] <-  array(var_s,dim=c(reh$M,reh$D,1)) # dimensions [M x D]
                    }
                    rm(var_s) 
                    env$statisticsREH$stats_names[[exos_names[s]]] <- exos_names[s]             
                }
                data$dyadic$variables[[exos_names[s]]] <- NULL
            }
            # ... Processing actor variable
            else if(length(find_var_actor)>0){ 
                if(is.logical(data$actor$variables[[exos_names[s]]][,1])){ # logical variable (TRUE/FALSE as 1/0)
                    var_s <- NULL
                    # time-varying variable with M time points as columns
                    var_s <- cbind(data$actor$variables[[exos_names[s]]][rep(1:reh$N,each = (reh$N-1)*1),]) # expanding the variable from reh$N to reh$D , here 1 = reh$C (one event type)
                    var_s <- t(var_s[rearranged_ids_actor[[attr(data$actor$variables[[exos_names[s]]],"which")]],]) # dimensions [M x D]
                    if(dim(data$actor$variables[[exos_names[s]]])[2]>1){
                        env$stats[[exos_names[s]]] <- array(var_s*1,dim=c(reh$M,reh$D,1))
                    }
                    else{
                        var_s <- array(rep(as.vector(var_s),each=reh$M),dim=c(reh$M,reh$D,1)) # creating a matrix of dimensions [M x D]
                        env$stats[[exos_names[s]]] <- (var_s*1)
                    }
                    env$statisticsREH$stats_names[[exos_names[s]]] <- exos_names[s]
                    rm(var_s)
                }
                else if(is.factor(data$actor$variables[[exos_names[s]]][,1]) | is.character(data$actor$variables[[exos_names[s]]][,1])){ # factor variable
                    factor_s <- factor(unique(as.vector(as.matrix(data$actor$variables[[exos_names[s]]])))) # because there can be levels observed at some time points and not in others 
                    ctr_s <- stats::contrasts(factor_s)
                    names_s <- paste(exos_names[s],".",colnames(ctr_s),sep="")
                    var_s <- data.frame(data$actor$variables[[exos_names[s]]][rep(1:reh$N,each = (reh$N-1)*1),])# expanding the variable from reh$N to reh$D, here 1 = reh$C (one event type)
                    var_s <- t(var_s[rearranged_ids_actor[[attr(data$actor$variables[[exos_names[s]]],"which")]],]) # dimensions [M x D]
                    if(dim(data$actor$variables[[exos_names[s]]])[2]>1){ # time-varying: this step might require lots of memory
                        env$stats[[exos_names[s]]] <- array(0,dim=c(reh$M,r/reh$C,(length(levels(factor_s))-1))) # dimensions [M x D x (levels-1)]
                        dimnames(env$stats[[exos_names[s]]])[[3]] <- names_s 
                        for(i in 1:(reh$D)){
                            for(j in 1:reh$M){
                                k <- which(rownames(ctr_s)==var_s[i,j])
                                env$stats[[exos_names[s]]][j,i,] <- as.numeric(ctr_s[k,])
                            }
                        }
                    }
                    else{
                        env$stats[[exos_names[s]]] <- array(0,dim=c(reh$M,reh$D,(length(levels(factor_s))-1))) # dimensions [M x D x (levels-1)]
                        dimnames(env$stats[[exos_names[s]]])[[3]] <- names_s
                        for(i in 1:dim(var_s)[2]){ # dim(var_s)[2] is the same as (reh$D)
                            k <- which(rownames(ctr_s)==var_s[1,i])
                            env$stats[[exos_names[s]]][,i,] <- matrix(rep(as.numeric(ctr_s[k,]),reh$M),nrow=reh$M,ncol=length(levels(factor_s))-1,byrow=TRUE)
                        }
                    }
                    env$statisticsREH$stats_names[[exos_names[s]]] <- names_s
                    rm(factor_s,var_s,ctr_s,names_s)
                }
                else if(is.numeric(data$actor$variables[[exos_names[s]]][,1])){ # numeric variable
                    var_s <- cbind(data$actor$variables[[exos_names[s]]][rep(1:reh$N,each = (reh$N-1)*1),]) # expanding the variable from reh$N to reh$D, here 1 = reh$C (one event type)
                    var_s <- t(var_s[rearranged_ids_actor[[attr(data$actor$variables[[exos_names[s]]],"which")]],]) # dimensions [M x D]
                    if(dim(data$actor$variables[[exos_names[s]]])[2]>1){ # time-varying: this step might require lots of memory
                        env$stats[[exos_names[s]]] <-  array(var_s,dim=c(reh$M,reh$D,1))# dimensions [M x D]
                    }
                    else{
                        env$stats[[exos_names[s]]] <- array(rep(as.vector(var_s),each = reh$M),dim=c(reh$M,reh$D,1))
                    }   
                    env$statisticsREH$stats_names[[exos_names[s]]] <- exos_names[s]
                    rm(var_s)
                }
                # ... Free-ing space (for now that data is provided as input object. When data will be on Hardware, we will minimize the memory usage)
                data$actor$variables[[exos_names[s]]] <- NULL
            }
            else{
                #stop("user should not reach this condition (!!)")
            }
        }
    }
    # ...  if 'data' is not supplied (NULL), do nothing
}





