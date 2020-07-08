#' initializeREH (preprocessing raw data in a suitable format to the 'bremory' package, I used this as starting function for the reh() in remstimate)
#'
#' @param edgelist a matrix of dimensions (number of events * 3). By column (time,sender,receiver) of relational events. The first column is the time: if the timing is "ordinal" then a vector of increasing integer (from 1 to the number of events) is provided. Sender and receivers can be provided as either char or integer type.
#' @param riskset a matrix of dimensions [number of dyads in the risk set * 2]. By column (sender,receiver) of the dyad in the risk set. Sender and receivers can be provided as either char or integer type. The risk set must contain only possible interactions between actors, possbile event types must be excluded in the definition of the risk set. 
#' @param directed it is set TRUE by default. If FALSE, statistics where directionality of dyads matters will be summed up: e.g., inertia(i,j)+inertia(j,i).
#' @param eventType vector of event types.
#' @param env environment where the user wants to store the preprocessed output (it is set to .GlobalEnv by default). 
#'
#' @return  The environment in which the user prefers to work is updated by including vectors, matrices and other data structures necessary to pursue the analysis.
#'
#' @export
initializeREH <- function(edgelist,
                          riskset = NULL,
                          directed = TRUE,
                          env = globalenv()){
                                   
        if(missing(edgelist) | !is.data.frame(edgelist)){stop("User must provide an edgelist as a dataframe object. See documentation for further information about the edgelist structure")}
        env$initializeREH <- new.env()            
        env$initializeREH$actors <- if(is.null(riskset)){unique(c(edgelist$sender,edgelist$receiver))}
                      else{unique(c(as.matrix(riskset)[,1],as.matrix(riskset)[,2]))}
        env$initializeREH$N <- length(env$initializeREH$actors)
        env$initializeREH$M <- dim(edgelist)[1]

 
        env$initializeREH$actors <- data.frame(actors = env$initializeREH$actors, integer_id = c(0:(env$initializeREH$N-1)))
        env$initializeREH$riskset <- if(is.null(riskset)){matrix(na.omit(getRiskset(actors_id = env$initializeREH$actors$integer_id, N = env$initializeREH$N, selfedges = FALSE)),ncol=2)}
                       else{riskset}
        env$initializeREH$n_dyads <- dim(env$initializeREH$riskset)[1]

        env$initializeREH$dyad_position_array <- numeric(env$initializeREH$n_dyads)  
        env$initializeREH$old_riskset <- t(apply(env$initializeREH$riskset,1, function(x) c(env$initializeREH$actors$actors[which(env$initializeREH$actors$integer_id==x[1])],env$initializeREH$actors$actors[which(env$initializeREH$actors$integer_id==x[2])])))
        
        for(i in 1:env$initializeREH$n_dyads) env$initializeREH$dyad_position_array[i] <- which((env$initializeREH$riskset[,1]==env$initializeREH$old_riskset[i,1]) & (env$initializeREH$riskset[,2]==env$initializeREH$old_riskset[i,2]))

        env$initializeREH$riskset_matrix <- getRisksetMatrix(riskset = env$initializeREH$riskset, N = env$initializeREH$N)               
        
        edgelist <- convertEdgelist(edgelist = edgelist, riskset = env$initializeREH$riskset,  actors = env$initializeREH$actors, rem = FALSE)               
        env$initializeREH$edgelist <- data.matrix(edgelist) # we save the edgelist as a matrix
        env$initializeREH$directed <- directed               
        env$initializeREH$t <- env$initializeREH$edgelist[,1] #time vector

        cat('\n Environment "initializeREH" has been created successfully !\n')
}
