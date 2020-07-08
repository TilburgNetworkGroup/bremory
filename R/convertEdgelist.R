#' convertEdgelist (edgelist) 
#'
#' A function which return the dimension of a specific environment. [remove export]
#' @param edgelist matrix of dimensions M*3  (by row[])
#' @param riskset matrix of dimensions N*(N-1) * (sender,receiver)
#' @param actors dataframe of actors names and corresponding integer id's.
#' @param rem boolean to understand whether the format of the edgelist has to be the 'rem' one (dyad,time) or the 'default' one (time,sender,receiver)
#'
#' @return  converted edgelist according to 'rem' argument
#' @export 
convertEdgelist <- function(edgelist, riskset, actors, rem = FALSE)
{
    if(rem){
        edgelist_out <- data.frame(dyad = rep(NA,dim(edgelist)[1]), time = edgelist[,1])
        for(i in 1:dim(edgelist)[1]){
            edgelist_out[i,1] <- which((riskset[,1]==edgelist[i,2]) & (riskset[,2]==edgelist[i,3])) 
        }
    }
    else{
        integer_id <- actors$integer_id
        actors <- actors$actors
        edgelist_out <- data.frame(time = edgelist$time, sender = rep(NA,dim(edgelist)[1]),  receiver = rep(NA,dim(edgelist)[1]))
        for(i in 1:dim(edgelist)[1]){
            edgelist_out$sender[i] <- integer_id[which(actors == edgelist$sender[i])]
            edgelist_out$receiver[i] <- integer_id[which(actors == edgelist$receiver[i])]
        }   
    }
    return(edgelist_out)
}