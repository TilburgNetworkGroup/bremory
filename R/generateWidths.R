#' generateWidths
#'
#' A function which return an array with all the statistics, both exogenous and endogenous.
#'
#' @param env environment where the user is currently working
#' @param K_range intervals dimensions
#' @param nsim_per_K number of simulations per K
#' @param max_width maximum time length
#'
#' @return  the function updates the object width in the getStepwiseModels environment.
#' @export
generateWidths <- function(env =  globalenv(), K_range = c(2:6), nsim_per_K = 2e02, max_width = NULL){

    Q <- length(K_range)*nsim_per_K
    min_diff <- max_width*0.02 #it was set equal to median(diff(c(0,env$initializeREH$t)))
    matrix_out <- matrix(NA, nrow = Q, ncol = (max(K_range)+1))

    for(i in 1:length(K_range))
    {
        dim_loc <- K_range[i]+1
        matrix_loc <- matrix(NA, nrow = nsim_per_K, ncol = dim_loc) 
        j <- 1
        while(j<=nsim_per_K)
        {
            diff_loc <- 0
            while(any(diff_loc<min_diff))
            {
                widths_loc <- sort(c(0,runif((K_range[i]-1),0,max_width),max_width))
                diff_loc <- diff(widths_loc)        
            }
            matrix_loc[j,] <- widths_loc #c(widths_loc,Inf)
            j <- j+1
        }
        matrix_out[c((i+(nsim_per_K-1)*(i-1)):(i*nsim_per_K)),c(1:dim_loc)] <- matrix_loc
        #print(i)
    }
    env$stepwiseModelsREH$widths <- matrix_out
    
}

