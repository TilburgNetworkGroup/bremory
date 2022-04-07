#' intervals
#'
#' A function which return an array with all the statistics, both exogenous and endogenous.
#'
#' @param K intervals dimensions
#' @param nsimK number of simulations per K
#' @param maxWidth maximum time length
#' @param minDiff minimum width (diff will be changed to width)
#' @param intervals type of intervals to be generated: "all" it will return an object of class "intervals" containing nsimK intervals increasing, nsimK intervals decreasing and 1 interval equal per each K; if "increasing" or "decreasing" only nsimK increasing intervals per each K will be generated; if "equal" only one equal size interval will be returned per each K.
#'
#' @return  the function updates the object width in the getStepwiseModels environment.
#' @export
intervals <- function(K = NULL, nsimK = 1e02, maxWidth = NULL, minDiff = NULL, intervals = c("all", "increasing", "decreasing","equal")){

    K_input <- K
    if(intervals == "all"){

        nsim <- nsimK/2  
        widths_matrix <- matrix(NA, nrow = (1+nsimK)*length(K), ncol = (max(K)+1) ) 
        widths_type <- NULL
        iter_K <- 1
        for(K in min(K):max(K)){
            out_K <- NULL
             # for the storing of widths
            matrix_loc <- matrix(NA, nrow = (2*nsim+1), ncol = (K+1)) 

            # (0) equal widths
            matrix_loc[1,] <- seq(0,maxWidth,length=(K+1))
            
            # (1) increasing widths
            out_K <- matrix(unlist(lapply(1:nsim,function(s){
            diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
                                while(min(diffs)<minDiff){
                                    diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
                                }
                                c(0,cumsum(diffs))
                                })),byrow=TRUE,ncol=K+1)
            out_K <- out_K * maxWidth
            matrix_loc[2:(nsim+1),] <- out_K
 

            # ...(2) decreasing widths (maybe here we should generate again increasing and reverse them, instead of using the same increasing)
            span_interval <- t(apply(matrix_loc[c(2:(nsim+1)),],1,diff))
            span_interval <- span_interval[,dim(span_interval)[2]:1]
            widths_decr <- t(apply(span_interval,1,function(x) cumsum(c(0,x))))
            matrix_loc[c((nsim+2):((nsim*2)+1)),] <- widths_decr
            
            # ...(3) store widths
            widths_matrix[c((1+(nsimK+1)*(iter_K-1)):(iter_K*(nsimK+1))),c(1:(K+1))] <- matrix_loc
            widths_type <- c(widths_type,c("equal",rep("increasing",nsim),rep("decreasing",nsim))) #,rep("equal",nsim)
            iter_K <- iter_K + 1

        }
        out_str <- structure(list(
                            widths=widths_matrix,
                            K=rep(K_input,each=(nsimK+1)),
                            type = widths_type), class="intervals")
        return(out_str)
    }

    if(intervals == "increasing"){
        diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
        while(min(diffs)<minDiff){
                diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
        }
        return( c(0,cumsum(diffs)) * maxWidth )
    }


    if(intervals == "decreasing"){
        diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
        while(min(diffs)<minDiff){
                diffs <- sort(c(extraDistr::rdirichlet(1,alpha=rep(1,K))))
        }
        diffs <- sort(diffs, decreasing = TRUE)
        return( c(0,cumsum(diffs)) * maxWidth )
    }


    if(intervals == "equal"){
        return(seq(0,maxWidth,length=(K+1)))
    }

    ##################
    # OLD CODE START #
    ##################


    #Q <- length(K_range)*nsimK
    #minDiff <- 0.1 #maxWidth*0.02 #it was set equal to median(diff(c(0,env$initializeREH$t)))
    #matrix_out <- matrix(NA, nrow = Q, ncol = (max(K_range)+1))

    #if(intervals == "default"){
    #for(i in 1:length(K_range))
    #{
    #    dim_loc <- K_range[i]+1
    #    matrix_loc <- matrix(NA, nrow = nsimK, ncol = dim_loc) 
    #    j <- 1
    #    while(j<=nsimK)
    #    {
    #        diff_loc <- 0
    #        while(any(diff_loc<minDiff))
    #        {
    #            widths_loc <- sort(c(0,runif((K_range[i]-1),0,maxWidth),maxWidth))
    #            diff_loc <- diff(widths_loc)        
    #        }
    #        matrix_loc[j,] <- widths_loc #c(widths_loc,Inf)
    #        j <- j+1
    #    }
    #    matrix_out[c((i+(nsimK-1)*(i-1)):(i*nsimK)),c(1:dim_loc)] <- matrix_loc
    #    #print(i)
    #}
    #env$stepwiseModelsREH$widths <- matrix_out
    #}


    # increasing
    #if(intervals == "increasing"){
    #    for(i in 1:length(K_range))
    #    {
    #        print(i+1)
    #        dim_loc <- K_range[i]+1
    #        matrix_loc <- matrix(NA, nrow = nsimK, ncol = dim_loc) 
    #        j <- 1
           
    #        while(j<=nsimK)
    #        {
    #            cond <- TRUE
    #            while(cond) 
    #            {
    #                widths_loc <- sort(c(0,runif((K_range[i]-1),minDiff,maxWidth-minDiff),maxWidth))
    #                diff_loc <- diff(widths_loc)   
    #                is_unsorted <- is.unsorted(diff_loc)
    #                cond <- (any(diff_loc<minDiff) | is_unsorted)

    #            }
    #            matrix_loc[j,] <- widths_loc #c(widths_loc,Inf)
    #            j <- j+1
    #        }
    #        matrix_out[c((i+(nsimK-1)*(i-1)):(i*nsimK)),c(1:dim_loc)] <- matrix_loc
    #        #print(i)
    #    }
    #    env$stepwiseModelsREH$widths <- matrix_out
    #}

    
   # if(intervals == "incr_eq_decr"){
   #     nsim <- (nsimK-1)/2 # old division with equal-ish intervals /3
   #     env$stepwiseModelsREH$widths <- matrix(NA, nrow = nsimK*length(K_range), ncol = (max(K_range)+1) )
   #     widths_type <- NULL
   #     iter_K <- 1
   #     for(K in min(K_range):max(K_range)){
   #     out_K <- NULL
   #     count <- 0
   #     matrix_loc <- matrix(NA, nrow = (2*nsim+1), ncol = (K+1)) # old row size with equal-ish (3*nsim+1) 
   #     matrix_loc[1,] <- seq(0,maxWidth,length=(K+1))
        
         # ...(1) increasing widths
   #     while(count!=nsim){
   #         widths_eq <- seq(0,maxWidth,length=(K+2))
   #         widths_eq[2] <- runif(1,minDiff,widths_eq[2]) #minDiff+0.1
   #         for(k in 3:(K+1)){  
   #             widths_eq[k] <- runif(1,widths_eq[k-1]+(widths_eq[k-1]-widths_eq[k-2]),widths_eq[k]) # maxWidth
   #         }
            
   #         widths_eq <- widths_eq[-(length(widths_eq)-1)]
            
   #         if(!any(is.na(widths_eq)) ){ # & (!is.unsorted(diff(widths_eq)))
   #             out_K <- rbind(out_K,widths_eq)
   #             count <- count+1
   #         }
   #     }
   #     matrix_loc[2:(nsim+1),] <- out_K

        # ...(2) decreasing widths
   #     span_interval <- t(apply(matrix_loc[c(2:(nsim+1)),],1,diff))
   #     span_interval <- span_interval[,dim(span_interval)[2]:1]
   #     widths_decr <- t(apply(span_interval,1,function(x) cumsum(c(0,x))))
   #     matrix_loc[c((nsim+2):((nsim*2)+1)),] <- widths_decr
        
        ########################
        ####### OLD CODE #######
        ########################
        # ...(3) small_large/large_small_widths (equal-ish widths)
        #deltas_K <- runif(nsim/2,0,minDiff)
        #seq_K <- seq(0,maxWidth,length=(K+1))
        #ncol_K <- K-1
        #ones_K <- 1:ncol_K
        #deltas_matrix_K <- matrix(deltas_K,nrow=nsim/2,ncol=ncol_K)
        
        # small - large
        #signed_vec_K <- sapply(1:length(ones_K),function(x) if(x%%2 ==0) 1 else{-1})
        #variations_small_large_K <- (deltas_matrix_K * matrix(signed_vec_K,nrow=nsim/2,ncol=ncol_K,byrow=TRUE)) + matrix(seq_K[-c(1,length(seq_K))],nrow=nsim/2,ncol=(K-1),byrow=TRUE)
        #small_large_K <- cbind(0,variations_small_large_K,maxWidth)
        
        # large - small
        #signed_vec_K <- sapply(1:length(ones_K),function(x) if(x%%2 ==0) -1 else{1})
        #variations_large_small_K <- (deltas_matrix_K * matrix(signed_vec_K,nrow=nsim/2,ncol=ncol_K,byrow=TRUE)) + matrix(seq_K[-c(1,length(seq_K))],nrow=nsim/2,ncol=(K-1),byrow=TRUE)
        #large_small_K <- cbind(0,variations_large_small_K,maxWidth)
        
        # save intervals
        #matrix_loc[c((2*nsim+2):((nsim*3)+1)),] <- rbind(small_large_K,large_small_K)
        ########################
        ####### OLD CODE #######
        ########################
        # ...(4) store widths
  #      env$stepwiseModelsREH$widths[c((1+nsimK*(iter_K-1)):(iter_K*nsimK)),c(1:(K+1))] <- matrix_loc
  #      iter_K <- iter_K + 1
  #      widths_type <- c(widths_type,c("equal",rep("increasing",nsim),rep("decreasing",nsim))) #,rep("equal",nsim)
  #      }
  #      env$stepwiseModelsREH$widths_type <- widths_type
  #  }
    
  #   if(intervals == "incr_eq_decr_new"){
  #     nsim <- (nsimK-1)/2 # old division with equal-ish intervals /3
  #      env$stepwiseModelsREH$widths <- matrix(NA, nrow = nsimK*length(K_range), ncol = (max(K_range)+1) )
  #      widths_type <- NULL
  #      iter_K <- 1
  #      for(K in min(K_range):max(K_range)){
  #      out_K <- NULL
  #      count <- 0
  #      matrix_loc <- matrix(NA, nrow = (2*nsim+1), ncol = (K+1)) # old row size with equal-ish (3*nsim+1) 
  #      matrix_loc[1,] <- seq(0,maxWidth,length=(K+1))

  #          # ...(1) increasing widths
  #      while(count!=nsim){
  #          widths_eq <- seq(0,maxWidth,length=(K+2))
  #          widths_eq[2] <- runif(1,minDiff,widths_eq[2]) #minDiff+0.1
  #          for(k in 3:(K+1)){  
  #              widths_eq[k] <- runif(1,widths_eq[k-1]+(widths_eq[k-1]-widths_eq[k-2]),maxWidth) 
  #          }
            
  #          widths_eq <- widths_eq[-(length(widths_eq)-1)]
            
  #          if(!any(is.na(widths_eq)) & (!is.unsorted(diff(widths_eq)))){ 
  #              out_K <- rbind(out_K,widths_eq)
  #              count <- count+1
  #          }
  #      }
  #      matrix_loc[2:(nsim+1),] <- out_K

  #      # ...(2) decreasing widths
  #      span_interval <- t(apply(matrix_loc[c(2:(nsim+1)),],1,diff))
  #      span_interval <- span_interval[,dim(span_interval)[2]:1]
  #      widths_decr <- t(apply(span_interval,1,function(x) cumsum(c(0,x))))
  #      matrix_loc[c((nsim+2):((nsim*2)+1)),] <- widths_decr
        
  #      # ...(3) store widths
  #      env$stepwiseModelsREH$widths[c((1+nsimK*(iter_K-1)):(iter_K*nsimK)),c(1:(K+1))] <- matrix_loc
  #      widths_type <- c(widths_type,c("equal",rep("increasing",nsim),rep("decreasing",nsim))) #,rep("equal",nsim)
  #      iter_K <- iter_K + 1
  #      }
  #      env$stepwiseModelsREH$widths_type <- widths_type
  #  }
    
    #if(intervals == "incr_eq_decr"){
    #    nsim <- (nsimK-1)/2
    #    env$stepwiseModelsREH$widths <- matrix(NA, nrow = nsimK*length(K_range), ncol = (max(K_range)+1) )

    #    for(k in 1:length(K_range)){
    #        print(K)
    #        K <- K_range[k]
    #        matrix_loc <- matrix(NA, nrow = (2*nsim+1), ncol = (K+1))
    #        matrix_loc[1,] <- seq(0,maxWidth,length=(K+1))

    #        # ...(1) increasing widths
    #        for(j in 2:(nsim+1)){
    #            cond <- TRUE
    #            while(cond){
    #                par <- runif(1)
    #                widths_loc <- c(cumsum(c(0,diff(sapply(0:(K-1),function(x) maxWidth*exp(par*x)/exp(par*maxWidth))))),maxWidth)
    #                diff_loc <- diff(widths_loc)
    #                cond <- any(diff_loc<minDiff)
    #                print("cond")
    #            }
    #            matrix_loc[j,] <- widths_loc
    #        }

    #        # ...(2) decreasing widths
    #        span_interval <- t(apply(matrix_loc[c(2:(nsim+1)),],1,diff))
    #        span_interval <- span_interval[,dim(span_interval)[2]:1]
    #        widths_decr <- t(apply(span_interval,1,function(x) cumsum(c(0,x))))
    #        matrix_loc[c((nsim+2):((nsim*2)+1)),] <- widths_decr

    #        # ...(3) store widths
    #        env$stepwiseModelsREH$widths[c((k+(nsimK-1)*(k-1)):(k*nsimK)),c(1:(K+1))] <- matrix_loc
    #        rm(matrix_loc)
    #    }       
    #}

    #if(intervals == "find_maxWidth"){
    #    gamma_1 <- seq(minDiff,maxWidth,length=nsimK)
    #    env$stepwiseModelsREH$widths <- cbind(0,gamma_1,gamma_1*2,gamma_1*3,gamma_1*4) #[1] cbind(0,gamma_1,max(env$initializeREH$edgelist[,1])) [2]cbind(0,minDiff,gamma_1)
    #}
    ##################
    # OLD CODE END   #
    ##################
    
}

