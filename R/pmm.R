#' Parametric Memory Model
#'
#' A function for the estimation of one or more step-wise models (at the moment working without remstats package)
#'
#' @param formula formula specifing the linear predictor (list of formulas if actor-oriented model)
#' @param reh either a 'reh' object (see remify::remify() function) or the edgelist to process
#' @param memory memory shape to use, if parameter is specified then no optimization is carried out
#' @param data [\emph{optional}] list of exogenous variables (time varying or not)
#' @param ncores [\emph{optional}] number of ncores to use for the parallelization
#' @param ... [\emph{optional}] other arguments referring to those in bremory::elpd()
#'
#' @return  object of class 'pmm' with attributes and methods explained in the vignette(topic = "pmm", package = "bremory")
#'
#' @export
#'
#' @examples 
#' 
#' # examples here
pmm <- function(formula,
                reh,
                memory,
                data = NULL,
                ncores = 1L,
                ...){
                    # ... parametric memory model here ...
    
    # ... ncores
    if(is.null(ncores)) ncores <- 1L
    else{
        if((parallel::detectCores() == 2L) & (ncores > 1L))
            stop("'ncores' is recommended to be set at most to 1.")
        else if((parallel::detectCores() > 2L) & (ncores > floor(parallel::detectCores()-2L)))
            stop("'ncores' is recommended to be set at most to: floor(parallel::detectCores()-2L)")
    }

                }