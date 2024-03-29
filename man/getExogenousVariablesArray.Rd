% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getExogenousVariablesArray}
\alias{getExogenousVariablesArray}
\title{getExogenousVariablesArray}
\usage{
getExogenousVariablesArray(exos, data, reh)
}
\arguments{
\item{exos}{vector of exogenous variables in the linear predictor (we avoid to process statistics that are not specified in the linear predictor but provided in 'data')}

\item{data}{list of exogenous variables (dyadic and actor variables), with actor/type names converted to ID's according to reh object}

\item{reh}{reh object obtained from remify::reh(). It is one of the input arguments needed by bremory::smm()}
}
\value{
list of arrays of 3 dimensions: reh$M rows (number of events in the sequence), reh$D columns (number of dyads at risk in the network), reh$S slices (number of exogenous statistics). No interactions are accounted for at this stage.
}
\description{
getExogenousVariablesArray
}
