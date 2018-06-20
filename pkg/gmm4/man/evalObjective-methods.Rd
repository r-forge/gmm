\name{evalObjective-methods}
\docType{methods}
\alias{evalObjective}
\alias{evalObjective-methods}
\alias{evalObjective,gmmModels,numeric,gmmWeights-method}
\alias{evalObjective,sysGmmModels,list,sysGmmWeights-method}
\title{ ~~ Methods for Function \code{evalObjective} in Package \pkg{Gmm} ~~}
\description{
 ~~ Methods to compute the GMM objective function. ~~
}
\usage{
\S4method{evalObjective}{gmmModels,numeric,gmmWeights}(object, theta,
wObj, \dots)

\S4method{evalObjective}{sysGmmModels,list,sysGmmWeights}(object, theta,
wObj, \dots)
}
\arguments{
  \item{object}{An object of class \code{"gmmModels"} or \code{"sysGmmModels"}.}
  \item{theta}{The vector for coefficients for single equation, or a
    list of vector for system of equations.}
  \item{wObj}{An object of class \code{"gmmWeights"} or \code{"sysGmmWeights"}.}
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "gmmModels", theta = "numeric",
    wObj = "gmmWeights")}}{
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)
w <- evalWeights(model1, theta)
evalObjective(model1, theta, w)
}

\keyword{methods}
