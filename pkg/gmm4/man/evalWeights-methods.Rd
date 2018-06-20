\name{evalWeights-methods}
\docType{methods}
\alias{evalWeights}
\alias{evalWeights-methods}
\alias{evalWeights,gmmModels-method}
\alias{evalWeights,sysGmmModels-method}
\alias{evalWeights,rslinearGmm-method}
\title{Methods for Function \code{evalWeights} in Package \pkg{Gmm}}
\description{
 This is a constructor for objects of class \code{gmmWeights} 
}
\usage{
\S4method{evalWeights}{gmmModels}(object, theta=NULL, w="optimal",
\dots)

\S4method{evalWeights}{sysGmmModels}(object, theta = NULL, w="optimal",
wObj=NULL)

\S4method{evalWeights}{rslinearGmm}(object, theta = NULL, w="optimal",
wObj=NULL)
}
\arguments{
  \item{object}{Object of class \code{gmmModels}}
  \item{theta}{The vector of coefficients to compute the optimal
    weights. If \code{NULL}, \code{theta0} for the object is used.}
  \item{w}{A matrix for fixed weights, one of \code{"optimal"} or
    \code{"ident"}}
  \item{wObj}{An object of class \code{"sysGmmWeights"}. Providing it
    avoid having to recompute Z'Z.}
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{
\item{\code{signature(object = "gmmModels")}}{
}

\item{\code{signature(object = "sysGmmModels")}}{
}

\item{\code{signature(object = "rslinearGmm")}}{
System of equations with restrictions on the coefficients. It only
affects the computation of the weights when there are cross-equation
restrictions.
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)

## Identity weights object
wObj1 <- evalWeights(model1, w="ident")

## Identity weights object (an alternative way less efficient)
wObj1 <- evalWeights(model1, w=diag(3))

## Optimal weights 
wObj1 <- evalWeights(model1, theta, w="optimal")


}
\keyword{methods}
