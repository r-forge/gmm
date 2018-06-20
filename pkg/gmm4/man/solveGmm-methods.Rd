\name{solveGmm-methods}
\docType{methods}
\alias{solveGmm}
\alias{solveGmm-methods}
\alias{solveGmm,allNLGmm,gmmWeights-method}
\alias{solveGmm,linearGmm,gmmWeights-method}
\alias{solveGmm,slinearGmm,sysGmmWeights-method}
\alias{solveGmm,rslinearGmm,sysGmmWeights-method}
\alias{solveGmm,snonlinearGmm,sysGmmWeights-method}
\alias{solveGmm,rnonlinearGmm,gmmWeights-method}
\title{ ~~ Methods for Function \code{solveGmm} in Package \pkg{Gmm} ~~}
\description{
 ~~ Methods for function \code{solveGmm} in package \pkg{Gmm} ~~
}
\usage{
\S4method{solveGmm}{linearGmm,gmmWeights}(object, wObj, theta0=NULL,
...)

\S4method{solveGmm}{allNLGmm,gmmWeights}(object, wObj, theta0=NULL,
...)

\S4method{solveGmm}{rnonlinearGmm,gmmWeights}(object, wObj, theta0=NULL,
...)

\S4method{solveGmm}{slinearGmm,sysGmmWeights}(object, wObj, theta0=NULL)

\S4method{solveGmm}{rslinearGmm,sysGmmWeights}(object, wObj, theta0=NULL)

\S4method{solveGmm}{snonlinearGmm,sysGmmWeights}(object, wObj,
theta0=NULL, \dots)
}
\arguments{
  \item{object}{An object of class \code{"linearGmm"} or a subclass of
    the class \code{"allNLGmm"}, \code{"sysGmmModels"} or \code{"rsysGmmModels".}}
  \item{theta0}{The vector of coefficients for the starting values used
  in \code{\link{optim}}. If NULL, the starting values in the object if
  used. For system of equations, it is a list of vectors.}
  \item{wObj}{An object of class \code{"gmmWeights"} or \code{"sysGmmWeights"}.}
  \item{\dots}{Arguments to pass to \code{\link{optim}}.}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "allNLGmm", wObj = "gmmWeights")}}{
Method to solve either nonlinear regressions or models in which moments
are computed with a function. The objective is minimized using optim.
}

\item{\code{signature(object = "rnonlinearGmm", wObj = "gmmWeights")}}{
Method to solve restricted nonlinear models. It computes the analytical
solution.
}

\item{\code{signature(object = "linearGmm", wObj = "gmmWeights")}}{
Method to solve linear models. It computes the analytical solution.
}

\item{\code{signature(object = "slinearGmm", wObj = "sysGmmWeights")}}{
Method to solve system of linear models. It computes the analytical solution.
}

\item{\code{signature(object = "rslinearGmm", wObj = "sysGmmWeights")}}{
Method to solve system of linear models in which restrictions have been
imposed on the coefficients. It computes the analytical solution.
}

\item{\code{signature(object = "slinearGmm", wObj = "sysGmmWeights")}}{
Method to solve system of nonlinear models. The solution is obtained
with optim using the analytical derivatives.
}

}}
\value{
  A list with the following:
  \item{theta}{The vector of solution}
  \item{convergence}{convergence code. 0 means normal convergence. For
    higher numbers, see \code{\link{optim}}}
  }
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)

## A manual two-step GMM
w0 <- evalWeights(model1, w="ident")
theta0 <- solveGmm(model1, w0)$theta
w <- evalWeights(model1, theta0)
theta1 <- solveGmm(model1, w)$theta 
}

\keyword{methods}

