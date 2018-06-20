\name{momentStrength-methods}
\docType{methods}
\alias{momentStrength}
\alias{momentStrength-methods}
\alias{momentStrength,functionGmm-method}
\alias{momentStrength,linearGmm-method}
\alias{momentStrength,rlinearGmm-method}
\alias{momentStrength,nonlinearGmm-method}
\title{ ~~ Methods for Function \code{momentStrength} in Package \pkg{gmm4} ~~}
\description{
It produces measures of the strength of the moment conditons. 
}
\usage{
\S4method{momentStrength}{linearGmm}(object, theta,
vcovType=c("OLS","HC","HAC"))
}
\arguments{
  \item{object}{An object of class \code{"linearGmm"}}
  \item{theta}{Coefficient vector at which the strength must be
    measured. It does not impact the measure for objects of class
    \code{linearGmm}.}
  \item{vcovType}{Type of covariance matrix used to
    compute the F-test of the first-stage regression. For \code{HC},
    the function \code{\link{vcovHC}} is used with "HC1", and
    \code{\link{vcovHAC}} is used with the default setup is "HAC" is
    chosen. In \code{summary} for \code{gmmfit} objects, it is adjusted
    to the type of covariance that is set in the object} 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "functionGmm")}}{
Not implemented yet. In that case, we want some measure of the rank of
the matrix of derivatives. 
}

\item{\code{signature(object = "linearGmm")}}{
It returns the F-test of the first stage regression. It is a measure of
the strength of the instruments.
}

\item{\code{signature(object = "rlinearGmm")}}{
Returns nothing for now.
}

\item{\code{signature(object = "nonlinearGmm")}}{
Not implemented yet. 
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)
momentStrength(model1)
}
\keyword{methods}

