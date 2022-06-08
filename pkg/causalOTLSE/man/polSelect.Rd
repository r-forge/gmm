\name{polSelect}
\alias{selASY}
\alias{selIC}
\alias{selCV}
\title{
Polynomial selection methods.
}
\description{
This is a collection of methods to select the piecewise polynomials.
}
\usage{
selASY(X, Y, Z, pFact=0.3, splineMet=c("manual","bs"))
selIC(X, Y, Z, pFact=0.3, type=c("AIC", "BIC"), splineMet=c("manual","bs"))
selCV(X, Y, Z, pFact=0.3, splineMet=c("manual","bs"))
}
\arguments{
  \item{X}{A vector of covariates}
  \item{Y}{A vector of observed outcomes}
  \item{Z}{A vector of treatment indicators}
  \item{pFact}{The maximum number of knots when the argument \code{knots} is set to
    \code{NA} if \code{n^pFact}, where n is the length of \code{X}.}    
  \item{type}{The type of informatrion criterion}
  \item{splineMet}{Should the method be homemade (manual) of based on the
    \code{bs} function from the splines package?}  
}

\keyword{selection, polynomial}
