\name{otlse}
\alias{otlse}
\title{
Optimal Thresholding Least Squares
}
\description{
This is the main function to estimate the causal effects using the
optimal thresholding least squares method.
}
\usage{
otlse(X, Y, Z, crit = c("ASY", "AIC", "BIC", "CV"),
      pFact=0.3, splineMet=c("manual","bs"))
}
\arguments{
  \item{X}{A vector of covariates}
  \item{Y}{A vector of observed outcomes}
  \item{Z}{A vector of treatment indicators}
  \item{crit}{The method to select the piecewise polynomial knots.}
  \item{pFact}{The maximum number of knots when the argument \code{knots} is set to
    \code{NA} if \code{n^pFact}, where n is the length of \code{X}.}
  \item{splineMet}{Should the method be homemade (manual) of based on the
    \code{bs} function from the splines package?}
}

\examples{
data(simData)
fit <- otlse(simData$X, simData$Y, simData$Z) 
fit
}

\keyword{causal effects}
