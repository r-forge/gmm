\name{splineMatrix}
\alias{splineMatrix}
\title{
Spline design matrix.
}
\description{
The function builds a matrix for piecewise polynomials fit. 
}
\usage{
splineMatrix(X, knots = NA, pFact=0.3, deg=1, method=c("manual","bs"))
}
\arguments{
  \item{X}{A numeric vector}
  \item{knots}{The piecewise polynomial knots. If set to \code{NA}, the
    knots are set to the \code{p} equally spaced quantiles of \code{X}. If
    \code{NULL}, the function returns \code{X}.}
  \item{pFact}{The number of knots when the argument \code{knots} is set to
    \code{NA} if \code{n^pFact}, where n is the length of \code{X}.}
  \item{deg}{The degree of the piecewise polynomials.}
  \item{method}{Should the method be homemade (manual) of based on the
    \code{bs} function from the splines package?}
}

\examples{
data(simData)
X <- splineMatrix(simData$X, quantile(simData$X, c(.25,.5,.75)))
}

\keyword{causal effects}
