\name{quadra-methods}
\docType{methods}
\alias{quadra}
\alias{quadra-methods}
\alias{quadra,gmmWeights,matrixORnumeric,matrixORnumeric-method}
\alias{quadra,gmmWeights,matrixORnumeric,missing-method}
\alias{quadra,gmmWeights,missing,missing-method}
\alias{quadra,sysGmmWeights,matrixORnumeric,matrixORnumeric-method}
\alias{quadra,sysGmmWeights,matrixORnumeric,missing-method}
\alias{quadra,sysGmmWeights,missing,missing-method}

\title{ ~~ Methods for Function \code{quadra} in Package \pkg{Gmm} ~~}
\description{
 ~~ Computes the quadratic form, where the center matrix is a class
 \code{gmmWeights} object ~~
}
\usage{
\S4method{quadra}{gmmWeights,missing,missing}(w, x, y)

\S4method{quadra}{gmmWeights,matrixORnumeric,missing}(w, x, y)

\S4method{quadra}{gmmWeights,matrixORnumeric,matrixORnumeric}(w, x, y)
}
\arguments{
  \item{w}{An object of class \code{"gmmWeights"}}
  \item{x}{A matrix or numeric vector}
  \item{y}{A matrix or numeric vector}
}
\section{Methods}{
\describe{
\item{\code{signature(w = "gmmWeights", x = "matrixORnumeric",  y =
    "matrixORnumeric")}}{
It computes \eqn{x'Wy}, where \eqn{W} is the weighting matrix.
}
\item{\code{signature(w = "gmmWeights", x = "matrixORnumeric",  y =
    "missing")}}{
It computes \eqn{x'Wx}, where \eqn{W} is the weighting matrix.  
}
\item{\code{signature(w = "gmmWeights", x = "missing",  y =
    "missing")}}{
It computes \eqn{W}, where \eqn{W} is the weighting matrix.  When
\eqn{W} is the inverse of the covariance matrix of the moment
conditions, it is saved as either a QR decompisition, a Cholesky
decomposition or a covariance matrix into the \code{\link{gmmWeights}}
object. The \code{quadra} method with no \code{y} and \code{x} is
therefore a way to invert it. The same applies to system of equations

}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)

gbar <- evalMoment(model1, theta)
gbar <- colMeans(gbar)

### Onjective function of GMM with identity matrix
wObj <- evalWeights(model1, w="ident")
quadra(wObj, gbar)

### Onjective function of GMM with efficient weights
wObj <- evalWeights(model1, theta)
quadra(wObj, gbar)

}
\keyword{methods}
