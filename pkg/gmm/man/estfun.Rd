\name{estfun}
\alias{estfun.gmmFct}
\alias{estfun.gmm}
\alias{estfun.gel}
\title{Extracts the empirical moment function}
\description{
It extracts the matrix of empirical moments so that it can be used by the \code{\link{kernHAC}} function. 
}
\usage{
\method{estfun}{gmmFct}(x, y = NULL, theta = NULL, ...)
\method{estfun}{gmm}(x, ...)
\method{estfun}{gel}(x, ...)
}
\arguments{
\item{x}{A function of the form \eqn{g(\theta,y)} or a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,y_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n} or an object of class \code{gmm}. See \code{\link{gmm}} for more details.}
\item{y}{The matrix or vector of data from which the function \eqn{g(\theta,y)} is computed if \code{g} is a function.}
\item{theta}{Vector of parameters if \code{g} is a function.}
\item{...}{Other arguments when \code{estfun} is applied to another class object}
}

\details{
For \code{estfun.gmmFct}, it returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,y_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. It is only used by \code{gmm} to obtain the estimates.

For \code{estfun.gmm}, it returns the matrix of first order conditions of \eqn{\min_\theta \bar{g}'W\bar{g}/2}, which is a \eqn{n \times q} matrix with the \eqn{t^{th}} row being \eqn{g(\theta, y_t)W G}, where \eqn{G} is \eqn{d\bar{g}/d\theta}. It allows to compute the sandwich covariance matrix using \code{link{kernHAC}} or \code{link{vcovHAC}} when \eqn{W} is not the optimal matrix.

The method if not yet available for \code{gel} objects.
}

\value{
A \eqn{n \times q} matrix (see details). 
}

\references{
   Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.

}


