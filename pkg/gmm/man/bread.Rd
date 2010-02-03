\name{bread}
\alias{bread.gmm}
\alias{bread.gel}
\title{Bread for sandwiches}
\description{
Computes the bread of the sandwich covariance matrix
}
\usage{
\method{bread}{gmm}(x, ...)
\method{bread}{gel}(x, ...)
}
\arguments{
\item{x}{A fitted model of class \code{gmm} or \code{gel}.}
\item{...}{Other arguments when \code{bread} is applied to another class object}
}

\details{
When the weighting matrix is not the optimal one, the covariance matrix of the estimated coefficients is:
\eqn{(G'WG)^{-1} G'W V W G(G'WG)^{-1}}, 
where \eqn{G=d\bar{g}/d\theta}, \eqn{W} is the matrix of weights, and \eqn{V} is the covariance matrix of the moment function. Therefore, the bread is \eqn{(G'WG)^{-1}}, which is the second derivative of the objective function. 

The method if not yet available for \code{gel} objects.
}

\value{
A \eqn{k \times k} matrix (see details). 
}

\references{
   Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.

}


