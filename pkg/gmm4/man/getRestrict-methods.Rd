\name{getRestrict-methods}
\docType{methods}
\alias{getRestrict}
\alias{getRestrict-methods}
\alias{getRestrict,gmmModels-method}
\alias{getRestrict,rlinearGmm-method}
\alias{getRestrict,rslinearGmm-method}
\alias{getRestrict,sysGmmModels-method}
\alias{getRestrict,rnonlinearGmm-method}
\alias{getRestrict,rfunctionGmm-method}
\title{ ~~ Methods for Function \code{getRestrict} in Package \pkg{gmm4} ~~}
\description{
It computes the matrices related to linear and nonlinear
contraints. Those matrices are used to perform hypothesis tests.
}
\usage{
\S4method{getRestrict}{rlinearGmm}(object, theta)

\S4method{getRestrict}{rslinearGmm}(object, theta)

\S4method{getRestrict}{rnonlinearGmm}(object, theta)

\S4method{getRestrict}{gmmModels}(object, theta, R, rhs=NULL)

\S4method{getRestrict}{sysGmmModels}(object, theta, R, rhs=NULL)

\S4method{getRestrict}{rfunctionGmm}(object, theta)

}
\arguments{
  \item{object}{Object of class included in \code{gmmModels},
    \code{rGmmModels}, and \code{rsysGmmModels}.}
  \item{theta}{A vector of coefficients for the unrestricted model (see
    examples).}
  \item{R}{A matrix, character or list of formulas that specifies the
    contraints to impose on the coefficients. See
    \code{\link{restGmmModel}} for more details.}
  \item{rhs}{The right hand side for the restriction on the
    coefficients. See \code{\link{restGmmModel}} for more details. It is
    ignored for objects of class \code{"nonlinearGmm"}.}
  }
\section{Methods}{
\describe{

 \item{\code{signature(object = "gmmModels")}}{
    A restricted model is created from the constraints, and the
    restriction matrices are returned. The methods is applied to linear
    and nonlinear models in a regression form.
}

\item{\code{signature(object = "sysGmmModels")}}{
    A restricted model is created from the constraints, and the
    restriction matrices are returned. The methods is applied to systems
    of linear and nonlinear models.
}

\item{\code{signature(object = "rlinearGmm")}}{
  The restriction matrices are evaluated at the coefficient vector
  \code{theta} of the unrestricted representation.
}

\item{\code{signature(object = "rslinearGmm")}}{
  The restriction matrices are evaluated at the coefficient vector
  \code{theta} of the unrestricted representation.
}

\item{\code{signature(object = "rnonlinearGmm")}}{
  The restriction matrices are evaluated at the coefficient vector
  \code{theta} of the unrestricted representation.  
}

\item{\code{signature(object = "rfunctionGmm")}}{
  The restriction matrices are evaluated at the coefficient vector
  \code{theta} of the unrestricted representation.  
}

}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## Unrestricted model
model1 <- gmmModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)

## The restricted model
R1 <- c("x1","2*x2+z1=2", "4+x3*5=3")
res <- gmmFit(model1)
rest <- getRestrict(model1, coef(res), R1)

## it allows to test the restriction
g <- rest$R-rest$q
v <- rest$dR\%*\%vcov(res)\%*\%t(rest$dR)
(test <- crossprod(g, solve(v, g)))
(pv <- 1-pchisq(test, length(rest$R)))


## Delta Method:
## To impose nonlinear restrictions, we need to convert
## the linear model into a nonlinear one
NLmodel <- as(model1, "nonlinearGmm")
R1 <- c("theta2=2", "theta3=theta4^2")
res <- gmmFit(NLmodel)
rest <- getRestrict(NLmodel, coef(res), R1)

g <- rest$R-rest$q
v <- rest$dR\%*\%vcov(res)\%*\%t(rest$dR)
(test <- crossprod(g, solve(v, g)))
(pv <- 1-pchisq(test, length(rest$R)))

## See hypothesisTest method for an easier approach.

}
\keyword{methods}
