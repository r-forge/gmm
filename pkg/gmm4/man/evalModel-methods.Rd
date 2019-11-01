\name{evalModel-methods}
\docType{methods}
\alias{evalModel}
\alias{evalModel-methods}
\alias{evalModel,gelModels-method}
\alias{evalModel,gmmModels-method}
\title{ ~~ Methods for Function \code{evalModel} in Package \pkg{gmm4} ~~}
\description{
Method to simply evaluate a model at a fixed coefficient vector. It
creates a \code{"gmmfit"} or \code{"gelfit"} object using that fixed vector. 
}
\usage{
\S4method{evalModel}{gmmModels}(object, theta, wObj=NULL, \dots)
\S4method{evalModel}{gelModels}(object, theta, lambda=NULL,
          gelType=NULL, rhoFct=NULL, lamSlv=NULL, lControl=list(),
          \dots)


}
\arguments{
  \item{object}{An object of class \code{"gmmModels"}.}
  \item{theta}{A vector of coefficients at which the model is estimated}
  \item{wObj}{An object of class \code{"gmmWeights"}. If not provided,
    the optimal weights based on the specification of the model
    evaluated at \code{theta} will be computed. }
  \item{lambda}{The Lagrange multiplier vector. If not provided, the
    optimal vector is obtained for the given theta}
  \item{gelType}{The type of GEL. This argument is only used if we want
    to fit the model with a different GEL method. see \code{\link{gelModel}}.}
  \item{rhoFct}{An alternative objective function for GEL. This argument
    is only used if we want to fit the model with a different GEL
    method. see \code{\link{gelModel}}.}
  \item{lamSlv}{An alternative solver for the Lagrange multiplier. By
    default, either \code{\link{Wu_lam}}, \code{\link{EEL_lam}},
    \code{\link{REEL_lam}} or \code{\link{getLambda}} is
    used.}
  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}
  \item{\dots}{Other arguments to pass. Not used for the moment.}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "gmmModels")}}{
%%  ~~describe this method here~~
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearGmm
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)
evalModel(model1, c(1,1))

## A nonlinearGmm
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- gmmModel(g, h, c(beta0=1, beta1=2), data=simData)
evalModel(model2, theta=c(beta1=2, beta0=0.5))

## A functionGmm
fct <- function(tet, x)
    {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
    }
dfct <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
			   -6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }
model3 <- gmmModel(fct, simData$x3, theta0=c(beta0=1, beta1=2), grad=dfct)
evalModel(model3, theta=c(beta1=.1, beta0=0.3))
}
\keyword{methods}
