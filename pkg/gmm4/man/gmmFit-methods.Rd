\name{gmmFit-methods}
\docType{methods}
\alias{gmmFit}
\alias{gmmFit-methods}
\alias{gmmFit,gmmModels-method}
\alias{gmmFit,sysGmmModels-method}
\alias{gmmFit,rlinearGmm-method}
\alias{gmmFit,rslinearGmm-method}
\alias{gmmFit,rnonlinearGmm-method}
\title{ ~~ Methods for Function \code{gmmFit} in Package \pkg{gmm4} ~~}
\description{
Method to fit a model using GMM, from an object of class \code{"gmmModels"}.
}
\usage{
\S4method{gmmFit}{gmmModels}(object, type=c("twostep", "iter","cue",
              "onestep"), itertol=1e-7, initW=c("ident", "tsls"),
              weights="optimal", itermaxit=100,
              efficientWeights=FALSE, start=NULL, \dots)

\S4method{gmmFit}{sysGmmModels}(object, type=c("twostep", "iter","cue",
              "onestep"), itertol=1e-7, initW=c("ident", "tsls", "EbyE"),
              weights="optimal", itermaxit=100,
              efficientWeights=FALSE, start=NULL, EbyE=FALSE, \dots)

\S4method{gmmFit}{rnonlinearGmm}(object, type=c("twostep", "iter","cue",
              "onestep"), itertol=1e-7, initW=c("ident", "tsls"),
              weights="optimal", itermaxit=100,
              efficientWeights=FALSE, start=NULL, \dots)

\S4method{gmmFit}{rlinearGmm}(object, type=c("twostep", "iter","cue",
              "onestep"), itertol=1e-7, initW=c("ident", "tsls"),
              weights="optimal", itermaxit=100,
              efficientWeights=FALSE, \dots)

\S4method{gmmFit}{rslinearGmm}(object, type=c("twostep", "iter","cue",
              "onestep"), itertol=1e-7, initW=c("ident", "tsls", "EbyE"),
              weights="optimal", itermaxit=100,
              efficientWeights=FALSE, start=NULL, EbyE=FALSE, \dots)
}
\arguments{
  \item{object}{An object of class \code{"gmmModels"}}
  \item{type}{What GMM methods should we use? for
    \code{type=="onestep"}, if \code{"weights"} is not a matrix, the
    model will be estimated with the weights equals to the identity
    matrix. For restricted }
  \item{itertol}{Tolance for the stopping rule in iterative GMM}
  \item{initW}{How should be compute the initial coefficient vector in
    the first. For single equation GMM, it only makes a difference for
    linear models for which the choice is GMM with identity matrix or
    two-stage least quares. For system of equations, \code{"tsls"},
    refers to equation by equation two-stage least squares. It is also
    possible to start at the equation by equation estimate using the
    same GMM type as specified by \code{"type"}.}    
  \item{weights}{What weighting matrix to use? The choices are
    \code{"optimal"}, in which case it is the inverse of the moment
    vovariance matrix, \code{"ident"} for the identity matrix, or a
    fixed matrix. It is also possible for weights to be an object of
    class \code{gmmWeights}.}
  \item{itermaxit}{Maximum iterations for iterative GMM}
  \item{efficientWeights}{If \code{weights} is a matrix or a
    \code{gmmWeights} class object, setting \code{efficientWeights} to
    \code{TRUE} implies that the resulting one-step GMM is
    efficient. As a result, the default covariance matrix for the
    coefficient estimates will not be a sandwich type.}
  \item{start}{An optional initial vector for \code{\link{optim}} when
    the model is nonlinear. By default, the theta0 argument of the model
    is used}
  \item{EbyE}{Should we estimate the system equation by equation?}
  \item{\dots}{Arguments to pass to other methods}
  }
\section{Methods}{
\describe{
\item{\code{signature(object = "gmmModels")}}{
The main method for all GMM-type models.
}

\item{\code{signature(object = "rnonlinearGmm")}}{
It makes a difference only if the number of contraints is equal to the
number of coefficients, in which case, the method \code{\link{evalGmm}}
is called at the contrained vector. If not, the next method is called.
}

\item{\code{signature(object = "rlinearGmm")}}{
It makes a difference only if the number of contraints is equal to the
number of coefficients, in which case, the method \code{\link{evalGmm}}
is called at the contrained vector. If not, the next method is called.
}

\item{\code{signature(object = "sysGmmModels")}}{
Method to estimate system of equations using GMM methods.
}

}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)

## Efficient GMM with HAC vcov and tsls as first step.
res1 <- gmmFit(model1, init="tsls")

## GMM with identity. Two ways.
res2 <- gmmFit(model1, type="onestep")
res3 <- gmmFit(model1, weights=diag(3))

## nonlinear regression with iterative GMM.
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- gmmModel(g, h, c(beta0=1, beta1=2), data=simData)
res4 <- gmmFit(model2, type="iter")

## GMM for with no endogenous vaiables is
## OLS with Robust standard error

library(lmtest)
model3 <- gmmModel(y~x1, ~x1, data=simData, vcov="MDS")
resGmm <- gmmFit(model3)
resLm <- lm(y~x1, simData)
summary(resGmm)
coeftest(resLm, vcov=vcovHC(resLm, "HC0"))
summary(resGmm, df.adj=TRUE)
coeftest(resLm, vcov=vcovHC(resLm, "HC1"))

### All constrained
R <- diag(2)
q <- c(1,2)
rmodel1 <- restGmmModel(model1, R, q)
gmmFit(rmodel1)

## Only one constraint
R <- matrix(c(0,1), ncol=2)
q <- 2
rmodel1 <- restGmmModel(model1, R, q)
gmmFit(rmodel1)

}
\keyword{methods}
