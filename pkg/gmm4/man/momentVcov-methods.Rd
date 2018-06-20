\name{momentVcov-methods}
\docType{methods}
\alias{momentVcov-methods}
\alias{momentVcov}
\alias{estfun.gmmFct}
\alias{momentVcov,gmmModels-method}
\alias{momentVcov,sysGmmModels-method}
\title{ ~~ Methods for Function \code{momentVcov} in Package \pkg{Gmm} ~~}
\description{
It computes the estimated covariance matrix of the moment
conditions based on the it is defined in the \code{"gmmModels"}
class object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "gmmModels")}}{
}

\item{\code{signature(object = "sysGmmModels")}}{
}
}}
\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)
momentVcov(model1, theta)[,]

model2 <- gmmModel(y~x1, ~z1+z2, data=simData, vcov="iid")
momentVcov(model2, theta)

model3 <- gmmModel(y~x1, ~z1+z2, data=simData, vcov="MDS")
momentVcov(model3, theta)
}

\keyword{methods}
