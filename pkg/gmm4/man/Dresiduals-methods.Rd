\name{Dresiduals-methods}
\docType{methods}
\alias{Dresiduals}
\alias{Dresiduals-methods}
\alias{Dresiduals,linearGmm-method}
\alias{Dresiduals,nonlinearGmm-method}
\alias{Dresiduals,sysGmmModels-method}
\title{ ~~ Methods for Function \code{Dresiduals} in Package \pkg{Gmm} ~~}
\description{
It returns the matrix of derivatives of the residuals with respect to
the coefficients.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearGmm")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "nonlinearGmm")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "sysGmmModels")}}{
%%  ~~describe this method here~~
}

}}

\examples{
data(simData)

theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)

Dresiduals(model1, theta)[1:3,]
}

\keyword{methods}