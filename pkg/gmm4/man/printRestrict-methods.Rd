\name{printRestrict-methods}
\docType{methods}
\alias{printRestrict}
\alias{printRestrict-methods}
\alias{printRestrict,rlinearGmm-method}
\alias{printRestrict,rslinearGmm-method}
\alias{printRestrict,rnonlinearGmm-method}
\alias{printRestrict,rfunctionGmm-method}
\title{ ~~ Methods for Function \code{printRestrict} in Package \pkg{gmm4} ~~}
\description{
It prints the detailed restrictions imposed on \code{gmmModels}.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "rlinearGmm")}}{
}

\item{\code{signature(object = "rnonlinearGmm")}}{
}

\item{\code{signature(object = "rfunctionGmm")}}{
}

\item{\code{signature(object = "rslinearGmm")}}{
}
}}
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## Unrestricted model
model1 <- gmmModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)

## restricted model
R <- matrix(c(1,1,0,0,0,0,0,2,0,0,0,0,0,1,-1),3,5, byrow=TRUE)
q <- c(0,1,3)
rmodel1 <- restGmmModel(model1, R, q)
printRestrict(rmodel1)
}
\keyword{methods}
