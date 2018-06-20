\name{modelResponse-methods}
\docType{methods}
\alias{modelResponse}
\alias{modelResponse-methods}
\alias{modelResponse,linearGmm-method}
\alias{modelResponse,slinearGmm-method}
\alias{modelResponse,rslinearGmm-method}
\alias{modelResponse,rlinearGmm-method}
\title{ ~~ Methods for Function \code{modelResponse} in Package \pkg{gmm4} ~~}
\description{
Return the response vector in models with and without restrictions
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearGmm")}}{
For linear models without restrictions on the coefficients.
}

\item{\code{signature(object = "slinearGmm")}}{
For system of linear models without restrictions on the coefficients.
}

\item{\code{signature(object = "rslinearGmm")}}{
For system of linear models with restrictions on the coefficients.
}

\item{\code{signature(object = "rlinearGmm")}}{
For linear models with restrictions on the coefficients.
}
}}
\examples{
data(simData)

## Unrestricted model
model1 <- gmmModel(y~x1+x2+x3, ~x2+x3+z1+z2, data=simData)
y <- modelResponse(model1)

## Restrictions change the response
R <- c("x2=2","x3=3")
rmodel1 <- restGmmModel(model1, R)
rmodel1
restY <- modelResponse(rmodel1)

}
\keyword{methods}
