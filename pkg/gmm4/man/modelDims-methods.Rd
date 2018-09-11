\name{modelDims-methods}
\docType{methods}
\alias{modelDims}
\alias{modelDims-methods}
\alias{modelDims,linearGmm-method}
\alias{modelDims,rslinearGmm-method}
\alias{modelDims,sysGmmModels-method}
\alias{modelDims,slinearGmm-method}
\alias{modelDims,snonlinearGmm-method}
\alias{modelDims,nonlinearGmm-method}
\alias{modelDims,functionGmm-method}
\alias{modelDims,formulaGmm-method}
\alias{modelDims,rformulaGmm-method}
\alias{modelDims,rlinearGmm-method}
\alias{modelDims,rfunctionGmm-method}
\alias{modelDims,rnonlinearGmm-method}
\title{Methods for Function \code{modelDims}}
\description{
It extracts important information from the model. It is mostly used by
other methods when a \code{gmmModels} has been modifed. An example is
when restrictions have been imposed on coefficients.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "rlinearGmm")}}{
}

\item{\code{signature(object = "rnonlinearGmm")}}{
}

\item{\code{signature(object = "rfunctionGmm")}}{
}

\item{\code{signature(object = "linearGmm")}}{
}

\item{\code{signature(object = "nonlinearGmm")}}{
}

\item{\code{signature(object = "functionGmm")}}{
}

\item{\code{signature(object = "formulaGmm")}}{
}

\item{\code{signature(object = "rformulaGmm")}}{
}

\item{\code{signature(object = "slinearGmm")}}{
}

\item{\code{signature(object = "rslinearGmm")}}{
}

\item{\code{signature(object = "snonlineatGmm")}}{
}
}}
\examples{
data(simData)

model1 <- gmmModel(y~x1+x2, ~x2+z1+z2, data=simData)
modelDims(model1)

## Unrestricted model

rmodel1 <- restGmmModel(model1, R=c("x1+x2=4"))
modelDims(rmodel1)


}

\keyword{methods}
