\name{coef-methods}
\docType{methods}
\alias{coef-methods}
\alias{coef,gmmfit-method}
\alias{coef,sgmmfit-method}
\alias{coef,rlinearGmm-method}
\alias{coef,rslinearGmm-method}
\alias{coef,rfunctionGmm-method}
\alias{coef,rformulaGmm-method}
\alias{coef,gmmModels-method}
\alias{coef,rnonlinearGmm-method}
\title{ ~~ Methods for Function \code{coef} in Package \pkg{stats} ~~}
\description{
It extract the coefficient estimates of some GMM objects.
}
\section{Methods}{
\describe{
 
\item{\code{signature(object = "gmmfit")}}{
}

\item{\code{signature(object = "sgmmfit")}}{
}

\item{\code{signature(object = "gmmModels")}}{
}

\item{\code{signature(object = "rlinearGmm")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rslinearGmm")}}{
It gives the unrestricted representation of a restricted model.  
}

\item{\code{signature(object = "rfunctionGmm")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rformulaGmm")}}{
It gives the unrestricted representation of a restricted model. See examples.  
}

\item{\code{signature(object = "rnonlinearGmm")}}{
It gives the unrestricted representation of a restricted nonlinear
model. 
}
}}
\examples{
data(simData)
model1 <- gmmModel(y~x1+x2+x3+z1, ~x1+x2+z1+z2+z3+z4, data=simData)
res1 <- modelFit(model1)
coef(res1)

### Restricted models
rmodel1 <- restGmmModel(model1, R=c("x1=1", "x2=2*x3"))
res2 <- modelFit(rmodel1)
res2
coef(rmodel1, coef(res2))
}
\keyword{methods}
