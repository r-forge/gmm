\name{model.matrix-methods}
\docType{methods}
\alias{model.matrix}
\alias{model.matrix-methods}
\alias{model.matrix,linearGmm-method}
\alias{model.matrix,nonlinearGmm-method}
\alias{model.matrix,rlinearGmm-method}
\alias{model.matrix,slinearGmm-method}
\alias{model.matrix,rslinearGmm-method}
\alias{model.matrix,snonlinearGmm-method}
\title{ ~~ Methods for Function \code{model.matrix} in Package \pkg{stats} ~~}
\description{
 Model matrix form \code{gmmModels}. It returns the matrix of regressors
 or the instruments. In restricted models, it returns the reduced matrix
 of regressors.
}
\usage{
\S4method{model.matrix}{linearGmm}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{rlinearGmm}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{nonlinearGmm}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{slinearGmm}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{rslinearGmm}(object,
type=c("regressors","instruments"))
\S4method{model.matrix}{snonlinearGmm}(object,
type=c("regressors","instruments"))
}
\arguments{
\item{object}{Object of class \code{linearGmm}, \code{rlinearGmm} or any
system of equations class.}
\item{type}{Should the function returns the matrix of instruments or the
matrix of regressors. For \code{nonlinearGmm} classes, type='regressors'
will produce an error message, because there is no such model matrix in
this case, at least not for now.}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "linearGmm")}}{
Linear models with not restrictions.
}

\item{\code{signature(object = "nonlinearGmm")}}{
Nonlinear models with not restrictions.
}

\item{\code{signature(object = "rlinearGmm")}}{
  linear models with restrictions.
}
  
 \item{\code{signature(object = "slinearGmm")}}{
   System of linear equations with no restrictions. 
 }

 \item{\code{signature(object = "rslinearGmm")}}{
   System of linear equations with restrictions. 
 }
 
 \item{\code{signature(object = "snonlinearGmm")}}{
   System of nonlinear equations with no restrictions. 
 }
  }}
  \examples{
data(simData)

## Unrestricted model
model1 <- gmmModel(y~x1+x2+x3, ~x2+x3+z1+z2, data=simData)
model.matrix(model1)[1:3,]

## Restrictions change the response
R <- c("x2=2","x3+x1=3")
rmodel1 <- restGmmModel(model1, R)
rmodel1
model.matrix(rmodel1)[1:3,]
}
\keyword{methods}
