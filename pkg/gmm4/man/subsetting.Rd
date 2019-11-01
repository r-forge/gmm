\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{subset}
\alias{subset,regGmm-method}
\alias{subset,functionGmm-method}
\alias{subset,formulaGmm-method}
\alias{subset,sysGmmModels-method}
\alias{[,gmmWeights,numeric,missing-method}
\alias{[,gmmWeights,missing,missing-method}
\alias{[,sysGmmWeights,numeric,list-method}
\alias{[,sysGmmWeights,numeric,missing-method}
\alias{[,sysGmmWeights,missing,list-method}
\alias{[,gmmModels,missing,missing-method}
\alias{[,regGmm,numeric,missing-method}
\alias{[,functionGmm,numeric,missing-method}
\alias{[,formulaGmm,numeric,missing-method}
\alias{[,rfunctionGmm,numeric,missing-method}
\alias{[,rfunctionGel,numeric,missing-method}
\alias{[,sysGmmModels,missing,list-method}
\alias{[,snonlinearGmm,numeric,missing-method}
\alias{[,slinearGmm,numeric,missing-method}
\alias{[,rslinearGmm,numeric,missing-method}
\alias{[,sysGmmModels,numeric,list-method}
\alias{[,sysGmmModels,missing,missing-method}

\title{Subsetting methods}
\description{
Different subsetting methods for S4 class objects of the package. The
subset method returns an new object with observations selected by the
second argument. See example.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "gmmWeights", i = "integer", j = "missing")}}{
  It creates a partition from the weighting matrix. 
}

\item{\code{signature(x = "gmmWeights", i = "missing", j = "missing")}}{  
  It generates the whole weighting matrix. 
}


\item{\code{signature(x = "sysGmmWeights", i = "missing", j = "list")}}{
  It creates a partition from the weighting matrix. j has no effect
  here.  It creates a partition from the weighting matrix in a systemof  
  equations. i selects the equation and the list j the moments in each
  equation. Missing i means all equations.
}

\item{\code{signature(x = "sysGmmWeights", i = "numeric", j =
    "missing")}}{
    It creates a partition from the weighting matrix. j has no effect
  here.  It creates a partition from the weighting matrix in a systemof  
  equations. i selects the equation and the list j the moments in each
  equation. Missing j means all moments.
}

\item{\code{signature(x = "sysGmmWeights", i = "missing", j = "missing")}}{
  No effect. It returns x.
}


\item{\code{signature(x = "snonlinearGmm", i = "numeric", j="missing")}}{
  It generates a system of equations with a subset of equations selected
  by i. If the number of remaining equations is one, it returns an
  object of class \code{"nonlinearGmm"}.
}


\item{\code{signature(x = "slinearGmm", i = "numeric", j="missing")}}{
  It generates a system of equations with a subset of equations selected
  by i. If the number of remaining equations is one, it returns an
  object of class \code{"linearGmm"}.
}

\item{\code{signature(x = "rslinearGmm", i = "numeric", j="missing")}}{
  It is only use to select one equation when no cross-equation
  restrictions are imposed. Only one equation can be selected. 
}

\item{\code{signature(x = "sysGmmModels", i = "numeric", j="list")}}{
  It generates a system of equations with a subset of equations selected
  by i and a subset of moment conditions selected by j. If the number of
  remaining equations is one, it returns an object of class \code{"linearGmm"}.
}

\item{\code{signature(x = "sysGmmModels", i = "missing", j="missing")}}{
  No effect. It returns x.
}


\item{\code{signature(x = "gmmModels", i = "missing", j = "missing")}}{
Returns the model without any change.
}

\item{\code{signature(x = "functionGmm", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions.
}

\item{\code{signature(x = "formulaGmm", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions.
}

\item{\code{signature(x = "rfuncionGmm", i = "numeric", j = "missing")}}{
  It generates the same model with a subset of the moment conditions. j
  has no effect here.
}

}}
\examples{
set.seed(1122)
x1 <- rchisq(50,5)
x2 <- rchisq(50,5)
x3 <- rnorm(50)
z1 <- .2*x1+rnorm(50)
z2 <- rnorm(50)
z3 <- rnorm(50)
y <- x1+rnorm(50)
dat <- data.frame(y=y,z1=z1,x1=x1,z2=z2,x2=x2,z3=z3,x3=x3)
model1 <- gmmModel(y~x1+x2, ~x2+x3+z1+z2+z3, data=dat, vcov="MDS")
w <- evalWeights(model1, theta=1:3)
w[]
w[1:3]

## A model with a subset of the instruments
model1[1:4]

## Selecting the observations:

subset(model1, x1<3)
subset(model1, 1:25)

}

\keyword{methods}
