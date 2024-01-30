\name{CDtest}
\alias{CDtest}
\alias{SYTables}
\title{
Cragg and Donald test for Weak Instruments
}
\description{
The Cragg and Donald test is based on the F-test from the first stage
regression. For models with one endogenous variable, it is the usual
F-test. For more than one endogenous variables, it is the smallest
eigenvalue of the matrix representation of the F statistic. Along with
the statistics, the critical values of Stock and Yogo (2005) for the
null hypothesis of weak instruments are provided. 
}
\usage{

CDtest(object, print=TRUE, ...)

SYTables(object, print=TRUE)

}

\arguments{
  \item{object}{
    A model of class \code{linearModel}
  }
  \item{print}{
    If code{TRUE}, the result is printed and nothing is returned. See
    below for details on what the functions return when it is set to
    \code{FALSE}
  }
  \item{...}{Arguments passed to \code{\link{formatC}} to print the statistic.
  }
}

\details{ The \code{CDtest} function computes the test for the model
  object and runs \code{SYTables} to extract the critical values that are
  specific to the model. }

\value{
  The function \code{CDtest} returns the Cragg anb Donald statistic when
  \code{print} is set to \code{FALSE}.

  The function \code{SYTable} returns a list with the following elements
  when \code{print} is set to \code{FALSE}:

  \item{biasTSLS, biasFuller}{Named vectors of critical values for TSLS
  or Fuller (see \code{\link{kclassfit}}). These critical values are
  used to achieve the maximum relative bias given by their names. The
  values are defined only for models with a number of overidentifying
  restrictions greater or equal to 2.} 

  \item{sizeTSLS, sizeLIML}{A named vector of critical values for TSLS
  or LIML (see \code{\link{kclassfit}}). These critical values are used
  to acheive a size that does not exceed the one given by their names.}

 Both functions return \code{NULL} invisibly when \code{print} is set to
 \code{TRUE}

}  
\examples{

data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- momentModel(y~x1, ~z1+z2, data=simData)
CDtest(model1)
}
 
\keyword{LIML}
\keyword{Fuller}
\keyword{Weak Instruments Test}
\keyword{Two Stage Least Squares}

