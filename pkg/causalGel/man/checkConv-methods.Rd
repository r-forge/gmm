\name{checkConv-methods}
\docType{methods}

\alias{checkConv}
\alias{checkConv-methods}
\alias{checkConv,causalGelfit-method}

\title{Methods that returns information about convergence.}
\description{
Returns convergence codes and compares moments for the treated and the
control groups.
}
\usage{

\S4method{checkConv}{causalGelfit}(object, tolConv=1e-4, verbose=TRUE, \dots)

}
\arguments{
  \item{object}{An object of class \code{"causalGelfit"}}
  \item{tolConv}{Tolerance to decide if balancing is achieved. }
  \item{verbose}{If FALSE nothinf is printed.}
   \item{\dots}{Arguments to pass to other methods.}
 }

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
           
model <- causalModel(g, balm, nsw)
fit <- gelFit(model, gelType="EL")
checkConv(fit)

}

  
\keyword{methods}
\keyword{Balancing moments}
