\name{print}
\alias{print.otlse}
\alias{print.summary.otlse}
\title{Print methods}
\description{
Print results of different methods applied to otlse Objects
}
\usage{
\method{print}{otlse}(x, ...)
\method{print}{summary.otlse}(x, digits = 4,
                              signif.stars = getOption("show.signif.stars"),
                              beta=FALSE, ...)
}
\arguments{
  \item{x}{Onject of class \code{otlse} or \code{summary.otlse}.}
  \item{digits}{The number of digits to print.}
  \item{signif.stars}{Should we print the significant stars?}
  \item{beta}{Should we print the coefficient matrix for the piecewise polynomials?}
\item{...}{Argument for other type of objects}
}

\examples{
data(simData)
fit <- otlse(simData$X, simData$Y, simData$Z) 
print(summary(fit), digits=5)
}












