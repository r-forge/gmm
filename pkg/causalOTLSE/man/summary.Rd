\name{summary}
\alias{summary.otlse}
\title{Summary for objects of class otlse}
\description{
The method builds the coefficient matrices. 
}
\usage{
\method{summary}{otlse}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{otlse}}
\item{...}{Argument for other type of objects}
}

\examples{
data(simData)
fit <- otlse(simData$X, simData$Y, simData$Z) 
summary(fit)
}

