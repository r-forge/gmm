\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,causalfit-method}
\title{Methods for Function \code{summary} in Package \pkg{base}}
\description{
Compute several results from a moment based model fit.
}
\usage{

\S4method{summary}{causalfit}(object, \dots)

}

\arguments{
  \item{object}{A fit object from the package}
  \item{\dots}{Other arguments to pass to other methods
    (not currently used)}
  }
\section{Methods}{
\describe{
\item{\code{signature(object = "causalfit")}}{
}
}}

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
ps <- treat~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
           
fit1 <- matching(g, balm, nsw, "ACE")
summary(fit1)

}

\keyword{methods}

