\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint,ANY-method}
\alias{confint,gelfit-method}
\alias{confint,gmmfit-method}
\alias{confint,numeric-method}
\alias{confint,data.frame-method}
\alias{confint,ANY-method}
\title{ ~~ Methods for Function \code{confint} in Package \pkg{stats} ~~}
\description{
Method to contruct confidence intervals for objects of class
\code{"gmmfit"} and \code{"gelfit"}.
}

\usage{
\S4method{confint}{gmmfit}(object, parm, level = 0.95, vcov=NULL, \dots)

\S4method{confint}{gelfit}(object, parm, level = 0.95, lambda = FALSE,
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, corr = NULL, vcov=NULL, \dots)

\S4method{confint}{data.frame}(object, parm, level = 0.95, gelType="EL", 
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, vcov="iid", n=10, 
                    cores=4, \dots) 

\S4method{confint}{ANY}(object, parm, level = 0.95, \dots)
}

\arguments{

  \item{object}{Object of class \code{"gmmfit"},
    \code{"gelfit"}, \code{"numeric"} or  \code{"data.frame"}.}

  \item{parm}{Vector of integers or characters for selecting the
    elements for which the intervals should be computed.}

  \item{level}{The confidence level.}

  \item{lambda}{Should be compute intervals for the Lagrange
    multipliers?}
  
  \item{type}{The type of confidence intervals. The default is the Wald
    interval, and the others are computed by inverting the LR, LM or J
    specification test.}

  \item{fact}{For the inversion of the specification tests,
    \code{\link{uniroot}} searches within \code{fact} standard error of
    the coefficient estimates}

  \item{corr}{Correction to apply to the specification tests}
  
  \item{vcov}{For Wald intervals, an optional covariance matrix can be
    provided}

  \item{cores}{The number of cores for \code{\link{mclapply}}. It is set
    to 1 for Windows OS.}

  \item{gelType}{Type of GEL confidence interval.}

  \item{n}{Number of lines to construct the confidence region.}

  \item{\dots}{Other arguments to pass to \code{\link{modelFit}}}
}

\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
The method from the \pkg{stats} in used in that case.
}

\item{\code{signature(object = "gelfit")}}{
Method for any GEL fit class.
}

\item{\code{signature(object = "gmmfit")}}{
Method for any GMM fit class.
}

\item{\code{signature(object = "numeric")}}{
It computes the GEL confidence interval for the mean.
}

\item{\code{signature(object = "data.frame")}}{
  It computes the 2D GEL confidence region for the means of two
  variables.
  
}



}}
\keyword{methods}
\keyword{Confidence Intervals}
