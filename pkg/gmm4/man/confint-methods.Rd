\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint,ANY-method}
\alias{confint,gelfit-method}
\alias{confint,gmmfit-method}
\title{ ~~ Methods for Function \code{confint} in Package \pkg{stats} ~~}
\description{
Method to contruct confidence intervals for objects of class
\code{"gmmfit"} and \code{"gelfit"}.
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
}}
\keyword{methods}
\keyword{Confidence Intervals}
