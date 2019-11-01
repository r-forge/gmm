\name{causalData-class}
\docType{class}
\alias{causalData-class}
\alias{causalData}
\alias{aceDmomFct,numeric,causalData-method}
\alias{aceMomFct,numeric,causalData-method}

\title{Class \code{"causalData"}}
\description{
On object used by the moment functions.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("causalData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{momType}:}{Object of class \code{"character"} ~~ }
    \item{\code{balCov}:}{Object of class \code{"character"} ~~ }    
    \item{\code{balMom}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{ACTmom}:}{Object of class \code{"integer"} ~~ }
    \item{\code{reg}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{bal}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{aceDmomFct}{\code{signature(theta = "numeric", object = "causalData")}: ... }
    \item{aceMomFct}{\code{signature(theta = "numeric", object = "causalData")}: ... }
	 }
}

\examples{
showClass("causalData")
}
\keyword{classes}
