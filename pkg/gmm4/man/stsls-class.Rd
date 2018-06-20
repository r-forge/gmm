\name{stsls-class}
\Rdversion{1.1}
\docType{class}
\alias{stsls-class}

\title{Class \code{"stsls"}}
\description{
The class for fitted system of equations by equation by equation
two-stage least squares.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stsls", ...)}.
The method \code{\link{tsls}} generates it.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"list"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"sysGmmWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"sysGmmModels"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sgmmfit}"}, directly.
}
\section{Methods}{
No methods defined with class "stsls" in the signature.
}
\examples{
showClass("stsls")
}
\keyword{classes}
