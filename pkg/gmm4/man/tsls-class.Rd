\name{tsls-class}
\docType{class}
\alias{tsls-class}

\title{Class \code{"tsls"}}
\description{
An object that contains a two-stage least squares fit to a linear model.
}
\section{Objects from the Class}{
Objects can be created using the \code{\link{tsls}} method.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"gmmWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"gmmModels"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gmmfit}"}, directly.
}
\section{Methods}{
No methods defined with class "tsls" in the signature.
}
\examples{
showClass("tsls")
}
\keyword{classes}
