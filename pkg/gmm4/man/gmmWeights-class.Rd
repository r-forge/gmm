\name{gmmWeights-class}
\docType{class}
\alias{gmmWeights}
\alias{gmmWeights-class}

\title{Class \code{"gmmWeights"}}
\description{
The \code{"gmmWeights"} is an object used to compute different values
like the objective function, covariance matrix, etc. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gmmWeights", ...)}.
The method \code{\link{evalWeights}} is the constructor.
}
\section{Slots}{
  \describe{
    \item{\code{w}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{HAC}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{evalObjective}{\code{signature(object = "gmmModels", theta = "numeric", wObj = "gmmWeights")}: ... }
    \item{quadra}{\code{signature(x = "matrixORnumeric", w = "gmmWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "allNLGmm", wObj = "gmmWeights")}: ... }
    \item{solveGmm}{\code{signature(object = "linearGmm", wObj = "gmmWeights")}: ... }
	 }
}
\examples{
showClass("gmmWeights")

}
\keyword{classes}
