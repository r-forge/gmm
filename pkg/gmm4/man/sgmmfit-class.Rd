\name{sgmmfit-class}
\docType{class}
\alias{sgmmfit-class}

\title{Class \code{"sgmmfit"}}
\description{
The class for fitted system of equations estimated by GMM.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sgmmfit", ...)}.
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
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "sgmmfit")}: ... }
    \item{print}{\code{signature(x = "sgmmfit")}: ... }
    \item{residuals}{\code{signature(object = "sgmmfit")}: ... }
    \item{show}{\code{signature(object = "sgmmfit")}: ... }
	 }
}
\examples{
showClass("sgmmfit")
}
\keyword{classes}
