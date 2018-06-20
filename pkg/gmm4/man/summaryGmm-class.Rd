\name{summaryGmm-class}
\docType{class}
\alias{summaryGmm-class}

\title{Class \code{"summaryGmm"}}
\description{
Objects containing the detailed results of a GMM fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summaryGmm", ...)}.
Mostly created my the method \code{"summary"}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{specTest}:}{Object of class \code{"specTest"} ~~ }
    \item{\code{strength}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"gmmModels"} ~~ }
    \item{\code{sandwich}:}{Object of class \code{"logical"}}
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{df.adj}:}{Object of class \code{"logical"} ~~} 
    \item{\code{breadOnly}}{Object of class \code{"logical"} ~~}
    \item{\code{niter}}{The number of iterations}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "summaryGmm")}: ... }
    \item{show}{\code{signature(object = "summaryGmm")}: ... }
	 }
}

\examples{
showClass("summaryGmm")
}
\keyword{classes}
