\name{summaryGel-class}
\docType{class}
\alias{summaryGel-class}

\title{Class \code{"summaryGel"}}
\description{
Objects containing the detailed results of a GEL fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summaryGel", ...)}.
Mostly created my the method \code{"summary"}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{specTest}:}{Object of class \code{"specTest"} ~~ }
    \item{\code{model}:}{Object of class \code{"gmmModels"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lconvergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{impProb}:}{Object of class \code{"list"} ~~ }
  }
}

\examples{
showClass("summaryGel")
}
\keyword{classes}
