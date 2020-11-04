\name{causalGelfit-class}
\docType{class}
\alias{causalGelfit-class}

\title{Class \code{"causalGelfit"}}
\description{
A class to store the fiited causality model obtained using the GEL method.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("causalGelfit",
  ...)}, but it mainly created using the \code{modelFit} method.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lconvergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"momentModel"} ~~ }
    \item{\code{restrictedLam}:}{Object of class \code{"integer"} ~~ }    
  }
}
\section{Extends}{
Class \code{"\linkS4class{gelfit}"}, directly.
}

\examples{
showClass("causalGelfit")
}
\keyword{classes}
