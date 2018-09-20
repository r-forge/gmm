\name{nonlinearGel-class}
\docType{class}
\alias{nonlinearGel-class}

\title{Class \code{"nonlinearGel"}}
\description{
GEL model class for nonlinear regression models and instruments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("nonlinearGel", ...)}.
It is however, recommended to use the constructor \code{\link{gelModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }    
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"expression"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"expressionORNULL"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }    
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{nonlinearGmm}"}, directly.
Class \code{"\linkS4class{gelModels}"}, directly.
Class \code{"regGmm"}, by class \code{"nonlinearGmm"}, distance 2.
Class \code{"allNLGmm"}, by class \code{"nonlinearGmm"}, distance 2.
Class \code{"gmmModels"}, by class \code{"nonlinearGmm"}, distance 2
}

\examples{
showClass("nonlinearGmm")
}
\keyword{classes}
