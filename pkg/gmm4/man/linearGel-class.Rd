\name{linearGel-class}
\docType{class}
\alias{linearGel-class}

\title{Class \code{"linearGel"}}
\description{
GEL model class characterized by a linar regression and a set of instruments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("linearGel", ...)}.
It is however, recommended to use the constructor \code{\link{gelModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }    
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{linearGmm}"}, directly.
Class \code{"\linkS4class{gelModels}"}, directly.
Class \code{"regGmm"}, by class \code{"linearGmm"}, distance 2.
Class \code{"gmmModels"}, by class \code{"linearGmm"}, distance 2
}

\examples{
showClass("linearGmm")
}
\keyword{classes}
