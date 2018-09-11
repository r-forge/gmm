\name{formulaGmm-class}
\docType{class}
\alias{formulaGmm-class}

\title{Class \code{"formulaGmm"}}
\description{
GMM model class for moment conditions defined by formulas. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("formulaGmm", ...)}.
It is however, recommended to use the constructor \code{\link{gmmModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{modelF}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{kernel}:}{Object of class \code{"character"} ~~ }
    \item{\code{bw}:}{Object of class \code{"numericORcharacter"} ~~ }
    \item{\code{prewhite}:}{Object of class \code{"integer"} ~~ }
    \item{\code{ar.method}:}{Object of class \code{"character"} ~~ }
    \item{\code{approx}:}{Object of class \code{"character"} ~~ }
    \item{\code{tol}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
    \item{\code{isMDE}:}{Object of class \code{"logical"} ~~ }
  }
}

\examples{
showClass("formulaGmm")
}
\keyword{classes}
