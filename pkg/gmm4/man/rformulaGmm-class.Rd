\name{rformulaGmm-class}
\docType{class}
\alias{rformulaGmm-class}

\title{Class \code{"rformulaGmm"}}
\description{
Class for GMM models in which the moment conditons are defined through a
formula and the coefficients are restricted.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rformulaGmm", ...)}.
The constructor is \code{\link{restGmmModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
\section{Extends}{
Class \code{"\linkS4class{functionGmm}"}, directly.
Class \code{"\linkS4class{rgmmModels}"}, directly.
Class \code{"\linkS4class{allNLGmm}"}, by class "functionGmm", distance 2.
Class \code{"\linkS4class{gmmModels}"}, by class "functionGmm", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "rfunctionGmm", i = "integer", j = "ANY")}: ... }
    \item{coef}{\code{signature(object = "rfunctionGmm")}: ... }
    \item{getRestrict}{\code{signature(object = "rfunctionGmm")}: ... }
    \item{modelDims}{\code{signature(object = "rfunctionGmm")}: ... }
    \item{print}{\code{signature(x = "rfunctionGmm")}: ... }
    \item{printRestrict}{\code{signature(object = "rfunctionGmm")}: ... }
	 }
}

\examples{
showClass("rfunctionGmm")
}
\keyword{classes}
