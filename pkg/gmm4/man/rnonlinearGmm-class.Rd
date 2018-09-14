\name{rnonlinearGmm-class}
\docType{class}
\alias{rnonlinearGmm-class}

\title{Class \code{"rnonlinearGmm"}}
\description{
It is a \code{nonlinearGmm} class with linear or nonlinear restrictions
imposed on the coefficients
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rnonlinearGmm",
  ...)}. But the \code{\link{restGmmModel}} method is the constructor.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
    \item{\code{kernel}:}{Object of class \code{"character"} ~~ }
    \item{\code{bw}:}{Object of class \code{"numericORcharacter"} ~~ }
    \item{\code{prewhite}:}{Object of class \code{"integer"} ~~ }
    \item{\code{ar.method}:}{Object of class \code{"character"} ~~ }
    \item{\code{approx}:}{Object of class \code{"character"} ~~ }
    \item{\code{tol}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{nonlinearGmm}"}, directly.
Class \code{"\linkS4class{regGmm}"}, by class "nonlinearGmm", distance 2.
Class \code{"\linkS4class{allNLGmm}"}, by class "nonlinearGmm", distance 2.
Class \code{"\linkS4class{gmmModels}"}, by class "nonlinearGmm", distance 2.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "rnonlinearGmm")}: ... }
    \item{Dresiduals}{\code{signature(object = "rnonlinearGmm")}: ... }
    \item{modelFit}{\code{signature(object = "rnonlinearGmm")}: ... }
    \item{modelDims}{\code{signature(object = "rnonlinearGmm")}: ... }
    \item{print}{\code{signature(x = "rnonlinearGmm")}: ... }
    \item{residuals}{\code{signature(object = "rnonlinearGmm")}: ... }
    \item{solveGmm}{\code{signature(object = "rnonlinearGmm", wObj = "gmmWeights")}: ... }
	 }
}
\examples{
showClass("rnonlinearGmm")
}
\keyword{classes}
