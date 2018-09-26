\name{nonlinearGmm-class}
\docType{class}
\alias{nonlinearGmm-class}

\title{Class \code{"nonlinearGmm"}}
\description{
GMM model class for nonlinear regression models and instruments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("nonlinearGmm", ...)}.
It is however, recommended to use the constructor \code{\link{gmmModel}}.
}
\section{Slots}{
  \describe{
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
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }        
  }
}
\section{Extends}{
Class \code{"\linkS4class{regGmm}"}, directly.
Class \code{"\linkS4class{allNLGmm}"}, directly.
}
\section{Methods}{
  \describe{
    \item{residuals}{\code{signature(object = "nonlinearGmm")}: ... }
	 }
}
\examples{
showClass("nonlinearGmm")
}
\keyword{classes}
