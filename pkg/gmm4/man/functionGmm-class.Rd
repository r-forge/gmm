\name{functionGmm-class}
\docType{class}
\alias{functionGmm-class}

\title{Class \code{"functionGmm"}}
\description{
GMM model class when moments matrix is defined by a function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("functionGmm", ...)}.
It is however, recommended to use the constructor \code{\link{gmmModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{X}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{fct}:}{Object of class \code{"function"} ~~ }
    \item{\code{dfct}:}{Object of class \code{"functionORNULL"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"numeric"} ~~ }
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
Class \code{"\linkS4class{allNLGmm}"}, directly.
Class \code{"\linkS4class{gmmModels}"}, directly.
}
\section{Methods}{
  \describe{
    \item{evalMoment}{\code{signature(object = "functionGmm")}: ... }
	 }
}

\examples{
showClass("functionGmm")
}
\keyword{classes}
