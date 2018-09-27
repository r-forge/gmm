\name{linearGmm-class}
\docType{class}
\alias{linearGmm-class}

\title{Class \code{"linearGmm"}}
\description{
GMM model class characterized by a linar regression and a set of instruments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("linearGmm", ...)}.
It is however, recommended to use the constructor \code{\link{gmmModel}}.
}
\section{Slots}{
  \describe{
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
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }    
  }
}
\section{Extends}{
Class \code{"\linkS4class{regGmm}"}, directly.
Class \code{"\linkS4class{gmmModels}"}, directly.
}
\section{Methods}{
  \describe{
    \item{residuals}{\code{signature(object = "linearGmm")}: ... }
	 }
}

\examples{
showClass("linearGmm")
}
\keyword{classes}
