\name{rlinearGmm-class}
\docType{class}
\alias{rlinearGmm-class}

\title{Class \code{"rlinearGmm"}}
\description{
It is a \code{linearGmm} class with linear restrictions imposed on the coefficients
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rlinearGmm", ...)}.
But the \code{\link{restGmmModel}} method is the constructor.
}
\section{Slots}{
  \describe{
    \item{\code{cstLHS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cstRHS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{linearGmm}"}, directly.
Class \code{"\linkS4class{regGmm}"}, by class "linearGmm", distance 2.
Class \code{"\linkS4class{gmmModels}"}, by class "linearGmm", distance 2.
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "rlinearGmm")}: ... }
	 }
}

\examples{
showClass("rlinearGmm")
}
\keyword{classes}
