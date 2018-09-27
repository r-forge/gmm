\name{rsnonlinearGmm-class}
\docType{class}
\alias{rsnonlinearGmm-class}

\title{Class \code{"rsnonlinearGmm"}}
\description{
A class for system of nonliear equations in which the coefficients are
restricted.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rsnonlinearGmm", ...)}.
The constructor is \code{\link{restGmmModel}}
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{instT}:}{Object of class \code{"list"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{theta0}:}{Object of class \code{"list"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{fRHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{fLHS}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqnNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }    
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sameMom}:}{Object of class \code{"logical"} ~~ }
    \item{\code{SUR}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"list"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }      	    
  }
}
\section{Extends}{
Class \code{"\linkS4class{snonlinearGmm}"}, directly.
Class \code{"\linkS4class{rsysGmmModels}"}, directly.
Class \code{"\linkS4class{sysGmmModels}"}, by class "snonlinearGmm", distance 2.
}
\section{Methods}{
No methods defined with class "rsnonlinearGmm" in the signature.
}
\examples{
showClass("rsnonlinearGmm")
}
\keyword{classes}
