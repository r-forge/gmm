\name{rslinearGmm-class}
\docType{class}
\alias{rslinearGmm-class}

\title{Class \code{"rslinearGmm"}}
\description{
A class for systems of linear equations with restrictions on the coefficients.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rslinearGmm", ...)}.
The constructor is \code{\link{restGmmModel}}.
}
\section{Slots}{
  \describe{
    \item{\code{cstLHS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cstRHS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{modelT}:}{Object of class \code{"list"} ~~ }
    \item{\code{instT}:}{Object of class \code{"list"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"character"} ~~ }
    \item{\code{n}:}{Object of class \code{"integer"} ~~ }
    \item{\code{q}:}{Object of class \code{"integer"} ~~ }
    \item{\code{k}:}{Object of class \code{"integer"} ~~ }
    \item{\code{parNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{momNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqnNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }    
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sameMom}:}{Object of class \code{"logical"} ~~ }
    \item{\code{SUR}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"list"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{slinearGmm}"}, directly.
Class \code{"\linkS4class{rsysGmmModels}"}, directly.
Class \code{"\linkS4class{sysGmmModels}"}, by class "slinearGmm", distance 2.
}
\section{Methods}{
  \describe{
    \item{getRestrict}{\code{signature(object = "rslinearGmm")}: ... }
    \item{print}{\code{signature(x = "rslinearGmm")}: ... }
    \item{printRestrict}{\code{signature(object = "rslinearGmm")}: ... }
	 }
}
\examples{
showClass("rslinearGmm")
}
\keyword{classes}
