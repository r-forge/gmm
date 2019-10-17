\name{causalGel-class}
\docType{class}
\alias{causalGel-class}
\alias{causalGel}

\title{Class \code{"causalGel"}}
\description{
This is the basic class for causality models. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("causalGel", ...)}.
It is however, recommended to use the constructor \code{aceModel}.
}
\section{Slots}{
  \describe{
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
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
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{functionGel}"}, directly.
Class \code{"\linkS4class{functionGmm}"}, by class "functionGel", distance 2.
Class \code{"\linkS4class{gelModels}"}, by class "functionGel", distance 2.
Class \code{"\linkS4class{allNLGmm}"}, by class "functionGel", distance 3.
Class \code{"\linkS4class{gmmModels}"}, by class "functionGel", distance 3.
}
\section{Methods}{
No methods defined with class "causalGel" in the signature.
}

\examples{
showClass("causalGel")
}
\keyword{classes}
