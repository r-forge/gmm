\name{rlinearGel-class}
\docType{class}
\alias{rlinearGel-class}

\title{Class \code{"rlinearGel"}}
\description{
Restricted GEL model class characterized by a linear regression and a set of instruments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rlinearGel", ...)}.
It is however, recommended to use the \code{\link{restModel} method}.
}
\section{Slots}{
  \describe{
    \item{\code{cstLHS}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cstRHS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }        
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }
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
Class \code{"\linkS4class{linearGel}"}, directly.
Class \code{"\linkS4class{rgelModels}"}, directly.
Class \code{"linearGmm"}, by class \code{"linearGel"}, distance 2.
Class \code{"gelModels"}, by class \code{"linearGel"}, distance 2.
Class \code{"regGmm"}, by class \code{"linearGel"}, distance 3.
Class \code{"gmmModels"}, by class \code{"linearGel"}, distance 3.
}

\examples{
showClass("rlinearGel")
}
\keyword{classes}
