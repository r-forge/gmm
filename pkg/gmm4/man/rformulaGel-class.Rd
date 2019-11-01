\name{rformulaGel-class}
\docType{class}
\alias{rformulaGel-class}

\title{Class \code{"rformulaGel"}}
\description{
Restricted GEL model class for moment conditions defined by formulas. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rformulaGel", ...)}.
It is however, recommended to use the \code{\link{restModel}} method.
}
\section{Slots}{
  \describe{
    \item{\code{R}:}{Object of class \code{"list"} ~~ }
    \item{\code{cstSpec}:}{Object of class \code{"list"} ~~ }            
    \item{\code{wSpec}:}{Object of class \code{"list"} ~~ }
    \item{\code{gelType}:}{Object of class \code{"list"} ~~ }    	
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
    \item{\code{vcovOptions}:}{Object of class \code{"list"} ~~ }    
    \item{\code{centeredVcov}:}{Object of class \code{"logical"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{isEndo}:}{Object of class \code{"logical"} ~~ }
    \item{\code{isMDE}:}{Object of class \code{"logical"} ~~ }
    \item{\code{omit}:}{Object of class \code{"integer"} ~~ }
    \item{\code{survOptions}:}{Object of class \code{"list"} ~~ }        
  }
}

\section{Extends}{
Class \code{"\linkS4class{formulaGel}"}, directly.
Class \code{"\linkS4class{rgelModels}"}, directly.
Class \code{"formulaGmm"}, by class \code{"formulaGel"}, distance 2.
Class \code{"gelModels"}, by class \code{"formulaGel"}, distance 2.
Class \code{"allNLGmm"}, by class \code{"formulaGel"}, distance 3.
Class \code{"gmmModels"}, by class \code{"formulaGel"}, distance 3.
}


\examples{
showClass("rformulaGel")
}
\keyword{classes}
