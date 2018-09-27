\name{functionGel-class}
\docType{class}
\alias{functionGel-class}

\title{Class \code{"functionGel"}}
\description{
GEL model class when moments matrix is defined by a function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("functionGel", ...)}.
It is however, recommended to use the constructor \code{\link{gelModel}}.
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
Class \code{"\linkS4class{functionGmm}"}, directly.
Class \code{"\linkS4class{gelModels}"}, directly.
Class \code{"allNLGmm"}, by class \code{"functionGmm"}, distance 2.
Class \code{"gmmModels"}, by class \code{"functionGmm"}, distance 2
}

\examples{
showClass("functionGmm")
}
\keyword{classes}
