\name{gmmToGel-methods}
\docType{methods}
\alias{gmmToGel}
\alias{gmmToGel-methods}
\alias{gmmToGel,gmmModels-method}
\alias{gmmToGel,gelModels-method}
\alias{gmmToGel,rgmmModels-method}
\title{ ~~ Methods for Function \code{gmmToGel} in Package \pkg{gmm4} ~~}
\description{
It converts \code{"gmmModels"} classes to \code{"gelModels"} class by
simply adding the specification about the type of GEL.
}

\usage{
\S4method{gmmToGel}{gmmModels}(object, gelType, rhoFct=NULL)

\S4method{gmmToGel}{rgmmModels}(object, gelType, rhoFct=NULL)
}

\arguments{
  \item{object}{An object of class \code{"gmmModels"}}
  \item{gelType}{The type pf GEL model}
  \item{rhoFct}{An alternative GEL function if the ones provided by the
  package do not satisfy}
}

\section{Methods}{
\describe{

\item{\code{signature(object = "gmmModels")}}{
The object it converted to \code{"gelModels"}.
}

\item{\code{signature(object = "gelModels")}}{  
It applies the transformation to the \code{"gmmModels"} contained in the
\code{"gelModels"} object.
}

\item{\code{signature(object = "rgmmModels")}}{
The object it converted to \code{"rgelModels"}.
}
}}
\keyword{methods}

