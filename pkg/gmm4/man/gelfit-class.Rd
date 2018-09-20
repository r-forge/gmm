\name{gelfit-class}
\docType{class}
\alias{gelfit-class}

\title{Class \code{"gelfit"}}
\description{
The main object that contains all information about the GEL fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gelfit", ...)}
but can be created with \code{\link{modelFit}}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lconvergence}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"gelModels"} ~~ }
  }
}

\seealso{
\code{\link{modelFit}}
}
\examples{
showClass("gelfit")
}
\keyword{classes}
