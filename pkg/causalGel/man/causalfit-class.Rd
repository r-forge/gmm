\name{causalfit-class}
\docType{class}
\alias{causalfit-class}

\title{Class \code{"causalfit"}}
\description{
This is the class for causal effect estimates other than GEL.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("causalfit", ...)}.
But, it is created by \code{\link{matching}}, \code{\link{LLmatching}},
or \code{\link{ipw}}.
}
\section{Slots}{
  \describe{
    \item{\code{estim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{se}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{coefNames}:}{Object of class \code{"character"} ~~ }        
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{form}:}{Object of class \code{"list"} ~~ }
    \item{\code{details}:}{Object of class \code{"list"} ~~ }
    \item{\code{info}:}{Object of class \code{"list"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{call}:}{Object of class \code{"callORNULL"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "causalfit" in the signature.
}

\examples{
showClass("causalfit")
}
\keyword{classes}
