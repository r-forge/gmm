\name{summaryCausalfit-class}
\docType{class}
\alias{summaryCausalfit-class}

\title{Class \code{"summaryCausalfit"}}
\description{
  This is the class for the summary of causal effect estimates other
  than GEL.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("summaryCausalfit", ...)}. But, it is created by
  \code{\link{summary-methods}}.
}
\section{Slots}{
  \describe{
    \item{\code{coef}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{form}:}{Object of class \code{"list"} ~~ }
    \item{\code{details}:}{Object of class \code{"list"} ~~ }
    \item{\code{info}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "summaryCausalfit" in the signature.
}

\examples{
showClass("summaryCausalfit")
}
\keyword{classes}
