\name{confint-class}
\docType{class}
\alias{confint-class}

\title{Class \code{"confint"}}
\description{
A class for confidence interval results.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("confint", ...)},
but it is recommended to use the \code{"confit"} method.

}
\section{Slots}{
  \describe{
    \item{\code{interval}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{level}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "confint")}: ... }
    \item{show}{\code{signature(object = "confint")}: ... }
	 }
}

\examples{
showClass("confint")
}
\keyword{classes}
