\name{specTest-class}
\docType{class}
\alias{specTest-class}

\title{Class \code{"specTest"}}
\description{
Object that contains information about specification tests.  
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("specTest", ...)}.
But it is mostly created by \code{\link{specTest}}
}

\section{Slots}{
  \describe{
    \item{\code{test}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{testname}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "specTest")}: ... }
    \item{show}{\code{signature(object = "specTest")}: ... }
	 }
}
\examples{
showClass("specTest")
}
\keyword{classes}
