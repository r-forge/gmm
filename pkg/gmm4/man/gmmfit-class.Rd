\name{gmmfit-class}
\docType{class}
\alias{gmmfit-class}

\title{Class \code{"gmmfit"}}
\description{
The main object that contains all information about the GMM fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gmmfit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{convergence}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{convIter}:}{Object of class \code{"numericORNULL"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{wObj}:}{Object of class \code{"gmmWeights"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{efficientGmm}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model}:}{Object of class \code{"gmmModels"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{bread}{\code{signature(x = "gmmfit")}: ... }
    \item{coef}{\code{signature(object = "gmmfit")}: ... }
    \item{meatGmm}{\code{signature(object = "gmmfit")}: ... }
    \item{print}{\code{signature(x = "gmmfit")}: ... }
    \item{residuals}{\code{signature(object = "gmmfit")}: ... }
    \item{show}{\code{signature(object = "gmmfit")}: ... }
    \item{specTest}{\code{signature(object = "gmmfit")}: ... }
    \item{summary}{\code{signature(object = "gmmfit")}: ... }
    \item{update}{\code{signature(object = "gmmfit")}: ... }
    \item{vcov}{\code{signature(object = "gmmfit")}: ... }
	 }
}

\seealso{
\code{\link{gmmFit}}, \code{\link{tsls}}, and \code{\link{gmm4}} to fit
a model with GMM and then create an object of class \code{"gmmfit"}..  
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("gmmfit")
}
\keyword{classes}
