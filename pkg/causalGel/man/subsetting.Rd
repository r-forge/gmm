\name{[-causalGel}
\docType{methods}
\alias{[,causalGel,missing,numeric-method}
\alias{[,causalGel,numeric,missing-method}
\alias{[,causalGel,numeric,numeric-method}

\title{Subsetting methods}
\description{
Different subsetting methods for S4 class objects of the package. The
subset method returns an new object with observations selected by the
second argument. See example.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "causalGel", i = "numeric", j = "missing")}}{
  Selects observations
}

\item{\code{signature(x = "causalGel", i = "missing", j = "numeric")}}{
  Selects balancing moments
}

\item{\code{signature(x = "causalGel", i = "numeric", j = "numeric")}}{
  \code{i} selects the observations and "j" selects the balancing
  moments.
}
}}

\keyword{methods}
\keyword{subsetting}
