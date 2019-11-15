\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,mconfint-method}
\title{ ~~ Methods for Function \code{plot} from package \pkg{graphics} ~~}
\description{
It plots the confidence region. 
}
\usage{

\S4method{plot}{ANY}(x, y, ...) 

\S4method{plot}{mconfint}(x, y, main=NULL, xlab=NULL, ylab=NULL, 
                          pch=21, bg=1, \dots)
}

\arguments{
  \item{x}{An object to plot}

  \item{y}{On used for \code{"ANY"}.}

  \item{main}{Optional title}

  \item{xlab}{Optional label for the x-axis.}

  \item{ylab}{Optional label for the y-axis.}

  \item{pch}{Type of points (see \code{\link{points}}).}

  \item{bg}{Background color for points.}

  \item{\dots}{Arguments to pass to \code{\link{polygon}}}
  }

\section{Methods}{
\describe{

  \item{\code{signature(object = "ANY")}}{
    It uses the \code{plot} from package \pkg{graphics}
  }

  \item{\code{signature(object = "mconfint")}}{
    Plot the 2D confidence region.
  }
}}

\keyword{methods}
