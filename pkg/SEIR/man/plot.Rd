\name{plot}
\alias{plot.seir}
\title{The SEIR model plot method}
\description{
 It is a plot method for \code{seir} objects.
}
\usage{
\method{plot}{seir}(x, which=c("acase","tcase", "S", "E", "I", "R",
                    "both","rzero"), add=FALSE, bothAdd=list(), ...)
}

\arguments{
  \item{x}{object of class \code{seir}.}
  \item{which}{A choice of lines to plot.}
  \item{add}{Should we add the lines to an existing plot?}
  \item{bothAdd}{A list if argument to pass to the active case plot when
  add is set to \code{TRUE}}
  \item{...}{Other arguments to pass to the plot method.}
}
