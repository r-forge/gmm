\name{SEIR}
\alias{SEIR}
\alias{solveSEIR}

\title{The SEIR model for virus spread}
\description{
The function that solves the SEIR system of differential equations.
}
\usage{

SEIR(t, y, y0, c=0, sigma=1/5.2, gamma=1/18,
     r0=matrix(c(0,20,70,84,90,3,2.6,1.9,1,.5),ncol=2),
     type=c("Lin", "Const"))    

solveSEIR(h=1e-2, T=180, c=0, y0=c(11e6, 40, 800, 0),
          sigma=1/5.2, gamma=1/18,
          r0=matrix(c(0,20,70,84,90,3,2.6,1.9,1,.5),ncol=2),
          type=c("Lin", "Const"))
}
\arguments{

  \item{t}{Time value to evaluate the derivatives.}
  
  \item{y}{The vector is variable of interest: S, E, I and R.}

  \item{y0}{Initial Values}

  \item{sigma}{Parameter.}

  \item{gamma}{Parameter.}

  \item{c}{Mutation Parameter}
  
  \item{h}{Step size for the RK4 method}
  
  \item{T}{Final t}

  \item{r0}{Matrix of breaking points for R-Zero. The first column is
  the time of the breaking points, and the second is the values of the
  R-zero.}

  \item{type}{Types of R-Zero function. The default is a linear
    interpolation and the second is a piecewise constant function.}
}

\value{
It returns an object of class \code{"seir"}  
}
\references{
Wang, H., Wang, Z., Dong, Y. et al (2020),  Phase-adjusted estimation of
the number of Coronavirus Disease 2019 cases in Wuhan, China.
\emph{Cell Discov}, \bold{6}, \bold{10}.
}


