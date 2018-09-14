\name{solveGel-methods}
\docType{methods}
\alias{solveGel}
\alias{solveGel-methods}
\alias{solveGel,gelModels-method}
\title{ ~~ Methods for Function \code{solveGel} in Package \pkg{gmm4} ~~}
\description{
 ~~ Methods for function \code{solveGel} in package \pkg{gmm4} ~~
}
\usage{
\S4method{solveGel}{gelModels}(object, theta0=NULL, lambda0=NULL,
                               lamSlv=NULL,
                               coefSlv=c("optim","nlminb","constrOptim"),
                               lControl=list(), tControl=list())
}
\arguments{
  \item{object}{An object of class \code{"gelModels"}}
  \item{theta0}{The vector of coefficients for the starting values used
    in minimization algorithm. If NULL, the starting values in the object is
    used. For linear models, it must be provided because
    \code{"linearGel"} objects do not have a theta0 slot.}
  \item{lambda0}{The \eqn{q \times 1} vector of starting values for the
    Lagrange multipliers. By default a zero vector is used.}
  \item{lamSlv}{An alternative solver for the Lagrange multiplier. By
    default, either \code{\link{EL.Wu}} or \code{\link{getLambda}} is
    used.}
  \item{coefSlv}{Minimization solver for the coefficient vector.}
  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}
  \item{tControl}{A list of controls for the coefficient algorithm.}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "gelModels")}}{
The method applies to all GEL classes.
}

}}
\value{
  A list with the following:
  \item{theta}{The vector of solution}
  \item{lambda}{The vector of Lagrange multiplier}
  \item{lconvergence}{convergence code for the Lagrange multiplier. 0
    means normal convergence.}
  \item{convergence}{convergence code for the coefficients. 0 means
    normal convergence. For higher numbers, see \code{\link{optim}},
    \code{\link{constrOptim}} or  \code{\link{nlminb}}}
  }
\examples{
data(simData)
theta <- c(beta0=1,beta1=2)
model1 <- gelModel(y~x1, ~z1+z2, gelType="EL", data=simData)

}

\keyword{methods}

