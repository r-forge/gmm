\name{evalDMoment-methods}
\docType{methods}
\alias{evalDMoment}
\alias{evalDMoment-methods}
\alias{evalDMoment,functionGmm-method}
\alias{evalDMoment,sysGmmModels-method}
\alias{evalDMoment,rslinearGmm-method}
\alias{evalDMoment,regGmm-method}
\title{ ~~ Methods for Function \code{evalDMoment} in Package \pkg{Gmm} ~~}
\description{
It computes the matrix of derivatives of the sample moments with respect
to the coefficients.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "functionGmm")}}{
}

\item{\code{signature(object = "regGmm")}}{
}

\item{\code{signature(object = "sysGmmModels")}}{
}

\item{\code{signature(object = "rslinearGmm")}}{
}
}}
\examples{
data(simData)
theta <- c(1,1)
model1 <- gmmModel(y~x1, ~z1+z2, data=simData)
G <- evalDMoment(model1, theta)

## A nonlinearGmm
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- gmmModel(g, h, c(beta0=1, beta1=2), data=simData)
G <- evalDMoment(model2, c(beta0=1, beta1=2))

## A functionGmm
fct <- function(tet, x)
    {
        m1 <- (tet[1] - x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3 - tet[1]*(tet[1]^2 + 3*tet[2]^2)
        f <- cbind(m1, m2, m3)
        return(f)
    }
dfct <- function(tet, x)
        {
        jacobian <- matrix(c( 1, 2*(-tet[1]+mean(x)), -3*tet[1]^2-3*tet[2]^2,0, 2*tet[2],
			   -6*tet[1]*tet[2]), nrow=3,ncol=2)
        return(jacobian)
        }
X <- rnorm(200)
model3 <- gmmModel(fct, X, tet0=c(beta0=1, beta1=2), grad=dfct)
G <- evalDMoment(model3, c(beta0=1, beta1=2))
}

\keyword{methods}
