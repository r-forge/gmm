\name{tsls-methods}
\docType{methods}
\alias{tsls}
\alias{tsls-methods}
\alias{tsls,linearGmm-method}
\alias{tsls,slinearGmm-method}
\title{ ~~ Methods for Function \code{tsls} in Package \pkg{gmm4} ~~}
\description{
It estimates a linear model using two-stage least squares.
}
\usage{
\S4method{tsls}{linearGmm}(object)
}
\arguments{
  \item{object}{An object of class \code{gmmModels}}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "linearGmm")}}{
}

\item{\code{signature(object = "slinearGmm")}}{
  2SLS for equation by equation estimation of a system of equations.
}
}}
\examples{
x <- rchisq(200,5)
z1 <- rnorm(200)
z2 <- .2*x+rnorm(200)
y <- x+rnorm(200)
dat <- data.frame(y=y,z1=z1,x=x,z2=z2)
theta <- c(beta0=1,beta1=2)
model1 <- gmmModel(y~x, ~z1+z2, data=dat)
res <- tsls(model1)
summary(res)

## Econometrics, Fumio Hayashi (2000)
## Empirical exercises (b) and (c)
data(Griliches)
Griliches$YEAR <- as.factor(Griliches$YEAR)
model1  <- gmmModel(LW~S+IQ+EXPR+TENURE+RNS+SMSA+YEAR-1,
                   ~S+EXPR+TENURE+RNS+SMSA+YEAR+MED+KWW+MRT+AGE-1,
                   data=Griliches, vcov="MDS")
res <- tsls(model1)
summary(res)
}
\keyword{methods}

