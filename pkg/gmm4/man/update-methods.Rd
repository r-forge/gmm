\name{update-methods}
\docType{methods}
\alias{update-methods}
\alias{update,ANY-method}
\alias{update,gmmfit-method}
\title{ ~~ Methods for Function \code{update} in Package \pkg{stats} ~~}
\description{
The method is used to refit the model by GMM with either a different
method or with modifications to the \code{gmmModels}
}
\usage{
\S4method{update}{gmmfit}(object, \dots)
}
\arguments{
  \item{object}{GMM fit object}
  \item{\dots}{Arguments to modify the model or the GMM method}
  }
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
}

\item{\code{signature(object = "gmmfit")}}{
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

(res <- gmmFit(model1))

## lets change to iterative
update(res, type="iter")

## Let change the HAC specification in the model1 object
## to MDS
update(res, vcov="MDS")


}
\keyword{methods}
