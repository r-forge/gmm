\name{smoothGel}

\alias{smoothGel}
	
\title{A kernel smoothing utility for \code{"gmmModels"} classes}

\description{
It either generates the optimal bandwidth and kernel weights when the
object is a GMM model, or the smoothed moments when the object is a GEL
model. 
}
\usage{
smoothGel(object, theta=NULL) 
}
\arguments{
  \item{object}{An object of class \code{"gmmModels"} or
    \code{"gelModels"}.}

  \item{theta}{An optional vector of coefficients. For
  \code{"gmmModels"}, it is the coefficient used to obtain the optimal
  bandwidth. If \code{NULL}, it is obtained using \code{\link{modelFit}}
  with the identity matrix as weights. For \code{"gelModels"}, the
  coefficient is require since the function returns the smoothed moments
  at a given vector of coefficients.}
}

\value{
  A list which contains:

 \item{k}{\eqn{2\times 1} vector of scaling factors used for GEL
   asymptotics. See Anatolyev (2005).}

 \item{w}{The kernel weights as an object of class "tskernel". See
   \code{\link{kernapply}}.}
 
 \item{bw}{A numeric bandwidth.}

 \item{kernel}{A character specifying th type of kernel used for smoothing}

 \item{smoothx}{Only when the object is of class \code{"gelModels"}, a
   matrix of smoothed moments}
}

\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Smith, R.J. (2011), GEL Criteria for Moment Condition Models.
\emph{Econometric Theory}, \bold{27}(6), 1192--1235.
}

\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearGmm
model1 <- gmmModel(y~x1, ~z1+z2, data=simData, kernel="Bartlett")

### get the bandwidth
### Notice that the kernel name is the not the same
### That's because a Truncated kernel for smoothing
### lead to a Bartlett kernel for the HAC of the moments
### See Smith (2011)
smoothGel(model1)


### The GEL model contains the info when it is created

model2 <- gmmToGel(model1, "EL")
model2@wSpec

smoothGel(model2, theta)$smoothx[1:5,]

}

