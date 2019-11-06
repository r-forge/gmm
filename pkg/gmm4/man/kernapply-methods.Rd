\name{kernapply-methods}
\docType{methods}
\alias{kernapply}
\alias{kernapply-methods}
\alias{kernapply,gelModels-method}
\alias{kernapply,gmmModels-method}

\title{A kernel smoothing utility for \code{"gelModels"} classes}

\description{
It either generates the optimal bandwidth and kernel
weights or the smoothed moments of GEL models.
}

\usage{
\S4method{kernapply}{gelModels}(x, theta=NULL, smooth=TRUE, \dots)

\S4method{kernapply}{gmmModels}(x, theta=NULL, \dots)
}

\arguments{
  \item{x}{An object of class \code{"gelModels"}.}
  
  \item{theta}{An optional vector of coefficients. For
    \code{smooth=FALSE}, it is used to obtain the optimal bandwidth. If
    \code{NULL}, the bandwidth is obtained using one step GMM with the
    identity matrix as weights. For \code{smooth=TRUE}, the coefficient is
    required since the function returns the smoothed moments at a given
    vector of coefficients.}
  \item{smooth}{By default, it returns the smoothed moment matrix. If
    \code{FALSE}, it computes the optimal bandwidth and kernel weights.}
  \item{\dots}{Other arguments to pass. Currently not used}
}
  
\value{
  A list which contains:

 \item{k}{\eqn{2\times 1} vector of scaling factors used for GEL
   asymptotics. See Anatolyev (2005).}

 \item{w}{The kernel weights as an object of class "tskernel". See
   \code{\link{kernel}}.}
 
 \item{bw}{A numeric bandwidth.}

 \item{kernel}{A character specifying th type of kernel used for smoothing}

 \item{smoothx}{Only when \code{smooth=TRUE}, a
   matrix of smoothed moments}
}

\references{ Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and
Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Smith, R.J. (2011), GEL Criteria for Moment Condition Models.
\emph{Econometric Theory}, \bold{27}(6), 1192--1235.
}

\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearGmm
model1 <- gmmModel(y~x1, ~z1+z2, data=simData,vcov="HAC",vcovOptions=list(kernel="Bartlett"))

### get the bandwidth
### Notice that the kernel name is the not the same
### That's because a Truncated kernel for smoothing
### lead to a Bartlett kernel for the HAC of the moments
### See Smith (2011)
kernapply(model1, smooth=FALSE)


### The GEL model contains the info when it is created

model2 <- gmmToGel(model1, "EL")
model2@wSpec

kernapply(model2, theta)$smoothx[1:5,]

}

