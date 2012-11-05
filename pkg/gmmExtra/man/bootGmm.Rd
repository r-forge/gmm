\name{bootGmm}
\alias{bootGmm}

\title{Bootstrap Gmm}
\description{It produces bootstrap estimates using a given DGP}
\usage{
bootGmm(obj, N, seed = NULL, niter=8, trace=TRUE)
}
\arguments{
 \item{obj}{Object of class "gmm" returned by \link{gmm}}
 \item{N}{The number of bootstrap estimates}
 \item{seed}{Seed for resampling. By default we let R sets the seed}
 \item{niter}{Number of simultaneous estimates sent to \link{mclapply}}
 \item{trace}{Should we keep track of the number of bootstraps done and the time left?}
}

\value{
A list of objects of class gmm.
}

\references{
 Inoue, A. and Shintani M. (2006),
  Bootstrapping GMM estimators for time series.
  \emph{Journal of Econometrics}, \bold{133},
  531-555,
}

\examples{
# Two-stage-least-squares (2SLS), or IV with iid errors.
# The model is:
# Y(t) = b[0] + b[1]C(t) + b[2]Y(t-1) + e(t)
# e(t) is an MA(1)
# The instruments are Z(t)={1 C(t) y(t-2) y(t-3) y(t-4)}
     
getdat <- function(n) 
     {
     e <- arima.sim(n,model=list(ma=.9))
     C <- runif(n,0,5)
     Y <- rep(0,n)
     Y[1] = 1 + 2*C[1] + e[1]
     for (i in 2:n){
     Y[i] = 1 + 2*C[i] + 0.9*Y[i-1] + e[i]
     }
     Yt <- Y[5:n]
     X <- cbind(1,C[5:n],Y[4:(n-1)])
     Z <- cbind(1,C[5:n],Y[3:(n-2)],Y[2:(n-3)],Y[1:(n-4)]) 
     return(list(Y=Yt,X=X,Z=Z))
     }
     
set.seed(123)
d <- getdat(500)
res <- gmm(d$Y~d$X-1,~d$Z-1,vcov="iid")
# Just resampling 25 time to save time
resB <- bootGmm(res, 25, seed = 123, niter = 2, trace=TRUE)
}
