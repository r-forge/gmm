\name{bootGmmMet}
\alias{plot.bootGmm}
\alias{summary.bootGmm}
\alias{print.bootGmm}

\title{Methods for Bootstrap Gmm}
\description{Some methods to analyze the properties of the bootstrap estimates}
\usage{
\method{plot}{bootGmm}(x, which = 1, type = c("points", "density"), ...)
\method{summary}{bootGmm}(object, ...)
\method{print}{bootGmm}(x, ...)
}
\arguments{
 \item{object}{Object of class "bootGmm" returned by \link{bootGmm}}
 \item{x}{Object of class "bootGmm" returned by \link{bootGmm}}
 \item{which}{which coefficients to plot. Enter either a number or the name of the coefficient}
 \item{type}{Type of graph}
 \item{...}{Other arguments to pass to plot}
}

\value{
"summary" returns a summary statistics of the matrix of bootstrap estimates
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
res <- gmm(Y~X-1,~Z-1, data=d)
# Just resampling 25 time to save time
resB <- bootGmm(res, 25, seed = 123, niter = 1, trace = TRUE)
resB
plot(resB, 2)
plot(resB, 2, "density")
summary(resB)
}
