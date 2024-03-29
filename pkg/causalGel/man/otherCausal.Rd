\name{otherCausal}
\alias{otherCausal}
\alias{matching}
\alias{LLmatching}
\alias{ipw}
	
\title{Popular Causal Effect estimation methods}

\description{ This documentation file presents a collection of popular
methods used to estimate the average causal effect, the causal effect on
the treated and the causal effect on the non-treated.
} 

\usage{

matching(form,  balm, data, type=c("ACE","ACT","ACC"), M=4,
         psForm=NULL, bcForm=NULL, vJ=4)

LLmatching(form, psForm, data, type=c("ACE","ACT","ACC"),
           kern=c("Gaussian","Epanechnikov"),tol=1e-4,
           h=NULL, from=.00001, to=5, ngrid=10, maxit=100,
           hMethod=c("Brent","Grid"))

ipw(form, psForm, data, type=c("ACE","ACT","ACC"),
    normalized=FALSE, ...)
}
\arguments{

  \item{form}{A formula that links the outcome variable to the treatment
    indicator. For the moment, only one treatment group is allowed.}
  
  \item{balm}{A formula or a matrix with balancing covariates to be
    matched.}

  \item{data}{A data.frame or a matrix with column names.}

  \item{type}{The type of causal effect to compute. \code{ACE} stands
    for average causal effect, \code{ACT} for causal effect on the treated
    and \code{ACC} for causal effect on the control or non-treated.}

  \item{M}{The minimum number of matches}
  
  \item{psForm}{It is the \code{formula} argument passed to
    \code{\link{glm}} to estimate the propensity scores by a logistic
    regression.}
  
  \item{bcForm}{A formula that represents the right hand side in the
    regression used for the bias correction.}

  \item{kern}{The type of kernel to use in the local linear regression
    method.}
    
  \item{tol}{The tolerance level for the stopping rule used to compute
    the optimal bandwidth.}

  \item{h}{A fixed bandwidth. By default, the optimal bandwidth is found
    by minimizing a cross-validation.}
  
  \item{from}{The lower bound for the search of the optimal bandwidth.}

  \item{to}{The upper bound for the search of the optimal bandwidth.}

  \item{ngrid}{The number of grid points if the optimal bandwidth is
    obtained by grid search.}

  \item{maxit}{The maximum number of iterations for the minimization of
    the cross-validation.}

  \item{hMethod}{The method used to find the optimal bandwidth.} 

  \item{normalized}{Should the weights be normalized. If set to
    \code{GPE}, the GPE method is used.}

  \item{vJ}{The minimum number of matches for the standard error estimation}
  
  \item{...}{Additional arguments to pass to \code{\link{optim}} when
    the \code{GPE} method is used.}
	    
}

\value{
  All methods return an object of classes \code{"causalfit"}.
 }

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
ps <- treat~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
           
fit1 <- matching(g, balm, nsw, "ACE")
fit1

fit2 <- LLmatching(g, ps, nsw, "ACE")
fit2

fit3 <- ipw(g, ps, nsw, "ACE", normalized=TRUE)
fit3

}

