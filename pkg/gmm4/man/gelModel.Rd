\name{gelModel}

\alias{gelModel}
\alias{gmmToGel}
	
\title{Constructor for \code{"gelModels"} classes}

\description{
It builds the object of either class \code{"linearGel"},
\code{"nonlinearGel"}, \code{"functionGel"} or
\code{"formulaGel"}. This is the first step before running any
estimation algorithm.
}
\usage{
gelModel(g, x=NULL, gelType, rhoFct=NULL, tet0=NULL,grad=NULL,
        vcov = c("HAC", "MDS", "iid"),
        kernel = c("Quadratic Spectral",  "Truncated", "Bartlett", "Parzen",
                   "Tukey-Hanning"), crit = 1e-06,
        bw = "Andrews", prewhite = 1L, ar.method = "ols", approx = "AR(1)", 
        tol = 1e-07, centeredVcov = TRUE, data=parent.frame())

gmmToGel(object, gelType, rhoFct=NULL)
}
\arguments{
\item{object}{An object of class \code{"gmmModels"}. The function
  convert the object to a \code{"gelModels"} class}
  
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a
  \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for
  \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build
  the q sample moment conditions. It can also be a formula if the model
  is linear (see detailsbelow).}

\item{x}{The matrix or vector of data from which the function
  \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n
    \times Nh} matrix of instruments or a formula (see details below).}

\item{gelType}{A character string specifying the type of GEL. The
  available types are \code{"EL"}, \code{"ET"}, \code{"EEL"} and
  \code{"HD"}.} 

\item{rhoFct}{An optional function that return \eqn{\rho(v)}. This is
  for users who want a GEL model that is not built in the package. The
  four arguments of the function must be \code{"gmat"}, the matrix of
  moments, \code{"lambda"}, the vector of Lagrange multipliers,
  \code{"derive"}, which specify the order of derivative to return, and
  \code{k} a numeric scale factor required for time series and kernel
  smoothed moments.}

\item{tet0}{A \eqn{k \times 1} vector of starting values. It is required
  only when "g" is a function because only then a numerical algorithm is
  used to minimize the objective function. If the dimension of
  \eqn{\theta} is one, see the argument "optfct".}

\item{grad}{A function of the form \eqn{G(\theta,x)} which returns a
  \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with
  respect to \eqn{\theta}. By default, the numerical algorithm
  \code{numericDeriv} is used. It is of course strongly suggested to
  provide this function when it is possible. This gradient is used to
  compute the asymptotic covariance matrix of \eqn{\hat{\theta}} and to
  obtain the analytical gradient of the objective function if the method
  is set to "CG" or "BFGS" in \code{\link{optim}} and if "type" is not
  set to "cue". If "g" is a formula, the gradiant is not required (see
  the details below).}

\item{vcov}{Assumption on the properties of the moment conditions. By
  default, they are weakly dependant processes. For \code{MDS}, we
  assume that the conditions are martingale difference sequences, which
  implies they are serially uncorrelated, but may be
  heteroscedastic. There is a difference between \code{iid} and
  \code{MDS} only when \code{g} is a formula. In that case, residuals
  are assumed homoscedastic as well as serially uncorrelated.}

\item{kernel}{type of kernel used to compute the covariance matrix of
  the vector of sample moment conditions (see \code{\link{kernHAC}} for
  more details)}

\item{crit}{The stopping rule for the iterative GMM. It can be reduce to
  increase the precision.}

\item{bw}{The method to compute the bandwidth parameter in the HAC
  weighting matrix. The default is \code{link{bwAndrews}} (as proposed in Andrews
  (1991)), which minimizes the MSE of the weighting matrix. Alternatives
  are \code{link{bwWilhelm}} (as proposed in Wilhelm
  (2015)), which minimizes the mean-square error (MSE) of the resulting
  GMM estimator, and \code{link{bwNeweyWest}} (as proposed in Newey-West(1994)).}

\item{prewhite}{logical or integer. Should the estimating functions be
  prewhitened? If \code{TRUE} or greater than 0 a VAR model of order
  \code{as.integer(prewhite)} is fitted via \code{ar} with method
  \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to
  \code{\link{ar}} for prewhitening.}

\item{approx}{A character specifying the approximation method if the
  bandwidth has to be chosen by \code{bwAndrews}.}

\item{tol}{Weights that exceed \code{tol} are used for computing the
  covariance matrix, all other weights are treated as 0.}

\item{centeredVcov}{Should the moment function be centered when
  computing its covariance matrix. Doing so may improve inference.}

\item{data}{A data.frame or a matrix with column names (Optional). }
}

\value{
'gmmModel' returns an object of one of the subclasses of \code{"gmmModels"}.
 }

\references{
Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Kitamura, Y. and Otsu, T. and Evdokimov, K. (2013), Robustness,
Infinitesimal Neighborhoods and Moment Restrictions.
\emph{Econometrica}, \bold{81}, 1185-1201.

Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2011), GEL Criteria for Moment Condition Models.
\emph{Econometric Theory}, \bold{27}(6), 1192--1235.
}

\examples{
data(simData)
theta <- c(beta0=1,beta1=2)

## A linearGmm
model1 <- gelModel(y~x1, ~z1+z2, gelType="EL", data=simData)
model1

## A nonlinearGmm
g <- y~beta0+x1^beta1
h <- ~z1+z2
model2 <- gelModel(g, h, tet0=c(beta0=1, beta1=2), data=simData, gelType="ET")
model2

}

