\name{sysGmmModel}

\alias{sysGmmModel}
	
\title{Constructor for \code{"sysGmmModels"} classes}

\description{
It builds the object of either class \code{"slinearGmm"} or
\code{"snonlinearGmm"}, which are system of equations to be estimated by
GMM.
}
\usage{
sysGmmModel(g, h, tet0=NULL, vcov = c("HAC", "MDS", "iid"),
            kernel = c("Quadratic Spectral",  "Truncated", "Bartlett",
            "Parzen", "Tukey-Hanning"), crit = 1e-06, bw = "Andrews",
            prewhite = 1L, ar.method = "ols", approx = "AR(1)",
            tol = 1e-07, centeredVcov = TRUE, data=parent.frame())
}
\arguments{
\item{g}{A list of linear or nonlinear regression formulas for each
  equation in the system.}

\item{h}{A list of linear formulas for the instruments in each
  equation in the system.}

\item{tet0}{A list of vectors of starting values. It is required
  only when the equations are nonlinear, in which case, it must be a
  list of named vector, with the names corresponding to the coefficient
  names in the regression formulas.}

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
'sysGmmModel' returns an object of one of the subclasses of
\code{"sysGmmModels"}.
 }

 \references{
 Hayashi, F. (2000). \emph{Econometrics}, New Jersey: Princeton
 University Press.
   
 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.
}
\examples{
set.seed(1122)
x1 <- rchisq(50,5)
x2 <- rchisq(50,5)
x3 <- rnorm(50)
x4 <- rnorm(50)
z1 <- .2*x1+rnorm(50)
z2 <- .2*x2+rnorm(50)
z3 <- rnorm(50)
z4 <- rnorm(50)
z5 <- rnorm(50)
y1 <- x1+rnorm(50)
y2 <- 2*x1+rnorm(50)
y3 <- 0.5*x2+rnorm(50)
dat <- data.frame(y1=y1,y3=y3,y2=y2, z1=z1,x1=x1,z2=z2,x2=x2,z3=z3,x3=x3,
                  x4=x4,z4=z4,z5=z5)

g1 <- y1~x1+x4; h1 <- ~z1+z2+z3+z4+x4
g2 <- y2~x1+x2+x3; h2 <- ~z1+z2+z3+z4+x3
g3 <- y3~x2+x3+x4; h3 <- ~z2+z3+z4+x3+x4
g <- list(g1,g2,g3)
h <- list(h1,h2,h3)

smodel <- sysGmmModel(g, h, data=dat)

## not really nonlinear
nlg <- list(y1~theta0+theta1*x1+theta2*x4,
            y2~alpha0+alpha1*x1+alpha2*x2+alpha3*x3,
            y3~beta0+beta1*x2+beta2*x3+beta3*x4)
tet0 <- list(c(theta0=1,theta1=2,theta2=3),
              c(alpha0=1,alpha1=2,alpha2=3, alpha3=4),
              c(beta0=1,beta1=2,beta2=3,beta3=4))
snmodel <- sysGmmModel(nlg, h, tet0, data=dat)


}

