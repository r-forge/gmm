\name{causalGEL}

\alias{causalGEL}
	
\title{Causal inference using balancing methods}

\description{
It fit a causality model using balancing based on GEL methods. It
creates  an object of class \code{"causalGelfit"}.
}
\usage{

causalGEL(g, balm, data, theta0=NULL,
          momType=c("ACE","ACT","ACC", "uncondBal","fixedMom"),
          popMom = NULL, rhoFct=NULL,ACTmom=1L, 
          gelType = c("EL", "ET", "EEL", "ETEL", "HD", "ETHD","REEL"),
          initTheta = c("gmm","theta0"), getVcov=FALSE,
          lambda0=NULL, 
          lamSlv=NULL, coefSlv= c("optim","nlminb","constrOptim"),
          lControl=list(), tControl=list())
}
\arguments{

  \item{g}{A formula that links the outcome variable to the treatment
    indicators.}

  \item{balm}{A formula or a matrix with balancing covariates}

  \item{data}{A data.frame or a matrix with column names.}

  \item{theta0}{A vector of starting values (optional). If not provided,
    the least squares method is use to generate them}

  \item{momType}{How the moments of the covariates should be
    balanced. By default, it is balanced using the sample mean of the
    covariates, which corresponds to the ACE. Alternatively, to the
    sample moments of the treated group (ACT), the control group (ACC),
    or to a known population mean. The option 'uncondBal' means that it
    is unconditionally balanced.}
        
  \item{popMom}{A vector of population moments to use for balancing. It
    can be used if those moments are available from a census, for
    example. When available, it greatly improves efficiency.}

  \item{rhoFct}{An optional function that return \eqn{\rho(v)}. This is
    for users who want a GEL model that is not built in the package. The
    four arguments of the function must be \code{"gmat"}, the matrix of
    moments, \code{"lambda"}, the vector of Lagrange multipliers,
    \code{"derive"}, which specify the order of derivative to return, and
    \code{k} a numeric scale factor required for time series and kernel
    smoothed moments.}

  \item{ACTmom}{When \code{momType} is set to 'ACT', that integer
    indicates which treated group to use to balance the covariates.} 
  
  \item{gelType}{"EL" for empirical likelihood, "ET" for exponential tilting,
    "EEL" for Euclidean empirical likelihood, "ETEL" for exponentially
    tilted empirical likelihood of Schennach(2007), "HD" for Hellinger
    Distance of Kitamura-Otsu-Evdokimov (2013), and "ETHD" for the
    exponentially tilted Hellinger distance of Antoine-Dovonon
    (2015). "REEL" is a restricted version of "EEL" in which the
    probabilities are bounded below by zero. In that case, an analytical
    Kuhn-Tucker method is used to find the solution.}

  \item{initTheta}{Method to obtain the starting values for the
    coefficient vector. By default the GMM estimate with identity matrix
    is used. The second argument means that the theta0 of the
    object, if any, should be used.}
  
  \item{lambda0}{Manual starting values for the Lagrange
    multiplier. By default, it is a vector of zeros.}

  \item{getVcov}{Should the method computes the covariance matrices of the
              coefficients and Lagrange multipliers.}

  \item{lamSlv}{An alternative solver for the Lagrange multiplier. By
    default, either \code{\link{Wu_lam}}, \code{\link{EEL_lam}},
    \code{\link{REEL_lam}} or \code{\link{getLambda}} is
    used.}
  
  \item{coefSlv}{Minimization solver for the coefficient vector.}
  
  \item{lControl}{A list of controls for the Lagrange multiplier
    algorithm.}
  
  \item{tControl}{A list of controls for the coefficient algorithm.}
	    
}

\value{
'causalGEL' returns an object of classesof \code{"causalGelfit"}.
 }

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
           
fit1 <-  causalGEL(g, balm, nsw, gelType="ET")
fit1

fit2 <-  causalGEL(g, balm, nsw, gelType="EL")
fit2

}

