\name{causalModel}

\alias{causalModel}
	
\title{Constructor for \code{"causalGel"} classes}

\description{
It builds the object of either class \code{"linearGmm"},
\code{"nonlinearGmm"} or \code{"functionGmm"}. This is the first step
before running any estimation algorithm.
}
\usage{
causalModel(g, balm, data,theta0=NULL,
            momType=c("ACE","ACT","ACC", "uncondBal","fixedMom"),
            popMom = NULL, rhoFct=NULL,ACTmom=1L, 
            gelType = c("EL", "ET", "EEL", "ETEL", "HD", "ETHD","REEL"))
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
  
}

\value{
'gmmModel' returns an object of one of the subclasses of \code{"gmmModels"}.
 }

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
           
model <- causalModel(g, balm, nsw)

}

