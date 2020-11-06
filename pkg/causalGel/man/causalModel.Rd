\name{causalModel}

\alias{causalModel}
	
\title{Constructor for \code{"causalModel"} classes}

\description{
It builds the object of class \code{"functionModel"}, but with a
specific function that depends on the type of causality. 
}
\usage{
causalModel(g, balm, data,theta0=NULL,
            momType=c("ACE","ACT","ACC", "uncondBal"),
            popMom = NULL, ACTmom=1L, orthoBases=FALSE) 
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

  \item{ACTmom}{When \code{momType} is set to 'ACT', that integer
    indicates which treated group to use to balance the covariates.}

  \item{orthoBases}{If \code{TRUE}, the matrix of balancing moments is
  replaced by the matrix of orthogonal bases that span the same
  space. It is likely to be more stable and to prevent some balancing
  moments to be collinear.}
   
}

\value{
It returns an object of class \code{"causalModel"}.
 }

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
           
model <- causalModel(g, balm, nsw)

}

