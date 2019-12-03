\name{[-causalGel}
\docType{methods}
\alias{[,causalGel,missing,numericORlogical-method}
\alias{[,causalGel,numeric,missing-method}
\alias{[,causalGel,numeric,numericORlogical-method}
\alias{[,causalGelfit,missing,numericORlogical-method}
\alias{[,causalGelfit,numeric,missing-method}
\alias{[,causalGelfit,numeric,numericORlogical-method}


\title{Subsetting methods}
\description{
Different subsetting methods for S4 class objects of the package. The
subset method returns an new object with observations selected by the
second argument. See example.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "causalGel", i = "missing", j = "numericORlogical")}}{
  Subsets observations.
}
  
\item{\code{signature(x = "causalGel", i = "numeric", j = "missing")}}{
  Selects balancing covatriates.
}

\item{\code{signature(x = "causalGel", i = "numeric", j = "numericORlogical")}}{
  Selects balancing covariates and observations.
}

\item{\code{signature(x = "causalGelfit", i = "missing", j = "numericORlogical")}}{
  Subsets observations and refit the model.
}
  
\item{\code{signature(x = "causalGelfit", i = "numeric", j = "missing")}}{
  Selects balancing covatriates and refit the model.
}

\item{\code{signature(x = "causalGelfit", i = "numeric", j = "numericORlogical")}}{
  Selects balancing covariates and observations and refit the model.
}



}}

\examples{
data(nsw)

balm <- ~age+ed+black+hisp:married+nodeg+re75+I(re75^2)
g <- re78~treat
         
model <- causalModel(g, balm, nsw)           
model[1:5, 1:500]

fit <-  modelFit(model)
fit[1:5,1:500]

}



\keyword{methods}
\keyword{subsetting}

