\name{KConfid}
\alias{KConfid}

\title{Confidence interval using the K statistics of Kleibergen}
\description{The confidence is either an interval or an ellipse.
}
\usage{
KConfid(obj, which, type = c("K", "KJ"), alpha = 0.05, alphaJ = 0.01, n = 4)
}
\arguments{
 \item{obj}{Object of class "gmm" returned by \link{gmm} (not restricted)}
 \item{type}{Should we base the confidence interval on the K or K-J statistics.}
 \item{which}{A 2x1 vector or a scalar. The interval is computed for \code{coef(obj)[which]}. }
 \item{alpha, alphaJ}{The overall size and the size for the J-test when type is "KS".}
 \item{n}{The number of points to compute the confidence region is 4(n-1). It must be greater than 2.}
}


\value{
Interval for \code{lenght(which)=1} or a series of points if \code{lenght(which)>1}.
}

\references{

 Kleibergen, F. (2005),
  Testing Parameters in GMM without assuming that they are identified.
  \emph{Econometrica}, \bold{73},
  1103-1123,
 
}

\examples{

data(Finance)
r <- Finance[1:300, 1]
rf <- Finance[1:300, "rf"]
z <- as.matrix(r-rf)
zm <-  Finance[1:300, "rm"]-rf
f1 <- zm
f2 <- Finance[1:300, "hml"] - rf
f3 <- Finance[1:300, "smb"] - rf
res <- gmm(z ~ f1 + f2 + f3, ~ f1 + f2 + f3)

KConfid(res,2)
sol <- KConfid(res,c(2,3))
plot(sol, main="Confidence Region")
polygon(sol,col="grey")
points(res$coef[2],res$coef[3],pch=21,bg=1)
text(res$coef[2],res$coef[3],expression(hat(theta)),pos=3)

}
