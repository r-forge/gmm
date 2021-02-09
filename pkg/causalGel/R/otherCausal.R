###################################
########  Matching methods
##################################

## units in x2 are matched to each element of x1
## minn is the minimum number of matches
## tol: if other individuals are within tol of 
## the highest distance in the first minn matches, then they are also included.
## y2 is y from group 2
## it returns:
##             list in indices for the match
##             E(y2 | z=1, x) (n1 x 1 vector of estimated missing Y(2))

getNN <- function(x,y,z,fromTo, minn, tol=1e-7)    
{
    if (length(fromTo)!=2)
        stop("fromTo must be a vector of 2 with 0s or 1s")
    if (!all(fromTo %in% c(0,1)))
        stop("fromTo must be a vector of 2 with 0s or 1s")
    x <- as.matrix(x)
    self <- z==fromTo[1]
    selt <- z==fromTo[2]
    x1 <- x[self,,drop=FALSE]
    x2 <- x[selt,,drop=FALSE]
    y2 <- y[selt]
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    k <- ncol(x1)
    type <- paste("Nearest neighbour of the ",
                  ifelse(fromTo[1]==1, "treated", "control"), " among the ",
                  ifelse(fromTo[2]==1, "treated", "control"), sep="")
    res <- .Fortran(F_findnn, as.double(x1), as.double(x2),as.double(y2),
                    as.double(tol),
                    as.integer(n1), as.integer(n2), as.integer(k),
                    as.integer(minn), K=double(n2), nn=integer(n1),
                    ind=integer(n1*n2), y2hat=double(n1))
    nn <- res$nn
    ind <- matrix(res$ind,nrow=n1)
    list(matches=lapply(1:n1, function(i) ind[i,1:nn[i]]), K=res$K,
         y2hat = res$y2hat, type=type)
}

## Function to get the missing counterfactual in y
## if which="y0", y[z==1] are replace by an estimate Y(0)|z=1
## if which="y1", y[z==0] are replace by an estimate Y(1)|z=0
## if type="match", only match method is used
## if type="reg", only regression is used
## if type="bc_match" the bias-corrected method is used.
## x is the matrix used for matching
## regMat is the regression matrix for E(Y|X,Z=1) or  E(Y|X,Z=0)

getMissingY <- function(y, z, x, which=c("y0", "y1"),
                        type=c("match","bc_match","reg","all"), minn=4, tol=1e-7,
                        regMat=cbind(1,x))
{
    which <- match.arg(which)
    type <- match.arg(type)
    ## regMat should have an intercept included
    x <- as.matrix(x)
    if (which == "y0")
    {
        sel <- z==1
        fromTo <- c(1,0)
    } else {
        sel <- z==0
        fromTo <- c(0,1)
    }
    res <- getNN(x,y, z, fromTo, minn, tol)
    if (type != "match" | type == "all")
        b <- lm.wfit(regMat[!sel,,drop=FALSE], y[!sel], res$K)$coefficients
    if (type == "match")
    {
        y[sel] <- res$y2hat
        y <- as.matrix(y)
        colnames(y) <- "match"
    } else if (type == "reg") {
        y[sel] <- c(regMat[sel,,drop=FALSE]%*%b)
        y <- as.matrix(y)
        colnames(y) <- "reg"
    } else if (type=="bc_match") {
        fity <- c(regMat%*%b)
        bc <- fity[sel] - sapply(res$matches, function(l) mean(fity[!sel][l]))
        y[sel] <- res$y2hat + bc
        y <- as.matrix(y)
        colnames(y) <- "bc"
    } else {
        ym <- y
        ym[sel] <- res$y2hat
        yreg <- y
        yreg[sel] <- c(regMat[sel,,drop=FALSE]%*%b)
        ybc <- y
        fity <- c(regMat%*%b)
        bc <- fity[sel] - sapply(res$matches, function(l) mean(fity[!sel][l]))
        ybc[sel] <- res$y2hat + bc
        y <- cbind(ym, ybc, yreg)
        colnames(y) <- c("match","bc","reg")
    }
    attr(y, "which") <- which
    y
}

## PSForm is the formula used to compute the propensity score.
## BCorForm is the bias correction regression model

matching <- function(form,  balm, data, type=c("ACE","ACT","ACC"), M=4,
                      adjust=FALSE, matchPS=FALSE, psForm=NULL,
                      bcForm=NULL)
{
    ## If bcForm does not have an intercept it will be added.
    type <- match.arg(type)
    Call <- match.call()
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    if (is.null(bcForm) & adjust)        
        bcForm <- balm
    if (is.null(psForm) & matchPS)
    {
        psForm <- attr(terms(balm), "term.labels")
        psForm <- reformulate(psForm, response=colnames(mf)[2])
    }
    BC <- NA
    if (!matchPS)
    {
        balm <- attr(terms(balm), "term.labels")
        balm <- reformulate(balm, intercept=FALSE)
        X <- model.matrix(balm, data)
    } else {
        res <- glm(psForm, data=data, family=binomial())
        X <- as.matrix(fitted(res))
        data$PScore <- X
        balm <- ~PScore-1
    }
    X <- scale(X)
    if (adjust)
    {
        bcForm <- reformulate(attr(terms(bcForm), "term.labels"))
        Z <- model.matrix(bcForm, data)
    } else {
        Z <- NULL
    }
    typeCor <- ifelse(adjust, "all", "match")
    if (type == "ACT")
    {
        Y0 <- getMissingY(Y, T, X, which="y0", type=typeCor, minn=M,
                          regMat=Z)
        NN <- mean(Y[T==1]-Y0[T==1,"match"])
        if (adjust)
            BC <- mean(Y[T==1]-Y0[T==1,"bc"])
    } else if (type == "ACC") {
        Y1 <- getMissingY(Y, T, X, which="y1", type=typeCor, minn=M,
                          regMat=Z)
        NN <- mean(Y1[T==0,"match"] - Y[T==0])
        if (adjust)                
            BC <- mean(Y1[T==0,"bc"] - Y[T==0])
    } else {
        Y0 <- getMissingY(Y, T, X, which="y0", type=typeCor, minn=M,
                          regMat=Z)
        Y1 <- getMissingY(Y, T, X, which="y1", type=typeCor, minn=M,
                          regMat=Z)
        NN <- mean(Y1[,"match"]-Y0[,"match"])
        if (adjust)                
            BC <- mean(Y1[,"bc"]-Y0[,"bc"])                           
    }
    estim <- c(NN=NN, NN.BiasCor=BC)
    form <- list(balForm=balm, psForm=psForm, bcForm=bcForm)
    info <- list()
    method <- ifelse(matchPS, "PS Matching Method", "Covariate Matching Method")
    details <- list(ifelse(matchPS,
                           paste("PS formula: ", deparse(psForm), sep=""), 
                           paste("Covariate formula: ",
                                 deparse(balm), sep="")),
                    paste("Min. number of matches: ", M, sep=""))
    if (adjust)
        details <- c(details, list(paste("Bias-cor. formula: ",
                                         deparse(bcForm),sep="")))
    shortm <- paste(c("","BiasCor_"),
                    paste(ifelse(matchPS, "ps","cov"),
                          "Matching", M, "_", type, sep=""),sep="")
    names(shortm) <- c("NN","NN.BiasCor")
    shortMet <- shortm
    new("causalfit", estim=estim, type=type, method=method,
        form=form, details=details, info=info, data=data, call=Call)
}



### ACE using local linear regression matching
### if h is set, it is used, if not, the optimal h is computed
### by minimizing the CV. A grid in seq(from, to, length=nh) is
### first obtain and the brent method is used to complete.
### In the output, if info=0, it converged, if info=1, the minimum
### is at from or to, so the brent merhod was not launch, and for info = 2
### the brent reached the minimum iteration. 
### the method requires propensity score so PSForm must be provided

LLmatching <- function(form, psForm, data, type=c("ACE","ACT","ACC"),
                    kern=c("Gaussian","Epanechnikov"),tol=1e-4,
                    h=NULL, from=.00001, to=5, ngrid=10, maxit=100,
                    hMethod=c("Brent","Grid"))
{
    type <- match.arg(type)
    Call <- match.call()
    kern <- match.arg(kern)
    hMethod <- match.arg(hMethod)
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    res <- glm(psForm, data=data, family=binomial())
    X <- as.matrix(fitted(res))
    if (type == "ACT")
    {
        Y0 <- getMissingY.LL(Y, T, X, "y0", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        ACE <- mean(Y[T==1]-Y0[T==1], na.rm=TRUE)
        info <- attr(Y0, "h")
    } else if (type == "ACC") {
        Y1 <- getMissingY.LL(Y, T, X, "y1", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        info <- attr(Y1, "h")
        ACE <- mean(Y1[T==0] - Y[T==0], na.rm=TRUE)
    } else {
        Y0 <- getMissingY.LL(Y, T, X, "y0", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        Y1 <- getMissingY.LL(Y, T, X, "y1", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        info <- list(Y0=attr(Y0, "h"), Y1=attr(Y1, "h"))
        ACE <- mean(Y1-Y0, na.rm=TRUE)
    }
    estim <- c(LL=ACE)
    info2 <- unlist(info)
    method <- "Local-Linear Regression Method"
    details <- list(paste("Kernel: ", kern, sep=""),
                    paste("PS formula: ", deparse(psForm), sep=""),
                    ifelse(is.null(h), paste("Bandwidth select: ",hMethod,sep=""),
                           "Bandwidth select: Fixed"))
    if (type == "ACE")
    {
        tmp <- paste("h (Y(0)) = ", round(info[[1]][[1]], 6), sep="")
        tmp2 <- paste("h (Y(1)) = ", round(info[[2]][[1]], 6), sep="")
        details <- c(details, tmp, tmp2)                                 
    } else {
        tmp <- paste("h = ", round(info[[1]], 6), sep="")
        details <- c(details, tmp)                
    }                        
    if (is.null(h))
    {
        if (type == "ACE")
        {
            tmp <- paste("Conv. code for bandwidth (Y(0)): ",
                         info[[1]][[2]], sep="")
            tmp2 <- paste("Conv. code for bandwidth (Y(1)): ",
                          info[[2]][[2]], sep="")
            details <- c(details, tmp, tmp2)
        } else {
            details <- c(details,
                         list(paste("Conv. code for  bandwidth: ",
                                    info[[2]],sep="")))
        }
    }
    if (type == "ACE")
        info2 <- list(meanh = mean(c(info[[1]][[1]], info[[2]][[1]])),
                      convergence = max(info[[1]][[2]], info[[2]][[2]]))
    shortm <- paste("LL", "_", strtrim(kern,3),"_",
                    ifelse(is.null(h), "autoh", "fixedh"),"_", type, sep="")
    shortMet <- shortm
    new("causalfit", estim=estim, type=type, method=method,
        form=list(psForm=psForm),data=data,
        details=details, info=info2, call=Call)
}

## p is the fitted propensity score
## The function estimate the counterfactual Y(0) or Y(1) 

getMissingY.LL <- function(y, z, p, which=c("y0", "y1"), h=NULL,
                           kern=c("Gaussian","Epanechnikov"),
                           from, to, ngrid, tol, maxit,
                           hMethod=c("Brent","Grid"))
{
    which <- match.arg(which)
    hMethod <- match.arg(hMethod)
    kern <- match.arg(kern)
    if (which == "y0")
        sel <- z==1
    else
        sel <- z==0
    missY <- llrF(p[sel], p[!sel], y[!sel], h,
                  kern, from, to, ngrid, tol, maxit, hMethod)
    y[sel] <- missY
    attr(y, "which") <- which
    attr(y, "h") <- attr(missY, "h")
    y
}

cvF <- function(p, y, h, kern=c("Gaussian","Epanechnikov"))
{
    kern <- match.arg(kern)
    nh <- length(h)
    kernF <- tolower(strtrim(kern,1))
    n <- length(p)
    res <- .Fortran(F_cvfct, as.double(p), as.double(y),
                    as.integer(n), as.double(h), as.character(kernF),
                    as.integer(nh), cv=double(nh))
    res$cv
}

llrF <- function(p1, p0, y0, h=NULL, kern=c("Gaussian","Epanechnikov"), 
                 from, to, ngrid, tol, maxit, hMethod=c("Brent","Grid"))
{
    hMethod <- match.arg(hMethod)
    kern <- match.arg(kern)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)
    n1 <- length(p1)
    n2 <- length(p0)
    convergence <- NA
    if (is.null(h))
    {
        resh <- optCVF(p0, y0, kern, from, to, ngrid,
                       tol, maxit, hMethod)
        h <- resh$h
        convergence <- resh$info
    }
    res <- .Fortran(F_llr, as.double(p1), as.double(p0), as.double(y0),
                    as.integer(n1), as.integer(n2), as.double(h),
                    as.integer(kernF), info = integer(n1),
                    y0hat = double(n1))
    y0hat <- res$y0hat
    y0hat[res$info == 1] <- NA
    attr(y0hat, "h") <- list(h=h,convergence=convergence)
    y0hat
}


gridcvF <- function(p, y, kern=c("Gaussian","Epanechnikov"), from, to, nh)
{
    kern <- match.arg(kern)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)
    n <- length(p)
    res <- .Fortran(F_gridcv, as.double(p), as.double(y),
                    as.integer(n), as.integer(kernF),
                    as.double(from), as.double(to), as.integer(nh),
                    info=integer(1), h=double(3), cv=double(3))
    cbind(res$h,res$cv)
}

### The gricv fortran subroutine returns, if info=0,  a 3x1 vector of h and cv, 
### such that cv[1]>cv[2] and cv[3]>cv[2]. The minimum is therefore the middle one.

optCVF <- function(p, y, kern=c("Gaussian","Epanechnikov"), from, to, nh,
                   tol=1e-4, maxit=100, method=c("Brent","Grid"))
{
    kern <- match.arg(kern)
    method <- match.arg(method)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)        
    n <- length(p)
    res <- .Fortran(F_gridcv, as.double(p), as.double(y),
                    as.integer(n), as.integer(kernF),
                    as.double(from), as.double(to), as.integer(nh),
                    info=integer(1), h=double(3), cv=double(3))
    if (res$info == 1)
    {
        return(list(cv=res$cv[1], h=res$h[1], info=res$info))
    }
    if (method == "Grid")
    {
        return(list(cv=res$cv[2], h=res$h[2], info=res$info))
    }
    res2 <- .Fortran(F_brentcv, as.double(p), as.double(y),
                     as.integer(n), as.integer(kernF),
                     as.double(tol), as.integer(maxit),
                     as.double(res$h), as.double(res$cv),
                     info=integer(1), h=double(1), cv=double(1))
    info <- ifelse(res2$info == 0, 0, 2)
    list(cv=res2$cv, h=res2$h, info=info)
}


### Weighting Methods
########################

ipw <- function(form, psForm, data, type=c("ACE","ACT","ACC"),
                normalized=FALSE, ...)
{
    type <- match.arg(type)
    Call <- match.call()
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    mnorm <- normalized
    info <- list()
    if (normalized == "GPE")
    {
        res <- getGPE.PS(form, psForm, data, ...)
        X <- res$phat
        info <- list(res$info)
        normalized <- TRUE
    } else {                        
        res <- glm(psForm, data=data, family=binomial())
        X <- as.matrix(fitted(res))
    }
    if (type == "ACT")
    {
        s1 <- sum(T)
        w1 <- rep(1/s1, s1)
        w0 <- X[T==0]/(1-X[T==0])/s1
        if (normalized)
            w0 <- w0/sum(w0)
    } else if (type == "ACC") {
        s0 <- sum(1-T)
        w1 <- (1-X[T==1])/X[T==1]/s0
        w0 <- rep(1/s0,s0)
        if (normalized)
            w1 <- w1/sum(w1)
    } else {
        n <- length(X)
        w1 <- 1/(n*X[T==1])
        w0 <- 1/(n*(1-X[T==0]))
        if (normalized)
        {
            w1 <- w1/sum(w1)
            w0 <- w0/sum(w0)
        }
    }
    ACE <- sum(w1*Y[T==1]) - sum(w0*Y[T==0])
    estim <- c(IPW=ACE)
    form <- list(psForm=psForm)
    method <- "Inverse Probability Weight Method"
    details <- list(paste("PS formula: ", deparse(psForm), sep=""))       
    if (mnorm == "GPE")
        details <- c(details,
                     list(paste("Method: normalized with GPE method")))
    else
        details <- c(details,
                     list(paste("Method: ", ifelse(mnorm, "normalized",
                                                   "non-normalized"), sep="")))
    if (mnorm == "GPE")
        details <- c(details,
                     paste("Convergence of optim for PS: ",
                           info[[1]][2], sep=""))
    shortm <- paste("IPW", "_",
                    ifelse(mnorm=="GPE", "GPE",
                    ifelse(mnorm, "normalized", "unnormalized")),
                    "_", type, sep="")
    shortMet <- shortm
    new("causalfit", estim=estim, type=type, method=method,
        form=list(psForm=psForm),
        details=details, info=info, call=Call)
}


getGPE.PS <- function(form, PSForm=z~kX, data, ...)
    {
        mf <- model.frame(form, data)
        Y <- c(mf[[1]])
        T <- c(mf[[2]])
        X <- model.matrix(PSForm, data)
        g <- function(theta, T, X)
            {
                kX <- c(X%*%theta)
                G <- plogis(kX)
                e <- (T-G)/(1-G)
                m <- as.matrix(e*X)
                mean(colMeans(m, na.rm=TRUE)^2)
            }
        res <- glm(PSForm, data=data, family=binomial())        
        res1 <- optim(coef(res), g, method="BFGS", T=T, X=X, ...)
        b <- res1$par
        fit <- c(X%*%b)
        phat <- plogis(fit)
        phat[phat>.999] <- .999
        list(phat=phat, info=c(obj=res1$value, convergence=res1$convergence))
    }


## print

setMethod("print", "gmmfit",
          function(x, model=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              ntype <- matrix(c("Two-Step GMM", "Iterated GMM", "CUE",
                                "One-Step GMM with fixed weights","Two-Stage Least Squares",
                                "Evaluated at a fixed Theta; No estimation",
                                "One-Step Efficient M.D.E.",
                                "twostep","iter","cue","onestep","tsls", "eval","mde"),
                              ncol=2)
              type <- ntype[match(x@type, ntype[,2]),1]
              spec <- modelDims(x@model)
              if (spec$q==spec$k && x@type != "eval")
                  type <- "One-Step, Just-Identified"
              cat("\nEstimation: ", type,"\n")
              if (!is.null(x@convergence))
                  cat("Convergence Optim: ", x@convergence, "\n")              
              if (!is.null(x@convIter))
                  cat("Convergence Iteration: ", x@convIter, "\n")
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
          })

## show

setMethod("show","causalfit", function(object) print(object))

## print

setMethod("print", "causalfit",
          function(x, digits=5, ...) {
              cat(x@method, "\n")
              cat("**************************\n")
              cat("Type of causal effect: ", x@type, "\n", sep="")
              if (length(x@estim)==1)
              {
                  cat("Causal estimate: ",
                      format(x@estim, digits=digits, ...), "\n")
              } else {
                  cat("Causal estimates: \n")
                  for (i in 1:length(x@estim))
                      cat("\t", names(x@estim)[i], " = ",
                          format(x@estim[i], digits=digits, ...),
                          "\n", sep="")
              }
              if (length(x@details))
              {
                  cat("Details:\n")
                  for (i in 1:length(x@details))
                      cat("\t", x@details[[i]], "\n")
              }
              invisible(x)
          })

