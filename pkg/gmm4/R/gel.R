gelModel <- function(g, x=NULL, gelType, rhoFct=NULL, tet0=NULL,grad=NULL,
                     vcov = c("HAC", "MDS", "iid"),
                     vcovOptions=list(), centeredVcov = TRUE, data=parent.frame())
    {
        vcov <- match.arg(vcov)
        model <- gmmModel(g=g, x=x, grad=grad, vcov=vcov, vcovOptions=vcovOptions,
                          centeredVcov=centeredVcov,
                          tet0=tet0, data=data)
        gmmToGel(model, gelType, rhoFct)
    }

gmmToGel <- function(object, gelType, rhoFct=NULL)
    {
        cls <- strsplit(class(object), "Gmm")[[1]][1]
        cls <- paste(cls, "Gel", sep="")
        if (object@vcov == "HAC")
            wSpec <- smoothGel(object)
        else
            wSpec <- list(k=c(1,1), w=kernel(1), bw=1, kernel="None")                
        new(cls, wSpec=wSpec, gelType=list(name=gelType, fct=rhoFct),
            object)
    }


rhoEL <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               log(1 - gml),
               -1/(1 - gml),
               -1/(1 - gml)^2)               
    }

rhoET <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -exp(gml)+1,
               -exp(gml),
               -exp(gml))               
    }

rhoEEL <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -gml - 0.5 * gml^2,
               -1 - gml,
               rep(-1, nrow(gmat)))               
    }

rhoREEL <- function(gmat, lambda, derive = 0, k = 1)
{
    rhoEEL(gmat, lambda, derive, k)
}

rhoHD <- function(gmat, lambda, derive = 0, k = 1) 
    {
        lambda <- c(lambda)*k
        gmat <- as.matrix(gmat)
        gml <- c(gmat %*% lambda)
        switch(derive+1,
               -1/(1 + gml)+1,
               1/((1 + gml)^2),
               -2/((1 + gml)^3))               
    }

Wu_lam <- function(gmat, tol=1e-8, maxiter=50, k=1)
    {
        gmat <- as.matrix(gmat)
        res <- .Fortran(F_wu, as.double(gmat), as.double(tol),
                        as.integer(maxiter), as.integer(nrow(gmat)),
                        as.integer(ncol(gmat)), as.double(k),
                        conv=integer(1), obj=double(1),
                        lambda=double(ncol(gmat)))
        list(lambda=res$lambda, convergence=list(convergence=res$conv),
             obj = res$obj)
    }

REEL_lam <- function(gmat, k=1, control=list())
{
        gmat <- as.matrix(gmat)
        n <- nrow(gmat)
        q <- ncol(gmat)
        maxit <- ifelse("maxit" %in% names(control),
                        control$maxit, 50)        
        res <- try(.Fortran(F_lamcuep, as.double(gmat),
                            as.integer(n), as.integer(q), as.double(k),
                            as.integer(maxit),conv=integer(1),
                            lam=double(q),pt=double(n),
                            obj=double(1)
                            ), silent=TRUE)
        if (class(res) == "try-error")
            return(list(lambda=rep(0,q), obj=0, pt=rep(1/n,n),
                        convergence=list(convergence=3)))
        list(lambda=res$lam, obj=res$obj, pt=res$pt,
             convergence=list(convergence=res$conv))
    }

EEL_lam <- function(gmat, k=1)
    {
        q <- qr(gmat)
        n <- nrow(gmat)
        l0 <- -qr.coef(q, rep(1,n))
        conv <- list(convergence=0)
        list(lambda = l0, convergence = conv, obj =
                 mean(rhoEEL(gmat,l0,0,k)))
    }

getLambda <- function (gmat, l0=NULL, gelType, rhoFct=NULL, 
                       tol = 1e-07, maxiter = 100, k = 1, method="BFGS", 
                       algo = c("nlminb", "optim", "Wu"), control = list()) 
    {
        algo <- match.arg(algo)
        gmat <- as.matrix(gmat)
        if (is.null(l0))
            l0 <- rep(0, ncol(gmat))
        if (is.null(rhoFct))
            rhoFct <- get(paste("rho",gelType,sep=""))    
        if (algo == "Wu" & gelType != "EL") 
            stop("Wu (2005) algo to compute Lambda is for EL only")
        if (algo == "Wu") 
            return(Wu_lam(gmat, tol, maxiter, k))
        if (gelType == "EEL")
            return(EEL_lam(gmat, k))
        if (gelType == "REEL")
            return(REEL_lam(gmat, k, control))
        
        fct <- function(l, X, rhoFct, k) {
            r0 <- rhoFct(X, l, derive = 0, k = k)
            -mean(r0)
        }
        Dfct <- function(l, X, rhoFct, k) {
            r1 <- rhoFct(X, l, derive = 1, k = k)
            -colMeans(r1 * X)
        }
        DDfct <- function(l, X, rhoFct, k) {
            r2 <- rhoFct(X, l, derive = 2, k = k)
            -crossprod(X * r2, X)/nrow(X)
        }
        if (algo == "optim") {
            if (gelType == "EL")
                {
                    ci <- -rep(1, nrow(gmat))
                    res <- constrOptim(l0, fct, Dfct, -gmat, ci, control = control,
                                       X = gmat, rhoFct = rhoFct, k = k)
                } else if (gelType == "HD") {
                    ci <- -rep(1, nrow(gmat))
                    res <- constrOptim(l0, fct, Dfct, -gmat, ci, control = control,
                                       X = gmat, rhoFct = rhoFct, k = k)
                } else {
                    res <- optim(l0, fct, gr = Dfct, X = gmat, rhoFct = rhoFct,
                                 k = k, method = method, control = control)
                }
        } else {
            res <- nlminb(l0, fct, gradient = Dfct, hessian = DDfct,
                          X = gmat, rhoFct = rhoFct, k = k, control = control)
        }
        l0 <- res$par
        if (algo == "optim") 
            conv <- list(convergence = res$convergence, counts = res$counts, 
                         message = res$message)
        else
            conv <- list(convergence = res$convergence, counts = res$evaluations, 
                         message = res$message)    
        return(list(lambda = l0, convergence = conv))
    }

smoothGel <- function (object, theta=NULL) 
{
    if (inherits(object, "gelModels"))
        {
            gt <- evalMoment(as(object, "gmmModels"), theta)
            x <- kernapply(gt, object@wSpec$w)
            sx <- list(smoothx = x, w = object@wSpec$w,
                       bw = object@wSpec$bw, k = object@wSpec$k)
            return(sx)
        }
    if (is.null(theta))        
        theta <- modelFit(as(object, "gmmModels"), weights="ident")@theta
    
    gt <- evalMoment(object, theta)
    gt <- scale(gt, scale=FALSE)
    class(gt) <- "gmmFct"
    vspec <- object@vcovOptions
    if (!(vspec$kernel%in%c("Bartlett","Parzen")))
        object@vcovOptions$kernel <- "Bartlett"
    kernel <- switch(object@vcovOptions$kernel,
                     Bartlett="Truncated",
                     Parzen="Bartlett")
    k <- switch(kernel,
                Truncated=c(2,2),
                Bartlett=c(1,2/3))
    if (is.character(vspec$bw))
        {
            bw <- get(paste("bw", vspec$bw, sep = ""))
            bw <- bw(gt, kernel = vspec$kernel, prewhite = vspec$prewhite,
                     ar.method = vspec$ar.method, approx = vspec$approx)
        } else {
            bw <- object@vcovOptions$bw
        } 
    w <- weightsAndrews(gt, bw = bw, kernel = kernel, prewhite = vspec$prewhite, 
                        ar.method = vspec$ar.method, tol = vspec$tol, verbose = FALSE, 
                        approx = vspec$approx)
    rt <- length(w)
    if (rt >= 2)
        {
            w <- c(w[rt:2], w)
            w <- w/sum(w)
            w <- kernel(w[rt:length(w)])
        } else {
            w <- kernel(1)
        }
    return(list(k=k, w=w, bw=bw, kernel=kernel))
}
