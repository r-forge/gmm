gelModel <- function(g, x=NULL, gelType, rhoFct=NULL, tet0=NULL,grad=NULL,
                     vcov = c("HAC", "MDS", "iid"),
                     kernel = c("Quadratic Spectral",  "Truncated", "Bartlett", "Parzen",
                          "Tukey-Hanning"), crit = 1e-06,
                     bw = "Andrews", prewhite = 1L, ar.method = "ols", approx = "AR(1)", 
                     tol = 1e-07, centeredVcov = TRUE, data=parent.frame())
    {
        vcov <- match.arg(vcov)
        kernel <- match.arg(kernel)
        args <- as.list(match.call())
        args$rhoFct <- NULL
        args$gelType <- NULL
        model <- do.call(gmmModel, args)
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

EL.Wu <- function (gmat, l0=NULL, tol = 1e-08, maxiter = 50, k=1) 
    {
        gmat <- as.matrix(gmat)*k
        if (is.null(l0))
            l0 <- rep(0, ncol(gmat))
        n = nrow(gmat)
        dif = 1
        j = 0
        while (dif > tol & j <= maxiter) {
            D1 = t(gmat) %*% ((1/(1 + gmat %*% l0)))
            DD = -t(gmat) %*% (c((1/(1 + gmat %*% l0)^2)) * gmat)
            D2 = solve(DD, D1, tol = 1e-40)
            dif = max(abs(D2))
            rule = 1
            while (rule > 0) {
                rule = 0
                if (min(1 + t(l0 - D2) %*% t(gmat)) <= 0) 
                    rule = rule + 1
                if (rule > 0) 
                    D2 = D2/2
            }
            l0 = l0 - D2
            j = j + 1
        }
        if (j >= maxiter) {
            l0 = rep(0, ncol(gmat))
            conv = list(convergence = 1)
        } else {
            conv = list(convergence = 0)
        }
        return(list(lambda = c(-l0), convergence = conv))
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
            return(EL.Wu(gmat, l0, tol, maxiter, k))
        
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
        theta <- gmmFit(as(object, "gmmModels"), weights="ident")@theta
    
    gt <- evalMoment(object, theta)
    gt <- scale(gt, scale=FALSE)
    class(gt) <- "gmmFct"
    if (!(object@kernel%in%c("Bartlett","Parzen")))
        object@kernel <- "Bartlett"
    kernel <- switch(object@kernel,
                     Bartlett="Truncated",
                     Parzen="Bartlett")
    k <- switch(kernel,
                Truncated=c(2,2),
                Bartlett=c(1,2/3))
    if (is.character(object@bw))
        {
            bw <- get(paste("bw", object@bw, sep = ""))
            bw <- bw(gt, kernel = object@kernel, prewhite = object@prewhite,
                     ar.method = object@ar.method, approx = object@approx)
        } else {
            bw <- object@bw
        } 
    w <- weightsAndrews(gt, bw = bw, kernel = kernel, prewhite = object@prewhite, 
                        ar.method = object@ar.method, tol = object@tol, verbose = FALSE, 
                        approx = object@approx)
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
