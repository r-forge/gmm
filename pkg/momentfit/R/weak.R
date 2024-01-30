## LSE for linearModel objects

setGeneric("lse", function(model, ...) standardGeneric("lse"))

setMethod("lse", "linearModel", function(model) {
    f <- formula(model@modelF)
    dat <- model@modelF
    attr(dat, "terms") <- NULL
    fit <- lm(f, dat)
    new("lsefit", fit, model=model)
})

## Methods for lsefit objects

setMethod("print", "lsefit",
          function(x, model=TRUE, ...) {
              if (model)
                  print(x@model)
              cat("\nEstimation: Least Squares\n")
              cat(paste(capture.output(print(as(x, "lm"), ...))[-(1:3)],
                        collapse="\n"))
              invisible()
          })

setMethod("show","lsefit", function(object) print(object))


## K-Class functions


getK <- function(object, alpha=1, returnRes=FALSE)
{
    if (!inherits(object, "linearModel"))
        stop("object must be a model of class linearModel")
    spec <- modelDims(object)
    Res <- NULL    
    if (spec$k == spec$q)
    {
        k1 <- 1
    } else {
        X <- model.matrix(object)
        endo <- cbind(model.response(object@modelF),
                      X[,object@isEndo])
        X <- X[,!object@isEndo, drop=FALSE]
        Z <- model.matrix(object, "instruments")
        e1 <- lm.fit(X, endo)$residuals
        e2 <- lm.fit(Z, endo)$residuals
        if (returnRes) Res <- e2[,-1,drop=FALSE]        
        e1 <- crossprod(e1)
        e2 <- crossprod(e2)                      
        k1 <- min(eigen(solve(e2,e1))$val)
    }
    k2 <- k1 - alpha/(spec$n-spec$q)
    if (!returnRes)
    {
        c(LIML=k1, Fuller=k2)
    } else {
        ans <- list(k=c(LIML=k1, Fuller=k2), Res=Res)
        class(ans) <- "kappa"
        ans
    }
}


kclassfit <- function(object, k, type=c("LIML", "Fuller", "BTSLS"), alpha = 1)
{
    if (!inherits(object, "linearModel"))
        stop("object must be a model of class linearModel")    
    type <- match.arg(type)
    spec <- modelDims(object)    
    if (missing(k))
    {
        if (type == "BTSLS")
        {
            Kres <- list()
            k2 <- spec$q-sum(!spec$isEndo)
            k <- spec$n/(spec$n-k2+2)
        } else {
            Kres <- getK(object, alpha, TRUE)
            k <- Kres$k[type]
        }
        method <- type
    } else {
        method="K-Class"
        if (is.list(k))
        {
            if (!inherits(k, "kappa"))
                stop("when k is a list, it must have been generared by getK")
            Kres <- k
            k <- Kres$k[type]
            method <- type
        } else {
            if (!is.numeric(k))
                stop("k must be a numeric vector")
            if (length(k)>1)
                stop("The length of k must be 1")
            Kres <- list()
        }
    }
    if (k==1)
        return(tsls(object))
    if (k==0)
        return(lse(object))
    EndoVars <- !(spec$parNames %in% spec$momNames)
    exoInst <- spec$momNames %in% spec$parNames
    if (all(!EndoVars))
    {
        warning("This model has no endogenous variables. Returning LSE.")
        return(lse(object))
    }
    g. <- formula(object@modelF)        
    X <- model.matrix(object)
    if (is.null(Kres$Res))
    {
        Z <- model.matrix(object, "instrument")
        Z <- X[, EndoVars, drop=FALSE] - as.matrix(k*lm.fit(Z, X[,EndoVars])$residuals)
    } else {
        Z <- X[, EndoVars, drop=FALSE] - as.matrix(k*Kres$Res)
    }
    colnames(Z) <- paste("KClass.", colnames(Z), sep="")
    parNames <- spec$parNames        
    parNames[EndoVars] <- colnames(Z)
    if (attr(terms(g.), "intercept"))
        parNames <- parNames[-1]
    h. <- reformulate(parNames, NULL, attr(terms(g.), "intercept"))
    dat <- object@modelF
    attr(dat, "terms") <- NULL
    dat <- cbind(dat, Z)
    object2 <- momentModel(g=g., x=h., data=dat,
                          vcov=object@vcov, vcovOptions=object@vcovOptions,
                          survOptions=object@survOptions,
                          centeredVcov=object@centeredVcov, smooth=object@smooth)
    fit <- gmmFit(object2)
    new("kclassfit", fit, method=method, kappa=k, origModel = object)
}


## Methods for kclassfit

setMethod("print", "kclassfit",
          function(x, model=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@origModel)
              type <- x@method
              k <- x@kappa
              spec <- modelDims(x@model)
              cat("\nEstimation: ", type, sep="")
              cat(" (k = ", k, ")\n", sep="")
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
          })
setMethod("show","kclassfit", function(object) print(object))


setMethod("specTest", c("kclassfit", "missing"),
          function(object, which, ...) {
              spec <- modelDims(object@origModel)
              AR <- log(object@kappa)*spec$n
              df <- spec$q-spec$k
              pv <- ifelse(df>0, pchisq(AR, df, lower.tail=FALSE), NA)
              AR <- cbind(AR, df, pv)
              dimnames(AR) <- list("Test E(g)=0:  ", c("Statistics", 
            "df", "pvalue"))
        ans <- new("specTest", test = AR, testname = "Anderson and Rubin")
        ans              
          })

setMethod("summary", "kclassfit",
          function(object, ...)
          {
              ans <- callNextMethod()
              new("summaryKclass", kappa=object@kappa, method=object@method,
                  origModel=object@origModel, ans)
    })


## Methods for summaryKclass

setMethod("print", "summaryKclass",
          function(x, ...) {
              p <- capture.output(print(as(x, "summaryGmm"), ...))
              type <- x@method
              k <- x@kappa
              method <- paste("Estimation: ", type,  " (k = ", k, ")", sep="")
              p[grepl("Estimation:", p)] <- method
              cat(paste(p, collapse="\n"))
              cat("\n")
              invisible()
          })
setMethod("show","summaryKclass", function(object) print(object))


## Stock and Yogo (2005)

CDtest <- function(object, print=TRUE, ...)
{
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    spec <- modelDims(object)
    Z <- model.matrix(object, "instrument")
    X <- model.matrix(object)
    z2n <- !(colnames(Z) %in% colnames(X))
    X1 <- X[,!spec$isEndo, drop=FALSE]
    X2 <- X[,spec$isEndo, drop=FALSE]
    secS <- lm.fit(Z, X2)        
    Z2 <- Z[,z2n,drop=FALSE]
    tmp <- if (ncol(X2)>1)
               Z[,z2n,drop=FALSE]%*%secS$coef[z2n,,drop=FALSE]
           else
               Z[,z2n,drop=FALSE]%*%secS$coef[z2n]
    e <- secS$residuals
    e2 <- lm.fit(X1, X2)$residuals
    e <- crossprod(e)/(spec$n-spec$q)
    e2 <- crossprod(e2, tmp)/ncol(Z2)
    test <- min(eigen(solve(e,e2))$val)
    if (!print)
        return(test)
    cat("Cragg and Donald Test for Weak Instruments\n",
        "******************************************\n", sep="")
    cat("Number of included Endogenous: ", ncol(X2), "\n", sep="")
    cat("Number of excluded Exogenous: ", ncol(Z2), "\n", sep="")
    cat("Statistics: ", formatC(test, ...), "\n\n", sep="")
    SYTables(object)
}

SYTables <- function(object, print=TRUE)
{
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    s <- modelDims(object)
    l <- sum(s$isEndo)
    k <- s$q - (s$k - l)       
    t <- sizeTSLS
    sizeTSLS <- t[attr(t, "excExo") == k, attr(t, "incEndo") == l]
    names(sizeTSLS) <- paste("size=",
                             attr(t, "size")[attr(t, "incEndo") == l], sep="")
    t <- biasFuller
    biasFuller <- t[attr(t, "excExo") == k, attr(t, "incEndo") == l]
    names(biasFuller) <- paste("bias=",
                               attr(t, "bias")[attr(t, "incEndo") == l], sep="")
    t <- sizeLIML
    sizeLIML <- t[attr(t, "excExo") == k, attr(t, "incEndo") == l]
    names(sizeLIML) <- paste("size=",
                             attr(t, "size")[attr(t, "incEndo") == l], sep="")
    if ((k-l)>=2)
    {
        t <- biasTSLS
        biasTSLS <- t[attr(t, "excExo") == k, attr(t, "incEndo") == l]
        names(biasTSLS) <- paste("bias=",
                                 attr(t, "bias")[attr(t, "incEndo") == l], sep="")        
    } else {
        biasTSLS <- NULL
    }

    crit <- list(biasTSLS=biasTSLS,
                 sizeTSLS=sizeTSLS,
                 biasFuller=biasFuller,
                 sizeLIML=sizeLIML)
    if (!print)
        return(crit)
    nCrit <- c("Target relative bias for TSLS:\n",
               "Target size for TSLS:\n",
               "Target relative bias for Fuller-K:\n",
               "Target size for LIML:\n")
    
    cat("Stock and Yogo (2005) critical values\n")
    cat("*************************************\n")
    for (i in 1:4)
    {
        if (!is.null(crit[[i]]))
        {
            cat(nCrit[i])
            print.default(format(crit[[i]]), print.gap = 2L, quote = FALSE)
            cat("\n")
        }
    }
    invisible()
}
