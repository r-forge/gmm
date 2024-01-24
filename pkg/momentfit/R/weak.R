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


getK <- function(object, alpha=1)
{
    if (!inherits(object, "linearModel"))
        stop("object must be a model of class linearModel")
    spec <- modelDims(object)
    if (spec$k == spec$q)
    {
        k1 <- 1
    } else {
        X <- model.matrix(object)
        endo <- cbind(model.response(object@modelF),
                      X[,object@isEndo])
        X <- X[,!object@isEndo, drop=FALSE]
        W <- model.matrix(object, "instruments")
        e1 <- lm.fit(X, endo)$residuals
        e2 <- lm.fit(W, endo)$residuals
        e1 <- crossprod(e1)
        e2 <- crossprod(e2)                      
        k1 <- min(eigen(solve(e2,e1))$val)
    }
    k2 <- k1 - alpha/(spec$n-spec$q)
    c(LIML=k1, Fuller=k2)
}


kclassfit <- function(object, k, type=c("LIML", "Fuller"))
{
    if (!inherits(object, "linearModel"))
        stop("object must be a model of class linearModel")    
    type <- match.arg(type)
    if (missing(k))
    {
        k <- getK(object)[type]
        method <- type
    } else {
        method="K-Class"
    }
    if (k==1)
        return(tsls(object))
    if (k==0)
        return(lse(object))
    spec <- modelDims(object)
    EndoVars <- !(spec$parNames %in% spec$momNames)
    exoInst <- spec$momNames %in% spec$parNames
    if (all(!EndoVars))
    {
        warning("This model has no endogenous variables. Returning LSE.")
        return(lse(object))
    }
    g. <- formula(object@modelF)        
    X <- model.matrix(object)
    Z <- model.matrix(object, "instrument")
    Z <- X[, EndoVars, drop=FALSE] - as.matrix(k*lm.fit(Z, X[,EndoVars])$residuals)
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
