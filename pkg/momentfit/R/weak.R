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
        X2 <- model.matrix(object, "includedEndo")
        endo <- cbind(modelResponse(object), X2)
        X <-  model.matrix(object, "includedExo")
        Z <- model.matrix(object, "instruments")
        if (!is.null(X))
        {
            e1 <- lm.fit(X, endo)$residuals
        } else {
            e1 <- endo
        }
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

CDtest <- function(object, print=TRUE, SWcrit=FALSE, ...)
{
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    spec <- modelDims(object)
    Z <- model.matrix(object, "instrument")
    Z2 <- model.matrix(object, "excludedExo")
    X1 <- model.matrix(object, "includedExo")
    X2 <- model.matrix(object, "includedEndo")
    z2n <- colnames(Z) %in% colnames(Z2)
    secS <- lm.fit(Z, X2)        
    df <- ncol(Z2)    
    if (ncol(Z2)==1 & SWcrit)
    {
        warning("The Sanderson and Windmeijer modification of CD-test and rejection rule only applies to models with more than 1 endogenous variables")
        SWcrit <- FALSE
    } else if (SWcrit) {       
        df <- ncol(Z2)-ncol(X2)+1
    }
    tmp <- if (ncol(X2)>1)
               Z[,z2n,drop=FALSE]%*%secS$coef[z2n,,drop=FALSE]
           else
               Z[,z2n,drop=FALSE]%*%secS$coef[z2n]
    e <- secS$residuals
    e2 <- if (!is.null(X1)) lm.fit(X1, X2)$residuals else X2
    e <- crossprod(e)/(spec$n-spec$q)
    e2 <- crossprod(e2, tmp)
    test <- min(eigen(solve(e,e2))$val)/df
    if (!print)
        return(test)
    cat("Cragg and Donald Test for Weak Instruments\n",
        "******************************************\n", sep="")
    if (SWcrit)
    {
        cat("Sanderson and Windmeijer specification\n")
        add1 <- "(-1 for the SW critical value)"
        add2 <- paste("(-", ncol(X2)-1, ")", sep="")
    } else {
        add1 <- add2 <- ""
    }
    cat("Number of included Endogenous: ", ncol(X2), add1, "\n", sep="")
    cat("Number of excluded Exogenous: ", ncol(Z2), add2, "\n", sep="")
    cat("The test is not robust to heteroskedasticity\n")
    cat("Statistics: ", formatC(test, ...), "\n\n", sep="")
    SYTables(object, TRUE, SWcrit)
}

SYTables <- function(object, print=TRUE, SWcrit=FALSE)
{
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    s <- modelDims(object)
    l <- sum(s$isEndo)
    k <- s$q - (s$k - l)
    if (l==1 & SWcrit)
        SWcrit <- FALSE
    if (SWcrit)
    {
        k <- k-l+1        
        l <- l-1
    }
    t <- sizeTSLS
    if (l>2)
    {
        cat("No critical values for models with more than 2 endogenous variables\n")
        return(invisible())
    }
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
    if (SWcrit)
        cat("Critical value adjusted to Sanderson and Windmeijer (2016) specification\n\n")
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

## Sanderson and Windmeijer (2016)

SWtest <- function(object, j=1, print=TRUE, ...)
{
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    spec <- modelDims(object)
    if (sum(spec$isEndo)<1)
    {
        warning("The model does not contain endogenous variables")
        return(NA)
    }
    if (sum(spec$isEndo)==1)
    {
        warning(paste("The number of endgenous variables is equal to 1\n",
                      "Returning the F-test", sep=""))
        return(CDtest(object, print))
    }
    Z2 <- model.matrix(object, "excludedExo")
    X1 <- model.matrix(object, "includedExo")
    X2 <- model.matrix(object, "includedEndo")
    if (!is.null(X1))
    {
        fitX1  <- lm.fit(X1, Z2)
        Z2 <- fitX1$residuals
        X2 <- qr.resid(fitX1$qr, X2)
    }
    Xj <- X2[,j]
    Xjm <- X2[,-j,drop=FALSE]
    fsReg <- lm.fit(Z2, Xjm)
    Xjmhat <- as.matrix(fsReg$fitted)
    fit <- lm.fit(Xjmhat, Xj)
    e <- Xj-c(Xjm%*%fit$coefficients)
    ehat <- qr.fitted(fsReg$qr, e)
    sig <- sum((e-ehat)^2)/(spec$n-spec$q)
    test <- sum(ehat^2)/sig/(ncol(Z2)-ncol(Xjm))
    if (!print)
        return(test)
    cat("Sanderson and Windmeijer Test for Weak Instruments\n")
    cat("***************************************************\n", sep="")
    add1 <- "(-1 for the critical values)"
    add2 <- paste("(-", ncol(X2)-1, " for the critical values)", sep="")
    cat("Number of included Endogenous: ", ncol(X2), add1, "\n", sep="")
    cat("Number of excluded Exogenous: ", ncol(Z2), add2, "\n", sep="")
    cat("The test is not robust to heteroskedasticity\n")    
    cat("Statistics: ", formatC(test, ...), "\n\n", sep="")
    SYTables(object, TRUE, TRUE)
}


## Montiel Olea and Pflueger (2013)

# computing x for the generalized test

getMOPx <- function(w, tau, type = c("TSLS", "LIML"), e=0.0001, nP = 10000,
                    maxi = 1000)
{
    type <- match.arg(type)
    fb <- function(beta, w, type)
    {
        s1 <- w$w1-2*beta*w$w12+beta^2*w$w2
        sig1 <- w$omega[1,1]-2*beta*w$omega[1,2]+beta^2*w$omega[2,2]
        s2 <- w$w2
        sig2 <- w$omega[2,2]
        s12 <- w$w12-beta*w$w2
        sig12 <- w$omega[1,2]-beta*w$omega[2,2]
        ts1 <- sum(diag(s1))
        ts2 <- sum(diag(s2))
        ts12 <- sum(diag(s12))
        if (type == "LIML")
        {
            tmp <-  2*s12-sig12*s1/sig1
            tmp <- 0.5*tmp + 0.5*t(tmp)
            lam <- eigen(tmp)$value[c(1,nrow(tmp))]
            num <- ts12 - ts1*sig12/sig1 - lam
            ne <- num/ts2
        } else {
            lam <- eigen(0.5*s12+0.5*t(s12))$value[c(1,nrow(s12))]
            num <- 1-2*lam/ts12
            ne <- num*ts12/ts2
        }
        BM <- sqrt(ts1/ts2)
        max(abs(ne))/BM
    }
    ew2 <- eigen(w$w2)$value
    maxf <- if (type == "LIML") max(ew2)/sum(ew2) else abs(1-2*min(ew2)/sum(ew2))
    b <- 1
    i <- 1
    while (TRUE)
    {
        crit <- min(abs(fb(b, w, type)/maxf - 1),
                    abs(fb(-b, w, type)/maxf - 1))
        if (crit <= e)
            break
        i <- i+1
        b <- b*2
        if (i>maxi)
        {
            warning("max iteration to find Brange reached")
            break
        }
    }
    res <- optimize(fb, c(-b,b), w=w, type=type, maximum=TRUE)
    Be <- res$objective
    Be/tau
}

MOPtest <- function(object, tau=0.10, size=0.05, print=TRUE,
                    estMethod = c("TSLS", "LIML"), simplified = TRUE,
                    digits = max(3L, getOption("digits") - 3L), ...)
{
    estMethod <- match.arg(estMethod)
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    spec <- modelDims(object)    
    if (sum(spec$isEndo)<1)
    {
        warning("The model does not contain endogenous variables")
        return(NA)
    }    
    if (sum(spec$isEndo)>1)
    {
        warning("The MOP test is defined for models with only one endogenous variable")
        return(NA)
    }
    Z2 <- model.matrix(object, "excludedExo")
    X1 <- model.matrix(object, "includedExo")
    X2 <- model.matrix(object, "includedEndo")
    y <- modelResponse(object)
    if (!is.null(X1))
    {
        fitX1  <- lm.fit(X1, Z2)
        Z2 <- as.matrix(fitX1$residuals)
        X2 <- qr.resid(fitX1$qr, X2)
        y <- qr.resid(fitX1$qr, y)
    }
    Z2 <- qr.Q(qr(Z2))*sqrt(nrow(Z2))
    colnames(Z2) <- paste("Z", 1:ncol(Z2), sep="")
    if ((ncol(Z2)-ncol(X2))<2)
        if (!simplified)
        {
            warning(paste("The generalized test is for models with 2 ",
                          "and more over-identifying restrictions:\n",
                          " simplified is changed to TRUE", sep=""))
            simplified <- TRUE
        }
    if (simplified)
    {
        if (estMethod == "LIML")
        {
            warning(paste("The simplified test is not defined for LIML\n",
                          "estMethod changed to TSLS", sep=""))
            estMethod <- "TSLS"
        }
        g <- reformulate(colnames(Z2), colnames(X2), FALSE)
        h <- reformulate(colnames(Z2), NULL, FALSE)
        dat <- data.frame(cbind(X2, Z2))
        m <- momentModel(g, h, data=dat, vcov=object@vcov,
                         vcovOptions=object@vcovOptions)
        b2 <- c(crossprod(Z2,X2))/nrow(Z2)
        w <- list(w2=vcov(m, b2))
        x <- 1/tau
    } else {
        b1 <- crossprod(Z2,y)/spec$n
        b2 <- crossprod(Z2,X2)/spec$n
        g <- list(reformulate(colnames(Z2), "y", FALSE),
                  reformulate(colnames(Z2), colnames(X2), FALSE))
        h <- reformulate(colnames(Z2), NULL, FALSE)
        dat <- as.data.frame(cbind(y=y,X2,Z2))
        m <- sysMomentModel(g=g, list(h), data = dat, vcov=object@vcov,
                            vcovOptions = object@vcovOptions)
        v <- cbind(y-Z2%*%b1, X2-Z2%*%b2)
        omega <- crossprod(v)/nrow(v)
        w <- vcov(m, list(b1,b2))
        w <- list(w1 = w[1:ncol(Z2), 1:ncol(Z2), drop=FALSE],
                  w2 = w[(ncol(Z2)+1):ncol(w), (ncol(Z2)+1):ncol(w), drop=FALSE],
                  w12 = w[1:ncol(Z2), (ncol(Z2)+1):ncol(w), drop=FALSE],
                  omega = crossprod(v)/nrow(v))
        x <- getMOPx(w, tau, estMethod, ...)
    }
    ev <- eigen(w$w2)$val
    se <- sum(ev)
    se2 <- sum(ev^2)
    me <- max(ev)
    ## Z'Y = Pi x n, so Y'Z'ZY = sum(Pi^2)*n^2
    ## there for Y'Z'ZY/se/n = sim(Pi^2)*n/se
    Feff <- sum(b2^2)/se*spec$n
    Keff <- se^2*(1+2*x)/(se2+2*x*se*me)
    crit <- qchisq(1-size, Keff, Keff*x)/Keff
    pv <- 1-pchisq(Feff*Keff, Keff, Keff*x)
    vcov <- object@vcov
    if (vcov=="MDS")
        vcov <- "HCCM"
    if (!print)
        return(c(Feff=Feff, Keff=Keff, x=x, critValue=crit, pvalue=pv))
    cat("Montiel Olea and Pflueger Test for Weak Instruments\n")
    cat("****************************************************\n", sep="")
    cat(ifelse(simplified, "Simplified Test", "Generalized Test"),
        " for ", estMethod, "\n", sep="")
    cat("Type of LS covariance matrix: ", vcov, "\n", sep="")
    cat("Number of included Endogenous: ", ncol(X2), "\n", sep="")
    cat("Effective degrees of freedom: ", Keff, "\n", sep="")
    cat("x: ", x, "\n", sep="")
    cat("Statistics: ", formatC(Feff, ...), "\n", sep="")
    cat(paste("Critical Value (size=",size,"): ", formatC(crit, digits=digits),
              "\n", sep=""))
    cat(paste("P-Value: ", formatC(pv, digits=digits), "\n\n", sep=""))    
    invisible()
}

getMOPw <- function(object)
{
    spec <- modelDims(object)
    if (sum(spec$isEndo)<1)
    {
        warning("The model does not contain endogenous variables")
        return(NA)
    }    
    if (sum(spec$isEndo)>1)
    {
        warning("The MOP test is defined for models with only one endogenous variable")
        return(NA)
    }
    Z2 <- model.matrix(object, "excludedExo")
    X1 <- model.matrix(object, "includedExo")
    X2 <- model.matrix(object, "includedEndo")
    y <- modelResponse(object)
    if (!is.null(X1))
    {
        fitX1  <- lm.fit(X1, Z2)
        Z2 <- as.matrix(fitX1$residuals)
        X2 <- qr.resid(fitX1$qr, X2)
        y <- qr.resid(fitX1$qr, y)
    }
    Z2 <- qr.Q(qr(Z2))*sqrt(nrow(Z2))
    colnames(Z2) = paste("Z", 1:ncol(Z2), sep="")    
    b <- c(b1 <- crossprod(Z2,y)/spec$n,
           b2 <- crossprod(Z2,X2)/spec$n)
    g <- list(reformulate(colnames(Z2), "y", FALSE),
              reformulate(colnames(Z2), colnames(X2), FALSE))
    h <- reformulate(colnames(Z2), NULL, FALSE)
    dat <- as.data.frame(cbind(y=y,X2,Z2))
    m <- sysMomentModel(g=g, list(h), data = dat, vcov=object@vcov,
                        vcovOptions = object@vcovOptions)
    v <- cbind(y-Z2%*%b1, X2-Z2%*%b2)
    omega <- crossprod(v)/nrow(v)
    w <- vcov(m, list(b1,b2))
    w1 <- w[1:ncol(Z2), 1:ncol(Z2), drop=FALSE]
    w2 <- w[(ncol(Z2)+1):ncol(w), (ncol(Z2)+1):ncol(w), drop=FALSE]
    w12 <- w[(ncol(Z2)+1):ncol(w), 1:ncol(Z2), drop=FALSE]
    list(w1=w1,w2=w2,w12=w12,omega=omega)
}


## Lewis and Mertens (2022)

phiMat <- function(w2, k2, l2)
{
    phi <- matrix(NA, k2,k2)
    if (length(w2) == 1)
        return(matrix(w2,1,1))
    for (i in 1:k2)
    {
        c0 <- 1+(i-1)*l2
        c1 <- ncol(w2)
        di <- diag(w2[,c0:c1])
        tmp <- colSums(matrix(di, nrow=l2))
        if (i == 1)
        {
            diag(phi) <- tmp
        } else if (i==k2) {
            phi[1,k2] <- phi[k2,1] <- tmp
        } else {
            j <- seq_len(length(tmp))
            j <- cbind(j,j)
            phi[-(1:(i-1)),][j] <- phi[,-(1:(i-1))][j] <- tmp
        }
    }
    phi
}

LewMertest <- function(object, tau=0.10, size=0.05, print=TRUE,
                       estMethod = c("TSLS", "LIML"), simplified = TRUE,
                       digits = max(3L, getOption("digits") - 3L),
                       dfCorr = TRUE, ...)
{
    estMethod <- match.arg(estMethod)
    if (!inherits(object, "linearModel"))
        stop("object must be of class linearModel")
    spec <- modelDims(object)    
    if (sum(spec$isEndo)<1)
    {
        warning("The model does not contain endogenous variables")
        return(NA)
    } 
    Z2 <- model.matrix(object, "excludedExo")
    X1 <- model.matrix(object, "includedExo")
    X2 <- model.matrix(object, "includedEndo")
    y <- modelResponse(object)
    if (!is.null(X1))
    {
        fitX1  <- lm.fit(X1, Z2)
        Z2 <- as.matrix(fitX1$residuals)
        X2 <- qr.resid(fitX1$qr, X2)
        y <- qr.resid(fitX1$qr, y)
    }
    Z2 <- qr.Q(qr(Z2))*sqrt(nrow(Z2))
    colnames(Z2) <- paste("Z", 1:ncol(Z2), sep="")
    colnames(X2) <- paste("endo",1:ncol(X2), sep="")
    b1 <- c(crossprod(Z2,y)/spec$n)
    Pi2 <- crossprod(Z2,X2)/spec$n
    g <- c(reformulate(colnames(Z2), "y", FALSE),
              lapply(colnames(X2), function(ei)
                  reformulate(colnames(Z2), ei, FALSE)))
    h <- reformulate(colnames(Z2), NULL, FALSE)
    dat <- as.data.frame(cbind(y=y,X2,Z2))
    m <- sysMomentModel(g=g, list(h), data = dat, vcov=object@vcov,
                        vcovOptions = object@vcovOptions)
    v <- cbind(y-Z2%*%b1, X2-Z2%*%Pi2)
    omega <- crossprod(v)/nrow(v)
    w <- vcov(m, c(list(b1),lapply(1:ncol(Pi2), function(i) Pi2[,i])))
    if (dfCorr)
        w <- w*nrow(Z2)/(nrow(Z2)-ncol(X1)-ncol(Z2))
    w2 <- w[(ncol(Z2)+1):ncol(w),(ncol(Z2)+1):ncol(w)]
    phi <- phiMat(w2, ncol(X2), ncol(Z2))
    if (dim(phi)[1] == 1)
    {
        gmin <- nrow(Z2)*c(crossprod(Pi2))/c(phi)
        sqPhi <- 1/sqrt(phi)
    } else {
        svPhi <- svd(phi)
        sqPhi <- svPhi$u%*%diag(1/sqrt(svPhi$d))%*%t(svPhi$v)
        gmin <- c(nrow(Z2)*min(eigen(sqPhi%*%crossprod(Pi2)%*%sqPhi)$value))
    }
    lam <- 1/tau
    cumul <- getCumul(w2, lam, sqPhi, ncol(Z2))    
    list(gmin=gmin, cumul=cumul,
         crit=crit(size, cumul$k1, cumul$k2, cumul$k3))
}

getCumul <- function(w2, lam, phi, l2)
{
    if (length(phi)==1)
    {
        sig <- w2/c(phi)*l2
    } else {
        tmp <- kronecker(phi, diag(l2))
        sig <- l2*(tmp%*%w2%*%tmp)
        svSig <- svd(sig)
        sig2 <- svSig$u%*%diag(svSig$d^2)%*%t(svSig$v)
        sig3 <- svSig$u%*%diag(svSig$d^3)%*%t(svSig$v)
    }
    k1 <- l2*(1+lam)
    if (length(sig) == 1)
    {
        k2 <- 2*(sig^2+2*lam*l2*sig)
        k3 <- 8*(sig^3+3*lam*l2*sig^2)
    } else {
    svSig <- svd(sig)
    sig2 <- svSig$u%*%diag(svSig$d^2)%*%t(svSig$v)
    sig3 <- svSig$u%*%diag(svSig$d^3)%*%t(svSig$v)
    tmp <- phiMat(sig2, dim(phi)[1], l2)
    k2 <- 2*(max(eigen(tmp)$val)+2*lam*l2*max(svSig$d))
    tmp <- phiMat(sig3, dim(phi)[1], l2)    
    k3 <- 8*(max(eigen(tmp)$val)+3*lam*l2*max(svSig$d)^2)
    }
    w <- k2/k3
    v <- 8*k2*w^2
    list(k1=k1, k2=k2, k3=k3, w=w, v=v)
}

crit <- function(alpha, k1, k2, k3)
{
    w <- k2/k3
    v <- 8*k2*w^2
    cr <- qchisq(1-alpha, v)
    (cr-v)/(4*w)+k1
}



