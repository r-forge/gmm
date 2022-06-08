splineMatrix <- function(X, knots = NA, pFact=0.3, deg=1, method=c("manual","bs"))
{
    method <- match.arg(method)
    n <- length(X)
    if (is.null(knots))
        return(as.matrix(X))
    if(any(is.na(knots)))
    {
        p <- floor(n^pFact)        
        prop.seq <- seq(from = 0, to = 1, length.out = p + 1)
        prop.seq <- prop.seq[-c(1, p + 1)]
        knots <- quantile(X, probs = prop.seq, type = 1)
    }
    if (method == "bs")
    {
        Xfi <- bs(x=X, knots=knots, degree=deg)
    } else {
        p <- length(knots) + 1
        Xfi <- matrix(nrow = n, ncol = p)
        Xfi[, 1] <- X * (X <= knots[1]) + knots[1] * (X > knots[1])
        Xfi[, p] <- (X - knots[p - 1]) * (X > knots[p - 1])
        if(p >= 3)
        {
            for(j in 2:(p - 1))
            {
                Xfi[, j] <- (X - knots[j - 1]) *
                    (X >= knots[j - 1]) * (X <= knots[j]) + 
                    (knots[j] - knots[j - 1]) * (X > knots[j])
            }
        }
        attr(Xfi, "knots") <- knots
    }
    Xfi
}

.getPval <- function(X, Y, Z, ppow, splineMet)
{
    n <- length(Y)
    id0 <- Z == 0
    id1 <- Z == 1
    Y0 <- Y[id0]
    Y1 <- Y[id1]
    X0 <- X[id0]
    X1 <- X[id1]

    Xfi0 <- splineMatrix(X=X0, pFact=ppow, method=splineMet)
    myknots0 <- attr(Xfi0, "knots")
    p0 <- ncol(Xfi0)
  
    Xfi1 <- splineMatrix(X=X1, pFact=ppow, method=splineMet)
    myknots1 <- attr(Xfi1, "knots")
    p1 <- ncol(Xfi1)

    p <- p0 + p1  

    Xf0 <- matrix(nrow = n, ncol = p0, 0)
    Xf1 <- matrix(nrow = n, ncol = p1, 0)
    Xf0[id0, ] <- Xfi0
    Xf1[id1, ] <- Xfi1

    Z <- factor(Z)
    pval0 <- rep(NA, p0 - 1)
    pval1 <- rep(NA, p1 - 1)
    lm.out0 <- lm(Y ~ 0 + factor(Z) + Xf0 + Xf1)
  
    for(j in 1 : (p0 - 1))
    {
        null_hyp <- paste("Xf0", j + 1, "-", "Xf0", j, "=0", sep = "")
        pval0[j] <- linearHypothesis(lm.out0, null_hyp, 
                                     vcov = vcovHC(lm.out0, type = "HC3"))[2, 4]
    }
  
    for(j in 1 : (p1 - 1))
    {
        null_hyp <- paste("Xf1", j + 1, "-", "Xf1", j, "=0", sep = "")
        pval1[j] <- linearHypothesis(lm.out0, null_hyp, 
                                     vcov = vcovHC(lm.out0, type = "HC3"))[2, 4]
    }  
    list(pval0=pval0, pval1=pval1, p0=p0, p1=p1, knots0=myknots0,
         knots1=myknots1)
}

selASY <- function(X, Y, Z, pFact=0.3, splineMet=c("manual","bs"))
{
    splineMet <- match.arg(splineMet)
    res <- .getPval(X, Y, Z, pFact, splineMet)
    pval=c(res$pval0, res$pval1)
    n <- length(X)
    q <- length(pval)
    p <- res$p0+res$p1
    id0 <- Z==0
    Jhat0 <- res$pval0 <= 1 / (p * log(p))
    Jhat1 <- res$pval1 <= 1 / (p * log(p))
    if(all(!Jhat0))
    {
        myknots0 <- NULL
    } else {
        myknots0 <- res$knots0[Jhat0]
    }
    Xfi0 <- splineMatrix(X=X[id0], knots=myknots0, method=splineMet)
    p00 <- ncol(Xfi0)
    Xf0 <- matrix(nrow = n, ncol = p00, 0)
    Xf0[id0, ] <- Xfi0
    
    if(all(!Jhat1))
    {
        myknots1 <- NULL
    } else {
        myknots1 <- res$knots1[Jhat1]
    }
    Xfi1 <- splineMatrix(X=X[!id0], knots=myknots1, method=splineMet)
    p10 <- ncol(Xfi1)
    Xf1 <- matrix(nrow = n, ncol = p10, 0)
    Xf1[!id0, ] <- Xfi1        
    list(Xf1=Xf1, Xf0=Xf0, knots0=myknots0, knots1=myknots1,
         pval=pval)
}

selIC <- function(X, Y, Z, pFact=0.3, type=c("AIC", "BIC"), splineMet=c("manual","bs"))
{
    type <- match.arg(type)
    splineMet <- match.arg(splineMet)    
    res <- .getPval(X, Y, Z, pFact, splineMet)
    pval=c(res$pval0, res$pval1)
    n <- length(X)    
    q <- length(pval)
    p <- res$p0+res$p1
    pval_sort <- sort(pval)
    Xf0 <- matrix(nrow = n, ncol = 1, 0)
    id0 <- Z==0
    Xf0[id0] <- X[id0]
    Xf1 <- matrix(nrow = n, ncol = 1, 0)
    Xf1[!id0] <- X[!id0]
    lm.out0 <- lm(Y ~ 0 + factor(Z) + Xf0 + Xf1)
    icV <- ic_seq0 <- get(type)(lm.out0)
    knots0 <- NULL
    knots1 <- NULL
    for(i in 1 : q)
    {    
        Jhat0i <- res$pval0 <= pval_sort[i]
        Jhat1i <- res$pval1 <= pval_sort[i]
        
        myknots0i <- if(all(!Jhat0i)) NULL else res$knots0[Jhat0i]
        Xfi0i <- splineMatrix(X=X[id0], knots=myknots0i, method=splineMet)
        p0i <- ncol(Xfi0i)
        Xf0i <- matrix(nrow = n, ncol = p0i, 0)
        Xf0i[id0, ] <- Xfi0i

        myknots1i <- if(all(!Jhat1i)) NULL else res$knots1[Jhat1i]
        Xfi1i <- splineMatrix(X=X[!id0], knots=myknots1i, method=splineMet)
        p1i <- ncol(Xfi1i)
        Xf1i <- matrix(nrow = n, ncol = p1i, 0)
        Xf1i[!id0, ] <- Xfi1i

        lm.out1 <- lm(Y ~ 0 + factor(Z) + Xf0i + Xf1i)
        ic_seq1 <- get(type)(lm.out1)
        icV <- c(icV, ic_seq1)
        if (ic_seq1<ic_seq0)
        {
            ic_seq0 <- ic_seq1
            Xf0 <- Xf0i
            Xf1 <- Xf1i
            knots0 <- myknots0i
            knots1 <- myknots1i
        } 
    }
    list(Xf1=Xf1, Xf0=Xf0, knots0=knots0, knots1=knots1, IC=icV,
         pval=pval)
}

.getCV <- function(Y, Z, X0, X1)
{
    n <- length(Y)
    X0 <- as.matrix(X0)
    X1 <- as.matrix(X1)
    myK <- floor(log(n))
    cv.outi <- cvFolds(n, K = myK)
    mspe_pred <- rep(NA, myK)    
    for(k in 1 : myK)
    {
        id.train <- cv.outi$subsets[cv.outi$which != k]
        id.valid <- cv.outi$subsets[cv.outi$which == k]
        train.datak <- list(Yk = Y[id.train], Zk = Z[id.train], 
                            Xf0ik = X0[id.train, , drop = FALSE],
                            Xf1ik = X1[id.train, , drop = FALSE])
        valid.datak <- list(Yk = Y[id.valid], Zk = Z[id.valid], 
                            Xf0ik = X0[id.valid, , drop = FALSE],
                            Xf1ik = X1[id.valid, , drop = FALSE])
        lm.outk <- lm(Yk ~ 0 + factor(Zk) + Xf0ik + Xf1ik, data = train.datak)
        pred.outk <- predict(lm.outk, newdata = valid.datak,
                             type = "response")
        mspe_pred[k] <- mean((valid.datak$Yk - pred.outk)^2,
                             na.rm = TRUE)
    }    
    mean(mspe_pred, na.rm = TRUE)
}

selCV <- function(X, Y, Z, pFact=0.3, splineMet=c("manual","bs"))
{
    splineMet <- match.arg(splineMet)    
    res <- .getPval(X, Y, Z, pFact, splineMet)
    pval=c(res$pval0, res$pval1)    
    n <- length(X)    
    q <- length(pval)
    p <- res$p0+res$p1
    myK <- floor(log(n))
    pval_sort <- sort(pval)
    mspe_seq <- rep(NA, q + 1)
    mspe_pred <- rep(NA, myK)
    cv.outi <- cvFolds(n, K = myK)
    id0 <- Z==0
    Xf0 <- matrix(nrow = n, ncol = 1, 0)
    Xf0[id0] <- X[id0]
    Xf1 <- matrix(nrow = n, ncol = 1, 0)
    Xf1[!id0] <- X[!id0]
    knots0 <- NULL
    knots1 <- NULL
    mspe0 <- mspe_seq[1] <- .getCV(Y, Z, Xf0, Xf1)
    for(i in 1 : q)
    {
        Jhat0i <- res$pval0 <= pval_sort[i]
        Jhat1i <- res$pval1 <= pval_sort[i]
        cv.outi <- cvFolds(n, K = myK)
        mspe_pred <- rep(NA, myK)

        myknots0i <- if(all(!Jhat0i)) NULL else res$knots0[Jhat0i]
        Xfi0i <- splineMatrix(X=X[id0], knots=myknots0i, method=splineMet)
        p0i <- ncol(Xfi0i)
        Xf0i <- matrix(nrow = n, ncol = p0i, 0)
        Xf0i[id0, ] <- Xfi0i

        myknots1i <- if(all(!Jhat1i)) NULL else res$knots1[Jhat1i]
        Xfi1i <- splineMatrix(X=X[!id0], knots=myknots1i, method=splineMet)
        p1i <- ncol(Xfi1i)
        Xf1i <- matrix(nrow = n, ncol = p1i, 0)
        Xf1i[!id0, ] <- Xfi1i
        mspe1 <- mspe_seq[i+1] <- .getCV(Y, Z, Xf0i, Xf1i)
        if (mspe1 < mspe0)
        {
            mspe0 <- mspe1
            knots0 <- myknots0i
            knots1 <- myknots1i
            Xf0 <- Xf0i
            Xf1 <- Xf1i
        }
    }
    list(Xf1=Xf1, Xf0=Xf0, knots0=knots0, knots1=knots1, IC=mspe_seq,
         pval=pval)
}

# currently implemented only for the case when X is univariate
otlse <- function(X, Y, Z, crit = c("ASY", "AIC", "BIC", "CV"),
                  pFact=0.3, splineMet=c("manual","bs"))
{
    crit <- match.arg(crit)
    splineMet <- match.arg(splineMet)
    optBasis <- switch(crit,
                       ASY = selASY(X, Y, Z, pFact, splineMet),
                       AIC = selIC(X, Y, Z, pFact, "AIC", splineMet),
                       BIC = selIC(X, Y, Z, pFact, "BIC", splineMet),
                       CV = selCV(X, Y, Z, pFact, splineMet))
    n <- length(Y)
    n1 <- sum(Z)
    n0 <- n-n1
    Xf0 <- optBasis$Xf0
    Xf1 <- optBasis$Xf1
    knots0 <- optBasis$knots0
    knots1 <- optBasis$knots1
    pval <- optBasis$pval
    p00 <- ncol(Xf0)
    p10 <- ncol(Xf1)
    lm.out <- lm(Y ~ 0 + factor(Z) + Xf0 + Xf1)
    vcov <- vcovHC(lm.out)
    idb0 <- 3 : (p00 + 2)
    idb1 <- (p00 + 3) : (p00 + p10 + 2)
    beta <- coef(lm.out)
    se.beta <- sqrt(diag(vcov))
    
    ## ACE
    
    X0 <- splineMatrix(X=X, knots=knots0, method=splineMet)
    X1 <- splineMatrix(X=X, knots=knots1, method=splineMet)
    Xbar0 <- apply(X0, 2, mean)
    Xbar1 <- apply(X1, 2, mean)
    vcovXf0 <- cov(X0)
    vcovXf1 <- cov(X1)
    ace <- c(beta[2] - beta[1] + crossprod(beta[idb1], Xbar1) - 
                                 crossprod(beta[idb0], Xbar0))
    
    Dvec <- rep(0, p00 + p10 + 2)
    Dvec[1] <- - 1
    Dvec[2] <- 1
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.ace <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta[idb0], crossprod(vcovXf0, beta[idb0])) / n + 
                 crossprod(beta[idb1], crossprod(vcovXf1, beta[idb1])) / n)^.5)

    ## ACT
    
    Xbar0 <- apply(X0[Z == 1, , drop = FALSE], 2, mean)
    Xbar1 <- apply(X1[Z == 1, , drop = FALSE], 2, mean)
    vcovXf0 <- cov(X0[Z == 1, , drop = FALSE])
    vcovXf1 <- cov(X1[Z == 1, , drop = FALSE])
    act <- c(beta[2] - beta[1] + crossprod(beta[idb1], Xbar1) - 
             crossprod(beta[idb0], Xbar0))
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.act <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta[idb0], crossprod(vcovXf0, beta[idb0])) / n1 + 
                 crossprod(beta[idb1], crossprod(vcovXf1, beta[idb1])) / n1)^.5)

    ## ACN
    
    Xbar0 <- apply(X0[Z == 0, , drop = FALSE], 2, mean)
    Xbar1 <- apply(X1[Z == 0, , drop = FALSE], 2, mean)
    vcovXf0 <- cov(X0[Z == 0, , drop = FALSE])
    vcovXf1 <- cov(X1[Z == 0, , drop = FALSE])
    acn <- c(beta[2] - beta[1] + crossprod(beta[idb1], Xbar1) - 
             crossprod(beta[idb0], Xbar0))
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.acn <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta[idb0], crossprod(vcovXf0, beta[idb0])) / n0 + 
                 crossprod(beta[idb1], crossprod(vcovXf1, beta[idb1])) / n0)^.5)
    
    ans <- list(beta = beta, se.beta = se.beta, 
                lm.out = lm.out, ace = ace, se.ace = se.ace, 
                act = act, se.act = se.act, acn = acn, se.acn = se.acn,
                myknots0 = knots0, myknots1 = knots1, pval = pval, crit=crit)
    class(ans) <- "otlse"
    ans
}

print.otlse <- function (x, ...)
{
    cat("Causal Effect using Optimal Thresholding Least Squares\n")
    cat("******************************************************\n")
    cat("Selection method: ", x$crit, "\n\n", sep="")
    cat("ACE = ", x$ace, "\nACT = ", x$act, "\nACN = ", x$acn,"\n")
}

summary.otlse <- function(object, ...)
{
    t <- c(object$ace, object$act, object$acn)/
        c(object$se.ace, object$se.act, object$se.acn)
    pv <- 2*pnorm(-abs(t))
    ace <- cbind(c(object$ace, object$act, object$acn),
                 c(object$se.ace, object$se.act, object$se.acn),
                 t, pv)
    dimnames(ace) <- list(c("ACE","ACT","ACN"),
                           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    nb0 <- ifelse(is.null(object$knots0), 1, length(object$knots0))
    nb1 <- ifelse(is.null(object$knots1), 1, length(object$knots1))
    t <- object$beta/object$se.beta
    pv <- 2*pnorm(-abs(t))
    beta <- cbind(object$beta, object$se.beta, t, pv)
    colnames(beta) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    ans <- list(causal=ace, beta=beta, crit=object$crit)
    class(ans) <- "summary.otlse"
    ans
}

print.summary.otlse <- function(x, digits = 4,
                                signif.stars = getOption("show.signif.stars"),
                                beta=FALSE, ...)
{
    cat("Causal Effect using Optimal Thresholding Least Squares\n")
    cat("******************************************************\n")
    cat("Selection method: ", x$crit, "\n\n", sep="")    
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    if (beta)
    {
        cat("Piecewise polynomials coefficients\n")
        cat("**********************************\n")
        printCoefmat(x$beta, digits = digits, signif.stars = signif.stars, 
                     na.print = "NA", ...)        
    }
}
        
