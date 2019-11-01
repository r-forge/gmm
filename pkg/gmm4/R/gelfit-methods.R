####  All methods for gelfit class
#####################################

### Hidden functions

.invTest <- function(object, which, level = 0.95, fact = 3,
                     type=c("LR", "LM", "J"), corr=NULL, ...)
{
    type <- match.arg(type)
    if (length(which) > 1)
        stop("tests are inverted only for one parameter")
    spec <- modelDims(object@model)
    df <- spec$q-spec$k
    if (df > 0)
    {
        test0 <- c(specTest(object, type=type, ...)@test)
        test0 <- test0[1]
    } else {
        test0 <- 0
    }
    v <- diag(vcov(object, ...)$vcov_par)
    sdcoef <- sqrt(v[which])
    coef <- coef(object)[which]
    int1 <- c(coef, coef + fact*sdcoef)
    int2 <- c(coef - fact*sdcoef, coef)
    if (length(coef(object)) == 2)
    {
        sd1 <- sqrt(v[-which])
        coef1 <- coef(object)[-which]
        rang <- c(coef1 - 2*fact*sd1, coef + 2*fact*sd1)
    } else {
        rang <- NULL
    }
    fct <- function(coef, which, type, fit, level, test0, corr=NULL, rang)
    {
        spec <- modelDims(fit@model)
        ncoef <- spec$parNames[which]
        R <- paste(ncoef, "=", coef)
        model <- restModel(fit@model, R)
        if (length(coef(fit))==2)
        {
            tControl <- list(method="Brent", lower=rang[1], upper=rang[2])
            fit2 <- suppressWarnings(update(fit, newModel=model, tControl=tControl,
                                            start.tet=coef(fit)[-which]))
        } else {
            fit2 <- suppressWarnings(update(fit, newModel=model,
                                            start.tet=coef(fit)[-which]))
        }
        test <- specTest(fit2, type=type, ...)@test[1] - test0
        if (is.null(corr))
            level - pchisq(test, 1)
        else
            level - pchisq(test/corr, 1)
    }
    res1 <- try(uniroot(fct, int1, which = which, type=type, level=level,
                        fit=object, test0=test0, corr=corr, rang=rang),
                silent=TRUE)
    res2 <- try(uniroot(fct, int2, which = which, type=type, level=level,
                        fit=object, test0=test0, corr=corr, rang=rang),
                silent=TRUE)
    if (any(c(class(res1), class(res2)) == "try-error"))
    {
        test <- c(NA,NA)
        mess <- "Could not compute the confidence interval because: \n"
        if (class(res1) == "try-error")
            mess <- paste(mess, "(1) ", res1[1], "\n", sep="")
        if (class(res2) == "try-error")
            mess <- paste(mess, "(2) ", res2[1], "\n", sep="")
        warning(mess)        
    } else {
        test <- sort(c(res1$root, res2$root))
    }
    test
}
        
## coef

setMethod("coef", "gelfit", function(object) object@theta)

## print

setMethod("print", "gelfit",
          function(x, model=TRUE, lambda=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              type <- x@type
              spec <- modelDims(x@model)
              if (spec$q==spec$k && x@type != "eval")
                  type <- paste("Just-Identified ", type, sep="")
              cat("\nEstimation: ", type,"\n")
              cat("Convergence Theta: ", x@convergence, "\n")
              cat("Convergence Lambda: ", x@lconvergence, "\n")              
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
              if (lambda)
                  {
                      cat("lambdas:\n")
                      print.default(format(x@lambda, ...), print.gap=2L, quote=FALSE)
                  }
          })

## show

setMethod("show","gelfit", function(object) print(object))

## residuals

setMethod("residuals", "gelfit", function(object) {
    residuals(object@model, object@theta)})

## getImpProb

setGeneric("getImpProb", function(object, ...) standardGeneric("getImpProb"))

setMethod("getImpProb", "gelfit",
          function(object, posProb=FALSE, normalize=TRUE) {
              rhoFct <- object@model@gelType
              if (is.null(rhoFct$fct))
                  rhoFct <- get(paste("rho", rhoFct$name, sep=""))
              else
                  rhoFct <- rhoFct$fct
              gt <- evalMoment(object@model, object@theta)
              k <- object@model@wSpec$k
              pt <- -rhoFct(gt, object@lambda, 1, k[1]/k[2])/nrow(gt)
              if (object@model@gelType$name == "EEL"  && posProb) {
                  eps <- -length(pt) * min(min(pt), 0)
                  pt <- (pt + eps/length(pt))/(1 + eps)
              }
              if (normalize)
                  pt <- pt/sum(pt)
              convMom <- colSums(pt * gt)
              convProb <- abs(sum(as.numeric(pt))-1)
              list(pt=pt, convMom=convMom, convProb=convProb)
          })


### To be removed once the above has need tested enough
#setMethod("getImpProb", "gelfit",
#          function(object) {
#              rhoFct <- object@model@gelType
#              if (is.null(rhoFct$fct))
#                  rhoFct <- get(paste("rho", rhoFct$name, sep=""))
#              else
#                  rhoFct <- rhoFct$fct
#              gt <- evalMoment(object@model, object@theta)
#              k <- object@model@wSpec$k
#              pt <- -rhoFct(gt, object@lambda, 1, k[1]/k[2])/nrow(gt)
#              if (object@model@gelType$name == "EEL") {
#                  eps <- -length(pt) * min(min(pt), 0)
#                  pt <- (pt + eps/length(pt))/(1 + eps)
#              }
#              convMom <- colSums(pt * gt)
#              convProb <- abs(sum(as.numeric(pt))-1)
#              pt <- pt/sum(pt)
#              list(pt=pt, convMom=convMom, convProb=convProb)
#          })

## vcov

setMethod("vcov", "gelfit",
          function(object, withImpProb=FALSE, tol=1e-10) {
              spec <- modelDims(object@model)
              q <- spec$q
              gt <- evalMoment(object@model, object@theta)
              n <- nrow(gt)
              bw <- object@model@wSpec$bw
              k <- object@model@wSpec$k
              if (withImpProb)
                  {
                      pt <- getImpProb(object)$pt
                      G <- evalDMoment(object@model, object@theta, pt)
                      G <- G/k[1]
                      gt <- gt * sqrt(pt * bw/k[2])
                  } else {
                      G <- evalDMoment(object@model, object@theta)
                      G <- G/k[1]
                      gt <- gt * sqrt(bw/k[2]/n)
                  }
              qrGt <- qr(gt)
              piv <- qrGt$pivot
              R <- qr.R(qrGt)
              X <- forwardsolve(t(R), G[piv,])
              Y <- forwardsolve(t(R), diag(q)[piv,])
              res <- lm.fit(as.matrix(X), Y)
              u <- res$residuals
              Sigma <- chol2inv(res$qr$qr)/n
              diag(Sigma)[diag(Sigma) < 0] <- tol
              if (q == ncol(G)) {
                  SigmaLam <- matrix(0, q, q)
              } else {
                  SigmaLam <- crossprod(Y, u)/n * bw^2
                  diag(SigmaLam)[diag(SigmaLam) < 0] <- tol
              }
              piv <- sort.int(piv, index.return = TRUE)$ix
              list(vcov_par = Sigma, vcov_lambda = SigmaLam)
          })

## Summary


setMethod("summary","gelfit",
          function (object, ...) 
              {
                  if (length(object@vcov) == 0)
                      v <- vcov(object, ...)
                  else
                      v <- object@vcov
                  se.t <- sqrt(diag(v$vcov_par))
                  se.l <- sqrt(diag(v$vcov_lambda))
                  theta <- object@theta
                  lambda <- object@lambda
                  tval.t <- theta/se.t
                  tval.l <- lambda/se.l
                  coef <- cbind(theta, se.t, tval.t,
                                2*pnorm(abs(tval.t), lower.tail = FALSE))
                  coefl <- cbind(lambda, se.l, tval.l,
                                 2*pnorm(abs(tval.l), lower.tail = FALSE))
                  stest <- specTest(object)
                  dimnames(coef) <- list(names(theta), c("Estimate", "Std. Error", 
                                                         "t value", "Pr(>|t|)"))
                  dimnames(coefl) <- list(names(lambda), c("Estimate", "Std. Error", 
                                                           "t value", "Pr(>|t|)"))
                  pt <- getImpProb(object)
                      
                  ans <- new("summaryGel", coef = coef, specTest = stest,
                             model = object@model, lambda=coefl,
                             convergence=object@convergence,
                             lconvergence=object@lconvergence, impProb=pt)
                  ans})

## confint

setMethod("confint", "gelfit",
          function (object, parm, level = 0.95, lambda = FALSE,
                    type = c("Wald", "invLR", "invLM", "invJ"),
                    fact = 3, corr = NULL, vcov=NULL, ...) 
          {
              type <- match.arg(type)
              spec <- modelDims(object@model)
              n <- spec$n
              theta <- coef(object)
              if (lambda)
              {
                  lam <- object@lambda
                  if (missing(parm))
                      parm <- 1:length(lam)                      
                  nlam <- spec$momNames
                  if (is.character(parm))
                      parm <- sort(unique(match(parm, nlam)))
                  nlam <- nlam[parm]
                  if (length(theta) == length(lam))
                  {
                      ntest <- paste("No confidence intervals for lambda",
                                     "when the model is just identified.")
                      ans <- matrix(NA, length(nlam), 2)
                  } else {
                      ntest <- "Wald confidence interval for Lambda"
                      if (is.null(vcov))
                          v <- vcov(object, ...)$vcov_lambda
                      se <- sqrt(diag(v))
                      if (missing(parm))
                          parm <- 1:length(lam)
                      se <- se[parm]
                      lam <- lam[parm]                  
                      zs <- qnorm((1 - level)/2, lower.tail = FALSE)              
                      ch <- zs * se
                      ans <- cbind(lam-ch, lam+ch)
                  }
                  dimnames(ans) <- list(nlam,
                                        c((1 - level)/2, 0.5 + level/2))
                  return(new("confint", interval=matrix,
                             type=ntest, level=level))                  
              }
              if (type == "Wald")
              {
                  if (is.null(vcov))
                      v <-  vcov(object, ...)
                  return(getMethod("confint", "gmmfit")(object, parm, level,
                      vcov=v$vcov_par))
              } else {
                  if (missing(parm)) 
                      parm <- 1:length(theta)
                  ntheta <- spec$parNames
                  if (is.character(parm))
                      parm <- sort(unique(match(parm, ntheta)))
                  ntheta <- ntheta[parm]
                  type <- strsplit(type, "v")[[1]][2]
                  ntest <- paste("Confidence interval based on the inversion of the ", 
                                 type, " test", sep = "")
                  ans <- lapply(parm, function(w)
                      .invTest(object, w, level = level, 
                               fact = fact, type = type, corr = corr, ...))
                  ans <- do.call(rbind, ans)
                  dimnames(ans) <- list(ntheta, c((1 - level)/2, 0.5 + level/2))
              }
              new("confint", interval=ans, type=ntest, level=level)
          })

## specTest

setMethod("specTest", signature("gelfit", "missing"),
          function(object, which, type=c("All", "LR", "LM", "J"))
          {
              type <- match.arg(type)
              spec <- modelDims(object@model)
              q <- spec$q
              n <- spec$n
              df <- q-spec$k              
              test <- numeric()
              if (type %in% c("All","LR"))
              {
                  LR <- evalObjective(object@model, object@theta, lambda=object@lambda)
                  test <- c(test, LR)
                  names(test) <- "LR: "
              }
              if (type %in% c("All","LM","J"))
                  gt <- evalMoment(object@model, object@theta)
              if (type %in% c("All","LM"))
              {
                  kHat <- crossprod(gt)/n
                  LM <- n * crossprod(object@lambda, crossprod(kHat, object@lambda))/
                      (object@model@wSpec$bw^2)
                  test <- c(test, LM)
                  names(test)[length(test)] <- "LM: "                  
              }
              if (type %in% c("All","J"))
              {
                  J <- sum(lm.fit(gt, rep(1,n))$fitted.values)
                  test <- c(test, J)
                  names(test)[length(test)] <- " J: "                  
              }
              if (df == 0)
                  pv <- NA
              else
                  pv <- 1-pchisq(test, df)
              test <- cbind(test, df, pv)
              colnames(test) <- c("Statistics", "df", "pvalue")
              new("specTest", test=test, testname="Test E(g)=0")
          })

setMethod("print", "summaryGel",
          function(x, digits=5, lambda=TRUE, ...)
          {
              print(x@model)
              cat("Convergence Theta: ", x@convergence, "\n", sep="")
              cat("Convergence Lambda: ", x@lconvergence, "\n", sep="")
              cat("Average |Sum of pt*gt()]|: ", format(mean(abs(x@impProb$convMom)),
                                                        digits=5), "\n", sep="")
              cat("|Sum of pt - 1|: ", format(mean(abs(x@impProb$convProb)),
                                              digits=5), "\n", sep="")
              
              cat("\ncoefficients:\n")
              printCoefmat(x@coef, digits=digits, ...)
              if (lambda)
                  {
                      cat("\nLambdas:\n")
                      printCoefmat(x@lambda, digits=digits, ...)
                  }
              print(x@specTest)
          })


## show
setMethod("show", "summaryGel", function(object) print(object)) 
    

## update    


setMethod("update", "gelfit",
          function(object, newModel=NULL, ..., evaluate=TRUE)
          {
              if (is.null(call <- getCall(object)))
                  stop("No call argument")
              arg <- list(...)
              model <- if(is.null(newModel))
                           object@model
                       else
                           newModel
              model <- update(model, ...)
              ev <- new.env(parent.frame())
              ev[["model"]] <- model
              call[["object"]] <- quote(model)
              arg <- arg[which(is.na(match(names(arg),
                                           c("rhoFct", slotNames(model)))))]
              if (length(arg) > 0) 
                  for (n in names(arg)) call[[n]] <- arg[[n]]
              if (evaluate)
                  eval(call, ev)
              else
                  call
          })
              


   
