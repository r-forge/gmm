
## Moment functions

setGeneric("causalMomFct", function(theta, object, ...) standardGeneric("causalMomFct"))

setMethod("causalMomFct", signature("numeric", "causalData"),
          function(theta, object) {
              Z <- model.matrix(terms(object@reg), object@reg)
              Y <- model.response(object@reg)
              X <- model.matrix(terms(object@bal), object@bal)
              k <- ncol(Z)
              e <- Y-c(Z%*%theta[1:k])
              m1 <- e*Z
              e <- t(t(Z[,-1,drop=FALSE])-theta[-(1:k)])
              m2 <- sapply(1:ncol(X), function(i) e*X[,i])
              if (object@momType == "uncondBal")
                  return(cbind(m1,m2))
              m3 <- sweep(X[,-1,drop=FALSE], 2, object@balMom, "-")
              cbind(m1,m2,m3)
          })

## evalDMoment functions

setMethod("evalDMoment", signature("causalModel"),
          function(object, theta, impProb=NULL, augmented=FALSE) {
              dat <- object@X
              Z <- model.matrix(terms(dat@reg), dat@reg)
              X <- model.matrix(terms(dat@bal), dat@bal)
              k <- ncol(Z)
              n <- nrow(Z)
              ntet <- length(theta)
              if (is.null(impProb))
                  impProb <- rep(1/n, n)
              ZT <- c(Z%*%theta[1:k])
              q <- 2*k + (k-1)*(ncol(X)-1) - 1
              G <- matrix(0, q, ntet)
              G11 <- lapply(1:k, function(i) -colSums(impProb*Z[,i]*Z))
              G[1:k, 1:k] <- do.call(rbind, G11)
              G[(k+1):ntet, (k+1):ntet] <- -sum(impProb)*diag(k-1)
              uK <- colSums(impProb*X[,-1,drop=FALSE])
              G[(2*k):q, (k+1):ntet] <- -kronecker(diag(k-1), uK)
              if (dat@momType != "uncondBal")
              {
                      G <- rbind(G, matrix(0, ncol(X)-1, ntet))
                      if (augmented & dat@momType != "fixedMom")
                      {
                          ncov <- length(object@X@balCov)
                          q <- nrow(G)- ncov
                          tmp <- rbind(matrix(0, q, ncov),
                                       -sum(impProb)*diag(ncov))
                          G <- cbind(G, tmp)
                      }
                  }
              G
          })


## Print

setMethod("print", "causalModel",
          function(x, printBalCov=FALSE, ...) {
              cat("Causal Model \n")
              cat("*************\n")
              momType <- switch(x@X@momType,
                                uncondBal = "Unconditional balancing",
                                ACT = "Causal effect on the treated",
                                ACE = "Average causal effect",
                                ACC = "Causal effect on the control",
                                fixedMom = "Balancing based on fixed Moments")
              if (x@X@momType == "ACT" & x@X@ACTmom > 1)
                  momType <- paste(momType, "(treatment group ",
                                   x@X@ACTmom, ")")
              cat("Model type: ", momType, "\n", sep="")
              d <- modelDims(x)
              cat("Number of treatments: ", (d$k-1)/2, "\n", sep="")
              cat("Number of moment conditions: ", d$q, "\n", sep="")
              cat("Number of balancing covariates: ", length(x@X@balCov), "\n", sep="")
              cat("Sample size: ", d$n, "\n")
              if (printBalCov)
              {
                  cat("Balancing covariates:\n ")
                  bal <- x@X@balCov
                  while (length(bal))
                  {
                      cat("\t", paste(head(bal,3), collapse=", "), "\n", sep="")
                      bal <- bal[-(1:min(3, length(bal)))]
                  }
              }
              invisible()
          })

## gelFit

setMethod("gelFit", signature("causalModel"), valueClass="causalGelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL              
              res <- callNextMethod()
              res@call <- Call
              obj <- new("causalGelfit", res)
              obj
          })

## model.matrix and modelResponse

setMethod("model.matrix", signature("causalModel"),
          function(object, type=c("regressors","balancingCov"))
          {
              type <- match.arg(type)
              if (type == "regressors")
              {
                  ti <- attr(object@X@reg, "terms")
                  mat <- as.matrix(model.matrix(ti, object@X@reg)[,])
              } else {
                  ti <- attr(object@X@bal, "terms")
                  mat <- as.matrix(model.matrix(ti, object@X@bal)[,-1])
              }
              mat
          })

setMethod("modelResponse", signature("causalModel"),
          function(object)
          {
              model.response(object@X@reg)
          })


## Residuals
# Not sure we will need it, but the residuals are well defined in this case

setMethod("residuals", signature("causalModel"), function(object, theta){
    X <- model.matrix(object)
    Y <- modelResponse(object)
    e <- Y-c(X%*%theta[1:ncol(X)])
    e
})

## Dresiduals 
# Same comment as for residuals

setMethod("Dresiduals", signature("causalModel"),
          function(object, theta) {
              -model.matrix(object)
          })

## modelDims

setMethod("modelDims", "causalModel",
          function(object) {
              res <- callNextMethod()
              res$balCov <- object@X@balCov
              res$momType <- object@X@momType
              res$balMom <- object@X@balMom
              res$ACTmom <- object@X@ACTmom
              res
          })

## subset for observations selection

setMethod("subset", "causalModel",
          function(x, i) {
              x@X@reg <- x@X@reg[i,,drop=FALSE]
              x@X@bal <- x@X@bal[i,,drop=FALSE]
              x@n <- nrow(x@X@reg)
              x})


## "["
## balancing moment selection

setMethod("[", c("causalModel", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              balCov <- spec$balCov
              nbal <- length(balCov)
              if (!all(abs(i) %in% (1:nbal))) 
                  stop(paste("Sub-balancing must be between 1 and ", nbal, sep=""))
              balCov <- balCov[i]
              if (length(balCov)<1)
                  stop("The number of balancing covariates cannot be 0")              
              momInd <- c(matrix((spec$k+1):spec$q, nrow=nbal)[i,])
              momNames <- x@momNames[c(1:spec$k, momInd)]
              q <- length(momNames)
              f <- reformulate(balCov, NULL, TRUE)
              x@X@bal <- model.frame(f, x@X@bal)
              x@q <- q
              x@momNames <- momNames
              x@X@balCov <- balCov
              if (!is.null(x@X@balMom))
                  x@X@balMom <- x@X@balMom[i]
              x
          })

setMethod("[", c("causalModel", "numeric", "numericORlogical"),
          function(x, i, j){
              x <- x[i]
              subset(x, j)
          })

setMethod("[", c("causalModel", "missing", "numericORlogical"),
          function(x, i, j){
              subset(x, j)
          })

setMethod("[", c("causalGelfit", "numeric", "missing"),
          function(x, i, j){
              mod <- x@model[i]
              update(x, newModel=mod)
          })

setMethod("[", c("causalGelfit", "numeric", "numericORlogical"),
          function(x, i, j){
              mod <- x@model[i,j]
              update(x, newModel=mod)
          })

setMethod("[", c("causalGelfit", "missing", "numericORlogical"),
          function(x, i, j){
              mod <- x@model[,j]
              update(x, newModel=mod)
          })

