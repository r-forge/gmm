####### All methods with sysGmmModels (and its subclasses) signature
#####################################################################

## print and show

setMethod("print", "sysGmmModels",
          function(x, ...)
              {
                  cat("System GMM Model\n")
                  cat("****************\n")
                  type <- gsub("s", "", is(x)[1])
                  cat("Moment type: ", strsplit(type, "G")[[1]][1], "\n", 
                      sep = "")
                  cat("Covariance matrix: ", x@vcov, sep = "")
                  if (x@vcov == "HAC") {
                      cat(" with ", x@kernel, " kernel and ")
                      if (is.numeric(x@bw)) 
                          cat("Fixed  bandwidth (", round(x@bw, 3), ")", sep = "")
                      else cat(x@bw, " bandwidth", sep = "")
                  }
                  cat("\n")
                  d <- modelDims(x)
                  for (i in 1:length(d$eqnNames))
                      cat(d$eqnNames[i], ": coefs=", d$k[i],
                          ", moments=", d$q[i], ", number of Endogenous: ",
                          sum(d$isEndo[[i]]), "\n", sep="")
                  cat("Sample size: ", d$n, "\n")
              })

setMethod("show", "sysGmmModels", function(object) print(object))


## modelDims

setMethod("modelDims", "slinearGmm",
          function(object) {
              list(q=object@q, k=object@k, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, eqnNames=object@eqnNames,
                   isEndo=object@isEndo)
          })


setMethod("modelDims", "snonlinearGmm",
          function(object) {
              list(k=object@k, q=object@q, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fRHS=object@fRHS, fLHS=object@fLHS, eqnNames=object@eqnNames,
                   isEndo=object@isEndo)
          })
## Subsetting '['

setMethod("[", c("sysGmmModels", "missing", "list"),
          function(x, i, j){
              if (length(j) != length(x@q))
                  stop("j must be a list with a length equals to the number of equations")
              spec <- modelDims(x)
              x@SUR <- FALSE
              if (x@sameMom)
              {
                  chk <- sapply(j[-1], function(ji) identical(ji, j[[1]]))
                  if (!all(chk))
                      x@sameMom <- FALSE
              }
              for (s in 1:length(j))
                  {
                      if (length(j[[s]]) > 0)
                      {
                          q <- spec$q[s]
                          if (!all(abs(j[[s]]) %in% (1:q))) 
                              stop("SubMoment must be between 1 and q")
                          momNames <- spec$momNames[[s]][j[[s]]]
                          if (length(momNames)<spec$k[s])
                          {
                              error <- paste("Equation", s, "is under-identified")
                              stop(error)
                          }
                              if (momNames[1] == "(Intercept)")
                              {
                                  f <- reformulate(momNames[-1], NULL, TRUE)
                              } else {
                                  f <- reformulate(momNames, NULL, FALSE)
                                  }
                          attr(f, ".Environment")<- .GlobalEnv
                          x@q[s] <- length(momNames)
                          x@instT[[s]] <- terms(f)
                          x@momNames[[s]] <- momNames                              
                      }    
                  }
              x
          })

setMethod("[", c("snonlinearGmm", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              neqn <- length(spec$k)              
              if (!all(abs(i) %in% (1:neqn)))
                  stop("Selected equations out of range")
              x@fLHS <- x@fLHS[i]
              x@fRHS <- x@fRHS[i]
              if (length(x@fLHS) == 0)
                  stop("Removed too many equations; the model is empty")
              x@instT <- x@instT[i]
              x@k=x@k[i]
              x@q <- x@q[i]
              x@parNames <- x@parNames[i]
              x@momNames <- x@momNames[i]
              x@eqnNames <- x@eqnNames[i]
              x@theta0 <- x@theta0[i]
              x@varNames <- x@varNames[i]
              x@isEndo <- x@isEndo[i]
              x@SUR <- FALSE
              if (length(x@q) > 1)
                  return(x)
              instF <- model.frame(x@instT[[1]], x@data)
              varN <- c(all.vars(x@fLHS[[1]]), all.vars(x@fRHS))
              varN <- varN[!(varN%in%names(x@theta0[[1]]))]
              modelF <- x@data[,varN, drop=FALSE]
              new("nonlinearGmm", instF=instF, modelF=modelF, q=x@q[[1]],
                  fLHS=x@fLHS[[1]], fRHS=x@fRHS[[1]], theta0=x@theta0[[1]],
                  k=x@k[[1]], parNames=x@parNames[[1]], momNames=x@momNames[[1]],
                  vcov=x@vcov, n=spec$n, kernel=x@kernel, bw=x@bw, prewhite=x@prewhite,
                  ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov=x@centeredVcov, varNames=x@varNames[[1]],
                  isEndo=x@isEndo[[1]])
          })


setMethod("[", c("slinearGmm", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              neqn <- length(spec$k)              
              if (!all(abs(i) %in% (1:neqn)))
                  stop("Selected equations out of range")
              x@modelT <- x@modelT[i]
              if (length(x@modelT) == 0)
                  stop("Removed too many equations; the model is empty")
              x@instT <- x@instT[i]
              x@k=x@k[i]
              x@q <- x@q[i]
              x@parNames <- x@parNames[i]
              x@momNames <- x@momNames[i]
              x@eqnNames <- x@eqnNames[i]
              x@varNames <- x@varNames[i]
              x@isEndo <- x@isEndo[i]
              x@SUR <- FALSE
              if (length(x@q) > 1)
                  return(x)
              instF <- model.frame(x@instT[[1]], x@data)
              modelF <- model.frame(x@modelT[[1]], x@data)
              new("linearGmm", instF=instF, modelF=modelF, q=x@q[[1]],
                  k=x@k[[1]], parNames=x@parNames[[1]], momNames=x@momNames[[1]],
                  vcov=x@vcov, n=spec$n, kernel=x@kernel, bw=x@bw, prewhite=x@prewhite,
                  ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov=x@centeredVcov, varNames=x@varNames[[1]],
                  isEndo=x@isEndo[[1]])
          })

setMethod("[", c("sysGmmModels", "numeric", "list"),
          function(x, i, j){
              x <- x[i]              
              if (!inherits(x, "sysGmmModels"))
                  {
                      if (length(j)>1)
                          warning("length(j)>1, only the first element used")
                      x[,j[[1]]]
                  } else {
                      x[,j]
                  }
          })

setMethod("[", c("sysGmmModels", "missing", "missing"),
          function(x, i, j) x)

## Observation subset

setMethod("subset", "sysGmmModels",
          function(x, i) {
              x@data <- x@data[i,,drop=FALSE]
              x@n <- nrow(x@data)
              x})


### merge

setGeneric("merge")
setMethod("merge", c("linearGmm", "linearGmm"),
          function(x, y, ...) {
              all <- c(list(x), y, list(...))              
              cl <- sapply(all, class)
              if (any(cl != "linearGmm"))
                  stop("Can only merge linearGmm with other linearGmm models")
              n <- sapply(all, function(s) modelDims(s)$n)
              if (any(n[-1L] != n[1L]))
                  stop("You can only merge models with the same number of observations")
              k <- sapply(all, function(s) modelDims(s)$k)
              q <- sapply(all, function(s) modelDims(s)$q)
              parNames <- lapply(all, function(s) modelDims(s)$parNames)
              momNames <- lapply(all, function(s) modelDims(s)$momNames)
              varNames <- lapply(all, function(s) s@varNames)
              isEndo <- lapply(all, function(s) s@isEndo)
              instT <- lapply(all, function(s) terms(s@instF))
              modelT <- lapply(all, function(s) terms(s@modelF))
              dat <- do.call(cbind, lapply(all, function(s) cbind(s@modelF, s@instF)))
              dat <- dat[,!duplicated(colnames(dat))]
              eqnNames <- paste("Eqn", 1:length(all), sep="")
              new("slinearGmm", data=dat, instT=instT, modelT=modelT,
                  eqnNames=eqnNames, vcov=x@vcov, kernel=x@kernel, bw=x@bw,
                  prewhite=x@prewhite, ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov = x@centeredVcov, k=k, q=q, n=n[1], parNames=parNames,
                  momNames=momNames, sameMom=FALSE, isEndo=isEndo, varNames=varNames,
                  SUR=FALSE)              
          })

setMethod("merge", c("nonlinearGmm", "nonlinearGmm"),
          function(x, y, ...) {
              all <- c(list(x), y, list(...))
              cl <- sapply(all, class)
              if (any(cl != "nonlinearGmm"))
                  stop("Can only merge nonlinearGmm with oter nonlinearGmm models")
              n <- sapply(all, function(s) modelDims(s)$n)
              if (any(n[-1L] != n[1L]))
                  stop("You can only merge models with the same number of observations")
              fRHS <- lapply(all, function(s) s@fRHS)
              fLHS <- lapply(all, function(s) s@fLHS)
              k <- sapply(all, function(s) modelDims(s)$k)
              q <- sapply(all, function(s) modelDims(s)$q)
              parNames <- lapply(all, function(s) modelDims(s)$parNames)
              momNames <- lapply(all, function(s) modelDims(s)$momNames)
              varNames <- lapply(all, function(s) s@varNames)
              isEndo <- lapply(all, function(s) s@isEndo)              
              instT <- lapply(all, function(s) terms(s@instF))
              theta0 <- lapply(all, function(s) s@theta0)
              eqnNames <- paste("Eqn", 1:length(all), sep="")
              dat <- do.call(cbind, lapply(all, function(s) cbind(s@modelF, s@instF)))
              dat <- dat[,!duplicated(colnames(dat))]
              new("snonlinearGmm", data=dat, instT=instT,
                  theta0=theta0,fRHS=fRHS,eqnNames=eqnNames,
                  fLHS=fLHS, vcov=x@vcov, kernel=x@kernel, bw=x@bw,
                  prewhite=x@prewhite,
                  ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov = x@centeredVcov, k=k, q=q,
                  n=n[1], parNames=parNames, isEndo=isEndo, varNames=varNames, 
                  momNames=momNames, sameMom=FALSE, SUR=FALSE)
          })


setMethod("merge", c("snonlinearGmm", "nonlinearGmm"),
          function(x, y, ...) {
              all <- c(list(y), list(...))
              cl <- sapply(all, class)
              if (any(cl != "nonlinearGmm"))
                  stop("Can only merge nonlinearGmm with oter nonlinearGmm models")
              n <- sapply(all, function(s) modelDims(s)$n)
              spec <- modelDims(x)
              if (any(n != spec$n))
                  stop("You can only merge models with the same number of observations")
              fRHS <- c(spec$fRHS, lapply(all, function(s) modelDims(s)$fRHS))
              fLHS <- c(spec$fLHS, lapply(all, function(s) modelDims(s)$fLHS))
              k <- c(spec$k, sapply(all, function(s) modelDims(s)$k))
              q <- c(spec$q, sapply(all, function(s) modelDims(s)$q))
              parNames <- c(spec$parNames, lapply(all, function(s) modelDims(s)$parNames))
              momNames <- c(spec$momNames, lapply(all, function(s) modelDims(s)$momNames))
              varNames <- c(x@varNames, lapply(all, function(s) s@varNames))
              isEndo <- c(x@isEndo, lapply(all, function(s) s@isEndo))
              instT <- c(x@instT, lapply(all, function(s) terms(s@instF)))
              theta0 <- c(spec$theta0, lapply(all, function(s) modelDims(s)$theta0))
              eqNames <- x@eqnNames
              eqnNames <- c(eqNames, paste("Eqn",
                                           (length(eqNames)+1):length(fRHS), sep=""))
              dat <- do.call(cbind, lapply(all, function(s) cbind(s@modelF, s@instF)))
              dat <- dat[,!duplicated(colnames(dat))]
              new("snonlinearGmm", data=dat, instT=instT,
                  theta0=theta0,fRHS=fRHS,eqnNames=eqnNames,
                  fLHS=fLHS, vcov=x@vcov, kernel=x@kernel, bw=x@bw,
                  prewhite=x@prewhite,
                  ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov = x@centeredVcov, k=k, q=q,
                  n=n[1], parNames=parNames, isEndo=isEndo, varNames=varNames,
                  momNames=momNames, sameMom=FALSE, SUR=FALSE)
          })

setMethod("merge", c("slinearGmm", "linearGmm"),
          function(x, y, ...) {
              all <- c(list(y), list(...))              
              cl <- sapply(all, class)
              if (any(cl != "linearGmm"))
                  stop("Can only merge linearGmm with other linearGmm models")
              n <- sapply(all, function(s) modelDims(s)$n)
              spec <- modelDims(x)
              if (any(n != spec$n))
                  stop("You can only merge models with the same number of observations")
              k <- c(spec$k, sapply(all, function(s) modelDims(s)$k))
              q <- c(spec$q, sapply(all, function(s) modelDims(s)$q))
              parNames <- c(spec$parNames, lapply(all, function(s) modelDims(s)$parNames))
              momNames <- c(spec$momNames, lapply(all, function(s) modelDims(s)$momNames))
              varNames <- c(x@varNames, lapply(all, function(s) s@varNames))
              isEndo <- c(x@isEndo, lapply(all, function(s) s@isEndo))
              instT <- c(x@instT, lapply(all, function(s) terms(s@instF)))
              modelT <- c(x@modelT, lapply(all, function(s) terms(s@modelF)))
              dat <- do.call(cbind, lapply(all, function(s) cbind(s@modelF, s@instF)))
              dat <- dat[,!duplicated(colnames(dat))]
              eqNames <- x@eqnNames
              eqnNames <- c(eqNames, paste("Eqn",
                                           (length(eqNames)+1):length(instT), sep=""))
              new("slinearGmm", data=dat, instT=instT, modelT=modelT,
                  eqnNames=eqnNames, vcov=x@vcov, kernel=x@kernel, bw=x@bw,
                  prewhite=x@prewhite, ar.method=x@ar.method, approx=x@approx, tol=x@tol,
                  centeredVcov = x@centeredVcov, k=k, q=q, n=n[1], parNames=parNames,
                  momNames=momNames, sameMom=FALSE, isEndo=isEndo, varNames=varNames,
                  SUR=FALSE)              
          })

## residuals

setMethod("residuals", "sysGmmModels",
          function(object, theta) {
              neqn <- length(object@eqnNames)
              r <- sapply(1:neqn, function(i) residuals(object[i], theta[[i]]))
              colnames(r) <- object@eqnNames
              r
          })


## Dresiduals

setMethod("Dresiduals", "sysGmmModels",
          function(object, theta) {
              neqn <- length(object@eqnNames)
              if (missing(theta))
                  r <- sapply(1:neqn, function(i) Dresiduals(object[i]))
              else
                  r <- sapply(1:neqn, function(i) Dresiduals(object[i], theta[[i]]))
              names(r) <- object@eqnNames
              r
          })


## evalMoment

setMethod("evalMoment", "sysGmmModels",
          function(object, theta) {
              neqn <- length(object@eqnNames)
              gt <- lapply(1:neqn, function(i) evalMoment(object[i], theta[[i]]))
              names(gt) <- object@eqnNames
              gt
          })

## evalDmoments

setMethod("evalDMoment", "sysGmmModels",
          function(object, theta) {
              neqn <- length(object@eqnNames)
              if (missing(theta))
                  dgt <- lapply(1:neqn, function(i) evalDMoment(object[i]))
              else
                  dgt <- lapply(1:neqn, function(i) evalDMoment(object[i], theta[[i]]))
              names(dgt) <- object@eqnNames
              dgt
          })

## model.matrix


setMethod("model.matrix", "slinearGmm",
          function(object, type =  c("regressors", "instruments")) {
              type <- match.arg(type)
              mm <- lapply(1:length(object@eqnNames), function(i)
                  model.matrix(object[i], type))
              names(mm) <- object@eqnNames
              mm
          })


setMethod("model.matrix", "snonlinearGmm",
          function(object, type =  c("regressors", "instruments"))
          {
              type <- match.arg(type)
              if (type == "regressors") 
                  stop("no model.matrix of type regressors for nonlinear Gmm. set type to 'instruments' to get the matrix of instruments")
              mm <- lapply(1:length(object@eqnNames), function(i)
                  model.matrix(object[i], type))
              names(mm) <- object@eqnNames
              mm
          })



## evalWeights

### The following multiplies each block matrix of ZZ, when the dimensions
### are defined by dimr and dimc (dim row and dim column) by each element
### of Sigma. Sigma must be length(dimr) by length(dimc)

.SigmaZZ <- function(ZZ, Sigma, dimr, dimc=NULL, lowerTri=FALSE, isSym=TRUE)
    {
        r1 <- 1        
        c1 <- 1
        if (is.null(dimc))
            dimc <- dimr
        for (i in 1:length(dimr))
            {
                r2 <- sum(dimr[1:i])
                start <- ifelse(isSym, i, 1)
                for (j in start:length(dimc))
                    {
                        c2 <- sum(dimc[1:j])
                        ZZ[r1:r2, c1:c2] <- ZZ[r1:r2, c1:c2]*Sigma[i,j]
                        c1 <- c1+dimc[j]
                    }
                r1 <- r1+dimr[i]
                c1 <- ifelse(isSym, sum(dimc[1:i])+1, 1)
            }                
        if (lowerTri && isSym)
            ZZ[lower.tri(ZZ)] <- t(ZZ)[lower.tri(ZZ)]
        ZZ
    }


setMethod("evalWeights", "sysGmmModels",
          function(object, theta = NULL, w="optimal", wObj=NULL)
          {
              spec <- modelDims(object)
              sameMom <- object@sameMom
              n <- object@n
              neqn <- length(spec$eqnNames)
              if (is.matrix(w))
              {
                  if (!all(dim(w) == sum(spec$q)))
                      stop("The weights matrix has the wrong dimension")
              }
              if (is.matrix(w) || (w == "ident"))
                  {
                      return(new("sysGmmWeights", type="weights",momNames=object@momNames,
                                 HAC=list(), w=w, Sigma=NULL, sameMom=sameMom,
                                 eqnNames=object@eqnNames))
              }
              if (w != "optimal")
                  stop("w is either 'ident', 'optimal' or a matrix")
              if (object@vcov == "iid")
              {
                  e <- residuals(object, theta)
                  type <- "iid"
                  Sigma <- chol(crossprod(e)/n)
                  if (!is.null(wObj))
                      {
                          if (wObj@type != "iid")
                              stop("wObj must come from a model with iid errors")
                          if (ncol(wObj@w$qr) != spec$q[1])
                              stop("The qr decomposition has the wrong dimension")
                          w <- wObj@w
                      } else {
                          Z <- model.matrix(object, type="instruments")
                          if (sameMom)
                              {
                                  w <- qr(Z[[1]]/sqrt(n))
                              } else {
                                  w <- crossprod(do.call(cbind,Z))/n
                              }
                      }
              } else if (object@vcov == "MDS") {
                  type <- "MDS"
                  gt <- evalMoment(object, theta)
                  w <- qr(do.call(cbind, gt)/sqrt(n))
                  Sigma <- NULL
              } else {
                  stop("Only identity, iid and MDS if allowed for now")
              }
              return(new("sysGmmWeights", type=type,momNames=object@momNames,
                         HAC=list(), w=w, Sigma=Sigma, sameMom=sameMom,
                         eqnNames=object@eqnNames))
          })

## evalObjective

setMethod("evalObjective", signature("sysGmmModels", "list", "sysGmmWeights"),
          function(object, theta, wObj, ...)
              {
                  gt <- evalMoment(object, theta)
                  gt <- lapply(gt, function(g) colMeans(g))
                  gt <- do.call("c", gt)
                  n <- object@n
                  obj <- quadra(wObj, gt)
                  n*obj
              })

## modelResponse

setMethod("modelResponse", signature("slinearGmm"),
          function(object)
          {
              neqn <- length(object@eqnNames)
              Y <- lapply(1:neqn, function(i) modelResponse(object[i]))
              Y
              })

## solveGMM

.GListToMat <- function(G)
{
    dimG <- sapply(G, dim)
    Gmat <- matrix(0, sum(dimG[1,]), sum(dimG[2,]))
    r1 <- 1
    c1 <- 1
    for (i in 1:length(G)) {
        r2 <- sum(dimG[1,1:i])
        c2 <- sum(dimG[2,1:i])
        Gmat[r1:r2, c1:c2] <- as.matrix(G[[i]])
        r1 <- r1 + dimG[1,i]
        c1 <- sum(dimG[2,1:i]) + 1
    }
    Gmat
}

.tetReshape <- function(theta, eqnNames, parNames)
{
    if (is.list(theta))
    {
        theta2 <- do.call("c", theta)
        k <- sapply(parNames, length)
        tn <- paste(rep(eqnNames, k), ".", do.call("c", parNames), 
                    sep = "")
        names(theta2) <- tn
    } else {
        k <- cumsum(sapply(parNames, length))
        names(theta) <- do.call("c", parNames)
        theta2 <- list(theta[1:k[1]])
        if (length(k)>1)
            theta2[2:length(k)] <- lapply(2:length(k),
                                          function(i) theta[(k[i-1]+1):k[i]])
        names(theta2) <- eqnNames
    }
    theta2
}


setMethod("solveGmm", c("slinearGmm", "sysGmmWeights"),
          function(object, wObj, theta0 = NULL) {
              if (wObj@type=="iid" && object@sameMom)
                  return(ThreeSLS(object, Sigma=wObj@Sigma, qrZ=wObj@w, coefOnly=TRUE))
              spec <- modelDims(object)              
              Y <- modelResponse(object)
              Z <- model.matrix(object, type="instruments")
              Syz <- lapply(1:length(Y),
                            function(i) colMeans(Y[[i]]*Z[[i]]))
              Syz <- do.call("c", Syz)
              G <- evalDMoment(object)
              G <- .GListToMat(G)
              T1 <- quadra(wObj, G)
              T2 <- quadra(wObj, G, Syz)
              theta <- -solve(T1, T2)
              theta <- .tetReshape(theta, object@eqnNames, object@parNames)
              list(theta=theta, convergence=NULL)
          })


setMethod("solveGmm", signature("snonlinearGmm", "sysGmmWeights"),
          function (object, wObj, theta0 = NULL, ...) 
          {
              if (is.null(theta0)) 
                  theta0 <- modelDims(object)$theta0
              if (!is.list(theta0))
                  stop("theta0 must be a list of vectors")
              g <- function(theta, wObj, object){
                  spec <- modelDims(object)
                  theta <- .tetReshape(theta, object@eqnNames, spec$parNames)
                  evalObjective(object, theta, wObj)
              }
              dg <- function(theta, wObj, object) {
                  spec <- modelDims(object)
                  theta <- .tetReshape(theta, object@eqnNames, spec$parNames)
                  gt <- evalMoment(object, theta)
                  gt <- do.call(cbind, gt)
                  n <- nrow(gt)
                  gt <- colMeans(gt)
                  G <- evalDMoment(object, theta)
                  G <- .GListToMat(G)
                  obj <- 2 * n * quadra(wObj, G, gt)
                  obj
              }
              spec <- modelDims(object)
              theta0 <- .tetReshape(theta0, object@eqnNames, spec$parNames)
              res <- optim(par = theta0, fn = g, gr = dg, method = "BFGS", 
                           object = object, wObj = wObj, ...)
              theta <- .tetReshape(res$par, object@eqnNames, object@parNames)
              list(theta = theta, convergence = res$convergence)
    })

## momentVcov

setMethod("momentVcov", signature("sysGmmModels"),
          function(object, theta, ...){
              spec <- modelDims(object)
              q <- spec$q
              if (object@vcov == "MDS")
                  {
                      gt <- evalMoment(object, theta)
                      gt <- do.call(cbind, gt)
                      if (object@centeredVcov)
                          gt <- scale(gt, scale=FALSE)
                      w <- crossprod(gt)/nrow(gt)
                  } else if (object@vcov == "iid") {
                      e <- residuals(object, theta)
                      Sigma <- crossprod(e)/nrow(e)
                      Z <- model.matrix(object, "instrument")
                      if (object@sameMom)
                          {
                              w <- kronecker(Sigma, crossprod(Z[[1]])/nrow(e))
                          } else {
                              Z <- crossprod(Z)/nrow(e)
                              w <- .SigmaZZ(Z, Sigma, q)
                          }
                  } else {
                      stop("not yet implemented for HAC")
                  }
              wn <- paste(rep(spec$eqnNames, q), ".", do.call("c", spec$momNames), 
                          sep = "")
              dimnames(w) <- list(wn,wn)
              w
          })

### tsls

setMethod("tsls", "slinearGmm",
          function(object)
          {
              call <- match.call()
              neqn <- length(object@eqnNames)
              res <- lapply(1:neqn, function(i) tsls(object[i]))
              w <- lapply(1:neqn, function(i) quadra(res[[i]]@wObj))
              w <- .GListToMat(w)
              wObj <- evalWeights(object, w=w)
              theta <- lapply(res, coef)
              names(theta) <- object@eqnNames
              new("stsls", theta=theta, convergence=NULL, convIter=NULL,
                  call=call, type="tsls", wObj=wObj, niter=1L,
                  efficientGmm=FALSE, model=object)
          })


### 3SLS and SUR

setGeneric("ThreeSLS", function(object, ...) standardGeneric("ThreeSLS"))

setMethod("ThreeSLS", "slinearGmm", 
          function(object, coefOnly=FALSE, qrZ=NULL, Sigma=NULL) {
              Call <- match.call()
              if (!inherits(object, "slinearGmm"))
                  stop("3SLS is for slinearGmm classes")
              if (!object@sameMom)
                  stop("For 3SLS, the instruments must be the same in each equation")
              efficientGmm <- object@vcov == "iid"
              spec <- modelDims(object)
              n <- spec$n
              neqn <- length(object@eqnNames)
              x <- model.matrix(object)
              y <- modelResponse(object)
              if (is.null(qrZ))
                  {
                      z <- model.matrix(object[1], "instruments")
                      qrZ <- qr(z/sqrt(n))
                  } else {
                      if (class(qrZ) != "qr")
                          stop("qrZ must be the qr decomposition of Z")
                      if (ncol(qrZ$qr) != length(spec$momNames[[1]]))
                          stop("The qr decomposition has the wrong dimension")
                  }
              if (is.null(Sigma))
                  {
                      theta <- lapply(1:neqn, function(i) lm.fit(qr.fitted(qrZ, x[[i]]),
                                                                 y[[i]])$coefficients)
                      e <- residuals(object, theta)
                      Sigma <- chol(crossprod(e)/n)
                  } else {
                      if (!all(Sigma[lower.tri(Sigma)] == 0))
                          stop("Sigma must be a cholesky and therefore upper tiangular")
                  }
              if (object@SUR)
                  {
                      xhat <- do.call(cbind, x)
                      type <- "SUR"
                  } else {
                      xhat <- qr.fitted(qrZ, do.call(cbind, x))
                      type <- "3SLS"
                  }
              iSigma <- chol2inv(Sigma)
              y <- do.call(cbind,y)
              A <- crossprod(xhat)
              A <- .SigmaZZ(A, iSigma, spec$k, spec$k, TRUE)
              C <- crossprod(xhat, y)
              C <- .SigmaZZ(C, iSigma, spec$k, rep(1, neqn), FALSE, FALSE)
              C <- rowSums(C)
              theta <- .tetReshape(solve(A, C), object@eqnNames, spec$parNames)
              if (coefOnly)
                  return(list(theta=theta, convergence=NULL))
              wObj <- new("sysGmmWeights", w=qrZ, Sigma=Sigma, type="iid",
                          momNames=spec$momNames,
                          HAC=list(), sameMom=TRUE, eqnNames=object@eqnNames)
              new("sgmmfit", theta=theta, convergence=NULL,
                  convIter=rep(NULL, neqn), call=Call, type=type, wObj=wObj,
                  niter=2L, efficientGmm=efficientGmm,  model=object)
          })


## modelFit

setMethod("modelFit", signature("sysGmmModels"), valueClass="sgmmfit", 
          function(object, type=c("twostep", "iter","cue", "onestep"),
                   itertol=1e-7, initW=c("ident", "tsls", "EbyE"), weights="optimal", 
                   itermaxit=100, efficientWeights=FALSE, start=NULL,
                   EbyE=FALSE, ...)
          {
              Call <- match.call()
              chk <- validObject(object)                  
              type <- match.arg(type)
              initW <- match.arg(initW)
              i <- 1L
              chk <- validObject(object, TRUE)
              if (!chk)
                  stop("object is not a valid gmmModels object")
              if (is.character(weights) && !(weights%in%c("optimal","ident")))
                  stop("weights is a matrix or one of 'optimal' or 'ident'")
              spec <- modelDims(object)
              if (all(spec$q==spec$k))
                  {
                      weights <- "ident"
                      type <- "onestep"
                      EbyE <- TRUE
                  } else if (type == "onestep" && !is.matrix(weights)) {
                      weights <- "ident"
                      EbyE <- TRUE
                  } else if (is.matrix(weights) || class(weights)=="sysGmmWeights") {
                      type <- "onestep"
                      EbyE <- FALSE
                  } else if (weights == "ident") {
                      type <- "onestep"
                      EbyE <- TRUE
                  }
              if (EbyE)
                  {
                      neqn <- length(object@eqnNames)
                      res <- lapply(1:neqn, function(i)
                          modelFit(object[i], type=type, weights=weights,itertol=itertol,
                                 initW=initW, itermaxit=itermaxit,
                                 efficientWeights=efficientWeights, start=start, ...))
                      theta <- lapply(res, coef)
                      convergence <- sapply(res, function(r) r@convergence)
                      if (is.list(convergence))
                          convergence <- do.call("c", convergence)
                      convIter <- sapply(res, function(r) r@convIter)
                      niter <- sapply(res, function(r) r@niter)
                      if (is.list(convIter))
                          convIter <- do.call("c", convIter)
                      type <- paste("EBE", type, sep="")
                      efficientGmm <- FALSE
                      if (is.character(weights) && weights=="ident")
                          {
                              wObj <- evalWeights(object, NULL, "ident")
                          } else {
                              wObj <- lapply(res, function(r) quadra(r@wObj))
                              wObj <- .GListToMat(wObj)
                              wObj <- evalWeights(object, w=wObj)
                          }
                      ans <- new("sgmmfit", theta=theta, convergence=convergence,
                                 convIter=convIter, call=Call, type=type, wObj=wObj,
                                 niter=niter, efficientGmm=efficientGmm, model=object)
                      return(ans)
                  }              
              if (type == "onestep")
                  {
                      if (class(weights)=="sysGmmWeights")
                          wObj <- weights
                      else
                          wObj <- evalWeights(object, w=weights)
                      res <- solveGmm(object, wObj, start, ...)
                      convergence <- res$convergence
                      efficientGmm <- efficientWeights
                      ans <- new("sgmmfit", theta=res$theta,
                                 convergence=convergence, convIter=NULL, type=type,
                                 wObj=wObj, model=object, call=Call, niter=i,
                                 efficientGmm=efficientGmm)
                      return(ans)
                  }
              if (type == "twostep")
                  {
                      itermaxit <- 1
                      if (inherits(object, "slinearGmm"))
                          {
                              if (object@vcov=="iid" && !object@sameMom && !object@SUR)
                                  type <- "FIVE"                              
                              if (initW=="tsls" && object@vcov=="iid" && object@sameMom)
                                  return(ThreeSLS(object))
                          }
                  }
              if (initW=="tsls")
                  {                          
                      theta0 <- try(coef(tsls(object)), silent=TRUE)
                      if (class(theta0)=="try-error")
                          stop("Cannot get the initial weights using 2SLS")
                  } else if (initW == "EbyE") {
                      neqn <- length(object@eqnNames)
                      res <- lapply(1:neqn, function(i)
                          modelFit(object[i], type=type, weights=weights,itertol=itertol,
                                 itermaxit=itermaxit,
                                 efficientWeights=efficientWeights, start=start, ...))
                      theta0 <- lapply(res, coef)
                  } else {
                      wObj <- evalWeights(object, NULL, "ident")
                      theta0 <- solveGmm(object, wObj, start, ...)$theta
                  }
              bw <- object@bw
              if (type != "cue")
                  {
                      while(TRUE)
                          {
                              if (i>1 && object@vcov=="iid")
                                  wObj0 <- wObj
                              else
                                  wObj0 <- NULL
                              wObj <- evalWeights(object, theta0, "optimal", wObj0)
                              if (is.character(object@bw) && object@vcov=="HAC")
                                  object@bw <- wObj@HAC$bw
                              res <- solveGmm(object, wObj, theta0, ...)
                              theta1 <- res$theta
                              convergence <- res$convergence
                              tet0 <- do.call("c", theta0)
                              dif1 <- do.call("c", theta1)-tet0
                              crit <- sqrt(sum(dif1^2))/(1+sqrt(sum(tet0^2)))
                              if (crit < itertol & type=="iter")
                                  {
                                      convIter <- 0
                                      break
                                  }
                              i <- i + 1L
                              theta0 <- theta1
                              if (i>itermaxit)
                                  {
                                      if (type %in% c("twostep", "FIVE"))
                                          convIter <- NULL
                                      else
                                          convIter <- 1
                                      break                                      
                                  }                              
                          }      
                  } else {
                      convIter <- NULL
                      if (is.character(object@bw) && object@vcov=="HAC")
                          {
                              w <- momentVcov(object, theta0)
                              object@bw <- bw
                          }
                      if (object@vcov == "iid")
                          wObj0 <- evalWeights(object, theta0)
                      else
                          wObj0 <- NULL
                      obj <- function(theta, object, wObj0, spec)
                          {
                              theta <- .tetReshape(theta, object@eqnNames,
                                                          spec$parNames)
                              wObj <- evalWeights(object, theta, "optimal", wObj0)
                              evalObjective(object, theta, wObj)
                          }
                      res <- optim(do.call("c",theta0), obj, object=object, wObj0=wObj0,
                                   spec=spec, ...)
                      theta1 <- .tetReshape(res$par, object@eqnNames,spec$parNames)
                      convergence <- res$convergence
                      wObj <- evalWeights(object, theta1, "optimal", wObj0)
                  }
              object@bw <- bw
              new("sgmmfit", theta=theta1, convergence=convergence,
                  convIter=convIter, call=Call, type=type, wObj=wObj,
                  niter=i, efficientGmm=TRUE, model=object)
          })


