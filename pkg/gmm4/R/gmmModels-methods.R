####### All methods with gmmModels (and its subclasses) signature
#################################################################

#######################  Print ########################
### The getGeneric for print is here only, so the file must be compiled
### before any other files containing print

setGeneric("print")
setMethod("print", "gmmModels",
          function(x, ...) {
              cat("GMM Model\n")
              cat("*********\n")
              cat("Moment type: ", strsplit(is(x)[1], "G")[[1]][1], "\n", sep="")
              cat("Covariance matrix: ", x@vcov, sep="")
              if (x@vcov == "HAC")
                  {
                      cat(" with ", x@kernel, " kernel and ")
                      if (is.numeric(x@bw))
                          cat("Fixed  bandwidth (", round(x@bw,3), ")",  sep="")
                      else
                          cat(x@bw, " bandwidth",  sep="")
                  }
              cat("\n")
              d <- modelDims(x)
              cat("Number of regressors: ", d$k, "\n", sep="")
              cat("Number of moment conditions: ", d$q, "\n", sep="")
              if (!inherits(x, "functionGmm"))
                  cat("Number of Endogenous Variables: ", sum(x@isEndo), "\n", sep="")
              cat("Sample size: ", d$n, "\n")})             

setMethod("show", "gmmModels", function(object) print(object))

##### coef  ########
### For this, it only attach the names to theta

setMethod("coef", "gmmModels",
          function(object, theta) {
              names(theta) <- object@parNames
              theta})

################## model.matrix and modelResponse #################
### I did not make model.response as generic because it is not
### a method in stats and I want different arguments

setGeneric("modelResponse", function(object, ...) standardGeneric("modelResponse"))

setMethod("modelResponse", signature("linearGmm"),
          function(object)
          {
              model.response(object@modelF)
          })

setGeneric("model.matrix")
setMethod("model.matrix", signature("linearGmm"),
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              if (type == "regressors")
              {
                  ti <- attr(object@modelF, "terms")
                  mat <- model.matrix(ti, object@modelF)[,]
              } else {
                  ti <- attr(object@instF, "terms")
                  mat <- model.matrix(ti, object@instF)[,]
              }
              mat
          })

setMethod("model.matrix", signature("nonlinearGmm"),
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              if (type == "regressors")
              {
                  stop("no model.matrix of type regressors for nonlinear Gmm. set type to 'instruments' to get the matrix of instruments")
              } else {
                  ti <- attr(object@instF, "terms")
                  mat <- model.matrix(ti, object@instF)[,]
              }
              mat
          })

#################### residuals ###########################
### The getGeneric for residuals is here only, so the file must be compiled
### before any other files containing print


### rlinearGmm inherits from linearGmm
setGeneric("residuals")
setMethod("residuals", signature("linearGmm"), function(object, theta){
    X <- model.matrix(object)
    Y <- modelResponse(object)
    e <- Y-c(X%*%theta)
    e
})

setMethod("residuals", signature("nonlinearGmm"), 
          function(object, theta)
              {
                  res <- modelDims(object)
                  nt <- names(theta)
                  nt0 <- names(res$theta0)
                  if (length(theta) != length(nt0))
                      stop("The length of theta is not equal to the number of parameters")
                  if (is.null(nt))
                      stop("theta must be a named vector")
                  if (!all(nt%in%nt0 & nt0%in%nt))
                      stop("names in theta dont match parameter names")
                  varList <- c(as.list(theta), as.list(object@modelF))
                  if (!is.null(res$fLHS))
                      {
                          lhs <- try(eval(res$fLHS, varList))
                          if (any(class(lhs)=="try-error"))
                              stop("Cannot evaluate the LHS")
                      } else {
                          lhs <- 0
                      }
                  rhs <- try(eval(res$fRHS, varList))
                  if (any(class(rhs)=="try-error"))
                      stop("Cannot evaluate the RHS")
                  c(lhs-rhs)
              })

################ evalMoment ##########################

setGeneric("evalMoment", function(object, theta, ...) standardGeneric("evalMoment"))

setMethod("evalMoment", signature("regGmm"),
          function(object, theta) {
              e <- residuals(object, theta)
              Z <- model.matrix(object, "instruments")
              as.matrix(Z)*e
          })

setMethod("evalMoment", signature("functionGmm"),
          function(object, theta) {
              spec <- modelDims(object)
              spec$fct(theta, object@X)
          })

setMethod("evalMoment", signature("formulaGmm"),
          function(object, theta) {
              res <- modelDims(object)
              nt <- names(theta)
              nt0 <- names(res$theta0)
              if (length(theta) != length(nt0))
                  stop("The length of theta is not equal to the number of parameters")
              if (is.null(nt))
                  stop("theta must be a named vector")
              if (!all(nt%in%nt0 & nt0%in%nt))
                  stop("names in theta dont match parameter names")
              varList <- c(as.list(theta), as.list(object@modelF))
              sapply(1:res$q, function(i) {
                  if (!is.null(res$fLHS[[i]]))
                      {
                          lhs <- try(eval(res$fLHS[[i]], varList))
                          if (any(class(lhs)=="try-error"))
                              stop("Cannot evaluate the LHS")
                      } else {
                          lhs <- 0
                      }
                  if (!is.null(res$fRHS[[i]]))
                      {
                          rhs <- try(eval(res$fRHS[[i]], varList))
                          if (any(class(lhs)=="try-error"))
                              stop("Cannot evaluate the RHS")
                      } else {
                          lhs <- 0
                      }
                  c(lhs-rhs)})
          })


################ evalDresiduals ##########################

setGeneric("Dresiduals", function(object, theta, ...) standardGeneric("Dresiduals"))

setMethod("Dresiduals", signature("linearGmm"),
          function(object, theta) {
              -model.matrix(object)
          })

setMethod("Dresiduals", signature("nonlinearGmm"),
          function(object, theta) {
              res <- modelDims(object)
              nt <- names(theta)
              nt0 <- names(res$theta0)
              if (length(theta) != length(nt0))
                  stop("The length of theta is not equal to the number of parameters")
              if (is.null(nt))
                  stop("theta must be a named vector")
              if (!all(nt%in%nt0 & nt0%in%nt))
                  stop("names in theta dont match parameter names")
              varList <- c(as.list(theta), as.list(object@modelF))
              De <- numeric()
              for (i in nt)
                  {
                      if (!is.null(res$fLHS))
                          d <- eval(D(res$fLHS, i), varList)      
                      else
                          d <- 0
                      De <-  cbind(De, d-matrix(eval(D(res$fRHS, i), varList),res$n,1))
                  }
              De
          })

############### modelDims #######################

setGeneric("modelDims", function(object, ...) standardGeneric("modelDims"))

setMethod("modelDims", "linearGmm",
          function(object) {
              list(k=object@k, q=object@q, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, isEndo=object@isEndo)
          })

setMethod("modelDims", "nonlinearGmm",
          function(object) {
              list(k=object@k, q=object@q, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fRHS=object@fRHS, fLHS=object@fLHS, isEndo=object@isEndo)
          })

setMethod("modelDims", "functionGmm",
          function(object) {
              list(k=object@k, q=object@q, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fct=object@fct, dfct=object@dfct, isEndo=object@isEndo)
          })

setMethod("modelDims", "formulaGmm",
          function(object) {
              list(k=object@k, q=object@q, n=object@n, parNames=object@parNames,
                   momNames=object@momNames, theta0=object@theta0,
                   fRHS=object@fRHS, fLHS=object@fLHS, isEndo=object@isEndo,
                   isMDE=object@isMDE)
          })

################ evalDMoment ##########################

setGeneric("evalDMoment", function(object, ...) standardGeneric("evalDMoment"))

setMethod("evalDMoment", signature("regGmm"),
          function(object, theta) {
              De <- Dresiduals(object, theta)
              Z <- model.matrix(object, "instrument")
              G <- apply(De,2, function(x) colMeans(Z*x))
              d <- modelDims(object)
              dimnames(G) <- list(d$momNames, d$parNames)
              G
          })

setMethod("evalDMoment", signature("functionGmm"),
          function(object, theta) {
              spec <- modelDims(object)
              if (is.null(spec$dfct))
              {
                  f <- function(theta, object)
                  {
                      gt <- evalMoment(object, theta)
                      colMeans(gt)
                  }
                  env <- new.env()
                  assign("theta", theta, envir=env)
                  assign("object", object, envir=env)
                  assign("f", f, envir=env)
                  G <- numericDeriv(quote(f(theta, object)), "theta", env)
                  G <- attr(G, "gradient")
              } else {
                  G <- spec$dfct(theta, object@X)
              }
              dimnames(G) <- list(spec$momNames, spec$parNames)
              G
              })

setMethod("evalDMoment", signature("formulaGmm"),
          function(object, theta) {
              res <- modelDims(object)
              nt <- names(theta)
              nt0 <- names(res$theta0)
              if (length(theta) != length(nt0))
                  stop("The length of theta is not equal to the number of parameters")
              if (is.null(nt))
                  stop("theta must be a named vector")
              if (!all(nt%in%nt0 & nt0%in%nt))
                  stop("names in theta dont match parameter names")
              varList <- c(as.list(theta), as.list(object@modelF))
              G <- numeric()
              for (i in nt)
                  {
                      lhs <- sapply(1:res$q, function(j) {
                          if (!is.null(res$fLHS[[j]]))
                              d <- mean(eval(D(res$fLHS[[j]], i), varList))
                          else
                              d <- 0
                          c(d)})
                      rhs <- sapply(1:res$q, function(j) {
                          if (!is.null(res$fRHS[[j]]))
                              d <- mean(eval(D(res$fRHS[[j]], i), varList))
                          else
                              d <- 0
                          c(d)})
                      G <- cbind(G, lhs-rhs)
                  }
              G
          })

###########   estfun :  Don't like it ###############

estfun.gmmFct <- function(x, ...) x


##########  vcovHAC from sandwich #################

setMethod("vcovHAC", "gmmModels",
          function (x, theta) { 
              gmat <- evalMoment(x, theta)
              if (x@centeredVcov) 
                  gmat <- scale(gmat, scale = FALSE)
              class(gmat) <- "gmmFct"
              if (is.character(x@bw))
              {
                  if (x@bw == "Wilhelm")
                  {
                      G <- evalDMoment(x, theta)
                      obj <- list(gt = gmat, G = G)
                      class(obj) <- "gmm"
                  } else {
                      obj <- gmat
                  }
                  bw <- get(paste("bw",x@bw,sep=""))(obj, order.by = NULL, 
                      kernel = x@kernel, prewhite = x@prewhite, 
                      ar.method = x@ar.method, approx = x@approx)
              } else {
                  bw <- x@bw
              }
              weights <- weightsAndrews(x = gmat, bw = bw, kernel = x@kernel, 
                                        prewhite = x@prewhite, tol = x@tol)
              w <- vcovHAC(x = gmat, order.by = NULL, weights = weights, 
                           prewhite = x@prewhite, sandwich = FALSE,
                           ar.method = x@ar.method)
              attr(w, "Spec") <- list(weights = weights, bw = bw, kernel = x@kernel)
              w
          })

################ momentVcov  ##########################

setGeneric("momentVcov", function(object, ...) standardGeneric("momentVcov"))

setMethod("momentVcov", signature("gmmModels"),
          function(object, theta, ...){
              if (class(object) == "functionGmm" & object@vcov == "iid")
                  object@vcov <- "MDS"
              if (object@vcov == "MDS")
                  {
                      gt <- evalMoment(object, theta)
                      if (object@centeredVcov)
                          gt <- scale(gt, scale=FALSE)
                      w <- crossprod(gt)/nrow(gt)
                  } else if (object@vcov == "iid") {
                      sig <- sd(residuals(object, theta))
                      Z <- model.matrix(object, "instrument")
                      w <- sig^2*crossprod(Z)/nrow(Z)
                  } else {
                      w <- vcovHAC(object, theta)
                      attr(w, "inv") <- TRUE
                  }
              w})

################### weights Object and methods: Is it too much??? #################


setGeneric("evalWeights", function(object, ...) standardGeneric("evalWeights"))

setMethod("evalWeights", signature("gmmModels"),valueClass="gmmWeights",
          function(object, theta=NULL, w="optimal", ...)
              {
                  hac <- list()
                  if (is.matrix(w))
                      {
                          type <- "weights"
                      } else {
                          if (w == "ident")
                              {
                                  type <- "weights"
                              } else {
                                  if (class(object) == "functionGmm" & object@vcov == "iid")
                                      object@vcov <- "MDS"
                                  if (object@vcov == "MDS")
                                      {
                                          gt <- evalMoment(object, theta)
                                          if (object@centeredVcov)
                                              gt <- scale(gt, scale=FALSE)
                                          w <- qr(gt/sqrt(nrow(gt)))
                                          if (w$rank < object@q)
                                              warning("The moment matrix is not full column rank")
                                          type <- "qr"
                                      } else if (object@vcov == "iid") {
                                          sig <- mean(residuals(object, theta)^2)
                                          sig <- sqrt(sig)
                                          ti <- attr(object@instF, "terms")
                                          Z <- model.matrix(ti, object@instF)[,]
                                          w <- qr(sig*Z/sqrt(nrow(Z)))
                                          if (w$rank < object@q)
                                              warning("The moment matrix is not full column rank")
                                          type <- "qr"
                                      } else {
                                          w <- vcovHAC(object, theta)
                                          hac <- attr(w,"Spec")
                                          w <- chol(w[,])
                                          type <- "chol"
                                      }
                              }
                      }
                  new("gmmWeights", type=type, w=w, HAC=hac)
              })

############ evalObjective #################################

setGeneric("evalObjective", function(object, theta, wObj, ...) standardGeneric("evalObjective"))

setMethod("evalObjective", signature("gmmModels", "numeric", "gmmWeights"),
          function(object, theta, wObj, ...)
              {
                  gt <- evalMoment(object, theta)
                  n <- modelDims(object)$n
                  gt <- colMeans(gt)
                  obj <- quadra(wObj, gt)
                  n*obj
              })

#########################  solveGmm  #########################

setGeneric("solveGmm", function(object, wObj, ...) standardGeneric("solveGmm"))

setMethod("solveGmm", signature("linearGmm", "gmmWeights"),
          function(object, wObj, theta0=NULL, ...)
          {
              X <- model.matrix(object)
              Z <- model.matrix(object, "instrument")
              Y <- modelResponse(object)
              d <- modelDims(object)
              n <- d$n
              Sig.zy <- crossprod(Z,Y)/n
              Sig.zx <- crossprod(Z,X)/n
              if (d$q == d$k)
              {
                  T1 <- Sig.zx
                  T2 <- Sig.zy
              } else {
                  T1 <- quadra(wObj, Sig.zx)
                  T2 <- quadra(wObj, Sig.zx, Sig.zy)
              }
              theta <- c(solve(T1, T2))
              names(theta) <- d$parNames
              list(theta=theta, convergence=NULL)
          })

setMethod("solveGmm", signature("allNLGmm", "gmmWeights"),
          function(object, wObj, theta0=NULL, algo=c("optim","nlminb"), ...)
              {
                  algo <- match.arg(algo)
                  if (is.null(theta0))
                      theta0 <- modelDims(object)$theta0
                  g <- function(theta, wObj, object)
                      evalObjective(object, theta, wObj)
                  dg <- function(theta, wObj, object)
                      {
                          gt <- evalMoment(object, theta)
                          n <- nrow(gt)
                          gt <- colMeans(gt)
                          G <- evalDMoment(object, theta)
                          obj <- 2*n*quadra(wObj, G, gt)
                          obj
                      }
                  if (algo == "optim")
                      {
                          if ("method" %in% names(list(...)))
                              res <- optim(par=theta0, fn=g, gr=dg, 
                                           object=object, wObj=wObj, ...)
                          else
                              res <- optim(par=theta0, fn=g, gr=dg, method="BFGS",
                                           object=object, wObj=wObj, ...)
                      } else {
                          res <- nlminb(start=theta0, objective=g, gradient=dg,
                                        object=object, wObj=wObj, ...)
                      }
                  theta <- res$par
                  names(theta) <- modelDims(object)$parNames
                  list(theta=theta, convergence=res$convergence)
              })


##################### momentStrength ####################

setGeneric("momentStrength", function(object, ...) standardGeneric("momentStrength"))

setMethod("momentStrength", signature("nonlinearGmm"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("functionGmm"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("formulaGmm"), 
          function(object, theta=NULL, ...)
          {
              list(strength=NULL, mess=NULL)
          })

setMethod("momentStrength", signature("linearGmm"), 
          function(object, theta, vcovType=c("OLS","HC","HAC")){
              spec <- modelDims(object)
              EndoVars <- !(spec$parNames %in% spec$momNames)
              exoInst <- spec$momNames %in% spec$parNames
              vcovType <- match.arg(vcovType)
              if (all(!EndoVars))
                  {
                      fstats <- NULL
                      mess <- "No endogenous variables: no strength measure"  
                  } else {
                      X <- model.matrix(object)
                      X <- X[,EndoVars,drop=FALSE]
                      Z <- model.matrix(object, "instrument")
                      fstats <- matrix(ncol=0, nrow=3)
                      df1 <- sum(!exoInst)                      
                      for (i in 1:ncol(X))
                          {
                              resu <- lm(X[,i   ]~Z-1)                              
                              v <- switch(vcovType,
                                         OLS=vcov(resu),
                                         HC=vcovHC(resu,"HC1"),
                                         HAC=vcovHAC(resu))[!exoInst,!exoInst]
                              b <- coef(resu)[!exoInst]
                              f <- b%*%solve(v, b)/df1
                              df2 <- resu$df.residual
                              fstats <- cbind(fstats, c(f, df1, df2))
                          }
                      fstats <- rbind(fstats, 1-pf(fstats[1,], fstats[2,],fstats[3,]))
                      colnames(fstats) <- colnames(X)
                      rownames(fstats) <- c("Stats","df1","df2","pv")
                      fstats <- t(fstats)
                      mess <- "Instrument strength based on the F-Statistics of the first stage OLS"
                  }
              list(strength=fstats, mess=mess)
          })

### Subsetting models

setGeneric("[")
setMethod("[", c("regGmm", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
              momNames <- x@momNames[i]
              if (length(momNames)<spec$k)
                  stop("The model is under-identified")
              if (momNames[1] == "(Intercept)") 
                  f <- reformulate(momNames[-1], NULL, TRUE)
              else 
                  f <- reformulate(momNames, NULL, FALSE)
              instF <- model.frame(f, x@instF)
              x@q <- length(momNames)
              x@instF <- instF
              x@momNames <- momNames
              x
          })


setMethod("[", c("functionGmm", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
               if (length(i)==q)
                  return(x)
               momNames <- x@momNames[i]
               if (length(momNames)<spec$k)
                   stop("The model is under-identified")
               attr(x@X, "subset") <- i
               attr(x@X, "fct") <- spec$fct
               x@fct <- function(theta, x)
                   attr(x, "fct")(theta,x)[,attr(x, "subset")]
               if (!is.null(x@dfct))
                   {
                       attr(x@X, "dfct") <- spec$dfct
                       x@dfct <- function(theta, x)
                           attr(x, "dfct")(theta,x)[,attr(x, "subset")]
                   }
               x@q <- length(momNames)
               x@momNames <- momNames
               x
          })

setMethod("[", c("formulaGmm", "numeric", "missing"),
          function(x, i, j){
              i <- unique(as.integer(i))
              spec <- modelDims(x)
              q <- spec$q
              if (!all(abs(i) %in% (1:q))) 
                  stop("SubMoment must be between 1 and q")
               if (length(i)==q)
                  return(x)
               momNames <- x@momNames[i]
               if (length(momNames)<spec$k)
                   stop("The model is under-identified")
              x@fRHS <- x@fRHS[i]
              x@fLHS <- x@fLHS[i]
              x@q <- length(momNames)
              x@momNames <- momNames
              x
           })

setMethod("[", c("gmmModels", "missing", "missing"),
          function(x, i, j) x)

### Observation subset

setGeneric("subset")
setMethod("subset", "regGmm",
          function(x, i) {
              x@modelF <- x@modelF[i,,drop=FALSE]
              x@instF <- x@instF[i,,drop=FALSE]
              x@n <- nrow(x@modelF)
              x})

setMethod("subset", "functionGmm",
          function(x, i) {
              if (is.matrix(x@X) || is.data.frame(x@X))
                  x@X <- x@X[,i,drop=FALSE]
              else if (is.numeric(x@X))
                  x@X <- x@X[i, drop=FALSE]
              else
                  stop("X is not subsetable")
              x@n <- nrow(x@X)
              x})

setMethod("subset", "formulaGmm",
          function(x, i) {
              x@modelF <- x@modelF[i,,drop=FALSE]
              x@n <- nrow(x@modelF)
              x})

## gmmFit

setGeneric("gmmFit", function(object, ...) standardGeneric("gmmFit"))

setMethod("gmmFit", signature("formulaGmm"), valueClass="gmmfit", 
          definition = function(object, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, start=NULL, ...)
              {
                  if (object@isMDE && object@centeredVcov)
                      {
                          if (is.character(weights) && weights == "optimal")
                              {
                                  spec <- modelDims(object)
                                  wObj <- evalWeights(object, spec$theta0, "optimal")
                                  met <- getMethod("gmmFit", "gmmModels")
                                  res <- met(object, weights=wObj, efficientWeights=TRUE,
                                             ...)
                                  res@type <- "mde"
                                  return(res)
                              } else {
                                  callNextMethod()
                              }
                      } else {
                          callNextMethod()
                      }
              })

setMethod("gmmFit", signature("gmmModels"), valueClass="gmmfit", 
          definition = function(object, type=c("twostep", "iter","cue", "onestep"),
              itertol=1e-7, initW=c("ident", "tsls"), weights="optimal", 
              itermaxit=100, efficientWeights=FALSE, start=NULL, ...)
              {
                  Call <- match.call()
                  chk <- validObject(object)                  
                  type <- match.arg(type)
                  initW <- match.arg(initW)
                  i <- 1L
                  chk <- validObject(object, TRUE)
                  if (!chk)
                      stop("object is not a valid gmmModels object")
                  if (initW == "tsls" && class(object) != "linearGmm")
                      stop("initW='tsls' is for linear models only")
                  if (is.character(weights) && !(weights%in%c("optimal","ident")))
                      stop("weights is a matrix or one of 'optimal' or 'ident'")
                  spec <- modelDims(object)
                  if (spec$q==spec$k)
                      {
                          # This allow to weight the moments in case of
                          # large scale difference.
                          if (!is.matrix(weights) && class(weights)!="gmmWeights")
                              weights <- "ident"
                          type <- "onestep"
                      } else if (type == "onestep" && !is.matrix(weights)) {
                          weights <- "ident"
                      } else if (is.matrix(weights) || class(weights)=="gmmWeights") {
                          type <- "onestep"
                      } else if (weights == "ident") {
                          type <- "onestep"
                      }
                  if (type == "onestep")
                      {
                          if (class(weights)=="gmmWeights")
                              wObj <- weights
                          else
                              wObj <- evalWeights(object, w=weights)
                          res <- solveGmm(object, wObj, start, ...)
                          convergence <- res$convergence
                          efficientGmm <- ifelse(is.character(weights), FALSE,
                                             efficientWeights)
                          ans <- new("gmmfit", theta=res$theta,
                                     convergence=convergence, convIter=NULL, type=type,
                                     wObj=wObj, model=object, call=Call, niter=i,
                                     efficientGmm=efficientGmm)
                          return(ans)
                      }
                  if (class(object) == "linearGmm")
                      {
                          if (object@vcov == "iid")
                              if (is.character(weights) && weights == "optimal")
                                  return(tsls(object))
                      }
                  if (type == "twostep")
                      {
                          itermaxit <- 1
                      }
                  if (initW=="tsls")
                      {                          
                          theta0 <- coef(tsls(object))
                      } else {
                          wObj <- evalWeights(object, NULL, "ident")
                          theta0 <- solveGmm(object, wObj, start, ...)$theta
                      }
                  bw <- object@bw
                  if (type != "cue")
                  {
                      while(TRUE)
                          {                              
                              wObj <- evalWeights(object, theta0, "optimal")
                              if (is.character(object@bw) && object@vcov=="HAC")
                                  object@bw <- wObj@HAC$bw
                              res <- solveGmm(object, wObj, theta0, ...)
                              theta1 <- res$theta
                              convergence <- res$convergence
                              crit <- sqrt( sum((theta1-theta0)^2)/(1+sqrt(sum(theta0^2))))
                              if (crit < itertol & type=="iter")
                                  {
                                      convIter <- 0
                                      break
                                  }
                              i <- i + 1L
                              theta0 <- theta1
                              if (i>itermaxit)
                                  {
                                      if (type=="twostep")
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
                      obj <- function(theta, object)
                          {
                              wObj <- evalWeights(object, theta, "optimal")
                              evalObjective(object, theta, wObj)
                          }
                      res <- optim(theta0, obj, object=object,
                                   ...)
                      theta1 <- res$par
                      convergence <- res$convergence
                      wObj <- evalWeights(object, theta1, "optimal")
                      
                  }
                  object@bw <- bw
                  names(theta1) <- spec$parNames
                  new("gmmfit", theta=theta1, convergence=convergence, type=type,
                      wObj=wObj, model=object, convIter=convIter, call=Call,
                      niter=i, efficientGmm=TRUE)
              })

## tsls

setGeneric("tsls", function(object, ...) standardGeneric("tsls"))

setMethod("tsls", signature("linearGmm"), valueClass="tsls", 
          function(object)
          {
              Call <- match.call()
              chk <- validObject(object)
              X <- model.matrix(object)
              Z <- model.matrix(object, "instrument")
              Y <- modelResponse(object)
              spec <- modelDims(object)
              EndoVars <- !(spec$parNames %in% spec$momNames)
              if (any(EndoVars))
                  {
                      res <- lm(X[,EndoVars]~Z)
                      X[,EndoVars] <- fitted(res)
                  }
              theta <- lm.fit(X, Y)$coefficients
              names(theta) <- spec$parNames
              vcov <- object@vcov
              object@vcov <- "iid"
              efficientGmm <- vcov == "iid"
              wObj <- evalWeights(object, theta, "optimal")
              object@vcov <- vcov
              obj <- new("tsls", theta=theta, convergence=NULL, type="tsls",
                         wObj=wObj, model=object, convIter=NULL, call=Call,
                         niter=1L, efficientGmm=efficientGmm)
              obj
          })


## evalGmm

setGeneric("evalGmm", function(object, ...) standardGeneric("evalGmm"))

setMethod("evalGmm", signature("gmmModels"),
          function(object, theta, wObj=NULL) {
              spec <- modelDims(object)
              if (!is.null(names(theta)))
                  {
                      if (!all(names(theta) %in% spec$parNames))
                          stop("You provided a named theta with wrong names")
                      theta <- theta[match(spec$parNames, names(theta))]
                  } else {
                      if (class(object) %in% c("formulaGmm","nonlinearGmm"))
                          stop("To evaluate nonlinear models, theta must be named")
                      names(theta) <- spec$parNames
                  }
              call <- match.call()
              if (is.null(wObj))
                  wObj <- evalWeights(object, theta)
              new("gmmfit", theta=theta, convergence=NULL, convIter=NULL,
                  call=call, type="eval", wObj=wObj, niter=0L, efficientGmm=FALSE,
                  model=object)
          })



