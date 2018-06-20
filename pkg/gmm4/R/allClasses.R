#####  All S4 classes of the package are defined here
######################################################


## gmmModel

setClassUnion("matrixORcharacter", c("matrix", "character"))
setClassUnion("matrixORnumeric", c("matrix", "numeric"))
setClassUnion("numericORcharacter", c("numeric", "character"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("numericORmatrixORNULL", c("matrix", "numeric", "NULL"))
setClassUnion("expressionORNULL", c("expression", "NULL"))
setClassUnion("functionORNULL", c("function", "NULL"))                                 
setClass("linearGmm", representation(modelF="data.frame", instF="data.frame",
                                     vcov="character",n="integer", q="integer", k="integer",
                                     parNames="character", momNames="character",
                                     kernel="character", bw="numericORcharacter",
                                     prewhite="integer", ar.method="character",
                                     approx="character", tol="numeric",
                                     centeredVcov="logical", varNames="character",
                                     isEndo="logical"),
         prototype(vcov="MDS", kernel="Quadratic Spectral", bw="Andrews", prewhite=1L,
                   ar.method="ols", approx="AR(1)", tol=1e-7))
setClass("nonlinearGmm", representation(modelF="data.frame", instF="data.frame",
                                        vcov="character",theta0="numeric",
                                        n="integer", q="integer",k="integer",
                                        parNames="character", momNames="character",
                                        fRHS="expression", fLHS="expressionORNULL",
                                        kernel="character", bw="numericORcharacter",
                                        prewhite="integer", ar.method="character",
                                        approx="character", tol="numeric",
                                        centeredVcov="logical", varNames="character",
                                        isEndo="logical"),
         prototype(vcov="MDS", kernel="Quadratic Spectral", bw="Andrews", prewhite=1L,
                   ar.method="ols", approx="AR(1)", tol=1e-7))
setClass("functionGmm", representation(X="ANY", fct="function",dfct="functionORNULL",
                                       vcov="character",theta0="numeric",
                                       n="integer", q="integer",k="integer",
                                       parNames="character", momNames="character",
                                       kernel="character", bw="numericORcharacter",
                                       prewhite="integer", ar.method="character",
                                       approx="character", tol="numeric",
                                       centeredVcov="logical", varNames="character",
                                       isEndo="logical"),
         prototype(vcov="MDS", kernel="Quadratic Spectral", bw="Andrews", prewhite=1L,
                   ar.method="ols", approx="AR(1)", tol=1e-7, dfct=NULL))
setClassUnion("regGmm", c("linearGmm", "nonlinearGmm"))
setClassUnion("allNLGmm", c("nonlinearGmm", "functionGmm"))
setClassUnion("gmmModels", c("linearGmm", "nonlinearGmm", "functionGmm"))

## gmmWeights

setClass("gmmWeights", representation(w="ANY", type="character", HAC="list"),
         prototype(HAC=list()))

## gmmfit

setClass("gmmfit", representation(theta = "numeric", convergence = "numericORNULL",
                                  convIter="numericORNULL",call="call",
                                  type="character", wObj="gmmWeights",niter="integer",
                                  efficientGmm="logical", model="gmmModels"))

setClass("tsls", contains="gmmfit")

## specTest

setClass("specTest", representation(test = "matrix", testname="character"))

## summaryGmm

setClass("summaryGmm", representation(coef="matrix", specTest = "specTest",
                                      strength="list", model="gmmModels",sandwich="logical",
                                      type="character", convergence = "numericORNULL",
                                      convIter="numericORNULL", wSpec="list",niter="integer",
                                      df.adj="logical", breadOnly="logical"))
## Restricted gmm Models

setClass("rlinearGmm", representation(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="linearGmm")

setClass("rnonlinearGmm", representation(R="list", cstSpec="list"),
         contains="nonlinearGmm")

setClass("rfunctionGmm", representation(R="list", cstSpec="list"),
         contains="functionGmm")

setClassUnion("rgmmModels", c("rlinearGmm", "rnonlinearGmm", "rfunctionGmm"))

## hypothesisTest

setClass("hypothesisTest", representation(test="numeric", hypothesis="character",
                                          dist="character", df="integer", pvalue="numeric",
                                          type="character"))
### System GMM

setClass("slinearGmm", representation(modelT="list", instT="list",data="data.frame",
                                      vcov="character",n="integer", q="integer",
                                      k="integer", parNames="list",
                                      momNames="list", eqnNames="character",
                                      kernel="character", bw="numericORcharacter",
                                      prewhite="integer", ar.method="character",
                                      approx="character", tol="numeric",
                                      centeredVcov="logical", sameMom="logical", SUR="logical",
                                      varNames="list", isEndo="list"),
         prototype(vcov="MDS", kernel="Quadratic Spectral", bw="Andrews", prewhite=1L,
                   ar.method="ols", approx="AR(1)", tol=1e-7))
setClass("snonlinearGmm", representation(data="data.frame", instT="list",
                                         vcov="character",theta0="list",
                                         n="integer", q="integer",k="integer",
                                         parNames="list", momNames="list",
                                         fRHS="list", fLHS="list", eqnNames="character",
                                         kernel="character", bw="numericORcharacter",
                                         prewhite="integer", ar.method="character",
                                         approx="character", tol="numeric",
                                         centeredVcov="logical", sameMom="logical",
                                         SUR="logical",
                                         varNames="list", isEndo="list"),
         prototype(vcov="MDS", kernel="Quadratic Spectral", bw="Andrews", prewhite=1L,
                   ar.method="ols", approx="AR(1)", tol=1e-7))
setClassUnion("sysGmmModels", c("slinearGmm", "snonlinearGmm"))

## Restricted System GMM

setClass("rslinearGmm", representation(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="slinearGmm")

setClass("rsnonlinearGmm", representation(R="list", cstSpec="list"),
         contains="snonlinearGmm")

setClassUnion("rsysGmmModels", c("rslinearGmm", "rsnonlinearGmm"))

### sysGmmWeights

setClass("sysGmmWeights", representation(w="ANY", type="character", HAC="list",
                                         Sigma="ANY", momNames="list", eqnNames="character",
                                         sameMom="logical"),
         prototype(w="ident", type="weights", momNames=list(), eqnNames=character(),
                   HAC=list(), sameMom=FALSE))

## summarySysGmm

setClass("summarySysGmm",
         representation(coef="list", specTest = "specTest",
                        strength="list", model="sysGmmModels",sandwich="logical",
                        type="character", convergence = "numericORNULL",
                        convIter="numericORNULL", wSpec="list",niter="integer",
                        df.adj="logical", breadOnly="logical"))

## Class converters

setAs("linearGmm", "nonlinearGmm",
      function(from) {
          spec <- modelDims(from)
          X <- model.matrix(from)
          theta0 <- rep(1,ncol(X))
          names(theta0) <- paste("theta", 1:ncol(X), sep="")         
          colnames(X) <- paste("X", 1:ncol(X), sep="")
          rhs <- paste(names(theta0), "*", colnames(X), sep="")
          rhs <- paste(rhs, collapse="+", sep="")
          rhs <- parse(text=rhs)
          X <- data.frame(Y=modelResponse(from), X)
          lhs <- expression(Y)
          new("nonlinearGmm", modelF=X, instF=from@instF, vcov=from@vcov,
              theta0=theta0, n=spec$n, q=spec$q, k=spec$k, parNames=names(theta0),
              momNames=spec$momNames, fRHS=rhs, fLHS=lhs, kernel=from@kernel,
              bw=from@bw, prewhite=from@prewhite, ar.method=from@ar.method,
              approx=from@approx, tol=from@tol, centeredVcov=from@centeredVcov,
              isEndo=from@isEndo, varNames=from@varNames)
      })

setAs("linearGmm", "functionGmm",
      function(from) {
          spec <- modelDims(from)          
          X <- model.matrix(from)
          theta0 <- rep(1,ncol(X))
          names(theta0) <- paste("theta", 1:ncol(X), sep="")         
          colnames(X) <- paste("X", 1:ncol(X), sep="")         
          Z <- model.matrix(from, "instruments")
          colnames(Z) <- paste("Z", 1:ncol(Z), sep="")         
          dat <- cbind(X, Z, Y=modelResponse(from))
          theta0 <- rep(0,ncol(X))
          names(theta0) <- paste("theta", 1:ncol(X), sep="")
          fct <- function(theta, x)
              {
                  wx <- which(strtrim(colnames(x),1) == "X")
                  wz <- which(strtrim(colnames(x),1) == "Z")
                  wy <- which(strtrim(colnames(x),1) == "Y")
                  e <- x[,wy]-c(x[,wx,drop=FALSE]%*%theta)
                  e*x[,wz]
              }
          dfct <- function(theta, x)
              {
                  wx <- which(strtrim(colnames(x),1) == "X")
                  wz <- which(strtrim(colnames(x),1) == "Z")
                  -crossprod(x[,wz],x[,wx])/nrow(x)
              }
          new("functionGmm", X=dat, fct=fct, dfct=dfct,  vcov=from@vcov,
              theta0=theta0, n=spec$n, q=spec$q, k=spec$k, parNames=names(theta0),
              momNames=colnames(Z), kernel=from@kernel,
              bw=from@bw, prewhite=from@prewhite, ar.method=from@ar.method,
              approx=from@approx, tol=from@tol, centeredVcov=from@centeredVcov)
      })

setAs("slinearGmm", "linearGmm",
      function(from) {
          spec <- modelDims(from)
          eqnNames <- from@eqnNames
          neqn <- length(eqnNames)
          datX <- lapply(1:neqn,
                         function(i) {
                             v <- from@varNames[[i]]
                             chk <- "(Intercept)" %in% v
                             v <- v[v!="(Intercept)"]
                             X <- from@data[,v]
                             colnames(X) <- paste(eqnNames[[i]],".", v, sep="")
                             if (chk)
                                 {
                                  X <- cbind(1, X)
                                  colnames(X)[1]<-paste(eqnNames[[i]], ".Intercept", sep="")
                                 }
                             X})
          datZ <- lapply(1:neqn,
                         function(i) {
                             v <- all.vars(from@instT[[i]])
                             chk <- attr(from@instT[[i]], "intercept")==1
                             Z <- from@data[,v]
                             colnames(Z) <- paste(eqnNames[[i]],".", v, sep="")
                             if (chk)
                                 {
                                  Z <- cbind(1, Z)
                                  colnames(Z)[1]<-paste(eqnNames[[i]], ".Intercept", sep="")
                                 }
                             Z})
          nZ <- do.call("c", lapply(datZ, colnames))
          nX <- do.call("c", lapply(datX, colnames))
          datZ <- .GListToMat(datZ)
          datX <- .GListToMat(datX)
          Y <- do.call("c", modelResponse(from))
          colnames(datZ) <- nZ
          colnames(datX) <- nX
          dat <- cbind(Y, datZ, datX)
          dat <- dat[,unique(colnames(dat))]
          dat <- data.frame(dat, row.names=1:nrow(datZ))
          g <- paste("Y~", paste(nX, collapse="+"), "-1")
          g <- formula(g, .GlobalEnv)
          h <- paste("~", paste(nZ, collapse="+"), "-1")
          h <- formula(h, .GlobalEnv)
          res <- gmmModel(g, h, vcov=from@vcov, kernel=from@kernel, bw=from@bw,
                          prewhite=from@prewhite, ar.method=from@ar.method,
                          approx=from@approx, tol=from@tol, centeredVcov=from@centeredVcov,
                          data=dat)
      })


setAs("rslinearGmm", "rlinearGmm",
      function(from) {
          m <- as(from, "slinearGmm")
          m <- as(m, "linearGmm")
          restGmmModel(m, from@cstLHS, from@cstRHS)
      })

setAs("sysGmmWeights", "gmmWeights",
      function(from) {
          w <- quadra(from)
          if (is.character(w))
              w <- "ident"
          new("gmmWeights", w=w, type="weights", HAC=list())
      })
          


### system GMM fit

setClass("sgmmfit", representation(theta = "list", convergence = "numericORNULL",
                                  convIter="numericORNULL",call="call",
                                  type="character", wObj="sysGmmWeights",niter="integer",
                                   efficientGmm="logical", model="sysGmmModels"))

setClass("stsls", contains="sgmmfit")
