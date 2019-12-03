#####  All S4 classes of the package are defined here
######################################################


## gmmModel

setClassUnion("matrixORcharacter", c("matrix", "character"))
setClassUnion("matrixORnumeric", c("matrix", "numeric"))
setClassUnion("numericORcharacter", c("numeric", "character"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("numericORlogical", c("numeric", "logical"))
setClassUnion("numericORmatrixORNULL", c("matrix", "numeric", "NULL"))
setClassUnion("expressionORNULL", c("expression", "NULL"))
setClassUnion("functionORNULL", c("function", "NULL"))
setClassUnion("callORNULL", c("call", "NULL"))
setClass("linearGmm", representation(modelF="data.frame", instF="data.frame",
                                     vcov="character",n="integer", q="integer", k="integer",
                                     parNames="character", momNames="character",
                                     vcovOptions="list", centeredVcov="logical",
                                     varNames="character", isEndo="logical",
                                     omit='integer', survOptions="list"))
setClass("nonlinearGmm", representation(modelF="data.frame", instF="data.frame",
                                        vcov="character",theta0="numeric",
                                        n="integer", q="integer",k="integer",
                                        parNames="character", momNames="character",
                                        fRHS="expression", fLHS="expressionORNULL",
                                        vcovOptions="list",
                                        centeredVcov="logical", varNames="character",
                                        isEndo="logical",omit='integer', survOptions="list"))
setClass("functionGmm", representation(X="ANY", fct="function",dfct="functionORNULL",
                                       vcov="character",theta0="numeric",
                                       n="integer", q="integer",k="integer",
                                       parNames="character", momNames="character",
                                       vcovOptions="list",
                                       centeredVcov="logical", varNames="character",
                                       isEndo="logical",omit='integer', survOptions="list"))
setClass("formulaGmm", representation(modelF="data.frame", 
                                      vcov="character",theta0="numeric",
                                      n="integer", q="integer",k="integer",
                                      parNames="character", momNames="character",
                                      fRHS="list", fLHS="list",
                                      vcovOptions="list",
                                      centeredVcov="logical", varNames="character",
                                      isEndo="logical", isMDE="logical",omit='integer',
                                      survOptions="list"))
setClassUnion("regGmm", c("linearGmm", "nonlinearGmm"))
setClassUnion("allNLGmm", c("nonlinearGmm", "functionGmm", "formulaGmm"))
setClassUnion("gmmModels", c("linearGmm", "nonlinearGmm", "functionGmm", "formulaGmm"))

## GEL models

setClass("linearGel", representation(wSpec="list", gelType="list"),
         contains="linearGmm")
setClass("nonlinearGel", representation(wSpec="list", gelType="list"),
         contains="nonlinearGmm")
setClass("functionGel", representation(wSpec="list", gelType="list"),
         contains="functionGmm")
setClass("formulaGel", representation(wSpec="list", gelType="list"),
         contains="formulaGmm")

setClassUnion("gelModels", c("linearGel", "nonlinearGel", "functionGel", "formulaGel"))

## gmmWeights

setClass("gmmWeights", representation(w="ANY", type="character", wSpec="list"))

## gmmfit

setClass("gmmfit", representation(theta = "numeric", convergence = "numericORNULL",
                                  convIter="numericORNULL",call="callORNULL",
                                  type="character", wObj="gmmWeights",niter="integer",
                                  efficientGmm="logical", model="gmmModels"))

setClass("tsls", contains="gmmfit")

## gelfit

setClass("gelfit", representation(theta = "numeric", convergence = "numeric",
                                  lambda = "numeric", lconvergence = "numeric",
                                  call="callORNULL", type="character", vcov="list",
                                  model="gelModels"))

## specTest

setClass("specTest", representation(test = "matrix", testname="character"))

## confint

setClass("confint", representation(interval = "matrix", type="character",
                                   level="numeric"))


setClass("mconfint", 
         representation(areaPoints="matrix", type="character", level="numeric"))

## summaryGmm

setClass("summaryGmm", representation(coef="matrix", specTest = "specTest",
                                      strength="list", model="gmmModels",sandwich="logical",
                                      type="character", convergence = "numericORNULL",
                                      convIter="numericORNULL", wSpec="list",niter="integer",
                                      df.adj="logical", breadOnly="logical"))

setClass("summaryGel", representation(coef="matrix", specTest = "specTest",
                                      model="gelModels", lambda="matrix",
                                      convergence="numeric",lconvergence="numeric",
                                      impProb="list"))

## Restricted gmm Models

setClass("rlinearGmm", representation(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="linearGmm")

setClass("rnonlinearGmm", representation(R="list", cstSpec="list"),
         contains="nonlinearGmm")

setClass("rfunctionGmm", representation(R="list", cstSpec="list"),
         contains="functionGmm")

setClass("rformulaGmm", representation(R="list", cstSpec="list"),
         contains="formulaGmm")

setClassUnion("rgmmModels", c("rlinearGmm", "rnonlinearGmm", "rfunctionGmm",
                              "rformulaGmm"))

## Restricted gel Models

setClass("rlinearGel", representation(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="linearGel")

setClass("rnonlinearGel", representation(R="list", cstSpec="list"),
         contains="nonlinearGel")

setClass("rfunctionGel", representation(R="list", cstSpec="list"),
         contains="functionGel")

setClass("rformulaGel", representation(R="list", cstSpec="list"),
         contains="formulaGel")

setClassUnion("rgelModels", c("rlinearGel", "rnonlinearGel", "rfunctionGel",
                              "rformulaGel"))

## hypothesisTest

setClass("hypothesisTest", representation(test="numeric", hypothesis="character",
                                          dist="character", df="integer", pvalue="numeric",
                                          type="character"))
### System GMM

setClass("slinearGmm", representation(modelT="list", instT="list",data="data.frame",
                                      vcov="character",n="integer", q="integer",
                                      k="integer", parNames="list",
                                      momNames="list", eqnNames="character",
                                      vcovOptions="list",
                                      centeredVcov="logical", sameMom="logical",
                                      SUR="logical", varNames="list", isEndo="list",
                                      omit='integer', survOptions="list"))

setClass("snonlinearGmm", representation(data="data.frame", instT="list",
                                         vcov="character",theta0="list",
                                         n="integer", q="integer",k="integer",
                                         parNames="list", momNames="list",
                                         fRHS="list", fLHS="list", eqnNames="character",
                                         vcovOptions="list",
                                         centeredVcov="logical", sameMom="logical",
                                         SUR="logical",
                                         varNames="list", isEndo="list",
                                         omit='integer', survOptions="list"))
setClassUnion("sysGmmModels", c("slinearGmm", "snonlinearGmm"))

## Restricted System GMM

setClass("rslinearGmm", representation(cstLHS="matrix", cstRHS="numeric", cstSpec="list"),
         contains="slinearGmm")

setClass("rsnonlinearGmm", representation(R="list", cstSpec="list"),
         contains="snonlinearGmm")

setClassUnion("rsysGmmModels", c("rslinearGmm", "rsnonlinearGmm"))

### sysGmmWeights

setClass("sysGmmWeights", representation(w="ANY", type="character", wSpec="list",
                                         Sigma="ANY", momNames="list", eqnNames="character",
                                         sameMom="logical"))
## summarySysGmm

setClass("summarySysGmm",
         representation(coef="list", specTest = "specTest",
                        strength="list", model="sysGmmModels",sandwich="logical",
                        type="character", convergence = "numericORNULL",
                        convIter="numericORNULL", wSpec="list",niter="integer",
                        df.adj="logical", breadOnly="logical"))

## Class converters

setAs("rgelModels", "rgmmModels",
      function(from) {
          obj <- as(from, "gmmModels")
          cls <- strsplit(class(from), "Gel")[[1]][1]
          cls <- paste(cls, "Gmm", sep="")
          if (grepl("linear", class(from)))
              new("rlinearGmm", cstLHS=from@cstLHS, cstRHS=from@cstRHS,
                  cstSpec=from@cstSpec, obj)
          else
              new(cls, R=from@R, cstSpec=from@cstSpec, obj)
      })

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
              momNames=spec$momNames, fRHS=rhs, fLHS=lhs,
              vcovOptions=from@vcovOptions, centeredVcov=from@centeredVcov,
              isEndo=from@isEndo, varNames=from@varNames,omit=from@omit,
              survOptions=from@survOptions)
      })

setAs("linearGel", "nonlinearGel",
      function(from) {
          model <- as(as(from, "linearGmm"), "nonlinearGmm")
          new("nonlinearGel", wSpec=from@wSpec, gelType=from@gelType, model)
      })


setAs("linearGmm", "functionGmm",
      function(from) {
          spec <- modelDims(from)          
          x <- from
          theta0 <- rep(0,spec$k)
          names(theta0) <- spec$parNames
          fct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalMoment(x, theta)
              }
          dfct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalDMoment(x, theta)
              }
          new("functionGmm", X=x, fct=fct, dfct=dfct,  vcov=from@vcov,
              theta0=theta0, n=spec$n, q=spec$q, k=spec$k, parNames=names(theta0),
              momNames=spec$momNames,vcovOptions=from@vcovOptions,
              centeredVcov=from@centeredVcov,omit=integer(),survOptions=from@survOptions)
      })

setAs("linearGel", "functionGel",
      function(from) {
          model <- as(as(from, "linearGmm"), "functionGmm")
          new("functionGel", wSpec=from@wSpec, gelType=from@gelType, model)
      })


setAs("allNLGmm", "functionGmm",
      function(from) {
          spec <- modelDims(from)          
          x <- from
          fct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalMoment(x, theta)
              }
          dfct <- function(theta, x)
              {
                  names(theta) <- modelDims(x)$parNames
                  gt <- evalDMoment(x, theta)
              }
          new("functionGmm", X=x, fct=fct, dfct=dfct,  vcov=from@vcov,
              theta0=from@theta0, n=spec$n, q=spec$q, k=spec$k,
              parNames=names(from@theta0),
              momNames=spec$momNames, vcovOptions=from@vcovOptions,
              centeredVcov=from@centeredVcov,omit=integer(), survOptions=from@survOptions)
      })

setAs("nonlinearGel", "functionGel",
      function(from) {
          model <- as(as(from, "nonlinearGmm"), "functionGmm")
          new("functionGel", wSpec=from@wSpec, gelType=from@gelType, model)
      })

setAs("formulaGel", "functionGel",
      function(from) {
          model <- as(as(from, "formulaGmm"), "functionGmm")
          new("functionGel", wSpec=from@wSpec, gelType=from@gelType, model)
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
                             X <- from@data[,v, drop=FALSE]
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
                             Z <- from@data[,v, drop=FALSE]
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
          res <- gmmModel(g, h, vcov=from@vcov, vcovOptions=from@vcovOptions,
                          centeredVcov=from@centeredVcov, data=dat)
      })

setAs("rslinearGmm", "rlinearGmm",
      function(from) {
          m <- as(from, "slinearGmm")
          m <- as(m, "linearGmm")
          restModel(m, from@cstLHS, from@cstRHS)
      })

setAs("sysGmmWeights", "gmmWeights",
      function(from) {
          w <- quadra(from)
          if (is.character(w))
              w <- "ident"
          new("gmmWeights", w=w, type="weights", wSpec=list())
      })
          


### system GMM fit

setClass("sgmmfit", representation(theta = "list", convergence = "numericORNULL",
                                   convIter="numericORNULL",call="callORNULL",
                                   type="character", wObj="sysGmmWeights",niter="integer",
                                   efficientGmm="logical", model="sysGmmModels"))

setClass("stsls", contains="sgmmfit")
