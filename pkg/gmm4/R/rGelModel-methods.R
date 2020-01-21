
setMethod("restModel", signature("linearGel"),
          function(object, R, rhs=NULL)
          {
              mod <- restModel(as(object, "gmmModels"), R, rhs)
              gmmToGel(mod, object@gelType$name, object@gelType$rhoFct)
          })

setMethod("restModel", signature("nonlinearGel"),
          function(object, R, rhs=NULL)
          {
              mod <- restModel(as(object, "gmmModels"), R, rhs)
              gmmToGel(mod, object@gelType$name, object@gelType$rhoFct)
          })

setMethod("restModel", signature("formulaGel"),
          function(object, R, rhs=NULL)
          {
              mod <- restModel(as(object, "gmmModels"), R, rhs)
              gmmToGel(mod, object@gelType$name, object@gelType$rhoFct)
          })

setMethod("restModel", signature("functionGel"),
          function(object, R, rhs=NULL)
          {
              mod <- restModel(as(object, "gmmModels"), R, rhs)
              gmmToGel(mod, object@gelType$name, object@gelType$rhoFct)
          })

## printRestrict

setMethod("printRestrict", signature("rgelModels"),
          function(object) printRestrict(as(object, "rgmmModels")))

## print

setMethod("print", "rgelModels",
          function(x)
          {
              print(as(x, "gelModels"))
              printRestrict(x)
          })

## modelDims

setMethod("modelDims", "rlinearGel",
          function(object) modelDims(as(object, "rgmmModels")))

setMethod("modelDims", "rnonlinearGel",
          function(object) modelDims(as(object, "rgmmModels")))

setMethod("modelDims", "rfunctionGel",
          function(object) modelDims(as(object, "rgmmModels")))

setMethod("modelDims", "rformulaGel",
          function(object) modelDims(as(object, "rgmmModels")))

## model.matrix and modelResponse


setMethod("model.matrix", "rlinearGel",
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              model.matrix(as(object, "rgmmModels"), type)
          })

setMethod("modelResponse", "rlinearGel",
          function(object) modelResponse(as(object, "rgmmModels")))
                                       

## getRestrict


setMethod("getRestrict", "rgelModels",
          function(object, theta) getRestrict(as(object,"rgmmModels"), theta))


setMethod("getRestrict", "gelModels",
          function(object, theta, R, rhs=NULL)
              getRestrict(as(object,"gmmModels"), theta, R, rhs))

## coef

setMethod("coef", "rgelModels",
          function(object, theta) coef(as(object, "rgmmModels"), theta))

## subset

 
setMethod("[", c("rfunctionGel", "numeric", "missing"),
          function(x, i, j){
               callNextMethod()
          })

## evalDMoment

setMethod("evalDMoment", "rgelModels",
          function(object, theta, impProb=NULL, lambda=NULL)
          {
              spec <- modelDims(object)
              if (object@vcov != "HAC")
              {
                  G <- evalDMoment(as(object, "rgmmModels"), theta, impProb, lambda)
              } else {
                  G <- getMethod("evalDMoment","gelModels")(object, theta, impProb, lambda)
              }
              G})

## modelFit

setMethod("modelFit", signature("rlinearGel"), valueClass="gelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              met <- getMethod("modelFit", "rgelModels")
              obj <- met(model, gelType, rhoFct, initTheta, theta0, lambda0, vcov, ...)
              obj@call <- Call
              obj
          })

setMethod("modelFit", signature("rnonlinearGel"), valueClass="gelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              met <- getMethod("modelFit", "rgelModels")
              obj <- met(model, gelType, rhoFct, initTheta, theta0, lambda0, vcov, ...)
              obj@call <- Call
              obj
          })

setMethod("modelFit", signature("rformulaGel"), valueClass="gelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              met <- getMethod("modelFit", "rgelModels")
              obj <- met(model, gelType, rhoFct, initTheta, theta0, lambda0, vcov, ...)
              obj@call <- Call
              obj
          })

setMethod("modelFit", signature("rfunctionGel"), valueClass="gelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              met <- getMethod("modelFit", "rgelModels")
              obj <- met(model, gelType, rhoFct, initTheta, theta0, lambda0, vcov, ...)
              obj@call <- Call
              obj
          })

 
setMethod("modelFit", signature("rgelModels"), valueClass="gelfit", 
          definition = function(model, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (inherits(Call,"try-error"))
                  Call <- NULL
              k <- modelDims(model)$k
              if (k == 0)
              {
                  if (!is.null(gelType))
                      model@gelType$name <- gelType
                  if (!is.null(rhoFct))
                      model@gelType$rhoFct <- rhoFct
                  return(evalModel(model, numeric(), ...))
              }
              initTheta <- match.arg(initTheta)

              if (is.null(theta0))
              {
                  if (initTheta == "gmm")
                  {
                      theta0 <- modelFit(as(model, "rgmmModels"))@theta
                  } else {
                      theta0 <- modelDims(model)$theta0
                  }
              }
              obj <- getMethod("modelFit", "gelModels")(model=model, gelType=gelType,
                  rhoFct=rhoFct, initTheta=initTheta, theta0=theta0,
                  lambda0=lambda0, vcov=vcov, ...)
              obj@call <- Call
              obj
              })


