### restModel

#setMethod("restModel", signature("linearGel"),
#          function(object, R, rhs=NULL)
#          {
#              mod <- getMethod("restModel", "linearGmm")(object, R, rhs)
#              new("rlinearGel",  cstLHS=mod@cstLHS, cstRHS=mod@cstRHS,
#                  cstSpec=mod@cstSpec, object)
#})
setMethod("restModel", signature("gelModels"),
          function(object, R, rhs=NULL)
          {
              mod <- callNextMethod()
              gmmToGel(mod, object@gelType$name, object@getType$rhoFct)
          })

##setMethod("restModel", signature("nonlinearGel"),
##         function(object, R, rhs=NULL)
##          {
##              mod <- getMethod("restModel", "nonlinearGmm")(object, R, rhs)
##              new("rnonlinearGel",  R=mod@R, cstSpec=mod@cstSpec, object)
##          })


## printRestrict

setMethod("printRestrict", signature("rgelModels"),
          function(object)
          {
              cl <- strsplit(class(object)[1],"Gel")[[1]][1]
              cl <- paste(cl, "Gmm", sep="")
              getMethod("printRestrict", cl)(object)
          })


## print

setMethod("print", "rgelModels",
          function(x)
          {
              cl <- class(x)[1]
              getMethod("print", "gelModels")(x)
              printRestrict(x)
          })

## modelDims


setMethod("modelDims", "rgelModels",
          function(object)
          {
              cl <- strsplit(class(object)[1],"Gel")[[1]][1]
              cl <- paste(cl, "Gmm", sep="")
              getMethod("modelDims", cl)(object)
              
          })

## model.matrix and modelResponse


setMethod("model.matrix", "rlinearGel",
          function(object, type=c("regressors","instruments"))
          {
              type <- match.arg(type)
              getMethod("model.matrix", "rlinearGmm")(object, type)
          })

setMethod("modelResponse", "rlinearGel",
          function(object)
          {
              getMethod("modelResponse", "rlinearGmm")(object)
          })
                                       

## getRestrict


setMethod("getRestrict", "rgelModels",
          function(object, theta)
          {
              cl <- strsplit(class(object)[1],"Gel")[[1]][1]
              cl <- paste(cl, "Gmm", sep="")
              getMethod("getRestrict", cl)(object, theta)
              
          })

setMethod("getRestrict", "gelModels",
          function(object, theta, R, rhs=NULL) {
              getMethod("getRestrict", "gmmModels")(object)
          })


## coef


setMethod("coef", "rgelModels",
          function(object, theta)
          {
              cl <- strsplit(class(object)[1],"Gel")[[1]][1]
              cl <- paste(cl, "Gmm", sep="")
              getMethod("coef", cl)(object, theta)
          })

## subset

 
setMethod("[", c("rfunctionGel", "numeric", "missing"),
          function(x, i, j){
               callNextMethod()
          })


## modelFit

 
setMethod("modelFit", signature("rgelModels"), valueClass="gelfit", 
          definition = function(object, gelType=NULL, rhoFct=NULL,
                                initTheta=c("gmm", "modelTheta0"), theta0=NULL,
                                lambda0=NULL, vcov=FALSE, ...)
          {
              Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
              if (class(Call)=="try-error")
                  Call <- NULL
              k <- modelDims(object)$k
              if (k == 0)
              {
                  if (!is.null(gelType))
                      object@gelType$name <- gelType
                  if (!is.null(rhoFct))
                      object@gelType$rhoFct <- rhoFct
                  return(evalModel(object, numeric(), ...))
              }
              initTheta <- match.arg(initTheta)
              if (is.null(theta0))
              {
                  if (initTheta == "gmm")
                      theta0 <- modelFit(as(object, "rgmmModels"))@theta
                  else
                      theta0 <- modelDims(object)$theta0
              }
              obj <- getMethod("modelFit", "gelModels")(object=object, gelType=gelType,
                  rhoFct=rhoFct, initTheta=initTheta, theta0=theta0,
                  lambda0=lambda0, vcov=vcov, ...)
              obj@call <- Call
              obj
              })


