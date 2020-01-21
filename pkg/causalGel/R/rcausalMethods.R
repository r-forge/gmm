## restricted model constructor

setMethod("restModel", signature("causalModel"),
          function(object, R, rhs=NULL)
          {
              mod <- restModel(as(object, "momentModel"), R, rhs)
              new("rcausalModel", mod)
          })


## print

setMethod("print", "rcausalModel",
          function(x)
          {
              print(as(x, "causalModel"))
              cat("Additional Specifications: Restricted model\n")
              printRestrict(x)
          })


## modelFit

setMethod("gelFit", signature("rcausalModel"), valueClass="causalGelfit", 
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


## modelDims

setMethod("modelDims", "rcausalModel",
          function(object) {
              res <- callNextMethod()
              res$balCov <- object@X@balCov
              res$momType <- object@X@momType
              res$balMom <- object@X@balMom
              res$ACTmom <- object@X@ACTmom
              res
          })


## subsetting

setMethod("subset", "rcausalModel",
          function(x, i) {
              x@X@reg <- x@X@reg[i,,drop=FALSE]
              x@X@bal <- x@X@bal[i,,drop=FALSE]
              x@n <- nrow(x@X@reg)
              x})

setMethod("[", c("rcausalModel", "numeric", "missing"),
          function(x, i, j){
              obj <- as(x, "causalModel")[i]
              mod <- new("rfunctionModel", R=x@R, cstSpec=x@cstSpec, obj)
              new("rcausalGel", mod)
          })

setMethod("[", c("rcausalModel", "numeric", "numericORlogical"),
          function(x, i, j){
              x <- x[i]
              subset(x, j)
          })

setMethod("[", c("rcausalModel", "missing", "numericORlogical"),
          function(x, i, j){
              subset(x, j)
          })

