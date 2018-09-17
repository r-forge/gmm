####### All methods with gelModels (and its subclasses) signature
#################################################################

#######################  Print ########################
### The getGeneric for print is here only, so the file must be compiled
### before any other files containing print

setMethod("print", "gelModels",
          function(x, ...) {
              cat("GEL Model: Type ", x@gelType$name, "\n")
              cat("*******************************\n")
              cat("Moment type: ", strsplit(is(x)[1], "G")[[1]][1], "\n", sep="")
              if (x@vcov == "HAC")
                  {
                      cat("Smoothing: ")
                      cat(x@wSpec$kernel, " kernel and ", sep="")
                      cat(x@bw, " bandwidth",  sep="")
                      cat(" (", round(x@wSpec$bw, 3), ")", sep="")
                  } else {
                      cat("No Smoothing required\n")
                  }
              cat("\n")
              d <- modelDims(x)
              cat("Number of regressors: ", d$k, "\n", sep="")
              cat("Number of moment conditions: ", d$q, "\n", sep="")
              if (!inherits(x, "functionGmm"))
                  cat("Number of Endogenous Variables: ", sum(x@isEndo), "\n", sep="")
              cat("Sample size: ", d$n, "\n")})             

################ evalMoment ##########################

setMethod("evalMoment", "gelModels", function(object, theta)
    {
        if (object@vcov != "HAC")
            {
                evalMoment(as(object, "gmmModels"), theta)
            } else {
                smoothGel(object, theta)$smoothx
            }
    })

################ evalDMoment ##########################

setMethod("evalDMoment", "gelModels", function(object, theta, impProb=NULL)
    {
        if (object@vcov != "HAC")
            {
                evalDMoment(as(object, "gmmModels"), theta, impProb)
            } else {
                f <- function(theta, object, impProb)
                    {
                        gt <- evalMoment(object, theta)
                        if (is.null(impProb))
                            colMeans(gt)
                        else
                            colSums(gt*impProb)
                    }
                env <- new.env()
                assign("theta", theta, envir = env)
                assign("object", object, envir = env)
                assign("f", f, envir = env)
                assign("impProb", impProb, envir=env)
                G <- numericDeriv(quote(f(theta, object, impProb)), "theta", 
                                  env)
                G <- attr(G, "gradient")
                spec <- modelDims(object)
                dimnames(G) <- list(spec$momNames, spec$parNames)
                G
            }
    })

################ momentVcov  ##########################

setMethod("momentVcov", signature("gelModels"),
          function(object, theta, ...){
              if (object@vcov != "HAC")
                  {
                      momentVcov(as(object, "gmmModels"), theta)
                  } else {
                      gt <- evalMoment(object, theta)
                      w <- crossprod(gt)/nrow(gt)
                      w
                  }
          })

############ evalObjective #################################

setMethod("evalObjective", signature("gelModels", "numeric", "missing"),
          function(object, theta, wObj, lambda, ...)
              {
                  gt <- evalMoment(object, theta)
                  k <- object@wSpec$k
                  if (is.null(object@gelType$fct))
                      rhoFct <- get(paste("rho",object@gelType$name,sep=""))
                  else
                      rhoFct <- object@gelType$fct
                  rho <- rhoFct(gmat=gt, lambda=lambda, derive = 0, k = k[1]/k[2])
                  n <- modelDims(object)$n
                  2*n*sum(rho)*k[2]/(k[1]^2*object@wSpec$bw)
              })

#########################  solveGel  #########################

setGeneric("solveGel", function(object, ...) standardGeneric("solveGel"))

setMethod("solveGel", signature("gelModels"),
          function(object, theta0=NULL, lambda0=NULL, lamSlv=NULL,
                   coefSlv=c("optim","nlminb","constrOptim"),
                   lControl=list(), tControl=list())
              {
                  coefSlv <- match.arg(coefSlv)
                  f <- function(theta, model, lambda0, slv, lcont,returnL=FALSE)
                      {
                          gt <- evalMoment(model, theta)
                          gelt <- model@gelType
                          k <- model@wSpec$k
                          args <- c(list(gmat=gt, l0=lambda0, gelType=gelt$name,
                                         rhoFct=gelt$fct), lcont, k=k[1]/k[2])
                          res <- do.call(slv, args)
                          if (returnL)
                              return(res)
                          evalObjective(model, theta, , res$lambda)
                      }
                  if (is.null(lambda0))
                      lambda0 <- rep(0, modelDims(object)$q)
                  if (is.null(theta0))
                      {
                          if (!("theta0"%in%slotNames(object)))
                              stop("Theta0 must be provided")
                          theta0 <- object@theta0
                      }
                  if (is.null(lamSlv))
                      lamSlv <- getLambda
                  if (coefSlv == "nlminb")
                      args <- c(list(start=theta0, objective=f,
                                     model=object, lambda0=lambda0,
                                     slv=lamSlv, lcont=lControl), tControl)
                  else
                      args <- c(list(par=theta0, fn=f, model=object, lambda0=lambda0,
                                     slv=lamSlv, lcont=lControl), tControl)
                  res <- do.call(get(coefSlv), args)
                  resl <- f(res$par,  object, lambda0, lamSlv, lControl, TRUE)
                  list(theta=res$par, convergence=res$convergence,
                       lambda=resl$lambda, lconvergence=resl$convergence)
          })


#########################  modelFit  #########################

setMethod("modelFit", signature("gelModels"), valueClass="gelfit", 
          definition = function(object, gelType=NULL, rhoFct=NULL,
              initTheta=c("gmm", "theta0"), start.tet=NULL,
              start.lam=NULL, ...)
              {
                  Call <- match.call()
                  initTheta = match.arg(initTheta)
                  if (!is.null(gelType))
                      object <- gmmToGel(as(object, "gmmModels"), gelType, rhoFct)
                  if (is.null(start.tet))
                      {
                          if (initTheta == "gmm")
                              start.tet <- modelFit(as(object, "gmmModels"))@theta
                          else if ("theta0" %in% slotNames(object))
                              start.tet <- object@theta0
                          else
                              stop("starting values is missing for the coefficient vector")
                      }
                  res <- solveGel(object, theta0=start.tet, lambda0=start.lam,
                                  ...)
                  
                  new("gelfit", theta=res$theta, convergence=res$convergence,
                      lconvergence=res$lconvergence$convergence,
                      lambda=res$lambda, call=Call, type=object@gelType$name,
                      model=object)
                  })


