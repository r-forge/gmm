################### the main gmm functions ###################
########## These functions ar to avoid having to builf model objects

gmm4 <- function (g, x, theta0 = NULL, grad = NULL, 
                  type = c("twostep", "iter", "cue", "onestep"),
                  vcov = c("iid", "HAC", "MDS", "TrueFixed", "CL"),
                  initW = c("ident", "tsls", "EbyE"), weights = "optimal", 
                  itermaxit = 50, cstLHS=NULL, cstRHS=NULL,
                  vcovOptions=list(),survOptions=list(),
                  itertol = 1e-07, centeredVcov = TRUE,
                  data = parent.frame(), ...) 
{
    Call <- match.call()
    vcov <- match.arg(vcov)
    type <- match.arg(type)
    initW <- match.arg(initW)
    if (vcov == "TrueFixed")
    {
        if (!is.matrix(weights) ||
            !(class(weights) %in% c("gmmWeights", "sysGmmWeigths")))
            stop("With TrueFixed vcov the weights must be provided")
        efficientWeights <- TRUE
        vcov2 <- "iid"
    } else {
        efficientWeights <- FALSE
        vcov2 <- vcov
    }
    if (is.list(g))
    {
        ## Formula or sysGMM? Need to find a better way.
        model <- NULL
        if (is.null(x) & !is.null(theta0))
            model <- try(gmmModel(g=g, x=x, theta0=theta0, grad=grad, vcov=vcov2,
                                  vcovOptions=vcovOptions,survOptions=survOptions,
                                  centeredVcov=centeredVcov, data=data), silent=TRUE)
        if (is.null(model) || class(model)=="try-error")
            model <- sysGmmModel(g=g, h=x, theta0=theta0, vcov=vcov2,
                                 vcovOptions=vcovOptions,survOptions=survOptions,
                                 centeredVcov=centeredVcov, data=data)
    } else {
        model <- gmmModel(g=g, x=x, theta0=theta0, grad=grad, vcov=vcov2,
                          vcovOptions=vcovOptions,survOptions=survOptions,
                          centeredVcov=centeredVcov, data=data)
        if (initW == "EbyE")
        {
            warning("initW cannot be EbyE for single equations, initW set to ident")
            initW="ident"
        }
    }
    if (!is.null(cstLHS))
        model <- restModel(model, cstLHS, cstRHS)
    fit <- modelFit(object=model, type=type, itertol=itertol, initW=initW,
                    weights=weights, itermaxit=itermaxit,
                    efficientWeights=efficientWeights, ...)
    fit@call <- Call
    fit
}


setMethod("tsls", "formula",
          function(object, x, vcov = c("iid", "HAC", "MDS", "CL"),                   
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))
              vcov <- match.arg(vcov)
              model <- gmmModel(g = object, x = x, vcov = vcov,
                                vcovOptions=vcovOptions,survOptions=survOptions,
                                centeredVcov = centeredVcov, data = data)
              obj <- tsls(model)
              obj@call <- Call
              obj
              })

setMethod("tsls", "list",
          function(object, x=NULL, vcov = c("iid", "HAC", "MDS", "CL"),
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))              
              vcov <- match.arg(vcov)
              model <- sysGmmModel(g = object, h = x, vcov = vcov,
                                   vcovOptions=vcovOptions,survOptions=survOptions,
                                   centeredVcov = centeredVcov, data = data)
              obj <- tsls(model)
              obj@call <- Call
              obj
              })


setMethod("ThreeSLS", "list",
          function(object, x=NULL, vcov = c("iid", "HAC", "MDS", "CL"),
                   vcovOptions=list(), survOptions=list(), centeredVcov = TRUE,
                   data = parent.frame())
          {
              Call <- match.call(call=sys.call(sys.parent()-1L))              
              vcov <- match.arg(vcov)
              model <- sysGmmModel(g = object, h = x, vcov = vcov,
                                   vcovOptions=vcovOptions,survOptions=survOptions,
                                   centeredVcov = centeredVcov, data = data)
              obj <- ThreeSLS(model)
              obj@call <- Call
              obj
              })
