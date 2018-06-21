################### the main gmm functions ###################
########## These functions ar to avoid having to builf model objects



gmm4 <- function (g, x, tet0 = NULL, grad = NULL, 
                  type = c("twostep", "iter", "cue", "onestep"),
                  vcov = c("MDS", "HAC", "iid", "TrueFixed"),
                  initW = c("ident", "tsls", "EbyE"), weights = "optimal", 
                  itermaxit = 50, cstLHS=NULL, cstRHS=NULL, 
                  kernel = c("Quadratic Spectral", "Truncated",
                      "Bartlett", "Parzen", "Tukey-Hanning"), crit = 1e-06, 
                  bw = "Andrews", prewhite = 1L, ar.method = "ols", approx = "AR(1)", 
                  kerntol = 1e-07, itertol = 1e-07, centeredVcov = TRUE,
                  data = parent.frame(), ...) 
{
    Call <- match.call()
    vcov <- match.arg(vcov)
    kernel <- match.arg(kernel)
    type <- match.arg(type)
    initW <- match.arg(initW)
    if (vcov == "TrueFixed")
        {
            if (!is.matrix(weights) ||
                !(class(weights) %in% c("gmmWeights", "sysGmmWeigths")))
                stop("With TrueFixed vcov the weights must be provided")
            efficientWeights <- TRUE
        } else {
            efficientWeights <- FALSE
        }
    if (is.list(g))
        {
            model <- sysGmmModel(g=g, h=x, tet0=tet0, vcov=vcov,
                                 kernel=kernel, crit=crit, bw=bw, prewhite=prewhite,
                                 ar.method=ar.method, approx=approx, tol=kerntol,
                                 centeredVcov=centeredVcov, data=data)
        } else {
            model <- gmmModel(g=g, x=x, tet0=tet0, grad=grad, vcov=vcov,
                              kernel=kernel, crit=crit, bw=bw, prewhite=prewhite,
                              ar.method=ar.method, approx=approx, tol=kerntol,
                              centeredVcov=centeredVcov, data=data)
            if (initW == "EbyE")
                {
                    warning("initW cannot be EbyE for single equations, initW set to ident")
                    initW="ident"
                }
        }
    if (!is.null(cstLHS))
        model <- restGmmModel(model, cstLHS, cstRHS)

    fit <- gmmFit(object=model, type=type, itertol=itertol, initW=initW,
                  weights=weights, itermaxit=itermaxit,
                  efficientWeights=efficientWeights, ...)
    fit@call <- Call
    fit
}


setMethod("tsls", "formula",
          function(object, x, vcov = c("iid", "HAC", "MDS"),
                   kernel = c("Quadratic Spectral", "Truncated", "Bartlett", 
                       "Parzen", "Tukey-Hanning"), crit = 1e-06, bw = "Andrews", 
                   prewhite = 1L, ar.method = "ols", approx = "AR(1)", kerntol = 1e-07, 
                   centeredVcov = TRUE, data = parent.frame())
              {
                  vcov <- match.arg(vcov)
                  kernel <- match.arg(kernel)
                  model <- gmmModel(g = object, x = x, vcov = vcov, 
                                    kernel = kernel, crit = crit, bw = bw,
                                    prewhite = prewhite, ar.method = ar.method,
                                    approx = approx, tol = kerntol, 
                                    centeredVcov = centeredVcov, data = data)
                  tsls(model)
              })


setMethod("tsls", "list",
          function(object, x=NULL, vcov = c("iid", "HAC", "MDS"),
                   kernel = c("Quadratic Spectral", "Truncated", "Bartlett", 
                       "Parzen", "Tukey-Hanning"), crit = 1e-06, bw = "Andrews", 
                   prewhite = 1L, ar.method = "ols", approx = "AR(1)", kerntol = 1e-07, 
                   centeredVcov = TRUE, data = parent.frame())
              {
                  vcov <- match.arg(vcov)
                  kernel <- match.arg(kernel)
                  model <- sysGmmModel(g = object, h = x, vcov = vcov, 
                                    kernel = kernel, crit = crit, bw = bw,
                                    prewhite = prewhite, ar.method = ar.method,
                                    approx = approx, tol = kerntol, 
                                    centeredVcov = centeredVcov, data = data)
                  tsls(model)
              })


setMethod("ThreeSLS", "list",
          function(object, x=NULL, vcov = c("iid", "HAC", "MDS"),
                   kernel = c("Quadratic Spectral", "Truncated", "Bartlett", 
                       "Parzen", "Tukey-Hanning"), crit = 1e-06, bw = "Andrews", 
                   prewhite = 1L, ar.method = "ols", approx = "AR(1)", kerntol = 1e-07, 
                   centeredVcov = TRUE, data = parent.frame())
              {
                  vcov <- match.arg(vcov)
                  kernel <- match.arg(kernel)
                  model <- sysGmmModel(g = object, h = x, vcov = vcov, 
                                       kernel = kernel, crit = crit, bw = bw,
                                       prewhite = prewhite, ar.method = ar.method,
                                       approx = approx, tol = kerntol, 
                                       centeredVcov = centeredVcov, data = data)
                  ThreeSLS(model)
              })
