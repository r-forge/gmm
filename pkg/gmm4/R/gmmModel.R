


##################  Constructor for the gmmModels Classes  #####################

gmmModel <- function(g, x, tet0=NULL,grad=NULL,
                     vcov = c("HAC", "MDS", "iid"),
                     kernel = c("Quadratic Spectral",  "Truncated", "Bartlett", "Parzen",
                         "Tukey-Hanning"), crit = 1e-06,
                     bw = "Andrews", prewhite = 1L, ar.method = "ols", approx = "AR(1)", 
                     tol = 1e-07, centeredVcov = TRUE, data=parent.frame())
    {
        vcov <- match.arg(vcov)
        kernel <- match.arg(kernel)
        if (is.numeric(bw))
            names(bw) <- "Fixed"
        if (!is.list(data) && !is.environment(data)) 
            stop("'data' must be a list or an environment")
        if (any(class(g)=="formula"))
            {
                chk <- names(tet0) %in% all.vars(g)
                if (length(chk) == 0 | all(!chk))
                    {
                        model <- .lGmmData(g,x,data)
                        gmodel <- new("linearGmm", modelF=model$modelF, 
                                      instF=model$instF,
                                      vcov=vcov, kernel=kernel, bw=bw,
                                      prewhite=as.integer(prewhite),
                                      ar.method=ar.method, approx=approx, tol=tol,
                                      centeredVcov = centeredVcov, k=model$k,
                                      q=model$q, n=model$n, parNames=model$parNames,
                                      momNames=model$momNames, varNames=model$varNames,
                                      isEndo=model$isEndo)
                    } else {
                        if (!all(chk))
                            stop("All parameters in tet0 must be in g for nl Gmm")
                        model <- .nlGmmData(g, x, tet0, data)
                        gmodel <- new("nonlinearGmm", modelF=model$modelF, 
                                      instF=model$instF,theta0=tet0,fRHS=model$fRHS,
                                      fLHS=model$fLHS, vcov=vcov, kernel=kernel, bw=bw,
                                      prewhite=as.integer(prewhite),
                                      ar.method=ar.method, approx=approx, tol=tol,
                                      centeredVcov = centeredVcov, k=model$k, q=model$q,
                                      n=model$n, parNames=model$parNames,
                                      momNames=model$momNames, varNames=model$varNames,
                                      isEndo=model$isEndo)
                    }
            } else {
                model <- .fGmmData(g, x, tet0)
                gmodel <- new("functionGmm", X=x, fct=g,
                              theta0=tet0, vcov=vcov, kernel=kernel, bw=bw,
                              prewhite=as.integer(prewhite),dfct=grad,
                              ar.method=ar.method, approx=approx, tol=tol,
                              centeredVcov = centeredVcov, k=model$k, q=model$q,
                              n=model$n, parNames=model$parNames,
                              momNames=model$momNames, varNames=model$varNames,
                              isEndo=model$isEndo)
            }
        gmodel
    }


