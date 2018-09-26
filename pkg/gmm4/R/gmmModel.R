#############  Options for covariance matrix

.getVcovOptions <- function(type, ...)
    {
        addO <- list(...)
        if (type == "HAC")
            {
                option <- list(kernel = "Quadratic Spectral",
                               crit = 1e-06,
                               bw = "Andrews", prewhite = 1L,
                               ar.method = "ols", approx = "AR(1)", 
                               tol = 1e-07)
                if (length(addO) > 0)
                    {
                        if (!all(names(addO) %in% names(option)))
                            stop(paste("Wrong options for vcov of type", type))
                        option[names(addO)] <- addO
                    }
                option$kernel <- match.arg(option$kernel,
                                           c("Quadratic Spectral", "Truncated", "Bartlett",
                                             "Parzen", "Tukey-Hanning"))
                if (!(option$ar.method %in% eval(as.list(args(ar))$method)))
                    stop("wrong value for ar.method")
                if (!(option$approx %in% eval(as.list(bwAndrews)$approx)))
                    stop("wrong value for approx")
                if (is.numeric(option$bw))
                    names(option$bw) <- "Fixed"
            } else if (type=="CL") {
                option <- list(cluster=NULL, type="HC0", cadjust=TRUE,
                               milti0=FALSE)
                if (length(addO) > 0)
                    {
                        if (!all(names(addO) %in% names(option)))
                            stop(paste("Wrong options for vcov of type", type))
                        option[names(addO)] <- addO
                    }
                if (option$type != "HC0")
                    stop("Only meatCL with type HC0 is allowed for GMM")
            } else {
                option <- list()
            }
        option
    }

##################  Constructor for the gmmModels Classes  #####################

gmmModel <- function(g, x=NULL, tet0=NULL,grad=NULL,
                     vcov = c("iid", "HAC", "MDS"),
                     vcovOptions=list(), centeredVcov = TRUE, data=parent.frame())
    {
        vcov <- match.arg(vcov)
        if (!is.list(vcovOptions))
            stop("vcovOptions must be a list")
        vcovOptions <- do.call(.getVcovOptions, c(vcovOptions, type=vcov))
        if (!is.list(data) && !is.environment(data)) 
            stop("'data' must be a list or an environment")    
        if (any(class(g)=="formula"))
            {
                chk <- names(tet0) %in% all.vars(g)
                if (length(chk) == 0 | all(!chk))
                    {
                        model <- .lGmmData(g,x,data)
                        if (!is.null(model$eqnNames))
                            gmodel <- new("slinearGmm", data = model$data,instT=model$instT, 
                                          modelT = model$modelT, vcov = vcov,
                                          vcovOptions=vcovOptions,centeredVcov=centeredVcov, 
                                          k = model$k, q = model$q, n = model$n,
                                          parNames = model$parNames, 
                                          momNames = model$momNames,eqnNames=model$eqnNames, 
                                          sameMom = TRUE, SUR = FALSE,
                                          varNames = model$varNames, 
                                          isEndo = model$isEndo, omit=model$na.action)
                        else
                            gmodel <- new("linearGmm", modelF=model$modelF, 
                                          instF=model$instF,
                                          vcov=vcov, vcovOptions=vcovOptions,
                                          centeredVcov = centeredVcov, k=model$k,
                                          q=model$q, n=model$n, parNames=model$parNames,
                                          momNames=model$momNames, varNames=model$varNames,
                                          isEndo=model$isEndo, omit=model$na.action)
                    } else {
                        if (!all(chk))
                            stop("All parameters in tet0 must be in g for nl Gmm")
                        model <- .nlGmmData(g, x, tet0, data)
                        gmodel <- new("nonlinearGmm", modelF=model$modelF, 
                                      instF=model$instF,theta0=tet0,fRHS=model$fRHS,
                                      fLHS=model$fLHS, vcov=vcov,vcovOptions=vcovOptions,
                                      centeredVcov = centeredVcov, k=model$k, q=model$q,
                                      n=model$n, parNames=model$parNames,
                                      momNames=model$momNames, varNames=model$varNames,
                                      isEndo=model$isEndo, omit=model$na.action)
                    }
            } else if (class(g)=="function") {
                model <- .fGmmData(g, x, tet0)
                gmodel <- new("functionGmm", X=x, fct=g,
                              theta0=tet0, vcov=vcov,vcovOptions=vcovOptions,
                              centeredVcov = centeredVcov, k=model$k, q=model$q,
                              n=model$n, parNames=model$parNames,
                              momNames=model$momNames, varNames=model$varNames,
                              isEndo=model$isEndo, omit=model$na.action)
            } else {
                if (!is.null(x))
                    stop("For formula GMM, x must be NULL. The moments are only defined as a list of formulas")
                if (class(g) != "list")
                    stop("For formula GMM, g must be a list of formulas")
                if (any(sapply(g, function(gi) class(gi)) != "formula"))
                    stop("For formula GMM, g must be a list of formulas")
                model <- .formGmmData(g, tet0, data)
                gmodel <- new("formulaGmm", modelF=model$modelF, 
                              vcov=vcov, theta0=tet0,fRHS=model$fRHS,
                              fLHS=model$fLHS,vcovOptions=vcovOptions,
                              centeredVcov = centeredVcov, k=model$k, q=model$q,
                              n=model$n, parNames=model$parNames,
                              momNames=model$momNames, varNames=model$varNames,
                              isEndo=model$isEndo, isMDE=model$isMDE, omit=model$na.action)
            }
        gmodel
    }

