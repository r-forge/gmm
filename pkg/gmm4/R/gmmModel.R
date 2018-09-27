#############  Options for covariance matrix

.getVcovOptions <- function(type, data, addO=list())
    {
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
                if (!is.null(option$cluster))
                    {
                        if (!inherits(option$cluster,
                                      c("vector","data.frame","formula")))
                            stop("cluster must be a data.frame, a vector or a formula")
                        if (inherits(option$cluster, "formula"))
                            {
                                fn <- all.vars(option$cluster[[length(option$cluster)]])
                                option$cluster <- try(data[fn], silent=TRUE)
                                if (class(option$cluster) == "try-error")
                                    stop("variables in the cluster formula are not in data")
                            }
                        option$cluster <- as.data.frame(option$cluster)
                        if (is.null(colnames(option$cluster)))
                            colnames(option$cluster) <- paste("CL", 1:ncol(option$cluster),
                                                              sep="")
                    }
                if (option$type != "HC0")
                    stop("Only meatCL with type HC0 is allowed for GMM")
            } else {
                option <- list()
            }
        option
    }

.getSurvOptions <- function(data, opt=list())
    {
        if (length(opt) == 0)
            return(list())
        type <- c("sampling", "frequency")
        if (length(opt)>2 || !(names(opt) %in% c("type","weights")))
            stop("survOptions list must contain only two arguments: weights and type")
        opt$type <- match.arg(opt$type, type)
        if (!inherits(opt$weights, c("integer", "numeric", "formula")))
            stop("survey weights must be a numeric vector or a formula")
        if (inherits(opt$weights, "formula"))
            {
                fn <- all.vars(opt$weights[[length(opt$weights)]])
                if (length(fn)>1)
                    stop("weights must be a single variable")
                opt$weights <- try(c(data[[fn]]), silent=TRUE)
                if (class(opt$weights) == "try-error")
                    stop("variable in the weights formula is not in data")
            }
        opt
    }

##################  Constructor for the gmmModels Classes  #####################

gmmModel <- function(g, x=NULL, tet0=NULL,grad=NULL,
                     vcov = c("iid", "HAC", "MDS", "CL"),
                     vcovOptions=list(), centeredVcov = TRUE, data=parent.frame(),
                     na.action="na.omit", survOptions=list())
    {
        vcov <- match.arg(vcov)
        if (!is.list(vcovOptions) | !is.list(survOptions))
            stop("vcovOptions and survOptions must be a list")
        vcovOptions <- .getVcovOptions(vcov, data, vcovOptions)
        survOptions <- .getSurvOptions(data, survOptions)
        if (!is.list(data) && !is.environment(data)) 
            stop("'data' must be a list or an environment")    
        if (any(class(g)=="formula"))
            {
                chk <- names(tet0) %in% all.vars(g)
                if (length(chk) == 0 | all(!chk))
                    {
                        model <- .lGmmData(g,x,data, survOptions, vcovOptions, na.action)
                        if (!is.null(model$eqnNames))
                            gmodel <- new("slinearGmm", data = model$data,instT=model$instT, 
                                          modelT = model$modelT, vcov = vcov,
                                          vcovOptions=model$vcovOptions,
                                          centeredVcov=centeredVcov, 
                                          k = model$k, q = model$q, n = model$n,
                                          parNames = model$parNames, 
                                          momNames = model$momNames,eqnNames=model$eqnNames, 
                                          sameMom = TRUE, SUR = FALSE,
                                          varNames = model$varNames, 
                                          isEndo = model$isEndo, omit=model$omit,
                                          survOptions=model$survOptions)
                        else
                            gmodel <- new("linearGmm", modelF=model$modelF, 
                                          instF=model$instF,
                                          vcov=vcov, vcovOptions=model$vcovOptions,
                                          centeredVcov = centeredVcov, k=model$k,
                                          q=model$q, n=model$n, parNames=model$parNames,
                                          momNames=model$momNames, varNames=model$varNames,
                                          isEndo=model$isEndo, omit=model$omit,
                                          survOptions=model$survOptions)
                    } else {
                        if (!all(chk))
                            stop("All parameters in tet0 must be in g for nl Gmm")
                        model <- .nlGmmData(g, x, tet0, data, survOptions, vcovOptions,
                                            na.action)
                        gmodel <- new("nonlinearGmm", modelF=model$modelF, 
                                      instF=model$instF,theta0=tet0,fRHS=model$fRHS,
                                      fLHS=model$fLHS, vcov=vcov,
                                      vcovOptions=model$vcovOptions,
                                      centeredVcov = centeredVcov, k=model$k, q=model$q,
                                      n=model$n, parNames=model$parNames,
                                      momNames=model$momNames, varNames=model$varNames,
                                      isEndo=model$isEndo, omit=model$omit,
                                      survOptions=model$survOptions)
                    }
            } else if (class(g)=="function") {
                model <- .fGmmData(g, x, tet0, survOptions, vcovOptions, na.action)
                gmodel <- new("functionGmm", X=x, fct=g,
                              theta0=tet0, vcov=vcov,vcovOptions=model$vcovOptions,
                              centeredVcov = centeredVcov, k=model$k, q=model$q,
                              n=model$n, parNames=model$parNames,
                              momNames=model$momNames, varNames=model$varNames,
                              isEndo=model$isEndo, omit=model$omit, 
                              survOptions=model$survOptions)
            } else {
                if (!is.null(x))
                    stop("For formula GMM, x must be NULL. The moments are only defined as a list of formulas")
                if (class(g) != "list")
                    stop("For formula GMM, g must be a list of formulas")
                if (any(sapply(g, function(gi) class(gi)) != "formula"))
                    stop("For formula GMM, g must be a list of formulas")
                model <- .formGmmData(g, tet0, data, survOptions, vcovOptions, na.action)
                gmodel <- new("formulaGmm", modelF=model$modelF, 
                              vcov=vcov, theta0=tet0,fRHS=model$fRHS,
                              fLHS=model$fLHS,vcovOptions=model$vcovOptions,
                              centeredVcov = centeredVcov, k=model$k, q=model$q,
                              n=model$n, parNames=model$parNames,
                              momNames=model$momNames, varNames=model$varNames,
                              isEndo=model$isEndo, isMDE=model$isMDE, omit=model$omit,
                              survOptions=model$survOptions)
            }
        gmodel
    }

