
##################  Constructor for the sysGmmModels classes   #####################


sysGmmModel <- function(g, h=NULL, tet0=NULL,
                        vcov = c("iid", "HAC", "MDS"),
                        vcovOptions=list(), centeredVcov = TRUE, data=parent.frame())
    {
        vcov <- match.arg(vcov)
        if (!is.list(vcovOptions))
            stop("vcovOptions must be a list")
        vcovOptions <- do.call(.getVcovOptions, c(vcovOptions, type=vcov))
        if (!is.list(data) && !is.environment(data)) 
            stop("'data' must be a list or an environment")
        if (!is.list(g))
            stop("For system GMM, g must be lists of formulas")
        clg <- sapply(g, class)
        if (!all(clg=="formula"))
            stop("For system GMM, g must be lists of formulas")
        if (length(g) == 1)
            stop("There is only one equation. Use gmmModel() instead")
        ## Try to see if parameters are in formulas
        if (!is.null(tet0))
            {
                if (!is.list(tet0))
                    stop("For system GMM, tet0 must be a list of named vectors")
                if (length(tet0) != length(g))
                    stop("The length of tet0 must match the length of g and h")
                chk <- sapply(1:length(tet0),
                              function(i) isTRUE(all(names(tet0[[i]]) %in% all.vars(g[[i]]))))
                if (all(chk))
                    {
                        nonlin <- TRUE
                        chk2 <- sapply(tet0[-1],
                                       function(l) any(names(l) %in% names(tet0[[1]])))
                        if (any(chk2))
                            stop("Coefficient names across equations must be different")
                    } else if (all(!chk)){
                        nonlin <- FALSE
                    } else {
                        stop("g must all be linear or all nonlinear")
                    }
            } else {
                nonlin <- FALSE
            }
        varg <- lapply(1:length(g), function(i) .formSpec(g[[i]], tet0[[i]]))
        if (is.null(h))
        {
            SUR <- TRUE
            sameMom <- TRUE
            anyEndo <- rep(FALSE, length(g))
            chk <- sapply(varg, function(v) length(v$lhs)==0)
            if (any(chk))
                stop("For SUR, there must be a LHS")
            reg <- lapply(1:length(g), function(i) varg[[i]]$rhs)
            reg <- unique(do.call("c", reg)) 
            if (!nonlin)
                intercept <- any(sapply(varg, function(v) v$intercept == 1))
            else
                intercept <- TRUE
            reg <- paste(reg, collapse="+", sep="")
            reg <- paste("~", reg, ifelse(intercept, "", "-1"), sep="")
            reg <- as.formula(reg, .GlobalEnv)
            h <- rep(list(reg), length(g))
        } else {
            SUR <- FALSE
            if (class(h) == "formula")
                {
                    h <- rep(list(h), length(g))
                    sameMom <- TRUE
                } else if (!is.list(h)) {
                    stop("h must be a list or a formula")
                } else if (length(h) == 1) {
                    h <- rep(h, length(g))
                    sameMom <- TRUE
                } else if (length(h) != length(g)) {
                    stop("With different instruments, we must have length(h)= length(g)")
                } else {
                    sameMom <- FALSE
                }
            varh <- lapply(h, .formSpec)
            anyEndo <- !sapply(1:length(g),
                               function(i) all(varg[[i]]$rhs %in% varh[[i]]$rhs) &
                                           (varg[[i]]$intercept==varh[[i]]$intercept))
        }
        if (!nonlin)
        {
            model <- .slGmmData(g,h,data)
            isEndo <- lapply(1:length(g),
                             function(i) model$parNames[[i]] %in% model$momNames[[i]])
            gmodel <- new("slinearGmm", data=model$data, 
                          instT=model$instT, modelT=model$modelT,
                          vcov=vcov, vcovOptions=vcovOptions,
                          centeredVcov = centeredVcov, k=model$k,
                          q=model$q, n=model$n, parNames=model$parNames,
                          momNames=model$momNames, eqnNames=model$eqnNames,
                          sameMom=sameMom, SUR=SUR, varNames=model$varNames,
                          isEndo=model$isEndo, omit=model$na.action)
        } else {
            model <- .snlGmmData(g, h, tet0, data)
            gmodel <- new("snonlinearGmm", data=model$data, instT=model$instT,
                          theta0=tet0,fRHS=model$fRHS,eqnNames=model$eqnNames,
                          fLHS=model$fLHS, vcov=vcov, vcovOptions=vcovOptions,
                          centeredVcov = centeredVcov, k=model$k, q=model$q,
                          n=model$n, parNames=model$parNames,
                          momNames=model$momNames, sameMom=sameMom, SUR=SUR,
                          varNames=model$varNames, isEndo=model$isEndo, omit=model$na.action)
        }
        gmodel
    }


.formSpec <- function(f, theta=NULL)
{
    if (length(f) == 3)
    {  
        rhs <- all.vars(f[[3]])
        if (!is.null(theta))
            rhs <- rhs[!(rhs%in%names(theta))]
        lhs <- all.vars(f[[2]])
        if (!is.null(theta))
            lhs <- lhs[!(lhs%in%names(theta))]
    } else {
        lhs <- character()
        rhs <- all.vars(f[[2]])
        if (!is.null(theta))
            rhs <- rhs[!(rhs%in%names(theta))]
    }
    list(rhs=rhs, lhs=lhs, intercept=attr(terms(f), "intercept"))
}
