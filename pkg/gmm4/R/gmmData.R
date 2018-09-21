######### Function to arrange the data for the gmmModel objects #################

.multiToSys <- function(formula, h, data, omit=TRUE)
{
    modelF <- model.frame(formula, data, na.action="na.pass",
                          drop.unused.levels=TRUE)
    Y <- model.response(modelF)
    modelF <- modelF[-1]
    Yn <- formula[[2]]
    Yn <- paste(Yn, ".", 1:ncol(Y), sep="")
    g <- lapply(1:length(Yn), function(i) {
        f <- formula
        f[[2]] <- as.symbol(Yn[i])
        f})
    colnames(Y) <- Yn
    modelF <- cbind(Y, modelF)
    if (any(class(h) == "formula"))
        {
            instF <- model.frame(h, data, na.action="na.pass",
                                 drop.unused.levels=TRUE)
        } else {
            h <- as.data.frame(h)
            chk <- apply(h, 2, function(x) all(x==x[1]))
            h <- h[, !chk]
            intercept <- any(chk)
            if (ncol(h) == 0)
                {                        
                    formula <- ~1
                } else {
                    if (is.null(colnames(h)))
                        colnames(h) <- paste("h", 1:ncol(h), sep="")
                    formh <- paste(colnames(h), collapse="+")
                    if (!intercept)
                        formh <- paste(formh, "-1", sep="")
                    formula <- as.formula(paste("~",formh))
                }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
        }
    h <- lapply(1:ncol(Y), function(i) formula(terms(instF), .GlobalEnv))
    data <- cbind(modelF, instF)
    data <- data[,!duplicated(colnames(data))]
    return(.slGmmData(g,h,data,omit))
}

.lGmmData <- function(formula, h, data, omit=TRUE)
    {
        modelF <- model.frame(formula, data, na.action="na.pass",
                              drop.unused.levels=TRUE)
        if (is.matrix(modelF[[1]]))
            return(.multiToSys(formula, h, data))
        parNames <- colnames(model.matrix(terms(modelF), modelF))
        k <- length(parNames)
        if (any(class(h) == "formula"))
            {
                instF <- model.frame(h, data, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        formula <- as.formula(paste("~",formh))
                    }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            }
        momNames <- colnames(model.matrix(terms(instF), instF))
        q <- length(momNames)
        isEndo <- !(parNames %in% momNames)
        na <- attr(na.omit(cbind(modelF, instF)), "na.action")
        if (!is.null(na) && omit)
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
        }
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, n=n, k=k, q=q, momNames=momNames,
             parNames=parNames, isEndo=isEndo, varNames=parNames, na.action=na)
    }



.formGmmData <- function(formula, tet0, data,omit=TRUE)
    {
        res <- lapply(formula, function(f) .nlGmmData(f, ~1, tet0, data))
        fRHS <- lapply(res, function(r) r$fRHS)
        fLHS <- lapply(res, function(r) r$fLHS)
        parNames <- res[[1]]$parNames
        varNames <- do.call("c", lapply(res, function(r) r$varNames))
        varNames <- unique(varNames)       
        chkLHS <- sapply(fLHS, function(r) any(all.vars(r) %in% names(tet0)))
        chkRHS <- sapply(fRHS, function(r) any(all.vars(r) %in% names(tet0)))
        isMDE <- all(chkLHS) |  all(chkRHS)        
        modelF <- sapply(varNames, function(n) data[[n]])
        modelF <- as.data.frame(modelF)
        na <- attr(na.omit(modelF), "na.action")
        if (!is.null(na) && omit)
            modelF <- modelF[-na,,drop=FALSE]
        k <- length(tet0)
        q <- length(formula)
        if (is.null(names(formula)))
            momNames <- paste("Mom_", 1:q, sep="")
        else
            momNames <- names(formula)
        isEndo <- rep(FALSE, length(varNames))
        n <- nrow(modelF)
        list(modelF=modelF,  fRHS=fRHS, fLHS=fLHS, n=n, k=k, q=q,
             momNames=momNames, parNames=parNames, varNames=varNames, isEndo=isEndo,
             isMDE=isMDE,na.action=na)
    }



.nlGmmData <- function(formula, h, tet0, data, omit=TRUE)
    {
        varNames <- all.vars(formula)
        parNames <- names(tet0)
        varNames <- varNames[!(varNames %in% parNames)]
        modelF <- try(sapply(varNames, function(n) data[[n]]), silent=TRUE)
        if (any(class(modelF)=="try-error"))
            stop("some variables are missing from data")
        modelF <- as.data.frame(modelF)        
        allVar <- c(as.list(modelF), as.list(tet0))
        k <- length(tet0)
        if (length(formula) == 3L)
        { 
            fLHS <- as.expression(formula[[2]])
            chk <- try(eval(fLHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the LHS")
            fRHS <- as.expression(formula[[3]])
            chk <- try(eval(fRHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the RHS")
        } else {
            fLHS <- NULL
            fRHS <- as.expression(formula[[2]])
            chk <- try(eval(fRHS, allVar))
            if (any(class(chk)=="try-error"))
                stop("Cannot evaluate the RHS")
        }
        if (any(class(h) == "formula"))
            {
                instF <- model.frame(h, data, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        formula <- as.formula(paste("~",formh))
                    }
                instF <- model.frame(formula, h, na.action="na.pass",
                                     drop.unused.levels=TRUE)
            }
        momNames <- colnames(model.matrix(terms(instF), instF))
        isEndo <- !(varNames %in% momNames)
        q <- length(momNames)
        na <- attr(na.omit(cbind(modelF, instF)), "na.action")
        if (!is.null(na) && omit)
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
        }
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, fRHS=fRHS, fLHS=fLHS, n=n, k=k, q=q,
             momNames=momNames, parNames=parNames, varNames=varNames, isEndo=isEndo,
             na.action=na)
    }

.fGmmData <- function(g, x, thet0, omit=NULL)
    {
        mom <- try(g(thet0, x))
        k <- length(thet0)        
        if (is.null(names(thet0)))
            parNames <- paste("tet", 1:k, sep="")
        else
            parNames <- names(thet0)
        if (any(class(mom)=="try-error"))
            {
                msg <- paste("Cannot evaluate the moments at thet0\n",
                             attr(mom,"conditon"))
                stop(msg)
            } else {
                q <-  ncol(mom)
                n <- nrow(mom)                
                if (!is.null(colnames(mom)))
                    momNames <- colnames(mom)
                else
                    momNames <- paste("h", 1:q, sep="")
            }
        list(q=q,n=n,k=k, momNames=momNames, parNames=parNames,
             varNames=character(), isEndo=logical())
    }

.slGmmData <- function(g,h,data,omit=TRUE)
    {
        res <- lapply(1:length(g), function(i) .lGmmData(g[[i]], h[[i]], data, FALSE))
        modelT <- lapply(res, function(x) terms(x$modelF))
        instT <-  lapply(res, function(x) terms(x$instF))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
        allDat <- na.omit(allDat)
        na <- attr(allDat, "na.action")
        if (omit && !is.null(na))
            allDat <- allDat[-na,]
        parNames <- lapply(1:length(g), function(i) res[[i]]$parNames)
        momNames <- lapply(1:length(g), function(i) res[[i]]$momNames)
        isEndo <- lapply(1:length(g), function(i) res[[i]]$isEndo)
        varNames <- lapply(1:length(g), function(i) res[[i]]$varNames)
        k <- sapply(parNames, length)
        q <- sapply(momNames, length)
        n <- nrow(allDat)
        if (!is.null(names(g)))
            eqnNames=names(g)
        else
            eqnNames <- paste("Eqn", 1:length(g), sep="")
        list(data=allDat, modelT=modelT, instT=instT, parNames=parNames,
             momNames=momNames, k=k,q=q,n=n, eqnNames=eqnNames,
             varNames=varNames, isEndo=isEndo, na.action=na)
    }

.snlGmmData <- function(g,h,tet0, data, omit=TRUE)
    {
        res <- lapply(1:length(g), function(i) .nlGmmData(g[[i]], h[[i]],
                                                          tet0[[i]], data, FALSE))
        fRHS <- lapply(res, function(x) x$fRHS)
        fLHS <- lapply(res, function(x) x$fLHS)
        instT <-  lapply(res, function(x) terms(x$instF))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
        allDat <- na.omit(allDat)
        na <- attr(allDat, "na.action")
        if (omit && !is.null(na))
            allDat <- allDat[-na,]
        parNames <- lapply(1:length(g), function(i) res[[i]]$parNames)
        momNames <- lapply(1:length(g), function(i) res[[i]]$momNames)
        isEndo <- lapply(1:length(g), function(i) res[[i]]$isEndo)
        varNames <- lapply(1:length(g), function(i) res[[i]]$varNames)
        k <- sapply(parNames, length)
        q <- sapply(momNames, length)
        n <- nrow(allDat)
        if (!is.null(names(g)))
            eqnNames=names(g)
        else
            eqnNames <- paste("Eqn", 1:length(g), sep="")
        list(data=allDat, fRHS=fRHS, fLHS=fLHS, parNames=parNames,
             momNames=momNames, k=k,q=q,n=n, eqnNames=eqnNames, instT=instT,
             varNames=varNames, isEndo=isEndo, na.action=na)
    }
