######### Function to arrange the data for the gmmModel objects #################

.lGmmData <- function(formula, h, data)
    {
        mf <- match.call()
        m <- match(c("formula", "data"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf$na.action <- "na.pass"
        mfh <- mf
        mf[[1L]] <- quote(stats::model.frame)
        modelF <- eval(mf, parent.frame())        
        parNames <- colnames(model.matrix(terms(modelF), modelF))
        k <- length(parNames)
        if (any(class(h) == "formula"))
            {
                mfh$formula <- h
                mfh[[1L]] <- quote(stats::model.frame)
                instF <- eval(mfh, parent.frame())
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        mfh$formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        mfh$formula <- as.formula(paste("~",formh))
                        mfh$data <- quote(h)
                    }
                mfh[[1L]] <- quote(stats::model.frame)
                instF <- eval(mfh)
            }
        momNames <- colnames(model.matrix(terms(instF), instF))
        q <- length(momNames)
        isEndo <- !(parNames %in% momNames)
        na <- attr(na.omit(cbind(modelF, instF)), "na.action")
        if (!is.null(na))
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
        }
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, n=n, k=k, q=q, momNames=momNames,
             parNames=parNames, isEndo=isEndo, varNames=parNames)
    }

.nlGmmData <- function(formula, h, tet0, data)
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
            mfh <- match.call()
            m <- match(c("formula", "data"), names(mfh), 0L)
            mfh <- mfh[c(1L, m)]
            mfh$drop.unused.levels <- TRUE
            mfh$na.action <- "na.pass"
            mfh$formula <- h
            mfh[[1L]] <- quote(stats::model.frame)
            instF <- eval(mfh, parent.frame())
            } else {
                h <- as.data.frame(h)
                chk <- apply(h, 2, function(x) all(x==x[1]))
                h <- h[, !chk]
                intercept <- any(chk)
                if (ncol(h) == 0)
                    {                        
                        mfh$formula <- ~1
                    } else {
                        if (is.null(colnames(h)))
                            colnames(h) <- paste("h", 1:ncol(h), sep="")
                        formh <- paste(colnames(h), collapse="+")
                        if (!intercept)
                            formh <- paste(formh, "-1", sep="")
                        mfh$formula <- as.formula(paste("~",formh))
                        mfh$data <- quote(h)
                    }
                mfh[[1L]] <- quote(stats::model.frame)
                instF <- eval(mfh)
            }
        momNames <- colnames(model.matrix(terms(instF), instF))
        isEndo <- !(varNames %in% momNames)
        q <- length(momNames)
        na <- attr(na.omit(cbind(modelF, instF)), "na.action")
        if (!is.null(na))
        {
            modelF <- modelF[-na,,drop=FALSE]
            instF <- instF[-na,,drop=FALSE]
        }
        n <- nrow(modelF)
        list(modelF=modelF,  instF=instF, fRHS=fRHS, fLHS=fLHS, n=n, k=k, q=q,
             momNames=momNames, parNames=parNames, varNames=varNames, isEndo=isEndo)
    }

.fGmmData <- function(g, x, theta0)
    {
        mom <- try(g(theta0, x))
        k <- length(theta0)        
        if (is.null(names(theta0)))
            parNames <- paste("tet", 1:k, sep="")
        else
            parNames <- names(theta0)
        if (any(class(mom)=="try-error"))
            {
                msg <- paste("Cannot evaluate the moments at theta0\n",
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

.slGmmData <- function(g,h,data)
    {
        res <- lapply(1:length(g), function(i) .lGmmData(g[[i]], h[[i]], data))
        modelT <- lapply(res, function(x) terms(x$modelF))
        instT <-  lapply(res, function(x) terms(x$instF))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
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
             varNames=varNames, isEndo=isEndo)
    }

.snlGmmData <- function(g,h,tet0, data)
    {
        res <- lapply(1:length(g), function(i) .nlGmmData(g[[i]], h[[i]],
                                                          tet0[[i]], data))
        fRHS <- lapply(res, function(x) x$fRHS)
        fLHS <- lapply(res, function(x) x$fLHS)
        instT <-  lapply(res, function(x) terms(x$instF))
        allDat <-  do.call(cbind, lapply(res, function(x) cbind(x$modelF, x$instF)))
        allDat <- allDat[,!duplicated(colnames(allDat))]
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
             varNames=varNames, isEndo=isEndo)
    }
