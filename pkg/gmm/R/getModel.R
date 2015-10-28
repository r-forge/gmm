#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

getModel <- function(object, ...)
    {
        UseMethod("getModel")
    }

getModel.constGmm <- function(object, ...)
    {
        class(object) <- "baseGmm"
        obj <- getModel(object)
        if (!is.null(object$t0))
            {
                if (!is.null(dim(object$eqConst)))
                    stop("When t0 is provided, eqConst must be a vector which indicates which parameters to fix")
                if (length(object$eqConst)>=length(object$t0))
                    stop("Too many constraints; use evalGmm() if all coefficients are fixed")
                if (is.character(object$eqConst))
                    {
                        if (is.null(names(object$t0)))
                            stop("t0 must be a named vector if you want eqConst to be names")
                        if (any(!(object$eqConst %in% names(object$t0))))
                            stop("Wrong coefficient names in eqConst")
                        object$eqConst <- sort(match(object$eqConst,names(object$t0)))
                    }
                restTet <- object$t0[object$eqConst]
                obj$t0 <- object$t0[-object$eqConst]
                object$eqConst <- cbind(object$eqConst,restTet)
            }  else {
                if (is.null(dim(object$eqConst)))
                    stop("When t0 is not provided, eqConst must be a 2xq matrix")
            }
        attr(obj$x, "eqConst") <- list(eqConst = object$eqConst)
        rownames(attr(obj$x, "eqConst")$eqConst) <- obj$namesCoef[object$eqConst[,1]]
        object$eqConst <- attr(obj$x, "eqConst")$eqConst
        if(is(object$g, "formula"))
            {
                if (obj$x$ny>1)
                    stop("Constrained GMM not implemented yet for system of equations")
                if (obj$x$k<=0)
                    stop("Nothing to estimate")
            }
        obj$eqConst <- object$eqConst
        attr(obj$x, "k") <- attr(obj$x, "k")-nrow(object$eqConst)
        obj$namesCoef <- obj$namesCoef[-object$eqConst[,1]]
        obj$type <- paste(obj$type,"(with equality constraints)",sep=" ")	
        mess <- paste(rownames(object$eqConst), " = " , object$eqConst[,2], "\n",collapse="")
        mess <- paste("#### Equality constraints ####\n",mess,"##############################\n\n",sep="")
        obj$specMod <- mess
        return(obj)
    }

getModel.baseGmm <- function(object, ...)
    {
        object$allArg <- c(object, list(...))
        if(is(object$g, "formula"))
            {
                object$gradv <- .DmomentFct
                object$gradvf <- FALSE
                dat <- getDat(object$g, object$x, data = object$data)
                if(is.null(object$weightsMatrix))
                    {
                        clname <- paste(class(object), ".", object$type, ".formula", sep = "")
                    } else {    
                        clname <- "fixedW.formula"
                        object$type <- "One step GMM with fixed W"
                    }
                object$x <- dat
                object$gform<-object$g
                namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k), drop=FALSE])
                nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh), drop=FALSE]) 
                if (dat$ny > 1)
                    {
                        namey <- colnames(dat$x[,1:dat$ny, drop=FALSE])
                        object$namesCoef <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
                        object$namesgt <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
                    } else {
                        object$namesCoef <- namex
                        object$namesgt <- nameh
                    }
                attr(object$x,"ModelType") <- "linear"
                attr(object$x, "k") <- object$x$k
                attr(object$x, "q") <- object$x$ny*object$x$nh
                attr(object$x, "n") <- NROW(object$x$x)
            } else {
                attr(object$x,"ModelType") <- "nonlinear"
                attr(object$x, "momentfct") <- object$g
                attr(object$x, "k") <- length(object$t0)
                attr(object$x, "q") <- NCOL(object$g(object$t0, object$x))
                attr(object$x, "n") <- NROW(object$x)
                if(is.null(names(object$t0)))
                    object$namesCoef <- paste("Theta[" ,1:attr(object$x, "k"), "]", sep = "")
                else
                    object$namesCoef <- names(object$t0)
                if(is.null(object$weightsMatrix))
                    {
                        clname <- paste(class(object), "." ,object$type, sep = "")
                    } else {
                        clname <- "fixedW"
                        object$type <- "One step GMM with fixed W"
                        attr(object$x, "weight")$w <- object$weightsMatrix
                    }
                if (!is.function(object$gradv))
                    { 
                        object$gradvf <- FALSE
                    } else {
                        attr(object$x, "gradv") <- object$gradv    
                        object$gradvf <- TRUE
                    }
                object$gradv <- .DmomentFct
            }   
        object$TypeGmm <- class(object)
        attr(object$x, "weight") <- list(w=object$weightsMatrix,
                                         centeredVcov=object$centeredVcov)
        attr(object$x, "weight")$WSpec <- list()
        attr(object$x, "weight")$WSpec$sandwich <- list(kernel = object$kernel, bw = object$bw,
                                                        prewhite = object$prewhite,
                                                        ar.method = object$ar.method,
                                                        approx = object$approx, tol = object$tol)
        attr(object$x, "weight")$vcov <- object$vcov
        object$g <- .momentFct
        class(object)  <- clname
        return(object)
    }

getModel.constGel <- function(object, ...)
    {
        class(object) <- "baseGel"
        obj <- getModel(object)
        if (!is.null(dim(object$eqConst)))
            stop("eqConst must be a vector which indicates which parameters to fix")
	if (length(object$eqConst)>=length(object$tet0))
		stop("Too many constraints; use evalGel() if all coefficients are fixed")
        if (is.character(object$eqConst))
            {
		if (is.null(names(object$tet0)))
                    stop("tet0 must be a named vector if you want eqConst to be names")
		if (any(!(object$eqConst %in% names(object$tet0))))
                    stop("Wrong coefficient names in eqConst")
		object$eqConst <- sort(match(object$eqConst,names(object$tet0)))
		}
	restTet <- object$tet0[object$eqConst]
	obj$tet0 <- object$tet0[-object$eqConst]
	object$eqConst <- cbind(object$eqConst,restTet)    
        attr(obj$x, "eqConst") <- list(eqConst = object$eqConst)
        rownames(attr(obj$x, "eqConst")$eqConst) <- obj$namesCoef[object$eqConst[,1]]
        object$eqConst <- attr(obj$x, "eqConst")$eqConst
        if(is(object$g, "formula"))
            {
                if (obj$x$ny>1)
                    stop("Constrained GMM not implemented yet for system of equations")
            }
        obj$eqConst <- object$eqConst
        attr(obj$x, "k") <- attr(obj$x, "k")-nrow(object$eqConst)
        obj$namesCoef <- obj$namesCoef[-object$eqConst[,1]]
        obj$type <- paste(obj$type,"(with equality constraints)",sep=" ")	
        mess <- paste(rownames(object$eqConst), " = " , object$eqConst[,2], "\n",collapse="")
        mess <- paste("#### Equality constraints ####\n",mess,"##############################\n\n",sep="")
        obj$specMod <- mess
        return(obj)
    }

getModel.baseGel <- function(object, ...)
    {
    object$allArg <- c(object, list(...))
    if(is(object$g, "formula"))
        {
            dat <- getDat(object$g, object$x, data = object$data)
            k <- dat$k            
            if (is.null(object$tet0))
                {
                    if (!is.null(object$eqConst))
                        stop("You have to provide tet0 with equality constrains")
                    if (object$optfct == "optimize")
                        stop("For optimize, you must provide the 2x1 vector tet0")
                    res0 <- gmm(object$g, object$x, data=object$data)
                    object$tet0 <- res0$coefficients
                    if (object$smooth)
                        gt <- res0$gt
                } else {
                    if (object$optfct == "optimize")
                        {
                            if (k != 1)
                                stop("optimize() is for univariate optimization")
                            if (length(object$tet0) != 2)
                                stop("For optimize(), tet0 must be a 2x1 vector")
                        } else {
                            if (k != length(object$tet0))
                                stop("The number of starting values does not correspond to the number of regressors")
                        }                
                    if (object$smooth)
                        gt <- gmm(object$g, object$x, data=object$data)$gt
                }                
            clname <- paste(class(object), ".modFormula", sep = "")
            object$gradv <- .DmomentFct
            object$gradvf <- FALSE
            object$x <- dat
            object$gform<-object$g
            namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k), drop=FALSE])
            nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh), drop=FALSE])
            if (dat$ny > 1)
                {
                    namey <- colnames(dat$x[,1:dat$ny])
                    namesCoef <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
                    object$namesgt <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
                } else {
                    namesCoef <- namex
                    object$namesgt <- nameh
                }
            if (is.null(names(object$tet0)))
                object$namesCoef <- namesCoef
            else
                object$namesCoef <- names(object$tet0)
            attr(object$x,"ModelType") <- "linear"
            attr(object$x, "k") <- k
            attr(object$x, "q") <- object$x$ny*object$x$nh
            attr(object$x, "n") <- NROW(object$x$x)
        } else {
            if (is.null(object$tet0))
                stop("You must provide the starting values tet0 for nonlinear moments")
            if(any(object$optfct == c("optim", "nlminb")))
                k <- length(object$tet0)
            else
                k <- 1                    
            attr(object$x,"ModelType") <- "nonlinear"
            attr(object$x, "momentfct") <- object$g
            attr(object$x, "k") <- k
            attr(object$x, "q") <- NCOL(object$g(object$tet0, object$x))
            attr(object$x, "n") <- NROW(object$x)
            if(is.null(names(object$tet0)))
                object$namesCoef <- paste("Theta[" ,1:attr(object$x, "k"), "]", sep = "")
            else
                object$namesCoef <- names(object$tet0)
            if (!is.function(object$gradv) | object$smooth)
                { 
                    object$gradvf <- FALSE
                } else {
                    attr(object$x, "gradv") <- object$gradv    
                    object$gradvf <- TRUE                            
                }
            object$gradv <- .DmomentFct
            if (object$smooth)
                gt <- gmm(object$g, object$x, object$tet0, wmatrix = "ident", ...)$gt
            clname <- paste(class(object), ".mod", sep = "")
        }
    if (object$smooth)
    {
        if (is.function(object$gradv))
            warning("The provided gradv is not used when smooth=TRUE",
                    call. = FALSE)		
        if(object$kernel == "Truncated")
            {
                object$wkernel <- "Bartlett"
                object$k1 <- 2
                object$k2 <- 2
            }
        if(object$kernel == "Bartlett")
            {
                object$wkernel <- "Parzen"
                object$k1 <- 1
                object$k2 <- 2/3
            }
        gt <- scale(gt, scale=FALSE)
        class(gt) <- "gmmFct"
        if (is.function(object$bw))
	    object$bwVal <- object$bw(gt, kernel = object$wkernel, prewhite = object$prewhite, 
	               ar.method = object$ar.method, approx = object$approx)
        else
	    object$bwVal <- object$bw
        object$w <- smoothG(gt, bw = object$bwVal, kernel = object$kernel, tol = object$tol_weights)$kern_weights
        attr(object$x,"smooth") <- list(bw=object$bwVal, w=object$w, kernel=object$kernel)
    } else {
        object$k1 <- 1
        object$k2 <- 1
        object$w <- kernel(1)
        object$bwVal <- 1
    }
    object$g <- .momentFct
    object$CGEL <- object$alpha
    class(object) <- clname
    return(object)
}

