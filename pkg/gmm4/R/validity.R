#############  Validity Methods  ##################
#### All validity Methods for all objects are here

.checkBased <- function(object)
    {
        error <- character()
        vcov <- c("HAC","MDS","iid")
        kernel <- c("Quadratic Spectral", "Truncated",
                    "Bartlett", "Parzen", "Tukey-Hanning")
        bw <- c("Andrews","NeweyWest","Wilhelm")
        method <- as.list(args(ar))$method
        approx <- as.list(bwAndrews)$approx        
        if ( !(object@vcov%in%vcov))
            {
                vcov <- paste(vcov, collapse=", ")
                msg <- paste("vcov must be one of ",
                             vcov, "/n", sep="")
                error <- c(error, msg)
            }        
        if ( !(object@kernel%in%kernel))
            {
                kernel <- paste(kernel, collapse=", ")
                msg <- paste("kernel must be one of ",
                             kernel, "/n", sep="")
                error <- c(error, msg)
            }
        if (!is.numeric(object@bw))
            if ( !(object@bw%in%bw))
                {
                bw <- paste(bw, collapse=", ")
                msg <- paste("bw must be one of ",
                             bw, "/n", sep="")
                error <- c(error, msg)
            }
        if ( !(object@ar.method%in%eval(method)))
            {
                method <- paste(method, collapse=", ")
                msg <- paste("ar.method must be one of ",
                             method, "/n", sep="")
                error <- c(error, msg)
            }
        if ( !(object@approx%in%eval(approx)))
            {
                approx <- paste(approx, collapse=", ")
                msg <- paste("approx must be one of ",
                             approx, "/n", sep="")
                error <- c(error, msg)
            }
                error
    }

.checkLinGmm <- function(object)
    {
        error <- .checkBased(object)
        ny <- NCOL(object@modelF[[1L]])
        nM <- nrow(object@modelF)
        nI <- nrow(object@instF)
        k <- ncol(model.matrix(object))
        q <- ncol(model.matrix(object, "instrument"))
        if (k != object@k)
            {
                msg <- "k does not correspond to the number of regressor in modelF"
                error <- c(error, msg)
            }
        if (q != object@q)
            {
                msg <- "q does not correspond to the number of instruments in intF"
                error <- c(error, msg)
            }        
        if (ny > 1)
            {
                msg <- "This class only accepts one response variable"
                error <- c(error, msg)
            }
        if (!all(object@n == c(nM, nI)))
            {
                msg <- "n must be equal to the number of observations in modelF and instF"
                error <- c(error, msg)
            }            
        if (q<k | nM<q)
            {
                msg <- "The model is under-identified"
                error <- c(error, msg)
            }
        if (is.null(attr(object@modelF, "terms")))
            {
                msg("modelF does not have a terms attribute")
                error <- c(error, msg)
            }
        if (is.null(attr(object@instF, "terms")))
            {
                msg("instF does not have a terms attribute")
                error <- c(error, msg)
            }
        if (q != length(object@momNames))
            {
                msg <- "the length of momNames is not equal to q"
                error <- c(error, msg)                
            }
        if (k != length(object@parNames))
            {
                msg <- "the length of parNames is not equal to k"
                error <- c(error, msg)                
            }
        if (length(error)>0)
            error
        else
            TRUE
    }

setValidity("linearGmm", .checkLinGmm)

.checkRestLGmm <- function(object)
    {
        error <- character()
        R <- object@cstLHS
        rhs <- object@cstRHS
        spec <- object@cstSpec
        if (!all(names(spec) %in% c("theta","minY","newParNames","parNames",
                                    "originParNames","k", "isEndo")))
            {
                msg <- "cstSpec have missing arguments"
                error <- c(error, msg)
            }    
        if (dim(R)[1] != length(rhs))
            {
                msg <- "Number of rows of cstLHS must equal the length of cstRHS"
                error <- c(error, msg)
            }
        if (dim(R)[2] != length(spec$originParNames))
            {
                msg <- paste("Number of columns of cstLHS must ",
                             "equal the number of coefficients ",
                             "in the unconstrained model", sep="")
                error <- c(error, msg)
            }
        if (dim(R)[1] > length(spec$originParNames))
            {
                msg <- "Too many constraints"
                error <- c(error, msg)
            }
        if (qr(R)$rank != dim(R)[1])
            {
                msg <- "The Hypothesis matrix if not full rank"
                error <- c(error, msg)
            }
        #res <- try(.imposeRestrict(R,rhs, spec$originParNames), silent=TRUE)
        #if (any(class(res)=="try-error"))
        #    {
        #        msg <- "Fail to generate constraint specifications"
        #        error <- c(error, msg)
        #    }
        #if (!identical(res[-which(names(res)=="isEndo")], spec))
        #    {
        #        msg <- "cstSpec does not seem to correspond to cstLHS and cstRHS"
        #        error <- c(error, msg)            
        #    }    
        if (length(error)==0)
            TRUE
        else
            error              
    }

setValidity("rlinearGmm", .checkRestLGmm)

.checkNLinGmm <- function(object)
    {
        error <- .checkBased(object)
        nM <- nrow(object@modelF)
        nI <- nrow(object@instF)
        k <-  length(object@theta0)
        if (k != object@k)
            {
                msg <- "k does not correspond to the number of resgressor in modelF"
                error <- c(error, msg)
            }
        if (!all(object@n == c(nM, nI)))
            {
                msg <- paste("n must be equal to the number of observations ",
                             "in modelF and instF",sep="")
                error <- c(error, msg)
            }            
        mom <- try(evalMoment(object, object@theta0))        
        if (any(class(mom)=="try-error"))
            {
                msg <- paste("Cannot evaluate the moments at theta0\n",
                             attr(mom,"conditon"))
                error <- c(error, msg)
            } else {
                q <-  ncol(mom)
                if (q != object@q)
                    {
                        msg <- paste("q does not correspond to the number of ",
                                     "instruments in intF",sep="")
                        error <- c(error, msg)
                    }
                if (q<k | object@n<q)
                    {
                        msg <- "The model is under-identified"
                        error <- c(error, msg)
                    }
            }                
        if (is.null(attr(object@instF, "terms")))
            {
                msg("instF does not have a terms attribute")
                error <- c(error, msg)
            }
        if (q != length(object@momNames))
            {
                msg <- "the length of momNames is not equal to q"
                error <- c(error, msg)                
            }
        if (k != length(object@parNames))
            {
                msg <- "the length of parNames is not equal to k"
                error <- c(error, msg)                
            }
        allVar <- c(all.vars(object@fRHS), all.vars(object@fLHS))
        if (!all(object@parNames %in% allVar))
            {
                msg <- "For nl Gmm all parameters (names) must be in the formula g"
                error <- c(error, msg)                
            }            
        if (length(error)==0)
            TRUE
        else
            error
    }

setValidity("nonlinearGmm", .checkNLinGmm)

.checkfGmm <- function(object)
    {
        error <- .checkBased(object)
        k <-  length(object@theta0)
        if (k != object@k)
            {
                msg <- "k does not correspond to the number of resgressor in modelF"
                error <- c(error, msg)
            }
        mom <- try(object@fct(object@theta0, object@X))
        if (!is.null(object@dfct))
            {
                dmom <- try(object@dfct(object@theta0, object@X))
                if (any(class(dmom)=="try-error"))
                    {
                        msg <- paste("Cannot evaluate the provided derivatives of",
                                     "moments at theta0\n",
                                     attr(mom,"conditon"))
                        error <- c(error, msg)
                    } else {
                        if (ncol(dmom) != object@k | nrow(dmom) != object@q)
                            {
                                msg <- paste("The dimention of the gradian of the ",
                                             "moments is not qxk",sep="")
                                error <- c(error, msg)
                            }
                    }
            }
        if (any(class(mom)=="try-error"))
            {
                msg <- paste("Cannot evaluate the moments at theta0\n",
                             attr(mom,"conditon"))
                error <- c(error, msg)
            } else {
                q <-  ncol(mom)
                n <- nrow(mom)
                if (q != object@q)
                    {
                        msg <- paste("q does not correspond to the number of ",
                                     "instruments in intF",sep="")
                        error <- c(error, msg)
                    }
                if (q<k | n<q)
                    {
                        msg <- "The model is under-identified"
                        error <- c(error, msg)
                    }
            }
        if (q != length(object@momNames))
            {
                msg <- "the length of momNames is not equal to q"
                error <- c(error, msg)                
            }
        if (k != length(object@parNames))
            {
                msg <- "the length of parNames is not equal to k"
                error <- c(error, msg)                
            }        
        if (length(error)==0)
            TRUE
        else
            error
    }

setValidity("functionGmm", .checkfGmm)


.checkGmmWeights <- function(object)
    {
        error <- character()
        if (is.character(object@w))
            {
                if (object@w != "ident")
                    {
                        msg <- "Only 'ident' is allowed when w is a character type"
                        error <- c(error, msg)
                    } else {
                        if (object@type != "weights")
                            {
                                msg <- "type must be 'weights' when w='ident'"
                                error <- c(error, msg)                                
                            }
                    }
            } else {
                if (object@type == "qr")
                    {
                        if (!(class(object@w) == "qr"))
                            {
                                msg <- paste("w must be a 'qr' type object when ",
                                             "type is a 'qr'",sep="")
                                error <- c(error, msg)                                
                            }
                    } else if (object@type == "weights" | object@type == "vcov") {
                        if (!is.matrix(object@w) |
                            dim(object@w)[1L] != dim(object@w)[2L])
                            {
                                msg <- paste("With type='weights' or 'vcov', w must ",
                                             "be a square matrix",sep="")
                                error <- c(error, msg)                                
                            }
                    } else if (object@type == "chol") {
                        if (!all(object@w[lower.tri(object@w)] == 0))
                            {
                                msg <- "Cholesky matrix must be upper triangular"
                                error <- c(error, msg)                                

                            } else {
                                if (!(all(diag(object@w)>0)))
                                    {
                                        msg <- paste("Cholesky matrix must be ",
                                                     "positive definite",sep="")
                                        error <- c(error, msg)
                                    }
                            }
                    } else {
                        msg <- "type must be one of 'weights', 'vcov', 'qr', or 'chol'"
                        error <- c(error, msg)
                    }
            }
        if (length(object@HAC)>0)
            {
                if (!all(names(object@HAC) %in% c("bw", "kernel","weights")))
                    {
                        msg <- "HAC must contain 'bw', 'kernel', and 'weights'"
                        error <- c(error, msg)
                    }
            }
        if (length(error)==0)
            TRUE
        else
            error
    }

setValidity("gmmWeights", .checkGmmWeights)



.checkSysGmmWeights <- function(object)
    {
        error <- character()
        if (is.character(object@w))
            {
                if (object@w != "ident")
                    {
                        msg <- "Only 'ident' is allowed when w is a character type"
                        error <- c(error, msg)
                    } else {
                        if (object@type != "weights")
                        {
                            msg <- "type must be 'weights' when w='ident'"
                            error <- c(error, msg)                                
                        }
                    }
            } else {
                q <- sapply(object@momNames, length)
                nEqn <- length(object@eqnNames)
                if (object@type == "iid")
                {
                    if (object@sameMom)
                    {
                        if (!(class(object@w) == "qr"))
                        {
                            msg <- "w must be a 'qr' type object when type is 'iid'"
                            error <- c(error, msg)                                
                        } else {
                            if (ncol(object@w$qr) != q[1])
                            {
                                msg <- "With the same instruments, w must be the qr decomposition of only one of the equations"
                                error <- c(error, msg)                                
                            }
                        }
                    } else {
                        if (!is.matrix(object@w))
                        {
                            msg <- "w must be a matrix for iid and different instrument"
                            error <- c(error, msg)                               
                        } else {
                            if (!all(dim(object@w) == sum(q)))
                            {
                                msg <- "w must be with nrow = total number of instruments"
                                error <- c(error, msg)                                
                            }
                        }
                    }
                    if (is.null(object@Sigma))
                    {
                        msg <- "With iid vcov, Sigma must be provided"
                        error <- c(error, msg)
                    } else {
                        if (!is.matrix(object@Sigma))
                        {
                            msg <- "Sigma must be a matrix (assumed chol)"
                            error <- c(error, msg)
                        } else {
                            if (!all(dim(object@Sigma) == nEqn))
                            {
                                msg <- "Sigma must be a square matrix with nrow = num. of Equ."
                                error <- c(error, msg)
                            }
                        }
                    }
                } else if (object@type == "weights" | object@type == "vcov") {
                    if (!is.matrix(object@w))
                    {
                        msg <- "With type='weights' or 'vcov', w must be a matrix"
                        error <- c(error, msg)                                
                    } else {
                        if (!all(dim(object@w)==sum(q)))
                        {
                            msg <- "w has wrong dimension"
                            error <- c(error, msg)
                        }
                    }
                } else if (object@type == "MDS") {
                    if (!(class(object@w) == "qr"))
                    {
                        msg <- "w must be a 'qr' type object when type is 'MDS'"
                        error <- c(error, msg)                                
                    }
                    if (ncol(object@w$qr) != sum(q))
                    {
                        msg <- "The qr dimension does not correspond to the number of instruments"
                        error <- c(error, msg)                                
                    }
                } else {
                    msg <- "type must be one of 'weights', 'vcov', 'iid', or 'vcov'"
                    error <- c(error, msg)
                }
            }
        if (length(object@HAC)>0)
        {
            if (!all(names(object@HAC) %in% c("bw", "kernel","weights")))
            {
                msg <- "HAC must contain 'bw', 'kernel', and 'weights'"
                error <- c(error, msg)
            }
        }
        if (length(error)==0)
            TRUE
        else
            error
    }

setValidity("sysGmmWeights", .checkSysGmmWeights)
