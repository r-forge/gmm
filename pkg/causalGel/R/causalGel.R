## Model builder

.orth <- function (M) 
{
    if (length(M) == 0) 
        return(c())
    if (!is.numeric(M)) 
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(M)) 
        M <- matrix(c(M), nrow = length(M), ncol = 1)
    svdM <- svd(M)
    U <- svdM$u
    s <- svdM$d
    tol <- max(dim(M)) * max(s) * .Machine$double.eps
    r <- sum(s > tol)
    U[, 1:r, drop = FALSE]
}

causalModel <- function(g, balm, data,theta0=NULL,
                      momType=c("ACE","ACT","ACC", "uncondBal"),
                      popMom = NULL, ACTmom=1L, orthoBases=FALSE) 
{
    momType <- match.arg(momType)
    if (!is.null(popMom))
        momType <- "fixedMom"
    if (orthoBases)
    {
        X <- model.matrix(balm, data)[,-1]
        X <- .orth(X)
        colnames(X) <- paste("Basis", 1:ncol(X), sep="")
        balm <- as.formula(paste("~", paste(colnames(X), collapse="+")))
        data <- cbind(data, X)
    }    
    tmp_model <- momentfit:::.lModelData(g, balm, data)
    if (attr(terms(tmp_model$modelF), "intercept") != 1)
        stop("You cannot remove the intercept from g")
    if (attr(terms(tmp_model$instF), "intercept") != 1)
        stop("You cannot remove the intercept from balm")
    k <- tmp_model$k
    ncoef <- 1+2*(k-1)
    if (k>2)
        treatInd <- 1:(k-1)
    else
        treatInd <- ""
    name_coef <- c("control",paste("causalEffect", treatInd, sep=""),
                   paste("probTreatment", treatInd, sep=""))
    if (!is.null(theta0))
    {
        if (length(theta0) != ncoef)
            stop(paste("The length of theta0 must be ", ncoef, sep=""))
        if (is.null(names(theta0)))
            names(theta0) <- name_coef
    } else {
        theta0 <- numeric(ncoef)
        names(theta0) <- name_coef
        theta0[1:k] <- coef(lm(g,data))
        theta0[-(1:k)] <- colMeans(tmp_model$modelF, na.rm=TRUE)[-1]
    }
    if (momType != "uncondBal")
        {
            X <- model.matrix(terms(tmp_model$instF), tmp_model$instF)
            Z <- model.matrix(terms(tmp_model$modelF), tmp_model$modelF)
            if (momType == "ACE") {
                popMom <- colMeans(X[,-1, drop=FALSE])
            } else if (momType == "ACT") {
                popMom <- colMeans(X[Z[,1+ACTmom]==1,-1, drop=FALSE])
            } else if (momType == "ACC") {
                popMom <- colMeans(X[rowSums(Z)==1,-1, drop=FALSE])
            }
        }    
    modData <- new("causalData", reg=tmp_model$modelF, bal=tmp_model$instF,
                   momType=momType, balMom=popMom, ACTmom=ACTmom,
                   balCov=tmp_model$momNames[-1])
    mod <- momentModel(g=causalMomFct, x=modData, 
                       theta0=theta0, grad=NULL, vcov="MDS", vcovOptions=list(),
                       centeredVcov=TRUE, data=NULL)
    momNames <- lapply(treatInd, function(i)
        paste("treat", i, "_", tmp_model$momNames[-1], sep=""))
    momNames <- do.call("c", momNames)
    if (momType == "uncondBal")
        mod@momNames <- c(names(theta0), momNames)
    else
        mod@momNames <- c(names(theta0), momNames, tmp_model$momNames[-1])
    new("causalModel", mod)
}

causalGEL <- function(g, balm, data, theta0=NULL,
                   momType=c("ACE","ACT","ACC", "uncondBal","fixedMom"),
                   popMom = NULL, rhoFct=NULL,ACTmom=1L, 
                   gelType = c("EL", "ET", "EEL", "ETEL", "HD", "ETHD","REEL"),
                   initTheta = c("gmm","theta0"), getVcov=FALSE,
                   lambda0=NULL, cstLHS=NULL, cstRHS=NULL,
                   lamSlv=NULL, coefSlv= c("optim","nlminb","constrOptim"),
                   lControl=list(), tControl=list(), restrictLam=FALSE,
                   orthoBases=FALSE)
{
    Call <- try(match.call(call=sys.call(sys.parent())), silent=TRUE)
    if (class(Call)=="try-error")
        Call <- NULL              
    momType <- match.arg(momType)
    initTheta <- match.arg(initTheta)
    coefSlv <- match.arg(coefSlv)
    gelType <- match.arg(gelType)
    if (initTheta=="theta0" & is.null(theta0))
        stop("theta0 is required when initTheta='theta0'")

    model <- causalModel(g, balm, data, theta0, momType, popMom, ACTmom,
                         orthoBases)
    
    if (initTheta == "theta0")
    {
        initTheta <- "modelTheta0"
        names(theta0) = model@parNames
    } else {
        theta0 <- NULL
    }
    if (!is.null(cstLHS)) {
        if (is.numeric(cstLHS))
        {
            parn <- model@parNames
            if (is.null(cstRHS))
            {
                cstLHS <- paste(parn[cstLHS], "=0", sep="")
            } else {
                if (!is.numeric(cstRHS))
                    stop("cstRHS is either NULL or numeric")
                if (length(cstRHS)!=length(cstLHS))
                    stop("cstRHS and csrLHS must have the can length")
                cstLHS <- paste(parn[cstLHS], "=", cstRHS, sep="")
                cstRHS <- NULL
            }
        }
        model <- restModel(model, cstLHS, cstRHS)
        spec <- modelDims(model)
        if (!is.null(theta0)) 
            theta0 <- theta0[(names(theta0) %in% spec$parNames)]
        if (inherits(model, "rcausalModel") & restrictLam)
        {
            warning(paste("restrictLam=TRUE is not recommended for restricted model.",
                          "\nIt is then set to FALSE. Use gelFit if you really what it",
                          sep=""))
            restrictLam <- FALSE
        }        
    }
    if (restrictLam)
        lControl$restrictedLam <- grep("control|causalEffect",
                                       modelDims(model)$momNames)
    fit <- gelFit(model=model, gelType=gelType, rhoFct=rhoFct,
                  initTheta=initTheta, theta0=theta0,
                  lambda0=lambda0, vcov=getVcov, coefSlv=coefSlv,
                  lamSlv=lamSlv, tControl=tControl, lControl=lControl)
    if (restrictLam)
    {
        spec <- modelDims(model)
        nz <- (spec$k-1)/2
        phi <- tail(coef(fit),nz)
        if (length(phi) == 1)
        {
            phi <- matrix(c(1,phi,phi,phi),2,2)
        } else {
            phi <- rbind(c(1,phi), cbind(phi, diag(phi)))
        }
        pt <- getImpProb(fit)$pt
        Z <- model.matrix(terms(model@X@reg), model@X@reg)
        Y <- model.response(model@X@reg)
        YZ <- colSums(pt*Y*Z) 
        theta <- solve(phi, YZ)
        fit@theta[1:(1+nz)] <- theta
    }
    fit@call <- Call
    fit    
}



