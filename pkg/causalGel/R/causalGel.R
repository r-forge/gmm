## Model builder

causalModel <- function(g, balm, data,theta0=NULL,
                      momType=c("ACE","ACT","ACC", "uncondBal","fixedMom"),
                      popMom = NULL, rhoFct=NULL,ACTmom=1L, 
                      gelType = c("EL", "ET", "EEL", "ETEL", "HD", "ETHD","REEL"))
{
    momType <- match.arg(momType)
    gelType <- match.arg(gelType)
    if (!is.null(popMom))
        {
            momType <- "fixedMom"
        } else {
            if (momType == "fixedMom")
                stop("With fixed moments, popMom must be provided")
        }    
    tmp_model <- gmm4:::.lGmmData(g, balm, data)
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
            stop(paste("The leangth of theta0 must be ", ncoef, sep=""))
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
    mod <- gelModel(g=causalMomFct, x=modData, gelType=gelType, rhoFct=rhoFct,
                    theta0=theta0, grad=NULL,vcov="MDS", vcovOptions=list(),
                    centeredVcov=TRUE, data=NULL)
    momNames <- lapply(treatInd, function(i)
        paste("treat", i, "_", tmp_model$momNames[-1], sep=""))
    momNames <- do.call("c", momNames)
    if (momType == "uncondBal")
        mod@momNames <- c(names(theta0), momNames)
    else
        mod@momNames <- c(names(theta0), momNames, tmp_model$momNames[-1])
    new("causalGel", mod)
}

causalGEL <- function(g, balm, data, theta0=NULL,
                   momType=c("ACE","ACT","ACC", "uncondBal","fixedMom"),
                   popMom = NULL, rhoFct=NULL,ACTmom=1L, 
                   gelType = c("EL", "ET", "EEL", "ETEL", "HD", "ETHD","REEL"),
                   initTheta = c("gmm","theta0"), getVcov=FALSE,
                   lambda0=NULL, cstLHS=NULL, cstRHS=NULL,
                   lamSlv=NULL, coefSlv= c("optim","nlminb","constrOptim"),
                   lControl=list(), tControl=list())
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

    model <- causalModel(g, balm, data, theta0, momType, popMom, rhoFct, ACTmom,
                       gelType)
    
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
    }    
    fit <- modelFit(model=model, initTheta=initTheta, theta0=theta0,
                    lambda0=lambda0, vcov=getVcov, coefSlv=coefSlv,
                    lamSlv=lamSlv, tControl=tControl, lControl=lControl)
    fit@call <- Call
    fit    
}



