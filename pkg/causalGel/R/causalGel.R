## Model builder

causalModel <- function(g, balm, data,theta0=NULL,
                      momType=c("ACE","ACT","uncondBal","fixedMom"),
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
    name_coef <- c("control",paste("treat", 1:(k-1), sep=""),
                   paste("ptreat", 1:(k-1), sep=""))
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
            }
        }    
    modData <- new("causalData", reg=tmp_model$modelF, bal=tmp_model$instF,
                   momType=momType, popMom=popMom, ACTmom=ACTmom)
    mod <- gelModel(g=causalMomFct, x=modData, gelType=gelType, rhoFct=rhoFct,
                    tet0=theta0, grad=causalDmomFct,vcov="MDS", vcovOptions=list(),
                    centeredVcov=TRUE, data=NULL)
    mod@momNames <- c(names(theta0), paste("Bal", 1:(mod@q-mod@k), sep=""))
    new("causalGel", mod)
}


