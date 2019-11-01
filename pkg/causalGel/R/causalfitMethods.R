### Hidden functions

### Helper for Covariance in the misspecified case

.psiGam <- function(object)
{
    spec <- modelDims(object@model)
    n <- spec$n
    q <- spec$q
    k <- spec$k
    ncov <- length(spec$balCov)
    Wk <- object@model@wSpec$k
    lam <- object@lambda
    theta <- coef(object)
    gt <- evalMoment(object@model, theta)
    rhoFct <- object@model@gelType
    if (is.null(rhoFct$fct))
    {
        rhoFct <- get(paste("rho", rhoFct$name, sep = ""))
    } else {
        rhoFct <- rhoFct$fct
    }
    rho1 <- rhoFct(gmat=gt, lambda=lam, derive=1, k=Wk[1]/Wk[2])
    rho2 <- rhoFct(gmat=gt, lambda=lam, derive=2, k=Wk[1]/Wk[2])
    Z <- model.matrix(object@model)
    l <- ncol(Z)
    ZT <- c(Z%*%theta[1:l])
    X <- model.matrix(object@model, "balancingCov")
    momType <- spec$momType
    balMom <-  spec$balMom
    lG1 <- sapply(1:l, function(i) -(Z[,i]*Z)%*%lam[1:l])
    q2 <- ncov*(l-1)+2*l-1
    lamM <- matrix(lam[(2*l):q2], ncol=(l-1))
    lG2 <- sapply(1:(l-1), function(i) -lam[l+i]-X%*%lamM[,i])
    lG <- cbind(lG1, lG2)
    G <- evalDMoment(object@model, theta, rho1, TRUE)
    G22 <- crossprod(rho2*gt, gt)/n
    if (momType %in% c("uncondBal", "fixedMom"))
    {
        Psi <- cbind(rho1*lG, rho1*gt)
        G11 <- crossprod(rho2*lG, lG)/n
        G12 <- t(G)/n + crossprod(rho2*lG, gt)/n
        Gamma <- rbind(cbind(G11, G12),
                       cbind(t(G12), G22))
        addPar <- 0
    } else {
        lG <- cbind(lG, matrix(-tail(lam, ncov), n, ncov, byrow=TRUE))
        G11 <- crossprod(rho2*lG, lG)/n
        G12 <- t(G)/n + crossprod(rho2*lG, gt)/n
        if (momType == "ACE")
        {
            Xi <- rep(1,n)
        } else if (momType == "ACT") {
            Xi <- Z[,spec$ACTmom+1]
        } else if (momType == "ACC") {
            Xi <- as.numeric(rowSums(Z)==1)
        } else {
            stop("Wrong balancing type")
        }
        nj <- sum(Xi)
        lam2 <- -sum(rho1)*tail(lam,ncov)/nj
        theta4 <- colSums(Xi*X)/nj
        G13 <- rbind(matrix(0, 2*l-1, ncov), -nj/n*diag(ncov))
        G23 <- matrix(0,q, ncov)
        G33 <- matrix(0, ncov, ncov)
        Psi <- cbind(rho1*lG, rho1*gt,
                     Xi*sweep(X, 2, theta4, "-"))
        Psi[,(2*l):(2*l+ncov-1)] <- Psi[,(2*l):(2*l+ncov-1)]-Xi%*%t(lam2)
        Gamma <- rbind(cbind(G11, G12, G13),
                       cbind(t(G12), G22, G23),
                       cbind(t(G13), t(G23), G33))
        addPar <- ncov
    }
    list(Psi=Psi, Gamma=Gamma, k=length(theta), q=q, addPar=addPar, n=n,
         qrGt= qr(gt/sqrt(n)))
}


####  Methods for causalGelfit class
####################################


## print

setMethod("print", "causalGelfit",
          function(x, model=TRUE, lambda=FALSE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              type <- x@type
              spec <- modelDims(x@model)
              cat("Convergence Theta: ", x@convergence, "\n")
              cat("Convergence Lambda: ", x@lconvergence, "\n")              
              cat("coefficients:\n")
              print.default(format(theta, ...), print.gap=2L, quote=FALSE)
              if (lambda)
                  {
                      cat("lambdas:\n")
                      print.default(format(x@lambda, ...), print.gap=2L, quote=FALSE)
                  }
          })

## vcov

setMethod("vcov", "causalGelfit",
          function(object, robToMiss = TRUE, withImpProb=FALSE, tol=1e-10) {
              if (!robToMiss)
                  {
                      allV <- getMethod("vcov","gelfit")(object, withImpProb, tol)
                      return(allV)
                  }
              res <- .psiGam(object)
              k <- res$k
              q <- res$q
              addPar <- res$addPar
              qrPsi <- qr(res$Psi/sqrt(res$n))
              piv <- sort.int(qrPsi$pivot, index.return=TRUE)$ix
              R <- qr.R(qrPsi)[,piv]
              T1 <- solve(res$Gamma, t(R))
              V <- T1%*%t(T1)/res$n
              allV <- list()
              allV$vcov_par <-  V[1:k, 1:k]
              allV$vcov_lambda <- V[(k+addPar+1):(k+addPar+q), (k+addPar+1):(k+addPar+q)]
              if (addPar > 0)
              {
                  allV$vcov_Allpar <-  V[1:(k+addPar), 1:(k+addPar)]
                  allV$vcov_Alllambda <- V[-(1:(k+addPar)), -(1:(k+addPar))]
              }
              allV$gtR <- qr.R(res$qrGt)
              allV              
          })
