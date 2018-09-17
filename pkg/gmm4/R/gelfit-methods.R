####  All methods for gmmfit class
#####################################
                          

## coef

setMethod("coef", "gelfit", function(object) object@theta)

## print

setMethod("print", "gelfit",
          function(x, model=TRUE, lambda=TRUE, ...) {
              theta <- coef(x)
              if (model)
                  print(x@model)
              type <- x@type
              spec <- modelDims(x@model)
              if (spec$q==spec$k && x@type != "eval")
                  type <- paste("Just-Identified ", type, sep="")
              cat("\nEstimation: ", type,"\n")
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

## show

setMethod("show","gelfit", function(object) print(object))

## residuals

setMethod("residuals", "gelfit", function(object) {
    residuals(object@model, object@theta)})

## getImpProb

setGeneric("getImpProb", function(object, ...) standardGeneric("getImpProb"))

setMethod("getImpProb", "gelfit",
          function(object) {
              rhoFct <- object@model@gelType
              if (is.null(rhoFct$fct))
                  rhoFct <- get(paste("rho", rhoFct$name, sep=""))
              else
                  rhoFct <- rhoFct$fct
              gt <- evalMoment(object@model, object@theta)
              k <- object@model@wSpec$k
              pt <- -rhoFct(gt, object@lambda, 1, k[1]/k[2])/nrow(gt)
              if (object@model@gelType$name == "EEL") {
                  eps <- -length(pt) * min(min(pt), 0)
                  pt <- (pt + eps/length(pt))/(1 + eps)
              }
              convMom <- colSums(pt * gt)
              convProb <- abs(sum(as.numeric(pt))-1)
              pt <- pt/sum(pt)
              list(pt=pt, convMom=convMom, convProb=convProb)
          })

## vcov

setMethod("vcov", "gelfit",
          function(object, withImpProb=FALSE, tol=1e-10) {
              spec <- modelDims(object@model)
              q <- spec$q
              gt <- evalMoment(object@model, object@theta)
              n <- nrow(gt)
              bw <- object@model@wSpec$bw
              k <- object@model@wSpec$k
              if (withImpProb)
                  {
                      pt <- getImpProb(object)$pt
                      G <- evalDMoment(object@model, object@theta, pt)
                      G <- G/k[1]
                      gt <- gt * sqrt(pt * bw/k[2])
                  } else {
                      G <- evalDMoment(object@model, object@theta)
                      G <- G/k[1]
                      gt <- gt * sqrt(bw/k[2]/n)
                  }
              qrGt <- qr(gt)
              piv <- sort.int(qrGt$pivot, index.return = TRUE)$ix
              R <- qr.R(qrGt)[, piv]
              X <- forwardsolve(t(R), G)
              Y <- forwardsolve(t(R), diag(q))
              res <- lm.fit(X, Y)
              u <- res$residuals
              Sigma <- chol2inv(res$qr$qr)/n
              diag(Sigma)[diag(Sigma) < 0] <- tol
              if (q == ncol(G)) {
                  SigmaLam <- matrix(0, q, q)
              } else {
                  SigmaLam <- backsolve(R, u)/n * bw^2
                  diag(SigmaLam)[diag(SigmaLam) < 0] <- tol
              }
              list(vcov_par = Sigma, vcov_lambda = SigmaLam)
          })


