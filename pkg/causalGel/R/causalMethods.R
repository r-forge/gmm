
## Moment functions

setGeneric("causalMomFct", function(theta, object, ...) standardGeneric("causalMomFct"))

setMethod("causalMomFct", signature("numeric", "causalData"),
          function(theta, object) {
              Z <- model.matrix(terms(object@reg), object@reg)
              Y <- model.response(object@reg)
              X <- model.matrix(terms(object@bal), object@bal)
              k <- ncol(Z)
              e <- Y-c(Z%*%theta[1:k])
              m1 <- e*Z
              e <- t(t(Z[,-1,drop=FALSE])-theta[-(1:k)])
              m2 <- sapply(1:ncol(X), function(i) e*X[,i])
              if (object@momType == "uncondBal")
                  return(cbind(m1,m2))
              m3 <- sweep(X[,-1,drop=FALSE], 2, object@popMom, "-")
              cbind(m1,m2,m3)
          })

## DMoment functions

setGeneric("causalDmomFct", function(theta, object, ...) standardGeneric("causalDmomFct"))

setMethod("causalDmomFct", signature("numeric", "causalData"),
          function(theta, object, pt=NULL) {
              Z <- model.matrix(terms(object@reg), object@reg)
              X <- model.matrix(terms(object@bal), object@bal)
              k <- ncol(Z)
              n <- nrow(Z)
              ntet <- length(theta)
              if (is.null(pt))
                  pt <- rep(1/n, n)
              ZT <- c(Z%*%theta[1:k])
              q <- 2*k + (k-1)*(ncol(X)-1) - 1
              G <- matrix(0, q, ntet)
              G11 <- lapply(1:k, function(i) -colSums(pt*Z[,i]*Z))
              G[1:k, 1:k] <- do.call(rbind, G11)
              G[(k+1):ntet, (k+1):ntet] <- -sum(pt)*diag(k-1)
              uK <- colSums(pt*X[,-1,drop=FALSE])
              G[(2*k):q, (k+1):ntet] <- -kronecker(diag(k-1), uK)
              if (object@momType != "uncondBal" |  object@momType=="fixedMon")
                  G <- rbind(G, matrix(0, ncol(X)-1, ntet))
              G
          })
