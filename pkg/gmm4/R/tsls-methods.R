## Method for tsls class objects


## bread

setMethod("bread", "tsls",
          function(x, ...) {
              sig <- sum(residuals(x)^2)
              X <- model.matrix(x@model)
              Xhat <- qr(qr.fitted(x@wObj@w, X))
              v <- matrix(ncol=ncol(X), nrow=ncol(X))
              v[Xhat$pivot, Xhat$pivot] <- chol2inv(qr.R(Xhat))
              sig*v
          })

## meat

setMethod("meatGmm", "tsls",
          function(object, robust=FALSE) {
              sig <- sum(residuals(object)^2)
              X <- model.matrix(object@model)
              if (object@model@vcov == "iid" || !robust)
                  {
                      Xhat <- qr.fitted(object@wObj@w, X)
                      meat <- crossprod(Xhat)/sig
                  } else {
                      T1 <- qr.coef(object@wObj@w, X)
                      v <- momentVcov(object@model, coef(object))
                      meat <- crossprod(T1, v)%*%T1/sig
                  }
              meat
          })


## vcov

setMethod("vcov", signature("tsls"),
          function(object, sandwich=TRUE, df.adj=FALSE)
              {
                  spec <- modelDims(object@model)
                  if (sandwich)
                      {
                          b <- bread(object)
                          meat <- meatGmm(object, TRUE)
                          vcov <- b %*% meat %*% b/spec$n
                      } else {
                          vcov <- bread(object)/spec$n
                      }
                  dimnames(vcov) <- list(spec$parNames, spec$parNames)
                  if (df.adj) 
                      vcov <- vcov * spec$n/(spec$n - spec$k)
                  attr(vcov, "type") <- list(sandwich = sandwich, df.adj = df.adj, 
                                             breadOnly = !sandwich)
                  vcov                  
              })

