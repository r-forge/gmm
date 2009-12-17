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

momentEstim <- function(object, ...)
  {
  UseMethod("momentEstim")
  }

momentEstim.baseGmm.twoStep <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    {
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    }
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	
  if (q == k2 | P$wmatrix == "ident")
    z = list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)	
  else
    {
    if (P$vcov == "iid")
      w <- P$iid(res$par, P$x, P$g)

    if (P$vcov == "HAC")
      w <- HAC(P$g(res$par, P$x), kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol)

    if (P$optfct == "optim")
      res2 <- optim(res$par, .obj1, x = P$x, w = w, gf = P$g, ...)

    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(res$par, .obj1, x = P$x, w = w, gf = P$g, ...)
      res2$value <- res2$objective
      }

    if (P$optfct == "optimize")
      {
      res2 <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }	

     z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df)	
    }

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

momentEstim.baseGmm.twoStep.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- getDat(P$gform, P$x)
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, P$g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    w=diag(rep(1, q))
    res1 <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, P$g)
    if (P$vcov == "iid")
      w <- P$iid(res1$par, x, g)
   if (P$vcov == "HAC")
      w <- HAC(g(res1$par, x), kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol)
     res2 <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df)	
    }
  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.iterative.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- getDat(P$gform, P$x)
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    w=diag(rep(1, q))
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "iid")
        w <- P$iid(tet, x, g)
      if (P$vcov == "HAC")
        w <- HAC(g(tet, x), kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, ar.method = P$ar.method, approx = P$approx, tol = P$tol)
      res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
      ch <- crossprod(abs(tet- res$par)/tet)^.5
      if (j>P$itermax)
        {
        cat("No convergence after ", P$itermax, " iterations")
        ch <- P$crit
        }
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df)	
   }
  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.iterative <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
    }	
  else
    {
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "iid")
        w <- P$iid(tet, P$x, P$g)
      if (P$vcov == "HAC")
        w <- HAC(P$g(tet, P$x), kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol)

      if (P$optfct == "optim")
        res <- optim(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
      if (P$optfct == "nlminb")
        {
        res <- nlminb(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
        res$value <- res$objective
        }
      if (P$optfct == "optimize")
        {
        res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
        res$par <- res$minimum
        res$value <- res$objective
        }	
        ch <- crossprod(abs(tet-res$par)/tet)^.5	
        if (j>P$itermax)
          {
          cat("No convergence after ", P$itermax, " iterations")
          ch <- P$crit
          }
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value,k=k, k2=k2, n=n, q=q, df=df)	
    }

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

momentEstim.baseGmm.cue.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- getDat(P$gform, P$x)
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    if (is.null(P$t0))
      P$t0 <- .tetlin(x,diag(q), dat$ny, dat$nh, dat$k, P$gradv, g)$par
    if (P$optfct == "optim")
      res2 <- optim(P$t0,.objCue, x = x, P = P, ...)
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0,.objCue, x = x, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = x, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }

  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.cue <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
    }	
  else
    {
    if (P$optfct == "optim")
      res2 <- optim(P$t0, .objCue, x = x, P = P, ...)
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0, .objCue, x = x, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = x, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients=res2$par,objective=res2$value, k=k, k2=k2, n=n, q=q, df=df)	
    }

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }



