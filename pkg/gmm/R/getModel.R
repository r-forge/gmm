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

getModel.baseGmm <- function(object, ...)
  {
  if(is(object$g, "formula"))
    {
    object$gradvf <- FALSE
    dat <- getDat(object$g, object$x)
    clname <- paste(class(object), ".", object$type, ".formula", sep = "")
    object$gform<-object$g
    g <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] - x[,(ny+1):(ny+k)] %*% t(tet)
      gt <- e * x[, ny+k+1]
      if(nh > 1)
	for (i in 2:nh)	  gt <- cbind(gt, e*x[, (ny+k+i)])
      return(gt)
      }
    gradv <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k, g = NULL)
      {
      a <- g
      tet <- NULL
      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)]) %*% x[,(ny+1):(ny+k)]) %x% diag(rep(1,ny))/nrow(x)
      return(dgb)
      }
    object$g <- g
    }
  else
    {
    clname<-paste(class(object), "." ,object$type, sep = "")
    if (!is.function(object$gradv))
      { 
      gradv <- .Gf
      object$gradvf <- FALSE
      }
    else
      {
      gradv <- object$gradv
      object$gradvf <- TRUE
      }
    }
	
  iid <- function(thet, x, g)
    {
    gt <- g(thet,x)
    n <- ifelse(is.null(nrow(x)), length(x), nrow(x))
    v <- crossprod(gt,gt)/n
    return(v)
    }
 
  object$iid<-iid
  object$TypeGmm <- class(object)
  object$gradv <- gradv	
 
  class(object)  <- clname
  return(object)
  }

getModel.baseGel <- function(object, ...)
  {
  P <- object
  if (P$type == "ETEL")
    {
    P$typel <- "ET"
    P$typet <- "EL"	
    }
  else
    {
    P$typel <- P$type
    P$typet <- P$type
    }
  if(P$optfct == "optim")
    P$k <- length(P$tet0)
  else
    P$k <- 1
  
  if (is(P$g, "formula"))
    {
    clname <- paste(class(P), ".modFormula", sep = "")
    dat <- getDat(P$g, P$x)
    x <- dat$x

    g <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] -  x[, (ny+1):(ny+k)]%*%t(tet)
      gt <- e*x[, ny+k+1]
      if (nh > 1)
        {	
        for (i in 2:nh)
          {
          gt <- cbind(gt, e*x[,(ny+k+i)])
          }
        }
      return(gt)
      }

    gradv <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1, ny))/nrow(x)
      return(dgb)
      }
    P$dat <- dat
    P$gform <- P$g
    P$g <- g
    P$gradv <- gradv
    }	
  else
    {
    clname <- paste(class(P), ".mod", sep = "")
    P$gform <- NULL
    }

  if (P$smooth)
    {
    P$g1 <- P$g
    rgmm <- gmm(P$g, x, P$tet0, wmatrix = "ident")

    if (is.function(P$weights))
      P$w <- P$weights(P$g(rgmm$coefficients, x), kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
               ar.method = P$ar.method, approx = P$approx, tol = P$tol_weights)
    else
      P$w <- P$weights

    P$sg <- function(thet, x, g1 = P$g1, w = P$w)
      {
      gf <- g1(thet, x)
      gt <- smoothG(gf, weights = w)$smoothx 
      return(gt)
      }
    }	
  class(P) <- clname
  return(P)
  }

