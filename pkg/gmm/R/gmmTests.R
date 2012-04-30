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


# This function compute what is needed for the K statistics of KleiBergen (2005)
#####################################################################################

.BigCov <- function(obj,theta0)
	{
	insertRC <- function(A,w,v)
		{
		NewA <- matrix(ncol=ncol(A)+length(w),nrow=nrow(A)+length(w))
		NewA[-w,-w] <- A
		NewA[w,] <- v
		NewA[,w] <- v
		NewA
		}
	dg <- function(obj)
		{
		dat <- obj$dat
		if (!is.null(attr(dat,"eqConst")))
			x <- attr(dat,"eqConst")$Xunc
		else
			x <- dat$x[,(dat$ny+1):(dat$ny+dat$k)]
		k <- ncol(x)
		h <- dat$x[,(dat$ny+dat$k+1):ncol(dat$x)]
		qt <- array(dim=c(dim(obj$gt),k))
		for (i in 1:k)
			qt[,,i] <- -x[,i]*h
		qt}

	if (attr(obj$dat,"ModelType") == "nonlinear")
		{
		Myenv <- new.env()
		assign("obj", obj, envir=Myenv)
		assign("theta", theta0, envir=Myenv)
		gFunct <- if (!is.null(attr(obj$dat,"eqConst")))
				attr(obj$dat,"eqConst")$unConstg
			  else
				obj$g
		assign("g",gFunct,envir=Myenv)
		res <- numericDeriv(quote(g(theta,obj$dat)),"theta",Myenv)
		qT <- attr(res,"gradient")
		} else {
		qT <- dg(obj)}

	qTmat <- apply(qT,3,colSums)
	qT <- matrix(qT,nrow=dim(qT)[1])
	gt <- obj$g(theta0,obj$dat)
	fT <- colSums(gt)
	n <- nrow(gt)
	q <- ncol(gt)
	All <- cbind(gt,qT)
	All <- sweep(All,2,colMeans(All),FUN="-")
	f <- function(x)
		all(abs(x)<1e-7)
	w <- which(apply(All,2,f))
	if (length(w) != 0)
		All <- All[,-w]
	if (dim(All)[2] >= dim(All)[1])
		stop("Too many moment conditions. Cannot estimate V")
	if (obj$WSpec$vcov == "iid") 
		{
		V <- crossprod(All)/nrow(All) 
	} else {
		class(All) <- "gmmFct"
		argSand <- obj$WSpec$sandwich
		argSand$x <- All
		argSand$sandwich <- FALSE
		V <- do.call(kernHAC,argSand)
		}
	if (length(w) != 0)
		V <- insertRC(V,w,0)
	Vff <- V[1:q,1:q]
	Vthetaf <- V[(q+1):nrow(V),1:q]
	list(Vff=Vff,Vthetaf=Vthetaf,qT=qTmat,fT=fT,n=n,q=q)
	}		

KTest <- function(obj, theta0=NULL, alphaK = 0.04, alphaJ = 0.01)
	{
	if (class(obj) != "gmm")
		stop("KTest is only for gmm type objects")

	if (!is.null(attr(obj$dat,"eqConst")))
		{
		if (!is.null(theta0))
			warning("setting a value for theta0 has no effect when the gmm is already constrained")
		resTet <- attr(obj$dat,"eqConst")$eqConst
		tet <- obj$coefficients
		theta0 <- vector(length=length(tet)+nrow(resTet))
		theta0[resTet[,1]] <- resTet[,2]
		theta0[-resTet[,1]] <- tet
		testName <- paste(rownames(resTet), " = ", resTet[,2], collapse="\n")
		if (is.list(obj$dat))		
			{
			x <- model.matrix(obj$dat$mt,obj$dat$mf)
			y <- model.response(obj$dat$mf)
			obj$dat$x <- cbind(y,x,obj$dat$x[,(obj$dat$ny+obj$dat$k+1):ncol(obj$dat$x)])
			obj$dat$k <- ncol(x)
			} else {
			obj$g <- attr(obj$dat,"eqConst")$unConstg
			}	
		dfK <- nrow(resTet)
		which <- resTet[,1]
	} else {
		if (is.null(theta0))
			stop("You must either estimate a restricted model first or set theta0 under H0")		
		if (length(theta0) != length(obj$coef))
			stop("theta0 is only for tests on the whole vector theta when obj is an unrestricted GMM")		
		dfK <- length(theta0)
		testName <- paste(names(obj$coef), " = ", theta0, collapse="\n")
		which <- 1:length(theta0)
		}
	V <- .BigCov(obj, theta0)
	Vff <- V$Vff
	Vtf <- V$Vthetaf
	qT <- V$qT
	fT <- V$fT
	dfJ <- V$q-length(theta0)
	# the following is vec(D)
	D <- c(qT)-Vtf%*%solve(Vff,fT)
	D <- matrix(D,ncol=length(theta0))
	meat <- t(D)%*%solve(Vff,D)
	bread <- t(fT)%*%solve(Vff,D)
	K <- bread%*%solve(meat,t(bread))/V$n
	pv <- 1-pchisq(K,dfK)
	J <- t(fT)%*%solve(Vff,fT)/V$n-K
	pvJ <- 1-pchisq(J,dfJ)
	type <- c("K statistics","J statistics")
	test <- c(K,J)
	test <- cbind(test,c(pv,pvJ),c(dfK,dfJ))
	dist <- paste("Chi_sq with ", c(dfK,dfJ), " degrees of freedom", sep="")
	if(dfJ>0)
		ans <- list(test=test,dist=dist,type=type,testName=testName)
	else
		ans <- list(test=matrix(test[1,],nrow=1),dist=dist[1],type=type[1],testName=testName)
	if (pvJ<alphaJ)	{
		message <- "reject"
	} else {
		if (pv < alphaK)
			message <- "reject"
		else
			message <- "do not reject"
	}
	ans$KJ <- (message == "do not reject")
	ans$Matrix <- list(D=D,bread=bread,meat=meat,qT=qT,fT=fT)	
	ans$message <- paste("KJ-test result: We ", message, " H0 (alphaJ = ", alphaJ, ", alphaK = ", alphaK, ")", sep="")	
	class(ans) <- "gmmTests"
	ans
	}	

print.gmmTests <- function(x, digits = 5, ...)
	{
	lab <- paste("%0.",digits,"f",sep="")
	cat("\nTest robust to weak identification\n")
	cat("**********************************\n\n")
	cat("The Null Hypothesis\n")
	cat(x$testName,"\n\n")
	for (i in 1:length(x$type))
		{
		cat(x$type[i],"\n")
		cat("Test: ", sprintf(lab,x$test[i,1]), "(",x$dist[i],")\n")
		cat("P-value: ", sprintf(lab,x$test[i,2]),"\n\n")
		}
	cat(x$message,"\n")
	}

gmmWithConst <- function(obj, which, value)
	{
	argCall <- obj$allArg
	argCall$call = NULL
	if (!is.null(attr(obj$w0,"Spec")))
		if (is.function(argCall$bw))
			argCall$bw <- attr(obj$w0,"Spec")$bw

	if (length(which)>=length(obj$coefficients))
		stop("Too many constraints")
        if (is.character(which))
		{
		if (any(!(which %in% names(obj$coefficients))))
		   stop("Wrong coefficient names in eqConst")
		if (attr(obj$dat,"ModelType") == "linear")
			which <- match(which,names(obj$coefficients))
		}  
	if (!is.null(argCall$t0))
		{
		argCall$t0 <- obj$coefficients
		argCall$t0[which] <- value
		}
		
	if (attr(obj$dat,"ModelType") == "nonlinear")
		{
		eqConst <- which
	} else {
		eqConst <- cbind(which,value)		
		}
	argCall$eqConst <- eqConst
	res <- do.call(gmm,argCall)
	res$call <- match.call()
	return(res)
	}

KConfid <- function(obj, which, type = c("K", "KJ"), alpha = 0.05, alphaJ = 0.01, n = 4)
	{
	type <- match.arg(type)
	if ( (obj$df == 0) & (type == "KJ"))
		stop("Only K type is available for just identified models")
		 
	if ( (alphaJ >= alpha) | (alphaJ <= 0) )
		stop("We must have 0 < alphaJ < alpha")
	if ( (alpha <= 0) | (alpha >=1) )
		stop("We must have 0 < alpha < 1")
	alphaK <- alpha-alphaJ
	
	if (!is.null(attr(obj$dat,"eqConst")))
		stop("Confidence intervals are constructed from unrestricted models")
	if (length(which) > 2)
		stop("This function computes confidence intervals for 1 or 2 coefficients only")	
	if (length(which)>length(obj$coefficients))
		stop("length(which) must not exceed length(coefficients)")
	if (is.character(which))
		{
		if (any(!(which %in% names(obj$coef))))
			stop("names in which do not match names in coefficients")
		else
			which <- sort(which(which %in% names(obj$coef)))
		} else {
		if (any(which>length(obj$coefficients)) | any(which<=0))
			stop("indices in which are not valid")
		}
	getUniInt <- function(tetx, obj, which, value=NULL, type , alpha, alphaJ, alphaK)
		{
		# value is the value of the other coefficient when we have 2 restrictions
		# tetx is coef(obj)[which[1]] when we have two restrictions

		if (is.null(value))
			value <- tetx
		else
			value <- c(tetx,value)
		res <- gmmWithConst(obj, which = which, value = value)
		test <- KTest(res)
		if (type == "K")
		   {
		      ifelse(alpha > test$test[1,2], 1, -1)
		   } else {
		      if (test$test[2,2]<alphaJ) 
		         {
			    return(1)
			 } else {
 			    ifelse(test$test[1,2]<alphaK, 1, -1) 
			 }
		   }	
	}
	gForMulti <- function(theta, obj, which, type , alpha, alphaJ, alphaK)
		getUniInt(theta[1], obj, which, theta[2], type , alpha, alphaJ, alphaK)

	if (length(which) == 1)
		{
		step <- 5*sqrt(diag(vcov(obj)))[which]
		X0 <- obj$coef[which]
		getBoth <- function(d)
			{
			res <- try(uniroot(getUniInt,c(X0,X0+d*step), obj = obj, which = which, type = type, 
					alpha = alpha, alphaJ = alphaJ, alphaK = alphaK, tol=1e-4), silent=TRUE)
			if (class(res) == "try-error")
				return(NA)
			else
				return(res$root)
			}
		res <- mclapply(c(-1,1), getBoth)
		c(res[[1]],res[[2]])
	} else {
		step <- 5*sqrt(diag(vcov(obj)))[which]
		sol <- .getCircle(obj$coef[which][1],obj$coef[which][2],gForMulti,n , step, obj=obj, 
					which=which, type=type, alpha=alpha, alphaJ=alphaJ, alphaK=alphaK)	
		colnames(sol) <- names(obj$coef)[which]
		sol
		}
	}

.getCircle <- function(x0,y0,g,n,b, trace=FALSE,  ...)
	{
	tol=1e-4
	if (any(b<=0))
		stop("b must be strictly positive")

	
	lambda <- seq(1,0,length=n)[-c(1,n)]
	f <- function(x, y0, x0 = NULL, xi = NULL, yi = NULL)
		{
		if (is.null(xi))
			y <- y0
		else
			y <- y0	+ (yi-y0)/(xi-x0)*(x-x0)
		g(c(x,y), ...)	
		}
	f2 <- function(y, y0)
		g(c(y0,y), ...)
	f3 <- function(x)
		ifelse(class(x)=="try-error",NA,x)

	# get the four points of the cross
	selectf <- list(f,f,f2,f2)
	selectd <- c(1,-1,1,-1)
	xy12 <- mclapply(1:4, function(i) try(uniroot(selectf[[i]],c(x0,x0+selectd[i]*b[1]),y0=y0,tol=tol)$root,silent=TRUE))

	x1 <- xy12[[1]]
	x2 <- xy12[[2]]
	y1 <- xy12[[3]]
	y2 <- xy12[[4]]

	getAll <- function(lambda, dir=1, x)
		{
		yi <- y1*lambda + y0*(1-lambda)
		xi <- x0*lambda + x*(1-lambda)
		xt <- try(uniroot(f,c(x0,x0+dir*b[1]),y0=y0, x0=x0, xi=xi, yi=yi,tol=tol)$root,silent=TRUE)
		xt <- f3(xt)
		yt <- y0 + (yi-y0)/(xi-x0)*(xt-x0)
		return(c(xt,yt))
		}

	res1 <- mclapply(lambda,getAll,dir=1,x=x1)
	res2 <- mclapply(lambda,getAll,dir=-1,x=x2)
	res1 <- t(simplify2array(res1))
	res2 <- t(simplify2array(res2))
	solU <- rbind(c(x2,y0),res2[length(lambda):1,],c(f3(x0), f3(y1)),res1,c(x1,y0))

	res1 <- mclapply(lambda,getAll,dir=-1,x=x1)
	res2 <- mclapply(lambda,getAll,dir=1,x=x2)
	res1 <- t(simplify2array(res1))
	res2 <- t(simplify2array(res2))
	solD <- rbind(res2[length(lambda):1,],c(f3(x0), f3(y2)),res1)

	sol <- rbind(solU,solD)
	return(sol)
	}

