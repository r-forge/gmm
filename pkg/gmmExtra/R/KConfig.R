KConfid <- function(obj, which, type = c("K", "KJ"), alpha = 0.05, alphaJ = 0.01, n = 4,
                    mc.cores=1)
	{
        theApply <- mclapply
        if (Sys.info()[["sysname"]] == "Windows")
            mc.cores <- 1

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
		      test$test[1,2]-alpha
		   } else {
		      if (test$test[2,2]<alphaJ) 
		         {
			    test$test[2,2] - alphaJ
			 } else {
 			    test$test[1,2] - alphaK
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
			res <- try(uniroot(getUniInt,c(X0,X0+d*step),
                                           obj = obj, which = which, type = type, 
                                           alpha = alpha, alphaJ = alphaJ, alphaK = alphaK,
                                           tol=1e-4), silent=TRUE)
			if (class(res) == "try-error")
                            return(NA)
			else
                            return(res$root)
                    }
		res <- theApply(c(-1,1), getBoth, mc.cores=mc.cores)
		c(res[[1]],res[[2]])
	} else {
		step <- 5*sqrt(diag(vcov(obj)))[which]
		sol <- .getCircle(obj$coef[which][1],obj$coef[which][2],
                                  gForMulti,n , step, obj=obj, 
                                  which=which, type=type, alpha=alpha,
                                  alphaJ=alphaJ, alphaK=alphaK)	
		colnames(sol) <- names(obj$coef)[which]
		sol
		}
	}

.getCircle <- function(x0,y0,g,n,b, trace=FALSE,  ...)
	{
	theApply <- mclapply
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
	xy12 <- theApply(1:2, function(i)
            try(uniroot(selectf[[i]],c(x0,x0+selectd[i]*b[1]),y0=y0,tol=tol)$root,
                silent=TRUE))
	xy34 <- theApply(3:4, function(i)
            try(uniroot(selectf[[i]],c(y0,y0+selectd[i]*b[2]),y0=x0,tol=tol)$root,
                silent=TRUE))

	x1 <- xy12[[1]]
	x2 <- xy12[[2]]
	y1 <- xy34[[1]]
	y2 <- xy34[[2]]

	getAll <- function(lambda, dir=1, x)
		{
		yi <- y1*lambda + y0*(1-lambda)
		xi <- x0*lambda + x*(1-lambda)
		xt <- try(uniroot(f,c(x0,x0+dir*b[1]),y0=y0, x0=x0,
                                  xi=xi, yi=yi,tol=tol)$root,silent=TRUE)
		xt <- f3(xt)
		yt <- y0 + (yi-y0)/(xi-x0)*(xt-x0)
		return(c(xt,yt))
		}

	res1 <- theApply(lambda,getAll,dir=1,x=x1)
	res2 <- theApply(lambda,getAll,dir=-1,x=x2)
	res1 <- t(simplify2array(res1))
	res2 <- t(simplify2array(res2))
	solU <- rbind(c(x2,y0),res2[length(lambda):1,],c(f3(x0), f3(y1)),res1,c(x1,y0))

	res1 <- theApply(lambda,getAll,dir=-1,x=x1)
	res2 <- theApply(lambda,getAll,dir=1,x=x2)
	res1 <- t(simplify2array(res1))
	res2 <- t(simplify2array(res2))
	solD <- rbind(res2[length(lambda):1,],c(f3(x0), f3(y2)),res1)

	sol <- rbind(solU,solD)
	return(sol)
	}

