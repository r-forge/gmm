.timeSimul <- function(iter,titer,num_proc=NULL,comment=NULL,start=1, every=1, env=.GlobalEnv)
	{

	# Function to keep track of the time remaining 
	# in repeating processes like simulations
	# usage: put it at the beginning of the processe inside the loop.

	# iter:  the iteration being done
	# titer: the total number of iterations
	# comment (optional) : Some comments to be add in string format. Useful if running several processes and want to know which one is running
	# num_proc( optional) : If simulation is made on more than one processes, num_proc[1] is the process being done and
	#            proc_num[2] is the total number of processes
 

	if(iter==start)
		{
		cat("starting the simulation process", "\n")
		cat("###############################",  "\n")
		cat("\n")
		assign("T_INITIAL_GLOBAL", proc.time()[3], envir = env)
		}
	else
		{
		T_INITIAL_GLOBAL <- get("T_INITIAL_GLOBAL", envir = env)
		t1 <- proc.time()[3]
		tm <- (t1 - T_INITIAL_GLOBAL)/(iter-start)	
		trest <- tm * (titer-iter+1)
		tresth <- floor(trest/3600)
		trestm <- floor((trest-tresth*3600)/60)
		trests <- round((trest-tresth*3600-trestm*60),0)
		if ((iter-1)/every==(iter-1)%/%every)
		{
		if(!is.null(comment))
				cat("## ", comment, " ##", "\n")
		if(is.null(num_proc))
			cat("iteration ",(iter-1)," from a total of ",titer,"; average time = ",round(tm,2)," sec.",  "\n")
		else
			{
			cat("Process ", num_proc[1], "from a total of ", num_proc[2], "\n")
			cat("iteration ",(iter-1)," from a total of ",titer,"; average time = ",round(tm,2)," sec.",  "\n")
			}
		cat("time left ",tresth," Hrs",trestm," Min",trests," sec.", "\n", "\n")
		}
		}
	}

.getBlock <- function(x, l)
	{
	aX <- attributes(x)
	x <- as.matrix(x)
	n <- nrow(x)
	b <- floor(n/l)
	reDist <- n-b*l
	if (reDist != 0)
		{
		N <- sample(1:(n-l), reDist, replace=TRUE)
		i <- lapply(N, function(j) j:(j+l))
		i <- simplify2array(i)
		} else {
		i <- vector()
		}		
	N <- sample(1:(n-l+1), (b-reDist), replace=TRUE)
	i2 <- lapply(N, function(j) j:(j+l-1))
	i2 <- simplify2array(i2)
	xF <- x[c(i,i2),]
	attributes(xF) <- aX
	attr(xF,"Blocks") <- c(i,i2)
	attr(xF,"l") <- l
	xF
	}

.getS <- function(gt, Blocks, l)
	{
	V <- matrix(0,ncol(gt),ncol(gt))
	T <- nrow(gt)-l+1
	j <- 1
	for (i in 1:length(Blocks))
		{
		l <- length(Blocks[[i]])
		gt0 <- gt[j:(j+l-1),,drop=FALSE]
		j <- j+l
		for (s in 1:l)
			for(t in 1:l)
				V <- V+outer(gt0[t,],gt0[s,])
		}
	V/T
	}

.SimByBlock <- function(f, niter, titer, trace=FALSE, ...)
	{

	# f is a function of i=1,....,titer
	# it returns the values as a list
	# if f(i) fails, the bad component of the list is set to TRUE
	# niter (number of simultaneous evaluations) is set to titer is niter>titer
	# X cannot be used in ...
	T_INITIAL_GLOBAL <- NULL

	res <- list()
	if (niter > titer)
		niter <- titer

	nfin <- 0
	i <- 1
	f2 <- function(i, ...)
		{
		ans <- list(ans=try(f(i, ...),silent=TRUE))
		if(class(ans$ans) == "try-error")
			return(list(bad=TRUE))
		else
			{	
			ans$bad <- FALSE
			return(ans)
			}
		}

	while (nfin != titer)
		{
		if (trace)
			.timeSimul(i,titer,comment=paste("(",niter," estimations per iteration)",sep=""), env=parent.frame())

		L1 <- list(...)
		L2 <- c(L1,list(X=c((nfin+1):(nfin+niter)), FUN=f2))
		res2 <- do.call(mclapply,L2)
		for (j in 1:length(res2))
			res[[(nfin+j)]] <- res2[[j]]

		nfin <- nfin + length(res2)

		if ( (titer-nfin) < niter)
			niter <- (titer-nfin)
		i <- nfin+1
		}

	return(res)
	}

bootGmm <- function(obj, N, seed = NULL, niter=8, trace=TRUE)
	{
	if(!is.null(seed))
		set.seed(seed)
	argCall <- obj$allArg
	argCall$call <- NULL
	argCall$t0 <- obj$coef
	bw <- attr(obj$w0, "Spec")$bw
	argCall$bw <- bw
	l <- attr(obj$w0, "Spec")$bw
	if (!is.null(l))
		{
		l <- floor(l)
		} else {
		l <- 1
		}
	l <- max(l,1)

	gtBlock <- na.omit(filter(obj$gt,rep(1/l,l)))
	mustar <- colMeans(gtBlock)
	
	if(attr(obj$dat,"ModelType") == "linear")
		{
		argCall$x <- obj$dat
		attr(argCall$x, "oldg") <- obj$g
		dat <- obj$dat$x
		} else {
		dat <- argCall$x
		attr(argCall$x, "ModelType") <- "nonlinear"
		attr(argCall$x, "oldg") <- argCall$g
		}
	attr(argCall$x, "mu") <- mustar
	
	Newdat <- list()
	for (i in 1:N)
		Newdat[[i]] <- .getBlock(dat, l)

	getAll <- function(i, argCall, Newdat)
		{
		g2 <- function(theta, x)
			{	
			mustar <- attr(x,"mu")
			gt <- attr(x,"oldg")(theta, x)
			sweep(gt, 2, mustar)
			}
		if (attr(argCall$x, "ModelType") == "nonlinear") {
			argCall$x <- Newdat[[i]]
		} else {
			argCall$x$x <- Newdat[[i]]}

		argCall$g <- g2
		argCall$wmatrix <- "ident"
			
		res0 <- do.call(gmm, argCall)
		argCall$wmatrix <- "optimal"
		argCall$vcov <- "TrueFixed"
		S <- .getS(res0$gt, attr(Newdat[[i]],"Blocks"), attr(Newdat[[i]],"l"))
		argCall$weightsMatrix <- try(solve(S),silent=TRUE)
		if (class(argCall$weightsMatrix) == "try-error")
			{
			stop("Singular S")
		} else {
			res1 <- do.call(gmm, argCall)
			res1$S <- S	
			res1
		}
		}
	res <- .SimByBlock(getAll, niter, N, trace=trace, argCall=argCall, Newdat=Newdat)
	chk <- sapply(1:N, function(i) !res[[i]]$bad)
	res <- res[chk]
	attr(res,"gmm") <- obj
	class(res) <- "bootGmm"
	return(res)
	}

summary.bootGmm <- function(object, ...)
	{
	n <- length(object)
	if (n == 0)
		stop("The bootGmm object is empty")
	coef <- sapply(1:n, function(i) object[[i]]$ans$coefficients)
	coef <- t(coef)
	conv <- sapply(1:n, function(i) object[[i]]$ans$algoInfo$convergence)
	cat("Summary Statistics of Boostrap Estimates (N=",n,")\n")
	cat("#######################################################\n")
	print(ans <- summary(coef))	
	if (!is.null(conv[1]))
		cat("\nThe number of estimation that did not converge is ", sum(conv),"\n")
	cat("The original estimates are\n")
	print(round(attr(object,"gmm")$coefficients,4))
	}


plot.bootGmm <- function(x, which = 1, type = c("points", "density"), ...)
	{
	if (length(x) == 0)
		stop("The bootGmm object is empty")
	type <- match.arg(type)
	check <- c("main", "sub", "xlab", "ylab") %in% names(list(...))
	coef <- sapply(1:length(x), function(i) x[[i]]$ans$coefficients[which])
	ncoef <- names(x[[1]]$ans$coefficients[which])
	message <- paste("The point estimate is ", ncoef, " = ", round(attr(x,"gmm")$coefficients[which],4), sep="")
	if (type == "points")
		{
		if (!check[1])
			main <- "Gmm Bootstrap Estimates"
		if (!check[2])
			sub <- message
		if (!check[3])
			xlab <- "Bootstrap index"
		if (!check[4])
			ylab <- ncoef
		plot(1:length(coef), coef, main = main, sub = sub, xlab=xlab, ylab=ylab, ...)
		}
	if (type == "density")
		{
		if (!check[1])
			main <- "Gmm Bootstrap Estimates (kernel density)"
		if (!check[2])
			sub <- message
		plot(density(coef), main = main, sub = sub, ...)
		}
	}

bootJ <- function(obj)
	{
	n <- length(obj)
	if (n == 0)
		stop("The bootGmm object is empty")
	J <- sapply(1:n, function(i) specTest(obj[[i]]$ans)$test[1])
	J0 <- specTest(attr(obj,"gmm"))$test[1]	
	F <- ecdf(J)
	pval <- 1-F(J0)
	list(test = c(Stats=J0, BootPVal=pval), JCDF = F)
	}	

linearHypothesis.bootGmm <- function(model, hypothesis.matrix, rhs = NULL, ...)
	{
	obj <- model
	n <- length(obj)
	if (n == 0)
		stop("The bootGmm object is empty")
	Test0 <- linearHypothesis(attr(obj,"gmm"), hypothesis.matrix, rhs=rhs, ...)[[2]][2]

	# The following is borrowed from the car::linearHypothesis.default
	if (is.character(hypothesis.matrix)) {
		L <- car:::makeHypothesis(names(attr(obj,"gmm")$coefficients), hypothesis.matrix, rhs = NULL)
		if (is.null(dim(L))) L <- t(L)
		rhs <- L[, NCOL(L)]
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	}
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
				else hypothesis.matrix
		if (is.null(rhs)) rhs <- rep(0,nrow(L))
	}
	hyp <- car:::printHypothesis(L, rhs, names(attr(obj,"gmm")$coefficients))
	
	##############
	rhsB <- L%*%attr(obj,"gmm")$coefficients

	Test <- sapply(1:n, function(i) linearHypothesis(obj[[i]]$ans, L, rhsB, ...)[[2]][2])
	F <- ecdf(Test)
	pval <- 1-F(Test0)
	list(hyp = hyp, test = c(Stats=Test0, BootPVal=pval), TCDF = F)
	}


