SEIR <- function(t, y, y0, c=0, sigma=1/5.2, gamma=1/18,
                 r0=matrix(c(0,20,70,84,90,3,2.6,1.9,1,.5),ncol=2),
                 type=c("Lin", "Const"))    
{
    type <- match.arg(type)
    type <- ifelse(type=="Lin", 0, 1)
    par <- c(c, sum(y0), sigma, gamma)
    res <- .Fortran(F_sysseir, as.double(t),
                    as.double(y), as.double(par), as.integer(nrow(r0)),
                    as.double(r0[,2]), as.double(r0[,1]), as.integer(type),
                    RZero=double(1), dy=double(4))$dy
}

solveSEIR <- function(h=1e-2, T=180,
                      c=0, y0=c(11e6, 40, 800, 0),
                      sigma=1/5.2, gamma=1/18,
                      r0=matrix(c(0,20,70,84,90,3,2.6,1.9,1,.5),ncol=2),
                      type=c("Lin", "Const"))
{
    type <- match.arg(type)
    typeF <- ifelse(type=="Lin", 0, 1)
    n <- floor(T/h)
    res <- .Fortran(F_solveseir, as.integer(n), as.double(0),
                    as.double(T), as.double(c), as.double(sigma),
                    as.double(gamma), as.double(y0),
                    as.integer(nrow(r0)),
                    as.double(r0[,2]), as.double(r0[,1]), as.integer(typeF),
                    RZero=double(n+1), t=double(n+1), y=double((n+1)*4),
                    h=double(1))
    ans <- list(y=matrix(res$y,ncol=4), t=res$t, h=res$h,
                RZero=res$RZero)
    colnames(ans$y) = c("S","E","I","R")
    class(ans) <- "seir"
    ans$Pop <- sum(y0)
    ans$h <= res$h
    ans
}

plot.seir <- function(x, which=c("acase","tcase", "S", "E", "I", "R", "both","rzero"),
                      add=FALSE, bothAdd=list(), ...)
{
    which <- match.arg(which)
    x$y <- x$y/1e3
    acase <- rowSums(x$y[,2:3])
    tcase <- rowSums(x$y[,2:4])
    if (which == "acase")
    {
        if (!add)
            {
                plot(x=x$t, y=acase, type="l", xlab="Days",
                     ylab = "Thousands",
                     lwd=2,
                     main="SEIR Model: Number of Active Cases", ...)
            } else {
                lines(x=x$t, y=acase, ...)
            }
    } else if (which == "tcase") {
        if (!add)
        {
            plot(x=x$t, y=tcase, type="l", xlab="Days",
                 ylab = "Thousands",
                 lwd=2, main="SEIR Model: Number of Cases", ...)
        } else {
            lines(x=x$t, y=tcase, ...)
        }
    } else if (which %in% c("S","E","I","R")) {
        w <- which(c("S","E","I","R") == which)
        wn <- c("Susceptibles", "Exposed", "Infectious", "Removed")[w]
        main <- paste("SEIR Model: Number of ", wn, sep="")
        y <- x$y[,w]
        if (!add)
        {
            plot(x=x$t, y=y, type="l", xlab="Days",
                 ylab =  "Thousands", lwd=2, main=main, ...)
        } else {
            lines(x=x$t, y=y, ...)
        }        
    } else if (which == "both") {
        ylim <- range(c(tcase,acase))
        if (!add)
        {
            plot(x=x$t, y=acase, type="l", xlab="Days",
                 ylab =  "Thousands", ylim=ylim, lwd=2,
                 main = "SEIR Model: Total and Active Cases", ...)
            lines(x=x$t, y=tcase, col=2, lwd=2) 
            legend("topleft", c("Active Cases", "Cumulative"),
                   col=1:2, lty=1, lwd=2, bty='n')
        } else {
            lines(x=x$t, y=acase, ...)
            bothAdd$x <- x$t
            bothAdd$y <- acase
            do.call(lines, bothAdd)
        }
    } else {
        plot(x=x$t, y=x$RZero,
             main="R-Zero function\n(number infected by one infectious)",
             xlab="Days", ylab="infected/infectious", type='l', lwd=2, ...)
    }
    invisible()
}

print.seir <- function(x, digits=2, ...)
{
    cat("SEIR model for the spread of a virus\n")
    cat("************************************\n")
    cat("Initial condition\n")
    cat("\tTotal Cases (thousands): ",
        formatC(round(sum(x$y[1,2:4])/1000,digits), format='f'),"\n",sep="")
    cat("\tActive Cases (thousands): ",
        formatC(round(sum(x$y[1,2:3])/1000,digits), format='f'),"\n",sep="")
    cat("Condition after ",
        formatC(round(tail(x$t,1),0), format='d')," Days\n",sep="")
    cat("\tTotal Cases (thousands): ",
        formatC(round(sum(tail(x$y[,2:4],1))/1000,digits), format='f'),"\n",sep="")
    cat("\tActive Cases (thousands): ",
        formatC(round(sum(tail(x$y[,2:3],1))/1000,digits), format='f'),"\n",sep="")
}

"[.seir" <- function(x, i)
{
    ans <- sapply(1:4, function(j) spline(x$t,x$y[,j],xout=i)$y)
    dimnames(ans) <- list(i, colnames(x$y))
    ans
}
