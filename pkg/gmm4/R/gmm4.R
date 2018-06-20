################### the main gmm function ###################

gmm4 <- function(g, x, tet0=NULL,grad=NULL, type=c("twostep", "iter","cue", "onestep"),
                vcov = c("HAC", "MDS", "iid", "TrueFixed"),
                initW=c("tsls","ident"), weights="optimal", itermaxit=50,
                kernel = c("Quadratic Spectral",  "Truncated", "Bartlett", "Parzen",
                    "Tukey-Hanning"), crit = 1e-06,
                bw = "Andrews", prewhite = 1L, ar.method = "ols", approx = "AR(1)", 
                kerntol = 1e-07, itertol=1e-7, centeredVcov = TRUE,
                data=parent.frame(), ...)
    {
        Call <- match.call()
        vcov <- match.arg(vcov)
        kernel <- match.arg(kernel)
        type <- match.arg(type)
        initW <- match.arg(initW)

        model <- gmmModel(g, x, tet0, grad, vcov, kernel, crit = 1e-06,
                          bw, prewhite, ar.method, approx, kerntol, centeredVcov, data)

        
        fit <- gmmFit(model, type, itertol, initW, weights, 
                      itermaxit=100, ...)
        fit@call <- Call
        fit
    }
