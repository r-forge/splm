sarem2REmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 3),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

    ## extensive function rewriting, Giovanni Millo 29/09/2010
    ## structure:
    ## a) specific part
    ## - set names, bounds and initial values for parms
    ## - define building blocks for likelihood and GLS as functions of parms
    ## - define likelihood
    ## b) generic part(independent from ll.c() and #parms)
    ## - fetch covariance parms from max lik
    ## - calc last GLS step
    ## - fetch betas
    ## - calc final covariances
    ## - make list of results

    # mark
    #print("uso versione 0") # done from saremsrREmod4.R

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi", "lambda", "psi")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999, -0.999)   # lag-specific line (3rd parm)
    upper.bounds <- c(1e+09, 0.999, 0.999)     # lag-specific line (idem)

    ## modules for likelihood
    B <- function(lambda, w) diag(1, ncol(w)) - lambda * w
    detB <- function(lambda, w) det(B(lambda, w))
    invSigma <- function(philambda, n, t, w) {
        Jt <- matrix(1, ncol = t, nrow = t)
        #In <- diag(1, n)
        It <- diag(1, t)
        Jbart <- Jt/t
        Et <- It - Jbart
        ## retrieve parms
        phi <- philambda[1]
        lambda <- philambda[2]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        BB <- crossprod(B(lambda, w))
        invSigma <- kronecker( (1/(t*phi+1)*Jbart + Et), BB )
        invSigma
    }
    detSigma <- function(phi, lambda, n, t, w) {
         Jt <- matrix(1, ncol = t, nrow = t)
        #In <- diag(1, n)
        It <- diag(1, t)
        Jbart <- Jt/t
        Et <- It - Jbart
        detSigma <- -n/2*log( det( (t*phi+1) * Jbart + Et) ) +
            t*log(detB(lambda, w))
        detSigma
    }

    ## likelihood function, both steps included
    ll.c <- function(philambda, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        phi <- philambda[1]
        lambda <- philambda[2]
        psi <- philambda[3]                       # lag-specific line
        ## calc inverse sigma
        sigma.1 <- invSigma(philambda, n, t, w)
        ## lag y
        Ay <- y - psi * wy                        # lag-specific line
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, Ay, sigma.1)         # lag-specific line (Ay for y)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- t*log(detB(psi, w2))              # lag-specific line (else zero <- 0)
        due <- detSigma(phi, lambda, n, t, w)
        tre <- -n * t/2 * log(s2e)
        quattro <- -1/(2 * s2e) * crossprod(e, sigma.1) %*% e
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + zero + due + tre + quattro
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here

    ## calc. Wy (spatial lag of y)
    Wy <- function(y, w, tind) {                  # lag-specific line
        wy<-list()                                # lag-specific line
        for (j in 1:length(unique(tind))) {       # lag-specific line
             yT<-y[tind==unique(tind)[j]]         # lag-specific line
             wy[[j]] <- w %*% yT                  # lag-specific line
             }                                    # lag-specific line
        return(unlist(wy))                        # lag-specific line
    }                                             # lag-specific line

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(crossprod(X, sigma.1) %*% X,
                       crossprod(X, sigma.1) %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (crossprod(ehat, sigma.1) %*% ehat)/(n * t)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y once for all
    wy <- Wy(y, w2, tind)                          # lag-specific line

    ## max likelihood
    optimum <- nlminb(start = myparms0, objective = ll.c,
                      gradient = NULL, hessian = NULL,
                      y = y, X = X, n = n, t = t, w = w, w2 = w2, wy = wy,
                      scale = 1, control = list(x.tol = x.tol,
                                 rel.tol = rel.tol, trace = trace),
                      lower = lower.bounds, upper = upper.bounds)

    ## log likelihood at optimum (notice inverted sign)
    myll <- -optimum$objective
    ## retrieve optimal parms
    myparms <- optimum$par

    ## one last GLS step at optimal vcov parms
    sigma.1 <- invSigma(myparms, n, t, w)
    Ay <- y - myparms[length(myparms)] * wy       # lag-specific line
    beta <- GLSstep(X, Ay, sigma.1)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, sigma.1) %*% X)

    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    nvcovpms <- length(nam.errcomp) - 1
    covAR <- covTheta[nvcovpms+1, nvcovpms+1, drop=FALSE]
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

    ## final parms
    betas <- as.vector(beta[[1]])
    arcoef <- myparms[which(nam.errcomp=="psi")]  # lag-specific line
    errcomp <- myparms[which(nam.errcomp!="psi")]
    names(betas) <- nam.beta
    names(arcoef) <- "psi"                        # lag-specific line
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="psi")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covAR) <- list(names(arcoef), names(arcoef))
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll)

    return(RES)
}