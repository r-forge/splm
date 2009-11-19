conclikpan <- function(lambda, eig, e0e0, e1e1,e0e1, N,T,NT,quiet=quiet)
{
	Nsig <- e0e0 - 2*lambda*e0e1 + lambda*lambda*e1e1
	sigma2 <- Nsig/NT
	if (is.complex(eig)) det <-  Re(prod(1 - lambda*eig)) 
	else det <-  prod(1 - lambda*eig)
	ret <-  T * log(det) - ( (NT)/2 )*log(2*pi*sigma2) - (1/(2*sigma2))*Nsig
	#ret <- - (NT/2)*log(Nsig)  + T * log(det)  
#print(T*log(det))
	  if (!quiet) 
        cat("(eigen) rho:\t", lambda, "\tfunction value:\t", ret, 
            "\n")

	ret
}

conclikpan.sp<-function (rho, W, I, e.a, e.b, e.c, n, quiet,NT,T) 
{
#stop(print((I - rho * W)))
    SSE <- e.a - 2 * rho * e.b + rho * rho * e.c
    s2 <- SSE/NT
    J1 <- try(determinant((I - rho * W), logarithm = TRUE)$modulus, 
        silent = TRUE)
        #print(class(J1))
    if (class(J1) == "try-error") {
        Jacobian <- NA
    }
    else {
        Jacobian <- J1
    }
    ret <- (T*Jacobian - ((NT/2) * log(2 * pi*s2))  - 
        (1/(2 * s2)) * SSE)
    if (!quiet) 
        cat("(spam) rho:\t", rho, "\tfunction value:\t", ret, 
            "\n")
    ret
}


conclikpan.M<-function (rho, W, I, e.a, e.b, e.c, n, nW, nChol, pChol, quiet,NT,T) 
{
	#print(interval)
    #    print(T)
    SSE <- e.a - 2 * rho * e.b + rho * rho * e.c
    s2 <- SSE/NT
    if (isTRUE(all.equal(rho, 0))) {
        Jacobian <- rho
    }
    else if (rho > 0) {
        detTRY <- try(Matrix:::ldetL2up(nChol, nW, 1/rho), silent = TRUE)
        if (class(detTRY) == "try-error") {
            Jacobian <- NaN
        }
        else {
            Jacobian <- n * log(rho) + T*detTRY
        }
    }
    else {
        detTRY <- try(Matrix:::ldetL2up(pChol, W, 1/(-rho)), 
            silent = TRUE)
        if (class(detTRY) == "try-error") {
            Jacobian <- NaN
        }
        else {
            Jacobian <- n * log(-(rho)) + T*detTRY
        }
    }
    #print(T)
#    print(Jacobian)
    ret <- (Jacobian - ((NT/2) * log(2 * pi*s2))  - 
        (1/(2 * s2)) * SSE)
    if (!quiet) 
        cat("(Matrix) rho:\t", rho, "\tfunction value:\t", ret, 
            "\n")
    ret
}
