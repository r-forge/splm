`conclikpan.M` <-
function (rho, W, I, e.a, e.b, e.c, n, nW, nChol, pChol, quiet,NT,T) 
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

