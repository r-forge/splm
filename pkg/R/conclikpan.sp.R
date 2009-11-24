`conclikpan.sp` <-
function (rho, W, I, e.a, e.b, e.c, n, quiet,NT,T) 
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

