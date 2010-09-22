`sarpanelerror.sp` <-
function (lambda, csrw, I, y, wy, x, WX, n,NT,T, quiet) 
{
    yl <- y - lambda * wy
    xl <- x - lambda * WX
    xl.q <- qr.Q(qr(xl))
    xl.q.yl <- t(xl.q) %*% yl
    SSE <- t(yl) %*% yl - t(xl.q.yl) %*% xl.q.yl
    s2 <- SSE/NT
    J1 <- try(determinant((I - lambda * csrw), logarithm = TRUE)$modulus, 
        silent = TRUE)
    if (class(J1) == "try-error") {
        Jacobian <- NA
    }
    else {
        Jacobian <- J1
    }
    ret <- (T*Jacobian - ((NT/2) * log(2 * pi)) - (NT/2) * log(s2) - 
        (1/(2 * (s2))) * SSE)
    if (!quiet) 
        cat("lambda:", lambda, " function:", ret, " Jacobian:", 
            Jacobian, " SSE:", SSE, "\n")
    ret
}

