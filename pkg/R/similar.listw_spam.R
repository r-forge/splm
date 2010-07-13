`similar.listw_spam` <-
function(listw) 
{
    nbsym <- attr(listw$neighbours, "sym")
    if (is.null(nbsym)) 
        nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
    if (!nbsym) 
        stop("Only symmetric nb can yield similar to symmetric weights")
    if (attr(listw$weights, "mode") == "general") 
        if (!attr(listw$weights, "glistsym")) 
            stop("General weights must be symmetric")
    n <- length(listw$neighbours)
    if (n < 1) 
        stop("non-positive number of entities")
    sww <- as.spam.listw(listw)
    if (listw$style == "W") {
        sd <- attr(listw$weights, "comp")$d
        sd1 <- 1/(sqrt(sd))
        sdd <- diag.spam(sd, n, n)
        sdd1 <- diag.spam(sd1, n, n)
        sww1 <- sdd %*% sww
        res <- sdd1 %*% sww1 %*% sdd1
    }
    else if (listw$style == "S") {
        q <- attr(listw$weights, "comp")$q
        Q <- attr(listw$weights, "comp")$Q
        eff.n <- attr(listw$weights, "comp")$eff.n
        q1 <- 1/(sqrt(q))
        qq <- diag.spam(q, n, n)
        qq1 <- diag.spam(q1, n, n)
        ww0 <- (Q/eff.n) * sww
        ww1 <- qq %*% ww0
        sim0 <- qq1 %*% ww1 %*% qq1
        res <- (eff.n/Q) * sim0
    }
    else stop("Conversion not suitable for this weights style")
    res
}

