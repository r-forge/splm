`similar.listw_Matrix` <-
function (listw) 
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
    ww <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
    if (listw$style == "W") {
        d <- attr(listw$weights, "comp")$d
        d1 <- 1/(sqrt(d))
        dd <- as(as(Diagonal(x = d), "symmetricMatrix"), "CsparseMatrix")
        dd1 <- as(as(Diagonal(x = d1), "symmetricMatrix"), "CsparseMatrix")
        ww1 <- dd %*% ww
        res <- dd1 %*% ww1 %*% dd1
    }
    else if (listw$style == "S") {
        q <- attr(listw$weights, "comp")$q
        Q <- attr(listw$weights, "comp")$Q
        eff.n <- attr(listw$weights, "comp")$eff.n
        q1 <- 1/(sqrt(q))
        qq <- as(as(Diagonal(x = q), "symmetricMatrix"), "CsparseMatrix")
        qq1 <- as(as(Diagonal(x = q1), "symmetricMatrix"), "CsparseMatrix")
        ww0 <- (Q/eff.n) * ww
        ww1 <- qq %*% ww0
        sim0 <- qq1 %*% ww1 %*% qq1
        res <- (eff.n/Q) * sim0
    }
    else stop("Conversion not suitable for this weights style")
    res
}

