	can.be.simmed <- function(listw) {
	res <- is.symmetric.nb(listw$neighbours, FALSE)
	if (res) {
		if (attr(listw$weights, "mode") == "general")
			res <- attr(listw$weights, "glistsym")
	} else return(res)
	res
}


similar.listw_Matrix <- function(listw) {
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	ww <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
	if (listw$style == "W") {
		d <- attr(listw$weights, "comp")$d
		d1 <- 1/(sqrt(d))
		dd <- as(as(Diagonal(x=d), "symmetricMatrix"), "CsparseMatrix")
		dd1 <- as(as(Diagonal(x=d1), "symmetricMatrix"),
		    "CsparseMatrix")
		ww1 <- dd %*% ww
		res <- dd1 %*% ww1 %*% dd1
	} else if (listw$style == "S") {
		q <- attr(listw$weights, "comp")$q
		Q <- attr(listw$weights, "comp")$Q
		eff.n <- attr(listw$weights, "comp")$eff.n
		q1 <- 1/(sqrt(q))
		qq <- as(as(Diagonal(x=q), "symmetricMatrix"), "CsparseMatrix")
		qq1 <- as(as(Diagonal(x=q1), "symmetricMatrix"),
		    "CsparseMatrix")
		ww0 <- (Q/eff.n) * ww
		ww1 <- qq %*% ww0
		sim0 <- qq1 %*% ww1 %*% qq1
		res <- (eff.n/Q) * sim0
	} else stop("Conversion not suitable for this weights style")
	res
}


similar.listw_spam <- function(listw) {
        if (!require(spam)) stop("spam not available")
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	sww <- as.spam.listw(listw)
	if (listw$style == "W") {
		sd <- attr(listw$weights, "comp")$d
		sd1 <- 1/(sqrt(sd))
		sdd <- diag.spam(sd, n, n)
		sdd1 <- diag.spam(sd1, n, n)
		sww1 <- sdd %*% sww
		res <- sdd1 %*% sww1 %*% sdd1
	} else if (listw$style == "S") {
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
	} else stop("Conversion not suitable for this weights style")
	res
}



similar.listw <- function(listw) {
	nbsym <- attr(listw$neighbours, "sym")
	if(is.null(nbsym)) nbsym <- is.symmetric.nb(listw$neighbours, FALSE)
	if (!nbsym) 
		stop("Only symmetric nb can yield similar to symmetric weights")
	if (attr(listw$weights, "mode") == "general")
		if (!attr(listw$weights, "glistsym"))
			stop("General weights must be symmetric")
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(listw$neighbours)
	if (listw$style == "W") {
		d <- attr(listw$weights, "comp")$d
		glist <- vector(mode="list", length=n)
		for (i in 1:n) glist[[i]] <- d[i] * listw$weights[[i]]
		sd1 <- 1/sqrt(d)
		for (i in 1:n) {
			inb <- listw$neighbours[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
				for (j in 1:icd) {
					glist[[i]][j] <- sd1[i] * 
						glist[[i]][j] * sd1[inb[j]]
				}
			}
		}
		res <- listw
		res$weights <- glist
		attr(res$weights, "mode") <- "sim"
		attr(res$weights, "W") <- TRUE
		attr(res$weights, "comp") <- attr(listw$weights, "comp")
		res$style <- "W:sim"
	} else if (listw$style == "S") {
		q <- attr(listw$weights, "comp")$q
		Q <- attr(listw$weights, "comp")$Q
		eff.n <- attr(listw$weights, "comp")$eff.n
		glist <- vector(mode="list", length=n)
		for (i in 1:n) {
			glist[[i]] <- (Q/eff.n) * listw$weights[[i]]
			glist[[i]] <- q[i] * glist[[i]]
		}
		sq1 <- 1/sqrt(q)
		for (i in 1:n) {
			inb <- listw$neighbours[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
				for (j in 1:icd) {
					glist[[i]][j] <- sq1[i] * 
						glist[[i]][j] * sq1[inb[j]]
				}
				glist[[i]] <- (eff.n/Q) * glist[[i]]
			}
		}
		res <- listw
		res$weights <- glist
		attr(res$weights, "mode") <- "sim"
		attr(res$weights, "S") <- TRUE
		attr(res$weights, "comp") <- attr(listw$weights, "comp")
		res$style <- "S:sim"
	} else stop("Conversion not suitable for this weights style")
	sym_out <- is.symmetric.glist(res$neighbours, res$weights)
	if (!sym_out) {
	    if (attr(sym_out, "d") < .Machine$double.eps ^ 0.5)
		res <- listw2U(res)
	    else stop("defective similarity")
	}
	res
}

