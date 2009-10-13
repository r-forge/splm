`dosparsepanel` <-
function (listw, y, x, wy, K, quiet=quiet, tol.opt, method, interval, 
    can.sim, zero.policy = FALSE,NT,T) 
{
    similar <- FALSE
    m <- ncol(x)
    n <- nrow(x)
    if (method == "spam") {
        if (listw$style %in% c("W", "S") & can.sim) {
            W <- listw2U_spam(similar.listw_spam(listw))
            similar <- TRUE
        }
        else W <- as.spam.listw(listw)
        I <- diag.spam(1, n/T, n/T)
        #print(I)
    }
    else if (method == "Matrix") {
        if (listw$style %in% c("W", "S") & can.sim) {
            W <- listw2U_Matrix(similar.listw_Matrix(listw))
            similar <- TRUE
        }
        else W <- as_dsTMatrix_listw(listw)
        W <- as(W, "CsparseMatrix")
        I <- as_dsCMatrix_I(n/T)
        #print(I)
        Imult <- 2
        if (listw$style == "B") {
            Imult <- ceiling((2/3) * max(apply(W, 1, sum)))
            interval <- c(-0.5, +0.25)
        }
        else interval <- c(-2, +1)
        nW <- -W
        pChol <- Cholesky(W, super = FALSE, Imult = Imult)
        nChol <- Cholesky(nW, super = FALSE, Imult = Imult)
        ns1 <- last <- 10
        prho1 <- seq(sqrt(.Machine$double.eps), interval[2], 
            length.out = ns1)
        while (last >= ns1) {
            pdet1 <- Matrix:::ldetL2up(nChol, nW, 1/prho1)
            wp1 <- which(is.finite(pdet1))
            last <- wp1[length(wp1)]
            if (last == ns1) 
                prho1 <- seq(interval[2], 1.5 * interval[2], 
                  length.out = ns1)
        }
        lwp1n <- prho1[last]
        lwp2n <- prho1[last + 1]
        prho2 <- seq(lwp2n, lwp1n, length.out = ns1)
        pdet2 <- Matrix:::ldetL2up(nChol, nW, 1/prho2)
        wp2 <- which(is.finite(pdet2))
        lwp2n <- prho2[wp2[length(wp2)]]
        nrho1 <- seq(interval[1], -sqrt(.Machine$double.eps), 
            length.out = ns1)
        first <- 1
        while (first == 1) {
            ndet1 <- Matrix:::ldetL2up(pChol, W, 1/(-nrho1))
            wn1 <- which(is.finite(ndet1))
            first <- wn1[1]
            if (first == 1) 
                prho1 <- seq(1.5 * interval[1], interval[1], 
                  length.out = ns1)
        }
        lwn1n <- nrho1[wn1[1]]
        lwn2n <- nrho1[wn1[1] - 1]
        nrho2 <- seq(lwn2n, lwn1n, length.out = ns1)
        ndet2 <- Matrix:::ldetL2up(pChol, W, 1/(-nrho2))
        wn2 <- which(is.finite(ndet2))
        lwn2n <- nrho2[wn2[1]]
        interval <- c(lwn2n, lwp2n)
        if (!quiet) 
            cat("using interval:", interval, "\n")
    }
    LLs <- NULL
    if (m > 1) {
        LLs <- vector(mode = "list", length = length(K:m))
        j <- 1
        for (i in K:m) {
            thisx <- x[, -i, drop = FALSE]
            lm.null <- lm.fit(thisx, y)
            lm.w <- lm.fit(thisx, wy)
            e.null <- lm.null$residuals
            e.w <- lm.w$residuals
            e.a <- t(e.null) %*% e.null
            e.b <- t(e.w) %*% e.null
            e.c <- t(e.w) %*% e.w
            if (method == "spam") {
                LLs[[j]] <- optimize(conclikpan.sp, interval = interval, 
                  maximum = TRUE, tol = tol.opt, W = W, I = I, 
                  e.a = e.a, e.b = e.b, e.c = e.c, n = n, T=T,NT=NT, quiet = quiet)$objective
            }
            else if (method == "Matrix") {
            	#print(interval)
        		#print(T)
                LLs[[j]] <- optimize(conclikpan.M, interval = interval, 
                  maximum = TRUE, tol = tol.opt, W = W, I = I, 
                  e.a = e.a, e.b = e.b, e.c = e.c, n = n, nW = nW, 
                  nChol = nChol, pChol = pChol, T=T, NT=NT, quiet = quiet)$objective
            }
            attr(LLs[[j]], "nall") <- n
            attr(LLs[[j]], "nobs") <- n
            attr(LLs[[j]], "df") <- (m + 2) - 1
            attr(LLs[[j]], "name") <- colnames(x)[i]
            class(LLs[[j]]) <- "logLik"
            j <- j + 1
        }
    }
    lm.null <- lm(y ~ x - 1)
    lm.w <- lm.fit(x, wy)
    e.null <- lm.null$residuals
    e.w <- lm.w$residuals
    e.a <- t(e.null) %*% e.null
    e.b <- t(e.w) %*% e.null
    e.c <- t(e.w) %*% e.w
    if (method == "spam") {
        opt <- optimize(conclikpan.sp, interval = interval, 
            maximum = TRUE, tol = tol.opt, W = W, I = I, e.a = e.a, 
            e.b = e.b, e.c = e.c, n = n,T=T, NT=NT, quiet = quiet)
    }
    else if (method == "Matrix") {
        opt <- optimize(conclikpan.M, interval = interval, 
            maximum = TRUE, tol = tol.opt, W = W, I = I, e.a = e.a, 
            e.b = e.b, e.c = e.c, n = n, nW = nW, nChol = nChol, 
            pChol = pChol, T=T, NT=NT, quiet = quiet)
    }
    maximum <- opt$maximum
    objective <- opt$objective
    res <- list(maximum = maximum, objective = objective, LLs = LLs, 
        lm.null = lm.null, similar = similar, opt = opt)
    res
}

