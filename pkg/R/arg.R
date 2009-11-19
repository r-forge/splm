`arg` <-
function (rhopar, v, verbose = FALSE) 
{
    vv <-  v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$smallg
    value <- sum(vv^2)
    if (verbose) 
        cat("function:", value, "rho:", rhopar[1], "sig2:", 
            rhopar[2], "\n")
    value
}

