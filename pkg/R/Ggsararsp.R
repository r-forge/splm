`Ggsararsp` <-
function (W, u, zero.policy = FALSE) 
{
      n <- length(u)
      tt<-matrix(0,n,1)
      tr<-sum(unlist(W$weights)^2)
      wu<-lag.listw(W,u)
      wwu<-lag.listw(W,wu)
    	uu <- crossprod(u, u)
    	uwu <- crossprod(u, wu)
 	uwpuw <- crossprod(wu, wu)
    	uwwu <- crossprod(u, wwu)
    	wwupwu <- crossprod(wwu, wu)
    	wwupwwu <- crossprod(wwu, wwu)
    	bigG <- matrix(0, 3, 3)
    	bigG[, 1] <- c(2 * uwu, 2 * wwupwu, (uwwu + uwpuw))/n
    	bigG[, 2] <-  -c(uwpuw, wwupwwu, wwupwu)/n
    	bigG[, 3] <- c(1, tr/n, 0)
    	litg <- c(uu, uwpuw, uwu)/n
    	list(bigG = bigG, litg = litg)
}

