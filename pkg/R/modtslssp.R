`modtslssp` <-
function(y,yend,X,Zinst) {
	Q <- Zinst
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
	QQ <- crossprod(Q,Q)
	Qye <- crossprod(Q,yend)
	bz <- solve(QQ,Qye)
	yendp <- Q %*% bz
	Zp <- cbind(yendp,X)
	ZpZp <- crossprod(Zp,Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	biv <- crossprod(ZpZpi,Zpy)
	#rownames(biv)<-c(colnames(yend), colnames(X))
	yp <- Z %*% biv
	e <- y - yp
    	s2 <- crossprod(e,e) / df
	    varb <- ZpZpi * s2[1,1]
	    sebiv <- sqrt(diag(varb))
	    tbiv <- biv / sebiv
	    pbiv <- pnorm(abs(tbiv),lower.tail=FALSE) * 2
	    result <- list(coefficients=biv,se=sebiv,t=tbiv,
	          p=pbiv,var=varb,s2=s2,
	          residuals=e,yhat=yp)
	    
	result
}

