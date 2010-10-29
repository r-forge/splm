`searg` <-
function (rhopar, v, verbose = FALSE) 
{
    vv <- v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$litg
    value <- sum(vv^2)
    if (verbose) 
        cat("function:", value, "rho:", rhopar[1], "sig2:", 
            rhopar[2], "\n")
    value
}



`arg` <-
function (rhopar, v, verbose = verbose) 
{
    vv <-  v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$smallg
    value <- sum(vv^2)
    if (verbose) 
        cat("function:", value, "rho:", rhopar[1], "sig2:", 
            rhopar[2], "\n")
    value
}


`arg1` <-
function (rhopar, v, ss,SS,T, verbose = verbose) 
{
	Ga<-cbind( (1/(T-1))*ss^2,0)
	Gb<-cbind( 0, SS^2)
	Gc<-rbind(Ga,Gb)
	Gamma<-kronecker(Gc,diag(3)) 
	Gammainv<-solve(Gamma)
    vv <-  v$GG %*% c(rhopar[1], rhopar[1]^2, rhopar[2], rhopar[3]) - v$gg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}



`arg2` <-
function (rhopar, v, ss,SS,T,TW, verbose = verbose) 
{
	Ga<-cbind( (1/(T-1))*ss^2,0)
	Gb<-cbind( 0, SS^2)
	Gc<-rbind(Ga,Gb)
	Gamma<-kronecker(Gc,TW) 
	Gammainv<-solve(Gamma)
    vv <-  v$GG %*% c(rhopar[1], rhopar[1]^2, rhopar[2], rhopar[3]) - v$gg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}




`arg3` <-
function (rhopar, v, ss,T,TW, verbose = verbose) 
{
	Ga<-(1/(T-1))*ss^2
	Gamma<-as.numeric(Ga)*TW 
	Gammainv<-solve(Gamma)
    vv <-  v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$smallg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}



`fs` <-
function(listw,u,N,T){
	ind<-seq(1,T)
	inde<-rep(ind,each=N)
	NT<-N*T
#ub<-matrix(,NT,1)
#for (i in 1:T) ub[(N*i-N+1):(N*i),]<-lag.listw(listw,u[(N*i-N+1):(N*i)])
#ubb<-matrix(,NT,1)
#for (i in 1:T) ubb[(N*i-N+1):(N*i),]<-lag.listw(listw,ub[(N*i-N+1):(N*i)])
#ubmt<-matrix(,N,1)
#for (i in 1:N) ubmt[i,]<-mean(ub[seq.int(from=i,to=NT,by=N)])
#ubbmt<-matrix(,N,1)		
#for (i in 1:N) ubbmt[i,]<-mean(ubb[seq.int(from=i,to=NT,by=N)])
#umt<-matrix(,N,1)		
#for (i in 1:N) umt[i,]<-mean(u[seq.int(from=i,to=NT,by=N)])
	ub<-unlist(tapply(u,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))###Wu
	ubb<-unlist(tapply(ub,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))###WWu
	ind1<-seq(1,N)
	inde1<-rep(ind1,T)
	umt<-tapply(u, inde1, mean) 
	ubmt<-tapply(ub, inde1, mean) ####mean over time 
	ubbmt<-tapply(ubb, inde1, mean)
	umNT<-rep(umt,T)
	ubmNT<-rep(ubmt,T)
	ubbmNT<-rep(ubbmt,T)
	Q0ub<-ub-ubmNT
	Q0ubb<-ubb-ubbmNT
	Q0u<-u-umNT
	uQ0ub<-crossprod(u,Q0ub)
	ubQ0ubb<-crossprod(ub,Q0ubb)
	ubbQ0ubb<-crossprod(ubb,Q0ubb)
	ubQ0ub<-crossprod(ub,Q0ub)
	uQ0u<-crossprod(u,Q0u)
	ubbQ0ub<-crossprod(ubb,Q0ub)
	uQ0ubb<-crossprod(u,Q0ubb)
	tr <- matrix(0, N, 1)
	for (i in 1:N) {
        tr[i] <- sum(listw$weights[[i]]^2)
	}      
	trwpw <- sum(tr)
	G1c<-(1/(N*(T-1)))*rbind(2*uQ0ub, 2*ubbQ0ub,(uQ0ubb+ubQ0ub) )
	G2c<- (-1/(N*(T-1)))* rbind(ubQ0ub,ubbQ0ubb,ubQ0ubb)
	G3c<- rbind(1,trwpw/N, 0)
	G<-cbind(G1c,G2c,G3c)	
	g<-(1/(N*(T-1)))*rbind(uQ0u,ubQ0ub,uQ0ub)
	output<-list(bigG=G, smallg=g, Q1u=umNT,Q1ub=ubmNT,Q1ubb=ubbmNT, ub=ub,ubb=ubb,TR=trwpw)
}

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


`tw` <-
function(W,N){
## depends on listw2dgCMatrix.R
	Ws<-listw2dgCMatrix(W)
	Wst<-t(Ws)
	WpW<-crossprod(Ws)
	WpWWpW<-WpW%*%WpW
	WppW<-Wst+Ws
	WpWWppW<-WpW%*%WppW
	WW<-Ws%*%Ws
	WWpWpW<-WW+WpW
	tr1<-sum(diag(WpW/N))
	tr2<-sum(diag(WpWWpW/N))
	tr3<-sum(diag(WpWWppW/N))
	tr4<-sum(diag(WWpWpW/N))
	TW1c<-rbind(2,2*tr1, 0)
	TW2c<-rbind(2*tr1, 2*tr2,tr3)
	TW3c<-rbind(0, tr3,tr4)
	TW<-cbind(TW1c,TW2c,TW3c)
	out<-list(TW=matrix(TW,3,3))
}

`pw` <-
function(bigG, smallg, Q1u,Q1ub,Q1ubb, u, ub,ubb,N, TR){
	uQ1u<-crossprod(u,Q1u)
	uQ1ub<-crossprod(u,Q1ub)
	ubbQ1ub<-crossprod(ubb,Q1ub)
	ubbQ1ubb<-crossprod(ubb,Q1ubb)
	uQ1ubb<-crossprod(u,Q1ubb)
	ubQ1ub<-crossprod(ub,Q1ub)
	ubQ1ubb<-crossprod(ub,Q1ubb)
	G1c1<-rbind(2*uQ1ub, 2*ubbQ1ub,  (uQ1ubb + ubQ1ub))/N
	G1c2<-rbind(ubQ1ub, ubbQ1ubb, ubQ1ubb)/(-N)
	G1c3<-rbind(1,TR/N,0)
	G1c<-cbind(G1c1,G1c2,rep(0,3),G1c3)
	g1<-rbind(uQ1u, ubQ1ub, uQ1ub)/N
	GG<-rbind(cbind(bigG,rep(0,3)),G1c)	
	gg<-rbind(smallg,g1)
	out<-list(GG=GG,gg=gg)
}



`tslssp` <-
function(y,yend,X,Zinst,robust=FALSE) {
	Q <- cbind(X,Zinst)
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
	#print(biv)
	rownames(biv)<-c(colnames(yend), colnames(X))
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

