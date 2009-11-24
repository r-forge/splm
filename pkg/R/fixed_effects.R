feerror<-function(y,x,ysms,xsms,ytms, xtms, beta,sige,yt,xt,N,T,NT,k,effects,lambda){
mx<-apply(x,2,mean)
intercept<-mean(y)-mx%*%beta
if (effects =="spfe"){
	sige<-as.numeric(sige)
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.sfe,T) + as.numeric(intercept)
	res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe<-sqrt(diag(res.var.sfe))
	res.t.sfe <- res.sfe / res.se.sfe 
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + N
	res.e <- y - xhat
	FE.out<-list(res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
	}
if (effects== "tpfe")	{
 	sige<-as.numeric(sige)
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.tfe,each=N) + as.numeric(intercept)
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe <-sqrt(diag(res.var.tfe))
	res.t.tfe <- res.tfe/res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + T
		res.e <- y - xhat
		FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects== "sptpfe"){
	sige<-as.numeric(sige)
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
		res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe <-sqrt(diag(res.var.sfe))
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe<-sqrt(diag(res.var.tfe))
	res.t.sfe <- res.sfe / res.se.sfe
	res.t.tfe <- res.tfe / res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	xhat<- x %*% as.matrix(beta) + rep(res.sfe,T) + rep(res.tfe,each=N) + as.numeric(intercept)
	N.vars <- k + N + T - 1
	res.e <- y - xhat
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects=="pooled") {
	xhat <-   x %*% as.matrix(beta)
	res.e <- y - xhat
	FE.out<-list(xhat=xhat,N.vars=k,res.e=res.e)
	}
	yhat <- xhat
	ywhat <-  xt %*% beta
	r1 <- as.matrix(yt - mean(yt))
	r2 <- as.matrix(ywhat - mean(ywhat))
	r1r2 <- crossprod(r1,r2)
	r1r1 <- crossprod(r1)
	r2r2 <- crossprod(r2)
	res.corr <- as.numeric(r1r2^2) / (as.numeric(r1r1)*as.numeric(r2r2))
FE.out <- list(FE.out, res.corr=res.corr)
	}
	
felag<-function(y,x,wy,ysms,xsms,ytms, xtms, wytms, wysms, beta,sige,yt,xt,N,T,NT,k,effects,method, rho,listw,inde){
		mx<-apply(x,2,mean)
		intercept <- mean(y)- mean(wy)*rho -  mx%*%beta
if (effects=="spfe"){
	res.sfe <- as.matrix(ysms) - as.matrix(wysms) *rho - xsms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.sfe,T) + as.numeric(intercept)
	res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe<-sqrt(diag(res.var.sfe))
	res.t.sfe <- res.sfe / res.se.sfe 
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + N
	res.e <- y - xhat - rho* wy
FE.out<-list(res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, 	res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
	}
if (effects== "tpfe")	{
	res.tfe <- as.matrix(ytms) - as.matrix(wytms)* rho - xtms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.tfe,each=N) + as.numeric(intercept)
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe <-sqrt(diag(res.var.tfe))
	res.t.tfe <- res.tfe/res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% 	solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	N.vars <- k + T
	res.e <- y - xhat - rho* wy
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, intercept=intercept, 	res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects== "sptpfe"){
	res.sfe <- as.matrix(ysms) - as.matrix(wysms) * rho - xsms %*% as.matrix(beta) - as.numeric(intercept)
	res.tfe <- as.matrix(ytms) - as.matrix(wytms) * rho - xtms %*% as.matrix(beta) - as.numeric(intercept)
	res.var.sfe<- (sige / T)  + (as.numeric(sige)*(xsms%*% solve(crossprod(xt)) %*% t(xsms)))
	res.se.sfe <-sqrt(diag(res.var.sfe))
	res.var.tfe <- (sige / N)  + (as.numeric(sige)*(xtms%*% solve(crossprod(xt)) %*% t(xtms)))
	res.se.tfe<-sqrt(diag(res.var.tfe))
	res.t.sfe <- res.sfe / res.se.sfe
	res.t.tfe <- res.tfe / res.se.tfe
	res.se.con<-sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% solve(crossprod(xt)) %*% as.matrix(mx))
	res.t.con <- as.numeric(intercept) / res.se.con
	xhat<- x %*% as.matrix(beta) + rep(res.sfe,T) + rep(res.tfe,each=N) + as.numeric(intercept)
	N.vars <- k + N + T - 1
	res.e <- y - xhat - rho* wy
FE.out<-list(res.tfe=res.tfe, res.se.tfe=res.se.tfe, res.sfe=res.sfe, res.se.sfe=res.se.sfe, intercept=intercept, res.se.con=res.se.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects=="pooled") {
	xhat <-   x %*% as.matrix(beta)
	res.e <- y - xhat - rho* wy
	FE.out<-list(xhat=xhat,N.vars=k,res.e=res.e)
	}
if (method=="eigen"){
	IrWi<-invIrW(listw,rho)
	xtb <- xt %*% beta
	yhat <- unlist(tapply(xhat,inde, function(u) IrWi %*% u))
	ywhat <- unlist(tapply(xtb,inde, function(u) IrWi %*% u))
	r1 <- as.matrix(yt - mean(yt))
	r2 <- as.matrix(ywhat - mean(ywhat))
	r1r2 <- crossprod(r1,r2)
	r1r1 <- crossprod(r1)
	r2r2 <- crossprod(r2)
	res.corr <- as.numeric(r1r2^2) / (as.numeric(r1r1)*as.numeric(r2r2))
}
else res.corr <- NULL 
FE.out <- list(FE.out, res.corr=res.corr)
FE.out
	}