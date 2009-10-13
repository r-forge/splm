`feerror` <-
function(y,x,ysms,xsms,ytms, xtms, beta,sige,yt,xt,N,T,NT,k,effects,lambda){
mx<-apply(x,2,mean)
intercept<-mean(y)-mx%*%beta
if (effects =="spfe"){
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.sfe,T) + as.numeric(intercept)
	res.t.sfe <- res.sfe / sqrt(sige / T* rep(1,N) + diag(as.numeric(sige) * xsms%*% crossprod(xt)%*% t(xsms) ))
	res.t.con <- as.numeric(intercept) / sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% crossprod(xt) %*% as.matrix(mx))
	N.vars <- k + N
	res.e <- y - xhat
	FE.out<-list(res.sfe=res.sfe, res.t.sfe=res.t.sfe, intercept=intercept, res.t.con=res.t.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
	}
if (effects== "tpfe")	{
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
	xhat <- x %*% as.matrix(beta) + rep(res.tfe,each=N) + as.numeric(intercept)
	res.t.tfe <- res.tfe / sqrt(sige / N* rep(1,T) + diag(as.numeric(sige) * xtms%*% crossprod(xt)%*% t(xtms) ))
	res.t.con <- as.numeric(intercept) / sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% crossprod(xt) %*% as.matrix(mx))
	N.vars <- k + T
		res.e <- y - xhat
		FE.out<-list(res.tfe=res.tfe, res.t.tfe=res.t.tfe, intercept=intercept, res.t.con=res.t.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
		}
if (effects== "sptpfe"){
	res.sfe <- as.matrix(ysms) - xsms %*% as.matrix(beta) - as.numeric(intercept)
	res.tfe <- as.matrix(ytms) - xtms %*% as.matrix(beta) - as.numeric(intercept)
	res.t.sfe <- res.sfe / sqrt(sige / T* rep(1,N) + diag(as.numeric(sige) * xsms%*% crossprod(xt)%*% t(xsms) ))
	res.t.tfe <- res.tfe / sqrt(sige / N* rep(1,T) + diag(as.numeric(sige) * xtms%*% crossprod(xt)%*% t(xtms) ))
	res.t.con <- as.numeric(intercept) / sqrt(as.numeric(sige) / NT + as.numeric(sige) * t(as.matrix(mx)) %*% crossprod(xt) %*% as.matrix(mx))
	xhat<- x %*% as.matrix(beta) + rep(res.sfe,T) + rep(res.tfe,each=N) + as.numeric(intercept)
	N.vars <- k + N + T - 1
	res.e <- y - xhat
FE.out<-list(res.tfe=res.tfe, res.t.tfe=res.t.tfe, res.sfe=res.sfe, res.t.sfe=res.t.sfe, intercept=intercept, res.t.con=res.t.con,xhat=xhat,N.vars=N.vars,res.e=res.e)
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

