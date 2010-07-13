`spregm` <-
function(formula, data=list(), index=NULL, w, 
                   method = c("init", "weigh", "fulweigh"), lag=FALSE){
 ## depends on:
#source('args.R')
#source('listw2dgCMatrix.R')
  ## reorder data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }
  
  index <- data[,1]
  tindex <- data[,2]

  ## record call
  cl <- match.call()

  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")

  ## reduce X,y to model matrix values (no NAs)
  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  ## reduce index accordingly
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  ## det. number of groups and df
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  ## det. max. group numerosity
  T<-max(tapply(x[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)

  ## check w
  if(is.matrix(w)) {
    if(dim(w)[[1]]!=N) stop("Non conformable spatial weights")
    require(spdep)
    listw <- mat2listw(w)
   } else {
    listw <- w
    rm(w)
   }

  ## check whether the panel is balanced
  balanced<-N*T==NT
  if(!balanced) stop("Estimation method unavailable for unbalanced panels")

  ## mostly unchanged from here:

  mt<-terms(formula,data=data)
  mf<-lm(formula,data,method="model.frame")#,na.action=na.fail
  na.act<-attr(mf,'na.action')
  xcolnames <- colnames(x)
  
  ind<-seq(1,T)
  inde<-rep(ind,each=N)


	wy <- unlist(tapply(y,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))    
	wy <- array(wy, c(length(y), 1))
   colnames(wy) <- ("lambda")

if(lag){
    K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[, 1] == 
        1), 2, 1)
    if (any(is.na(wy))) 
        stop("NAs in spatially lagged dependent variable")
    if (k > 1) {
        WX <- matrix(nrow = nrow(x), ncol = (k - (K - 1)))
        WWX <- matrix(nrow = nrow(x), ncol = (k - (K - 1)))
for (i in K:k) {
            wx <- unlist(tapply(x[,i],inde, function(TT) lag.listw(listw,TT), simplify=TRUE))    
            wwx <- unlist(tapply(wx,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))    
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WX[, (i - (K - 1))] <- wx
            WWX[, (i - (K - 1))] <- wwx
        }
    }
instr <- cbind(WX, WWX)

betaIV <- spdep:::tsls(y = y, yend = wy, X = x, Zinst = instr)
res <- residuals(betaIV)
NT<-N*T
df<-NT-dim(x)[2]
S<-sum(res^2)/(N-1)
 	}
 	
else{ 
  XpX<-crossprod(x)
  Xpy<-crossprod(x,y)
  betaOLS<-solve(XpX,Xpy)
  res<-y-x%*%betaOLS
  NT<-N*T
  df<-NT-dim(x)[2]
  S<-sum(res^2)/(N-1)
 }
 
  Gg<-fs(listw,res,N,T)
  pars<-c(0,0)
  estim1 <- optim(pars, arg, v = Gg, method='L-BFGS-B', verbose =FALSE)
  urub<-res- estim1$par[1]*Gg$ub
  Q1urQ1ub<-Gg$Q1u - estim1$par[1]*Gg$Q1ub
  S1<- crossprod(urub, Q1urQ1ub)/N

  method <- match.arg(method)

  switch(method, init = {
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=S1
    }, weigh = {
	Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      pars2<-c(estim1$par[1],estim1$par[2],S1)
      estim2 <- optim(pars2, arg1, v = Ggw,ss=estim1$par[2] ,SS=S1,T=T, method='L-BFGS-B',verbose =FALSE)
	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]
    }, fulweigh = {
	Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      weights<-tw(W=listw,N)
      pars2<-c(estim1$par[1],estim1$par[2],S1)
      estim3 <-optim(pars2, arg2, v = Ggw, ss=estim1$par[2] ,SS=S1,T=T,TW=weights$TW, method='L-BFGS-B',verbose =FALSE)
	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]
    })

  yt<-y-finrho*wy
  dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw,TT), simplify=TRUE))
  xl<-apply(x,2,dm)
  xt<- x-finrho*xl
if(lag){
  w2y<-unlist(tapply(wy,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))
  wyt<- wy-finrho*w2y   	
  	}

  theta<- 1-(sqrt(finsigmaV)/sqrt(finsigma1))
  ind1<-seq(1,N)
  inde1<-rep(ind1,T)

  ytmt<-tapply(yt, inde1, mean)
  ytNT<-rep(ytmt,T)
  yf<-(yt - theta*ytNT)

if(lag){
  dm1<- function(A) rep(unlist(tapply(A,inde1,mean,simplify=TRUE)),T)
  xt<-cbind(wyt,xt) 
  xtNT<-apply(xt,2,dm1)
  xf<-(xt - as.numeric(theta)*xtNT)
  colnames(xf)<-c("lambda", xcolnames)
  H<-cbind(xf[,-1], instr)  
  HH<-crossprod(H)
  Hye<-crossprod(H,xf[,1])  
  bfs<-solve(HH,Hye)
  yendf<-H %*% bfs
  Zf<-cbind(yendf, xf[,-1])
  xfpxf<-crossprod(Zf)
  xfpxfi<-solve(xfpxf)
  betaGLS<- xfpxfi %*% crossprod(Zf,yf)
	}

else{	
  dm1<- function(A) rep(unlist(tapply(A,inde1,mean,simplify=TRUE)),T)
  xtNT<-apply(xt,2,dm1)
  xf<-(xt - as.numeric(theta)*xtNT)
  xfpxf<-crossprod(xf)
  xfpxfi<-solve(xfpxf)
  betaGLS<-xfpxfi%*%crossprod(xf,yf)   
  }
  
  fv<-as.vector(xf%*%betaGLS)
  egls<-yf - fv
  #print(class(egls))
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*xfpxf/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
if(lag)  nam.beta <- c("lambda", dimnames(x)[[2]])
else nam.beta <- dimnames(x)[[2]]
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(betaGLS) <- nam.beta
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
  model.data <- data.frame(cbind(y,x[,-1]))
  sigma2 <- SGLS
  
  
  type <- "random effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=NULL,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL)
  class(spmod) <- "splm"
  return(spmod)
}

