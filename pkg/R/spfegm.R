`spfegm` <-
function(formula, data=list(), index=NULL, w, method = c("init","fulweigh"), effects=c('pooled','spfe','tpfe','sptpfe'), lag=FALSE){
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

  #make sure that the model has no intercept if effects !=pooled
  if (effects !="pooled" && colnames(x)[1]=="(Intercept)") {
  	x<-x[,-1]
  	cat('\n Warning: x may not contain an intercept if fixed effects are specified \n')
	}

if(is.vector(x)){
  k<-1 
  x<-matrix(x)
  colnames(x)<-clnames[-1]
  dimnames(x)[[1]]<-rwnames
  dimnames(x)[[2]]<-clnames[-1]
  }
  else   k<-dim(x)[[2]]



  ## mostly unchanged from here:

  mt<-terms(formula,data=data)
  mf<-lm(formula,data,method="model.frame")#,na.action=na.fail
  na.act<-attr(mf,'na.action')
  xcolnames <- colnames(x)
  
	ind<-seq(1,T)
	inde<-as.numeric(rep(ind,each=N)) ####takes the first n observations
	indic<-seq(1,N)
	indec<-rep(indic,T) ####takes observations 1,  n+1, 2n+1...

	wy <- unlist(tapply(y,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))    
	wy <- array(wy, c(length(y), 1))
   colnames(wy) <- ("lambda")



if (effects=="tpfe" | effects=="sptpfe"){
	ytms<-tapply(y,inde,mean) ####for each time period takes the mean for the cross section observations
	tpms<-function(q) tapply(q,inde,mean)
	xtms<-apply(x,2,tpms)   ###same thing for the X variable
	ytm<-rep(ytms,each=N) ##generate the NT variables
	xtm<-matrix(,NT,k)
	for (i in 1:k) xtm[,i]<-rep(xtms[,i],each=N)

		wytms<-tapply(wy,inde,mean) ###same thing for Wy
		wytm<-rep(wytms,each=N)

	}

if (effects=="spfe" | effects=="sptpfe"){
	ysms<-tapply(y,indec,mean) ###for each cross-sectional unit takes the mean over the time periods
	spms<-function(q) tapply(q,indec,mean)
	xsms<-apply(x,2,spms)
	ysm<-rep(ysms,T)
	xsm<-matrix(,NT,k)
	for (i in 1:k) xsm[,i]<-rep(xsms[,i],T)
			wysms<-tapply(wy,indec,mean)
			wysm<-rep(wysms,T)

	}

if (effects=='pooled'){
	yt<-y  	###keep the variables with no transformation
	xt<-x
	wyt<-wy
	}

if (effects=="tpfe"){ ####generate the demeaned variables for tpfe
	yt<-y-ytm
	xt<-x-xtm
wyt<-wy-wytm
						}
if(effects=="spfe"){ ####generate the demeaned variables for spfe
	yt<-y-ysm
	xt<-x-xsm
wyt<-wy-wysm
	 					}

if (effects=="sptpfe"){ ####generate the demeaned variables for both types of FE
	yt<-y - ysm - ytm + rep(mean(y),NT)
	xmm<-matrix(,NT,(k))
	for (i in 1:(k)) xmm[,i]<-rep(mean(x[,i]),NT)
	xt<-x - xsm - xtm + xmm
wyt<- wy - wysm - wytm + rep(mean(wy),NT)							}
	
	colnames(xt)<-dimnames(x)[[2]]

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

betaIV <- spdep:::tsls(y = yt, yend = wyt, X = xt, Zinst = instr)
res <- residuals(betaIV)
NT<-N*T
df<-NT-dim(x)[2]
S<-sum(res^2)/(N-1)
 	}
 	
else{ 
  XpX<-crossprod(xt)
  Xpy<-crossprod(xt,yt)
  betaOLS<-solve(XpX,Xpy)
  res<-yt-xt%*%betaOLS
  NT<-N*T
  df<-NT-dim(x)[2]
  S<-sum(res^2)/(N-1)
 }
 
##moments procedure

  Gg<-fs(listw,res,N,T)
  pars<-c(0,0)
  estim1 <- nlminb(pars, arg, v = Gg, control = list(), verbose =FALSE)
  urub<-res- estim1$par[1]*Gg$ub
  Q1urQ1ub<-Gg$Q1u - estim1$par[1]*Gg$Q1ub
  S1<- crossprod(urub, Q1urQ1ub)/N

  method <- match.arg(method)

  switch(method, init = {
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
    }, fulweigh = {
      weights<-tw(W=listw,N)
      pars2<-c(estim1$par[1],estim1$par[2])
      
      estim2 <- nlminb(pars2, arg3, v = Gg, ss=estim1$par[2], T=T,TW = weights$TW ,verbose =FALSE)

	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
    })

    

  yf<-yt-finrho*wyt
  dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw,TT), simplify=TRUE))
  xtl<-apply(xt,2,dm)
  xt2<- xt-finrho*xtl

if(lag){
  w2yt<-unlist(tapply(wyt,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))
  wyt2<- wyt - finrho*w2yt   	
  	}

  ind1<-seq(1,N)
  inde1<-rep(ind1,T)


if(lag){
  xf<-cbind(wyt2,xt2) 
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
  xf<-xt2
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
  errcomp<-rbind(finrho,finsigmaV)
if(lag)  nam.beta <- c("lambda", dimnames(x)[[2]])
else nam.beta <- dimnames(x)[[2]]
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(betaGLS) <- nam.beta
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
  model.data <- data.frame(cbind(y,x[,-1]))
  sigma2 <- SGLS
  
  
  type <- "fixed effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=NULL,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL)
  class(spmod) <- "splm"
  return(spmod)
}

