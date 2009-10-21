spfeml<-function(formula, data, index=NULL,listw, model=c("lag","error"),effects=c('pooled','spfe','tpfe','sptpfe'), method="eigen",na.action=na.fail,quiet=TRUE,zero.policy = FALSE,
    interval = c(-1, 0.999), tol.solve = 1e-10, tol.opt = .Machine$double.eps^0.5){
	###
	##model should be one between "lag"  or "error"
	##effects should be one of:
		#"pooled" (default) no space or time fixed effects
		#"spfe" : spatial fixed effects
		#"tpfe" : time period fixed effects
		#"sptpfe" : for both types of fixed effects
	##the method argument is taken from lagsarlm
##relates on the following functions:
	  ## reorder data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }

  index <- data[,1]
  tindex <- data[,2]

	  ## record call
  cl <- match.call()

#check the model
model<-match.arg(model)


#check the effects
effects<-match.arg(effects)


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

  #make sure that the model has no intercept if effects !=pooled
  if (effects !="pooled" && colnames(x)[1]=="(Intercept)") {
  	x<-x[,-1]
  	cat('\n Warning: x may not contain an intercept if fixed effects are specified \n')
	}
  ## det. number of groups and df
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  ## det. max. group numerosity
  T<-max(tapply(x[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)

##checks on listw
  if(is.matrix(listw)) {
    if(dim(listw)[[1]]!=N) stop("Non conformable spatial weights")
    require(spdep)
    listw <- mat2listw(listw)
   }
  if (!inherits(listw, "listw"))
        stop("No neighbourhood list")
    can.sim <- as.logical(NA)
    if (listw$style %in% c("W", "S"))
        can.sim <- spdep:::can.be.simmed(listw)
if (!quiet && model=="lag") cat(paste("\n Fixed Effects Spatial lag model\n", "Jacobian calculated using "))
if (!quiet && model=="error") cat(paste("\n Fixed Effects Spatial error model\n", "Jacobian calculated using "))
switch(method, eigen = if (!quiet)
        cat("neighbourhood matrix eigenvalues\n"),
Matrix = {
        if (listw$style %in% c("W", "S") && !can.sim)
            stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C") && !(is.symmetric.glist(listw$neighbours,
            listw$weights)) && model=='lag')
            stop("Matrix method requires symmetric weights")
            if (listw$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw$neighbours,
            listw$weights)) && model=='error')
            stop("Matrix method requires symmetric weights")
        if (listw$style == "U" && model=='lag')
            stop("U style not permitted with lag model, use C")
        if (!quiet)
            cat("sparse matrix techniques using Matrix\n")
}, spam = {
        if (listw$style %in% c("W", "S") && !can.sim)
            stop("spam method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw$neighbours,
            listw$weights)))
            stop("spam method requires symmetric weights")
        if (!quiet)
            cat("sparse matrix techniques using spam\n")
    }, stop("...\nUnknown method\n"))


  ## check whether the panel is balanced
  balanced<-N*T==NT
  if(!balanced) stop("Estimation method unavailable for unbalanced panels")


  mt<-terms(formula,data=data)
  mf<-lm(formula,data,method="model.frame")#,na.action=na.fail
  na.act<-attr(mf,'na.action')
  cl<-match.call()

	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####takes the first n observations
	indic1<-seq(1,N)
	inde1<-rep(indic1,T) ####takes observations 1,  n+1, 2n+1...
  ### generates Wy if model=='lag'
if (model=='lag') wy<-unlist(tapply(y,inde, function(u) lag.listw(listw,u), simplify=TRUE))

#demeaning of the y and x variables depending both on model and effects

if (effects=="tpfe" | effects=="sptpfe"){
	ytms<-tapply(y,inde,mean) ####for each time period takes the mean for the cross section observations
	tpms<-function(q) tapply(q,inde,mean)
	xtms<-apply(x,2,tpms)   ###same thing for the X variable
	ytm<-rep(ytms,each=N) ##generate the NT variables
	xtm<-matrix(,NT,k)
	for (i in 1:k) xtm[,i]<-rep(xtms[,i],each=N)
	if (model=="lag") {
		wytms<-tapply(wy,inde,mean) ###same thing for Wy
		wytm<-rep(wytms,each=N)
}
	}

if (effects=="spfe" | effects=="sptpfe"){
	ysms<-tapply(y,inde1,mean) ###for each cross-sectional unit takes the mean over the time periods
	spms<-function(q) tapply(q,inde1,mean)
	xsms<-apply(x,2,spms)
	ysm<-rep(ysms,T)
	xsm<-matrix(,NT,k)
	for (i in 1:k) xsm[,i]<-rep(xsms[,i],T)
		if (model=="lag") {
			wysms<-tapply(wy,inde1,mean)
			wysm<-rep(wysms,T)
			}
	}
if (effects=='pooled'){
	yt<-y  	###keep the variables with no transformation
	xt<-x
	if (model=="lag") wyt<-wy
	}
if (effects=="tpfe"){ ####generate the demeaned variables for tpfe
	yt<-y-ytm
	xt<-x-xtm
	if (model=="lag") wyt<-wy-wytm
						}
if(effects=="spfe"){ ####generate the demeaned variables for spfe
	yt<-y-ysm
	xt<-x-xsm
	if (model=="lag") wyt<-wy-wysm
	 					}
if (effects=="sptpfe"){ ####generate the demeaned variables for both types of FE
	yt<-y - ysm - ytm + rep(mean(y),NT)
	xmm<-matrix(,NT,(k))
	for (i in 1:(k)) xmm[,i]<-rep(mean(x[,i]),NT)
	xt<-x - xsm - xtm + xmm
	if (model=="lag") wyt<- wy - wysm - wytm + rep(mean(wy),NT)							}
if 	(model=='error'){
	wyt<-unlist(tapply(yt,inde, function(u) lag.listw(listw,u), simplify=TRUE))
	dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw,TT), simplify=TRUE))
   wxt<-apply(xt,2,dm)
   colnames(wxt)<-paste('W',colnames(x))
	}
	  ## uses the two estimation functions: splaglm & sperrorlm
if(model=='lag'){
   RES <- splaglm(xt,yt,wyt,listw,method,K=k,NT=NT,T=T,inde,zero.policy = FALSE,quiet=quiet,
    interval = c(-1, 0.999), tol.solve = 1e-10, tol.opt = .Machine$double.eps^0.5,can.sim=can.sim)
    res.eff<-felag(y,x,wy,ysms,xsms,ytms, xtms, wytms, wysms, beta=RES$coeff,sige=RES$s2,yt,xt,N,T,NT,k,effects,method,rho=RES$rho, listw,inde)
    }
    if (model=='error'){
    RES<- sperrorlm(xt=xt,yt=yt,wyt=wyt,WX=wxt,listw,method,K=k,NT=NT,T=T,zero.policy = FALSE,quiet=quiet,
    interval = c(-1, 0.999), tol.solve = 1e-10, tol.opt = .Machine$double.eps^0.5,can.sim)
    	res.eff<-feerror(y,x,ysms,xsms,ytms, xtms, beta=RES$coeff,sige=RES$s2,yt,xt,N,T,NT,k,effects, lambda=RES$lambda)
    }
    ##calculate the R-squared
    yme <- y-mean(y)
    rsqr2 <- crossprod(yme)
    rsqr1 <- crossprod(res.eff[[1]]$res.e)
    res.R2<- 1- rsqr1/rsqr2

	#generate fixed values (from fixed_effect)
	y.hat <- res.eff[[1]]$xhat
	res <- as.numeric(res.eff[[1]]$res.e)
	N.vars<-res.eff$N.vars

	nam.rows <- dimnames(x)[[1]]
   names(y.hat) <- nam.rows
   names(res) <- nam.rows

	## make model data
   model.data <- data.frame(cbind(y,x))
   dimnames(model.data)[[1]] <- nam.rows

if (model=="lag")   spat.coef<-RES$rho
else spat.coef<-RES$lambda

if (is.null(RES$lambda.se) && model=="error") Coeff<-RES$coeff
else  Coeff<-c(spat.coef,RES$coeff)


type <- paste("fixed effects", model)

  spmod <- list(coefficients=Coeff, errcomp=NULL,
                vcov=RES$asyvar1,spat.coef=spat.coef,
                vcov.errcomp=NULL,
                residuals=res, fitted.values=y.hat,
                sigma2=RES$s2, type=type, model=model.data,
                call=cl, logLik=RES$ll, method=method, effects=effects, res.eff=res.eff)
  class(spmod) <- "splm"
  return(spmod)

	}
