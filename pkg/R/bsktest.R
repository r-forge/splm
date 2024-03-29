`bsktest` <-
function(x,...){
  UseMethod("bsktest")
}

## define necessary function pchibar() (courtesy Jason Sinnwell)
## eliminates dependency on 'ibdreg'

pchibar <- function(x, df, wt){
  # compute P(X <= x), where P is the chi-bar distribution
  # df = vector of df
  # wt = vector of mixing proportions
  if(x<=0){
      return(0)
  }
  zed <- df==0
  cdf <- ifelse(any(zed), wt[zed], 0)
  cdf <- cdf + sum(pchisq(x, df[!zed])*wt[!zed])
  return(cdf)
}


#`bsktest.splm` <-
#function(x, listw, index=NULL, test=c("CLMlambda","CLMmu"), ...){
#	switch(match.arg(test), CLMlambda = {
#    bsk = clmltest.model(x,listw, index, ...)
#  }, CLMmu = {
#    bsk = clmmtest.model(x,listw, index, ... )
#  })
#  return(bsk)
#}

`bsktest.formula` <-
function(x, data, index=NULL, listw,
         test=c("LMH","LM1","LM2","CLMlambda","CLMmu"),
         standardize=FALSE, method = "eigen", ...){

  ## transform data if needed
  if(!("pdata.frame" %in% class(data))) {
    data <- pdata.frame(data, index)
  }

  
switch(match.arg(test), LM1 = {
    
    bsk = slm1test(x, data, index,  listw, standardize, ...)

  }, LM2 = {

    bsk = slm2test(x, data, index,  listw, standardize, ...)

  }, LMH = {

    bsk = LMHtest(x, data, index,  listw, ...)

  }, CLMlambda = {

    bsk = clmltest(x, data, index,  listw, ...)

  }, CLMmu = {

    bsk = clmmtest(x, data, index,  listw, method = method,...)

  })

  return(bsk)

}







`slm1test` <-
function(formula, data, index=NULL, listw, standardize, ...){
    ## temporary switch because of bad results in SLM1
    # if(standardize) stop("Standardized SLM1 test temporarily unavailable:
    # \n use 'standardize=FALSE' for LM1 test instead")
    
  #if(!is.null(index)) { ####can be deleted when using the wrapper
    #require(plm)
  #  data <- plm.data(data, index)
  #  }

tr<-function(R) sum(diag(R))
fun<-function(Q) tapply(Q,inde1,sum)

  index <- data[,1]
  tindex <- data[,2]

  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
   cl<-match.call()
  
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  time <-max(tapply(x[,1],ind,length))
  NT<-length(ind)
	
	ols<-lm(y~x-1)
	XpXi<-solve(crossprod(x))
    n<-dim(ols$model)[1]

	indic<-seq(1,time)
	inde<-as.numeric(rep(indic, each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,time)) ####indicator to get the time periods observations

	bOLS<-coefficients(ols)
	e<-as.matrix(residuals(ols))
	ee<-crossprod(e)



		JIe<-tapply(e,inde1,sum)
		JIe<-rep(JIe,time) 
		G <- (crossprod(e,JIe)/ee)-1 


		LM1<-sqrt((NT/(2*(time-1))))*as.numeric(G) 
		
		s<-NT-k 
		B<-XpXi%*%t(x)   
		
 
		JIx<-apply(x,2,fun)
		JIX<-matrix(,NT,k)
for (i in 1:k) JIX[,i]<-rep(JIx[,i], time) ## "NOTE ON THE TRACE.R"
		di <- as.numeric(NT)
		XpJIX<-crossprod(x,JIX)
		
		
		
		d1<- NT-tr(XpJIX%*%XpXi) 
		Ed1 <- d1/s 

		di2<-numeric(NT)
		JIJIx<-apply(JIX,2,fun)
		JIJIX<-matrix(,NT,k)
for (i in 1:k) JIJIX[,i]<-rep(JIJIx[,i],time)
		JIJIxxpx<-JIJIX%*%XpXi
		di1<- crossprod(x, JIJIxxpx)
		tr1<-tr(di1)
		XpIJX<-crossprod(x,JIX)
		fp<-XpIJX%*%B
		sp<-JIX%*%XpXi
		tr3<-tr(fp%*%sp)
		fintr<-NT*time-2*tr1+tr3 
		Vd1 <- 2*(s*fintr - (d1^2))/(s^2*(s+2))

SLM1 <- ((G+1)- Ed1)/sqrt(Vd1) 

  statistics <- if(standardize) SLM1 else LM1
  pval <- 2*pnorm(statistics, lower.tail=FALSE)
	
  names(statistics) <- if(standardize) "SLM1" else "LM1"
	method<- "Baltagi, Song and Koh SLM1 marginal test"
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Random effects")
  class(RVAL) <- "htest"
  return(RVAL)
}


`slm2test` <-
function(formula, data, index=NULL, listw, standardize, ...){
    # if(standardize) stop("Standardized SLM2 test temporarily unavailable: \n use 'standardize=FALSE' for LM2 test instead")                  


  #if(!is.null(index)) { 
    #require(plm)
  #  data <- plm.data(data, index)
  #  }

  index <- data[,1]
  tindex <- data[,2]

  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
   cl<-match.call()
	names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  N<-length(unique(ind))
  k<-dim(x)[[2]]
  time <-max(tapply(x[,1],ind,length))
  NT<-length(ind)
  ols<-lm(y~x-1)
	XpXi<-solve(crossprod(x))
   n<-dim(ols$model)[1]

	indic<-seq(1,time)
	inde<-as.numeric(rep(indic,each=N)) 
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,time)) 
	
	bOLS<-coefficients(ols)
	e<-as.matrix(residuals(ols))
	ee<-crossprod(e)


		Ws<-listw2dgCMatrix(listw) 
		Wst<-t(Ws)  
		WWp<-(Ws+Wst)/2 

yy<-function(q){ 
	wq<-WWp%*%q
	wq<-as.matrix(wq)
	}

		IWWpe<-unlist(tapply(e,inde,yy)) 
		H<-crossprod(e,IWWpe)/crossprod(e) 
		W2<-Ws%*%Ws 
		WW<-crossprod(Ws) 
    tr<-function(R) sum(diag(R))
	b<-tr(W2+WW) 
		LM2<-sqrt((N^2*time)/b)*as.numeric(H)
		s<-NT-k
lag<-function(QQ)lag.listw(listw,QQ)
fun2<-function(Q) unlist(tapply(Q,inde,lag))
	Wx<-apply(x,2,fun2)
	WX<-matrix(Wx,NT,k)
	XpWx<-crossprod(x,WX)
	D2M<-XpWx%*%XpXi 
	Ed2<- (time*sum(diag(Ws)) - tr(D2M))/s 

	WWx<-apply(WX,2,fun2)
	WWX<-matrix(WWx,NT,k)
	XpWWX<-crossprod(x,WWX)				
	spb<-XpWWX%*%XpXi
	spbb<-tr(spb)
	tpb<-XpWx%*%XpXi%*%XpWx%*%XpXi
	fintr2<-time*tr(W2) - 2* spbb + tr(tpb)
	Vd2<-2*(s*fintr2 - (sum(diag(D2M))^2))/(s^2*(s+2))
	We<-unlist(tapply(e,inde,function(W) lag.listw(listw,W)))
	d2<-crossprod(e,We)/ee
	
	SLM2<- (d2-Ed2)/sqrt(Vd2) 
  
  statistics <- if(standardize) SLM2 else LM2
  pval <- 2*pnorm(statistics, lower.tail=FALSE)
	
  names(statistics) <- if(standardize) "SLM2" else "LM2"
  
	method<- "Baltagi, Song and Koh LM2 marginal test"
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Spatial autocorrelation")
  class(RVAL) <- "htest"
  return(RVAL)
}


`LMHtest` <-
function(formula, data, index=NULL, listw, ...){
    ## depends on listw2dgCMatrix.R
  #require(ibdreg) # for mixed chisquare distribution
  # now imported

  #if(!is.null(index)) { ####can be deleted when using the wrapper
    #require(plm)
  #  data <- plm.data(data, index)
  #  }

  index <- data[,1]
  tindex <- data[,2]

  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  cl<-match.call()
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  N<-length(unique(ind))
  k<-dim(x)[[2]]
  time <-max(tapply(x[,1],ind,length))
  NT<-length(ind)
  ols<-lm(y~x-1)
	XpXi<-solve(crossprod(x))
   n<-dim(ols$model)[1]

	indic<-seq(1,time)
	inde<-as.numeric(rep(indic,each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,time)) ####indicator to get the time periods observations
	bOLS<-coefficients(ols)
	e<-as.matrix(residuals(ols))
	ee<-crossprod(e)


		JIe<-tapply(e,inde1,sum)
		JIe<-rep(JIe,time) 
		G<-(crossprod(e,JIe)/ee)-1 
      tr<-function(R) sum(diag(R))

		LM1<-sqrt((NT/(2*(time-1))))*as.numeric(G) 


####calculate the elements of LMj, LM1, SLM1
		Ws<-listw2dgCMatrix(listw) 
		Wst<-t(Ws)  
		WWp<-(Ws+Wst)/2 

yy<-function(q){
	wq<-WWp%*%q
	wq<-as.matrix(wq)
	}
	
		IWWpe<-unlist(tapply(e,inde,yy)) 
		H<-crossprod(e,IWWpe)/crossprod(e) 
		W2<-Ws%*%Ws 
		WW<-crossprod(Ws) 
      tr<-function(R) sum(diag(R))
	   b<-tr(W2+WW) 
		LM2<-sqrt((N^2*time)/b)*as.numeric(H)


if (LM1<=0){
		if (LM2<=0) JOINT<-0
		else JOINT<-LM2^2
		}		####this is chi-square_m in teh notation of the paper.

	else{
		if (LM2<=0) JOINT<-LM1^2
		else JOINT<-LM1^2 + LM2^2
		}
  statistics<-JOINT
  pval <- 1 - pchibar(statistics, df=0:2, wt=c(0.25,0.5,0.25))
    

  names(statistics)="LM-H"
	method<- "Baltagi, Song and Koh LM-H one-sided joint test"

  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Random Regional Effects and Spatial autocorrelation")
  class(RVAL) <- "htest"
  return(RVAL)
}


`clmmtest` <-
function(formula, data, index=NULL, listw, method, ...){

##ml <- spfeml(formula=formula, data=data, index=index, listw=listw, model="error", effects="pooling", method = method)

    ## modified for obtaining correct residuals
    ## (spfeml omitted the intercept)
    ml <- spreml(formula=formula, data=data, index=index, w=listw,
                 errors="sem", lag=F)
                                                        
    #if(!is.null(index)) {
    # data <- plm.data(data, index)
    #}

  index <- data[,1]
  tindex <- data[,2]
  X<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(X))]
  tind<-tindex[which(names(index)%in%row.names(X))]
  oo<-order(tind,ind)
  X<-X[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]
  N<-length(unique(ind))
  k<-dim(X)[[2]]
  time <-max(tapply(X[,1],ind,length))
  NT<-length(ind)


	indic<-seq(1,time)
	inde<-as.numeric(rep(indic,each=N))
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,time))

	lambda <- ml$errcomp["rho"] # was: ml$spat.coef,
                                    # modified for use with spreml
        eML<-residuals(ml)

 	Ws<-listw2dgCMatrix(listw)
	Wst<-t(Ws)
	B <- Diagonal(N)-lambda*Ws


	BpB<-crossprod(B)
	BpB2 <- BpB %*% BpB
	BpBi<- solve(BpB)
tr<-function(R) sum(diag(R))
	trBpB<-tr(BpB)

vc<-function(R) {
	BBu<-BpB %*% R
	BBu<-as.matrix(BBu)
	}

	eme<-unlist(tapply(eML,inde,vc))

	sigmav2<-crossprod(eML,eme)/(N*time)
	sigmav4<-sigmav2^2


yybis<-function(q){
	wq<-rep(q,time)
	tmp<-wq%*%eML
					}

	BBu<-apply(BpB2,1,yybis)
	BBu<-rep(BBu,time)
	upBBu<-crossprod(eML,BBu)

	Dmu<- -(time/(2*sigmav2))*trBpB + (1/(2*sigmav4))*upBBu

	WpB<-Wst%*%B
	BpW<-crossprod(B, Ws)
	WpBplBpW <-WpB + BpW
	bigG<-WpBplBpW %*% BpBi

	smalle<-tr(BpB2)
	smalld<-tr(WpBplBpW)
	smallh<-trBpB
	smallg<-tr(bigG)
	smallc<-tr(bigG%*%bigG)

	NUM<- ((2 * sigmav4)/time) * ((N*sigmav4*smallc)-(sigmav4*smallg^2))   ###equation 2.30 in the paper
	DENft<- NT*sigmav4* smalle * smallc
	DENst<- N*sigmav4* smalld^2
	DENtt<- time*sigmav4* smallg^2 * smalle
	DENfot<- 2* sigmav4 *smallg * smallh* smalld
	DENfit<- sigmav4 * smallh^2* smallc
	DEN<- DENft - DENst - DENtt + DENfot - DENfit
	LMmu <- Dmu^2*NUM / DEN

	LMmustar<- sqrt(LMmu)

  statistics<-LMmustar
  pval <- 2*pnorm(LMmustar, lower.tail=FALSE)

  names(statistics)="LM*-mu"
	method<- "Baltagi, Song and Koh LM*- mu conditional LM test (assuming lambda may or may not be = 0)"
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Random regional effects")
  class(RVAL) <- "htest"
  return(RVAL)

}


clmltest <- function (formula, data, index = NULL, listw)
{
   # ml <- spreml(formula, data = data, w = listw2mat(listw),
   #     errors = "re")
    #if (!is.null(index)) {
        #require(plm)
    #    data <- plm.data(data, index)
    #}
    index <- data[, 1]
    tindex <- data[, 2]
    data$tindex <- tindex

	#ml <- lme(formula, data, random=~1|tindex)

    ## modified to get rid of lme problems with "inner_perc_table" 
    ml <- spreml(formula=formula, data=data, index=index, w=listw,
                 errors="re", lag=F)
    
    X <- model.matrix(formula, data = data)
    y <- model.response(model.frame(formula, data = data))
    names(index) <- row.names(data)
    ind <- index[which(names(index) %in% row.names(X))]
    tind <- tindex[which(names(index) %in% row.names(X))]
    oo <- order(tind, ind)
    X <- X[oo, ]
    y <- y[oo]
    ind <- ind[oo]
    tind <- tind[oo]
    N <- length(unique(ind))
    k <- dim(X)[[2]]
    time <- max(tapply(X[, 1], ind, length))
    NT <- length(ind)
    eML <- residuals(ml)
    indic <- seq(1, time)
    inde <- as.numeric(rep(indic, each = N))
    ind1 <- seq(1, N)
    inde1 <- as.numeric(rep(ind1, time))
    eme <- tapply(eML, inde1, mean)
    emme <- eML - rep(eme, time)
    sigmav <- crossprod(eML, emme)/(N * (time - 1))
    sigma1 <- crossprod(eML, rep(eme, time))/N
    c1 <- sigmav/sigma1^2
    c2 <- 1/sigmav
    c1e <- as.numeric(c1) * eML
    Wst <- listw2dgCMatrix(listw)
    Ws <- t(Wst)
    WpsW <- Wst + Ws
    yybis <- function(q) {
        wq <- (WpsW) %*% q
        wq <- as.matrix(wq)
    }
    Wc1e <- unlist(tapply(eML, inde, yybis))
    sumWc1e <- tapply(Wc1e, inde1, sum)
    prod1 <- as.numeric(c1) * rep(sumWc1e, time)/time
    prod2 <- as.numeric(c2) * (Wc1e - rep(sumWc1e, time)/time)
    prod <- prod1 + prod2
    D <- 1/2 * crossprod(eML, prod)
    W2 <- Ws %*% Ws
    WW <- crossprod(Ws)
    tr <- function(R) sum(diag(R))
    b <- tr(W2 + WW)
    LMl1 <- D^2/(((time - 1) + as.numeric(sigmav)^2/as.numeric(sigma1)^2) *
        b)
    LMlstar <- sqrt(LMl1)
    statistics <- LMlstar
    pval <- 2*pnorm(LMlstar, lower.tail = FALSE)
    names(statistics) = "LM*-lambda"
    method <- "Baltagi, Song and Koh LM*-lambda conditional LM test (assuming sigma^2_mu >= 0)"
    dname <- deparse(formula)
    RVAL <- list(statistic = statistics, method = method, p.value = pval,
        data.name = deparse(formula), alternative = "Spatial autocorrelation")
    class(RVAL) <- "htest"
    return(RVAL)
}
