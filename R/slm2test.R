`slm2test` <-
function(formula, data, index=NULL, listw){
    ## depends on listw2dgCMatrix.R
  if(!is.null(index)) { ####can be deleted when using the wrapper
    require(plm)
    data <- plm.data(data, index)
    }

  index <- data[,1]
  tindex <- data[,2]

  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
   cl<-match.call()
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
	ols<-lm(y~x)
	XpXi<-solve(crossprod(x))
   n<-dim(ols$model)[1]

	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,T)) ####indicator to get the time periods observations
	bOLS<-coefficients(ols)
	e<-as.matrix(residuals(ols))
	ee<-crossprod(e)
####calculate the elements of LMj, LM1, SLM1
		Wst<-listw2dgCMatrix(listw) ###transform the listw object in a sparse matrix
		Ws<-t(Wst)  ### this is the real W since listw2dgCMatrix generate W'
		WWp<-(Ws+Wst)/2 ##generate (W+W')/2
yy<-function(q){ #### for very big dimension of the data this can be changed looping over the rows and columns of W or either the listw object
	wq<-WWp%*%q
	wq<-as.matrix(wq)
	}
		IWWpe<-unlist(tapply(e,inde,yy)) ####calculates (I_T kronecker (W+W')/2)*u
		H<-crossprod(e,IWWpe)/crossprod(e) #calculate H (same notation as in the paper)
		W2<-Ws%*%Ws ####generate W^2
		WW<-crossprod(Ws) ####generate W'*W
tr<-function(R) sum(diag(R))
		b<-tr(W2+WW) ###generates b (same notation as the paper)
#		LMj<-(NT/(2*(T-1)))*as.numeric(G)^2 + ((N^2*T)/b)*as.numeric(H)^2   ###LMj as in the paper
		LM2<-sqrt((N^2*T)/b)*as.numeric(H)^2 ###same notation as in Baltagi et al.
		s<-NT-k
lag<-function(QQ)lag.listw(listw,QQ)############################################
fun2<-function(Q) tapply(Q,inde,lag)
	Wx<-apply(x,2,fun2)
	WX<-matrix(unlist(Wx),NT,k)
#test<-matrix(,NT,k)   CAN BE DONE ALSO LIKE THIS. HOWEVER FOR BIG DATASET (WITH LARGE T, IT IS LIKELY TO BE SLOWER)
#for(i in 1:T) test[(N*i-N+1):(N*i),]<-lag.listw(listw,x[(N*i-N+1):(N*i),])
#all.equal(test,WX)
	XpWx<-crossprod(x,WX)##this part calculates the trace of D2M where D2=I_T kronecker W
	D2M<-XpWx%*%XpXi ############################################
	Ed2<- -tr(D2M)/s
	WWx<-apply(WX,2,fun2)
	WWX<-matrix(unlist(WWx),NT,k)
	XpWWX<-crossprod(x,WWX)				####This part calculates the tr(D2M)^2
	spb<-XpWWX%*%XpXi
	spbb<-tr(spb)
	tpb<-XpWx%*%XpXi%*%XpWx%*%XpXi
	fintr2<-T*tr(W2) - 2* spbb + tr(tpb) ############################################
	Vd2<-2*(s*fintr2 - (sum(diag(D2M))^2)/s^2*(s+2)) ####this is the variance of d2
	We<-unlist(tapply(e,inde,function(W) lag.listw(listw,W)))
	d2<-crossprod(e,We)/ee
	SLM2<- (d2-Ed2)/sqrt(Vd2) ###this is the expression for SLM2

STAT2<- qnorm(0.95,lower.tail=TRUE)
	statistics<-SLM2
  pval <- pnorm(SLM2, lower.tail=FALSE)

  names(statistics)="SLM2"
	method<- "Baltagi, Song and Koh SLM2 marginal test"
	#alt<-"serial corr. in error terms, sub RE and spatial dependence"
  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Spatial autocorrelation")
  class(RVAL) <- "htest"
  return(RVAL)
}

