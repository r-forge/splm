`clmmtest` <-
function(formula, data, index=NULL, listw){
    ## depends on listw2dgCMatrix.R, spfeml.R

	ml <- spfeml(formula, data=data, index=index, listw, model="error", effects="pooled")

	  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }
  index <- data[,1]
  tindex <- data[,2]
  X<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  ## reduce index accordingly
   names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(X))]
  tind<-tindex[which(names(index)%in%row.names(X))]

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  X<-X[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  ## det. number of groups and df
  N<-length(unique(ind))
  k<-dim(X)[[2]]
  ## det. max. group numerosity
  T<-max(tapply(X[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)

###maximum likelihood estimation under the null hypothesis that lambda is equal to zero. extract the residuals.
	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,T)) ####indicator to get the time periods observations

	lambda<-ml$spat.coef
	#print(lambda)
	eML<-residuals(ml)
#	print(length(eML))
 	Wst<-listw2dgCMatrix(listw) ###transform the listw object in a sparse matrix
	Ws<-t(Wst)  ### this is the real W since listw2dgCMatrix generate W'
	B<- -lambda*Ws
	diag(B)<- 1
	BpB<-crossprod(B)
	BpBi<- solve(BpB)
vc<-function(R) {
	BBu<-BpBi%*%R
	BBu<-as.matrix(BBu)
	}
	eme<-unlist(tapply(eML,inde,vc))
	sigmav2<-crossprod(eML,eme)/(N*(T-1)) ####estimate of the variance component sigma_v
	sigmav4<-sigmav2^2

tr<-function(R) sum(diag(R))
	trBpB<-tr(BpB)
	BpB2<-BpB%*%BpB
yybis<-function(q){ #### for very big dimension of the data this can be changed looping over the rows and columns of W or either the listw object
	wq<-rep(q,T)
	tmp<-wq%*%eML
					}
	BBu<-apply(BpB2,1,yybis)
	BBu<-rep(BBu,T)
	upBBu<-crossprod(eML,BBu)
	Dmu<--((T/(2*sigmav2))*trBpB) + ((1/(2*sigmav4))*upBBu)
	WpB<-Wst%*%B
	BpW<-t(B)%*%Ws
	WpBplBpW <-WpB + BpW
	G<-WpBplBpW %*% BpBi
	e<-tr(BpB2)
	d<-tr(WpBplBpW)
	h<-trBpB
	g<-tr(G)
	c<-tr(G%*%G)
	#print(c(e,d,h,g,c))
	NUM<- ((2*sigmav4)/T)*((N*sigmav4*c)-(sigmav4*g^2))   ###equation 2.30 in the paper
	DENft<- NT*sigmav4*e*c
	DENst<- N*sigmav4*d^2
	DENtt<- T*sigmav4*g^2 * e
	DENfot<- 2* sigmav4 *g*h*d
	DENfit<- sigmav4*h^2* c
	DEN<- DENft - DENst - DENtt + DENfot - DENfit
	LMmu <- Dmu^2*NUM / DEN
	LMmustar<- sqrt(LMmu)
	statistics<-LMmustar
  pval <- pnorm(LMmustar, lower.tail=FALSE)

  names(statistics)="LM*-mu"
	method<- "Baltagi, Song and Koh LM*- mu conditional LM test (assuming lambda may or may not be = 0)"
	#alt<-"serial corr. in error terms, sub RE and spatial dependence"
  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Random regional effects")
  class(RVAL) <- "htest"
  return(RVAL)

}

