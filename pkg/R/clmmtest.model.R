`clmmtest.model` <-
function(x, listw, index){
    ## depends on listw2dgCMatrix.R, spfeml.R

if(!inherits(x,"splm")) stop("argument should be an object of class splm")
frm<-x$call
if(x$type != "fixed effects error") stop("argument should be of type random effects ML")

  if(is.null(index))  stop("index should be specified to retrieve information on time and cross-sectional dimentions")

  if(!inherits(listw,"listw")) stop("object w should be of class listw")



  ind <- index[,1]
  tind <- index[,2]

if(names(x$coefficients)[1]=="(Intercept)")  X<-data.frame(cbind(rep(1,ncol(x$model)), x$model[,-1]))
else X<-x$model[,-1]
  y<-x$model[,1]
  	eML<-x$residuals

  ## reduce index accordingly

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  X<-X[oo,]
  y<-y[oo]

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

	lambda<-x$spat.coef
	#print(lambda)
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
               p.value = pval, data.name=deparse(frm), alternative="Random regional effects")
  class(RVAL) <- "htest"
  return(RVAL)

}

