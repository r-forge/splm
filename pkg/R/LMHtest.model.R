`LMHtest.model` <-
function(x, listw, index){
## depends on listw2dgCMatrix.R

if(!inherits(x,"lm")) stop("argument should be an object of class lm")

  if(is.null(index))  stop("index should be specified to retrieve information on time and cross-sectional dimentions")

  if(!inherits(listw,"listw")) stop("object w should be of class listw")

  ind <- index[,1]
  tind <- index[,2]

###extract objects from x
  y<-model.response(x$model)
  e<-as.matrix(residuals(x))
  	ee<-crossprod(e)
   n<-dim(x$model)[1]
  	bOLS<-coefficients(x)
  	  	  form<-x$call
  x<-model.matrix(eval(x$call),x$model)
  #print(x)
	XpXi<-solve(crossprod(x))

   cl<-match.call()
  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  e<-e[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  ## det. number of groups and df
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  ## det. max. group numerosity
  T<-max(tapply(x[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)
#print(c(N,k,T,NT))
#   print(ols$model)
#	k<-dim(ols$model)[2]-1
	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,T)) ####indicator to get the time periods observations
####calculate the elements of LMj, LM1, SLM1

		JIe<-tapply(e,inde1,sum)
		JIe<-rep(JIe,T) ####calculates (J_T kronecker I_N)*u
		G<-(crossprod(e,JIe)/crossprod(e))-1 ###calculate G in LMj (same notation as in the paper)
tr<-function(R) sum(diag(R))
		LM1<-sqrt((NT/(2*(T-1))))*as.numeric(G) ###same notation as in Baltagi et al.


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
		b<-tr(W2+WW) ###generates b (same notation as the paper)
#		LMj<-(NT/(2*(T-1)))*as.numeric(G)^2 + ((N^2*T)/b)*as.numeric(H)^2   ###LMj as in the paper
		LM2<-sqrt((N^2*T)/b)*as.numeric(H)^2 ###same notation as in Baltagi et al.
if (LM1<=0){
		if (LM2<=0) JOINT<-0
		else JOINT<-LM2^2
		}		####this is chi-square_m in teh notation of the paper.
	else{
		if (LM2<=0) JOINT<-LM1^2
		else JOINT<-LM1^2 + LM2^2
		}
STAT<- qchisq(0.05,1,lower.tail=FALSE)
STAT1<- qchisq(0.05,2,lower.tail=FALSE)
if (JOINT>=2.952) {
		if (JOINT<7.289 & JOINT>=4.321) pval<-0.05
		if (JOINT >= 7.289) pval<-0.01
		if (JOINT<= 4.321)	pval<-0.1
	}
else pval<-1

	statistics<-JOINT

  names(statistics)="LM-H"
	method<- "Baltagi, Song and Koh LM-H one-sided joint test"
	#alt<-"serial corr. in error terms, sub RE and spatial dependence"
  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(form), alternative="Random Regional Effects and Spatial autocorrelation")
  class(RVAL) <- "htest"
  return(RVAL)
}

