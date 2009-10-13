`slm1test.model` <-
function(x, listw, index){
### depends on listw2dgCMatrix.R

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
		G<-(crossprod(e,JIe)/ee)-1 ###calculate G in LMj (same notation as in the paper)
tr<-function(R) sum(diag(R))
		#LMj<-(NT/(2*(T-1)))*as.numeric(G)^2 + ((N^2*T)/b)*as.numeric(H)^2   ###LMj as in the paper
		LM1<-sqrt((NT/(2*(T-1))))*as.numeric(G) ###same notation as in Baltagi et al.
		s<-NT-k ###needed for SLM1
		B<-XpXi%*%t(x)   ## I call B XpXi %*% Xp
fun<-function(Q) tapply(Q,inde1,sum) ############################################
		JIx<-apply(x,2,fun)
		JIX<-matrix(,NT,k)
for (i in 1:k) JIX[,i]<-rep(JIx[,i],T) ##ALL this lines are needed to calculate the trace of D1M. You should work with the expression of D1and M=I - X(X'X)^{-1}X' to avoid the calculation of the trace   of a matrix of dimension NTxNT. Working with the properties of the trace (tr(A+B) = tr(A) + tr(B) and tr(ABC) = tr(BCA) = tr(CAB)). for more details see "NOTE ON THE TRACE.R"
		di<-numeric(NT)
		XpJIX<-crossprod(x,JIX)
		d1<-NT-tr(XpJIX%*%XpXi) ############################################
#		print(d1)
		Ed1<-d1/s ####expected value of d1 (same notation as in the paper)
		di2<-numeric(NT)
		JIJIx<-apply(JIX,2,fun)
		JIJIX<-matrix(,NT,k)
for (i in 1:k) JIJIX[,i]<-rep(JIJIx[,i],T)
		JIJIxxpx<-JIJIX%*%XpXi
		di1<- crossprod(x, JIJIxxpx)
		tr1<-tr(di1)
		XpIJX<-crossprod(x,JIX)
		fp<-XpIJX%*%B
		sp<-JIX%*%XpXi
		tr3<-tr(fp%*%sp)
#		print(tr3)
		fintr<-NT*T-2*tr1+tr3 #####trace of (D1M)^2 again see "NOTE ON THE TRACE.R"
		Vd1<-2*(s*fintr - (tr1^2)/s^2*(s+2)) ##this is therefore the variance of D1

SLM1<-((G+1)- Ed1)/sqrt(Vd1) ###and this is the expression for SLM1
#print(c(G,Ed1,Vd1,SLM1))
STAT2<- qnorm(0.95,lower.tail=TRUE)


	statistics<-SLM1
  pval <- pnorm(SLM1, lower.tail=FALSE)

  names(statistics)="SLM1"
	method<- "Baltagi, Song and Koh SLM1 marginal test"
	#alt<-"serial corr. in error terms, sub RE and spatial dependence"
  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval,  data.name= deparse(form),
               alternative="Random Regional Effects")
  class(RVAL) <- "htest"
  return(RVAL)
}

