`clmltest` <-
function(formula, data, index=NULL, listw){
    ## depends on listw2dgCMatrix.R, REmod.R, spreml.R

	ml <- spreml(formula, data=data, w=listw2mat(listw), errors="re")
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
	eML<-residuals(ml)
###maximum likelihood estimation under the null hypothesis that lambda is equal to zero. extract the residuals.
	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####indicator to get the cross-sectional observations
	ind1<-seq(1,N)
	inde1<-as.numeric(rep(ind1,T)) ####indicator to get the time periods observations

	eme<-tapply(eML,inde1,mean)
	emme<-eML - rep(eme,T)
	sigmav<-crossprod(eML,emme)/(N*(T-1)) ####estimate of the variance component sigma_v
	sigma1<-crossprod(eML,rep(eme,T))/N ####estimate of the variance component sigma_1
	c1<-sigmav/sigma1^2
	c2<-1/sigmav
	c1e<-as.numeric(c1)*eML
	Wst<-listw2dgCMatrix(listw) ###transform the listw object in a sparse matrix
	Ws<-t(Wst)  ### this is the real W since listw2dgCMatrix generate W'
	WpsW<-Wst+Ws
yybis<-function(q){ #### for very big dimension of the data this can be changed looping over the rows and columns of W or either the listw object
	wq<-(WpsW)%*%q
	wq<-as.matrix(wq)
	}
	Wc1e<-unlist(tapply(eML,inde,yybis)) #### (W'+W)*u
	sumWc1e<-tapply(Wc1e,inde1,sum)
	prod1<-as.numeric(c1)*rep(sumWc1e,T)/T
	prod2<-as.numeric(c2)* (Wc1e - rep(sumWc1e,T)/T)
	prod<-prod1+prod2
	D<-1/2*crossprod(eML,prod) ###calculates D in the notation of the paper.
	W2<-Ws%*%Ws ####generate W^2
	WW<-crossprod(Ws) ####generate W'*W
tr<-function(R) sum(diag(R))
	b<-tr(W2+WW) ###generates b (same notation as the paper)
	LMl1<-D^2/(((T-1)+as.numeric(sigmav)^2/as.numeric(sigma1)^2)*b) ###conditional LM test for lambda equal to zero
	LMlstar<-sqrt(LMl1) ###one-sided version
	statistics<-LMlstar
  pval <- pnorm(LMlstar, lower.tail=FALSE)

  names(statistics)="LM*-lambda"
	method<- "Baltagi, Song and Koh LM*-lambda conditional LM test (assuming sigma^2_mu >= 0)"
	#alt<-"serial corr. in error terms, sub RE and spatial dependence"
  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = statistics,
               method = method,
               p.value = pval, data.name=deparse(formula), alternative="Spatial autocorrelation")
  class(RVAL) <- "htest"
  return(RVAL)

}

