`spsegm` <-
function(formula,data=list(), spec='default', listw, method='spatialsim',zero.policy = FALSE){
	
	##make sure that the data are not plm data
if (inherits(data,"plm.dim")) stop("Method not available for plm.data")

if(attributes(terms(formula,"intercept"))$intercept==0) stop("The intercept shoudl be controlled through the argument spec")

##generates model terms and model frames

	mt<-terms(formula,data=data)
	mf<-lm(formula,data,na.action=na.fail,method="model.frame")
	na.act<-attr(mf,'na.action')
	
#generates x and y 
# N.B. y is not a vector but a matrix	
	cl<-match.call()
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)
	xcolnames <- colnames(x)
if (ncol(y)<=1) stop("spsegm requires more than one dependent variable")
if (ncol(x)<=1) stop("spsegm requires more explanatory variables than simply the intercept")

	
#number of equations is equal to the number of columns of y
	eq<-ncol(y)
	k<-ncol(x)
	obs<-nrow(x)
	type<-'spsegm'
	
##check if the model is well specified
if (spec!="default" && !is.list(spec)) stop("spec should either be default or a list type of object")

if (is.list(spec) && length(spec) != eq) stop("The model has not been specified correctly: length of spec and number of columns of y should match")

##check listw and its consistency with the data
if (!inherits(listw,"listw")) stop("No neighborhood list")

if (obs != length(listw$neighbours)) stop("listw objects and variables have different dimension")

#generates the spatial lag of all y's
wy<-matrix(,obs,eq)
for(i in 1:eq){ 
wy[,i] <- lag.listw(listw,y[,i]) 
}
colnames(wy)<-paste("W",colnames(y))


###generates the instruments 
##it takes care of the argument spec 
#In fact, if it is equal to default x may or may not contain an intercept, whereas if
# spec is a list, x always contains an intercept
if (!is.list(spec)){
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	Wx <- matrix(nrow = obs, ncol = (k - (K - 1)))
        for (i in K:k) {
            wx <- lag.listw(listw, x[, i], zero.policy = zero.policy)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            Wx[, (i - (K - 1))] <- wx
        }
colnames(Wx) <- xcolnames[K:k]
WWx<-matrix(,obs,ncol(Wx))
for(i in  1 : ncol(Wx))  WWx[,i]<-lag.listw(listw,Wx[,i])
	}
else {
Wx<-matrix(,obs,k-1)
xnc<-matrix(x[,-1],obs,k-1)
for(i in  1 : (k-1)){ 
Wx[,i]<-lag.listw(listw,xnc[,i])
}
WWx<-matrix(,obs,k-1)
for(i in  1 : (k-1)){ 
WWx[,i]<-lag.listw(listw,Wx[,i])
}
colnames(Wx)<-paste("W",colnames(x)[-1])
colnames(WWx)<-paste("WW",colnames(x)[-1])
}

####ESTIMATION OF THE FIRST STEP
# If spec is not a list (DEFAULT) all the instruments are used in all the equations
# as well as all the wy and y are contained in all equations

if (!is.list(spec)){ 
	inst=cbind(Wx,WWx)
	r<-matrix(,obs,eq)
	b<-matrix(,k+eq+(eq-1),eq)
	for (i in 1:eq){
		yend<-cbind(wy,y[,-i])
		est<-tslssp(y[,i],yend,x,inst)
			r[,i]<-est$residuals
			b[,i]<-est$coefficients

					}
							}
							
# otherwise, one should select different x's in different equations
else { 
		inst=cbind(x,Wx,WWx)
		#print(inst)
		r<-matrix(,obs,eq)
		#lg<-array(,c(eq,1))
		b<-vector(mode="list",eq)
#for (i in 1:eq) lg[i,]<-length(spec[[i]])
		#dm<-max(lg)
for (i in 1:eq){
		yend<-cbind(wy,y[,-i])
		sel<-spec[[i]]
		est<-modtslssp(y[,i],yend,x[,sel],inst)
			r[,i]<-est$residuals
			b[[i]]<-est$coefficients
}
	}
	
#once the first step coefficients are obtained, one can proceed to perform the GM estimator
rho<-matrix(,eq,1)
sigma<-matrix(,eq,1)
for (i in 1:eq){
	mom<-Ggsararsp(u=r[,i],W=listw)
	pars<-c(0,0)
	estim <- optim(pars, searg, v = mom, verbose =FALSE)
	rho[i,]<-estim$par[1]
	sigma[i,]<-estim$par[2]
}

##this part generates a series of objects to be used in 
#the transformation of the data to perform the GS2SLS estimator
##the object lg1 serves to establish the dimension of Z and PZ when spec is not default
yt<-matrix(,obs,eq)
Wyt<-matrix(,obs,eq)
xt<-matrix(,obs,k)
r2<-matrix(,obs,eq)
b2<-vector(mode="list",eq)
Yt<-matrix(,obs*eq,1)
lg<-matrix(,eq,1)
lg1<-matrix(0,eq+1,1)
ntc<-matrix(0,eq,1)
sel<-matrix(0,eq,1)
if (is.list(spec)) for (i in 1:eq)   sel[i,]<-length(spec[[i]])
Sel<-cumsum(sel)
Sel1<-c(0,Sel)
lg1[1,]<-1

##generates Z=[X,Y,Wy] and PZ whose dimensions depend on the specification
if (!is.list(spec)) {
	Z<-Matrix(0,eq*obs,((eq-1)+eq+k)*eq)
	PZ<-Matrix(0,eq*obs,((eq-1)+eq+k)*eq)
		}
else {
	for (i in 1:eq) lg[i,]<- eq+length(spec[[i]]) +(eq-1)
	Z<-Matrix(0,eq*obs, sum(lg)) 
	PZ<-Matrix(0,eq*obs, sum(lg))
	}	

## generates the matrix of intruments and the projections
H<-cbind(x,Wx,WWx)
HH<-crossprod(H,H)
HHinv<-solve(HH)
Hp<-t(H)
P<-H%*%HHinv%*%Hp

##transform the y's and Wy's
for (i in 1:eq){
for (j in 1:k)	xt[,j] <- x[,j] - as.numeric(rho[i,]) * lag.listw(listw,x[,j]) 
for (t in 1:eq)	{
	yt[,t]<- y[,t] - as.numeric(rho[i,]) * lag.listw(listw,y[,t])
	Wyt[,t] <- lag.listw(listw,yt[,t])
}
	yendt<-cbind(Wyt,yt[,-i])
	
	##on the transformed data perform S2SLS
	#note that this part also generates some of the objects (Z,PZ and Yt) needed for the system S3SLS estimator 
	#Again the procedure is different whether spec is the default or not  
if (!is.list(spec)){
		inst=cbind(Wx,WWx)
		est<-tslssp(yt[,i],yendt,xt,inst)
		b2[[i]]<-est$coefficients
		r2[,i]<-est$residuals
		tmp1<-cbind(yt[,-i],xt,Wyt)
		A<-(obs*i)-obs+1
		B<- (obs*i)
		E<- i*ncol(tmp1)-ncol(tmp1)+1
		F<-i*ncol(tmp1)
		Z[A:B, E:F]<-tmp1
		PZ[A:B, E:F]<-P%*%tmp1
		Yt[A:B,]<-yt[,i]
	}
else{ 
	inst=cbind(x,Wx,WWx)
	    sel<-spec[[i]]
	    ntc[i,]<- Sel1[i] + (i-1)*(eq+(eq-1)) +1
	    lg1[i,]<- Sel1[i+1] + i*(eq+(eq-1))  
		 est<-modtslssp(yt[,i],yendt,xt[,sel],inst)
    	 b2[[i]]<-est$coefficients
	    r2[,i]<-est$residuals	
	    A<-(obs*i)-obs+1
		B<- (obs*i)
		E<-ntc[i,] 
		F<-lg1[i,]
		tmp1<-cbind(yt[,-i],xt[,sel],Wyt)
		Z[A:B, E:F]<-tmp1
		PZ[A:B, E:F]<-P%*%tmp1
		Yt[A:B,]<-yt[,i]
				}
}


###this last part performs the full information estimation: GS3SLS
##here there is an option of  using Matrix (which implies a )
SIGMA<-crossprod(r2)
SIGMAinv<-solve(SIGMA)

AZ <- matrix(,obs*eq,ncol(Z))
for(i in 1:ncol(Z)){
	tmp <- matrix(Z[,i],obs,eq)
tmpt <- t(tmp)
tmp2 <- SIGMAinv%*%tmpt
tmp3 <- matrix(t(tmp2),nrow=obs*eq,ncol=1)
AZ[,i]<-tmp3
}

AZH <- matrix(,obs*eq,ncol(PZ))
for(i in 1:ncol(PZ)){
	tmp <- matrix(PZ[,i],obs,eq)
tmpt <- t(tmp)
tmp2 <- SIGMAinv%*%tmpt
tmp3 <- matrix(t(tmp2),nrow=obs*eq,ncol=1)
AZH[,i]<-tmp3
}

Ytm <- matrix(Yt,obs,eq)
Ytmt <- t(Ytm)
Ay1 <- SIGMAinv%*%Ytmt
Ay <- matrix(t(Ay1),nrow=obs*eq,ncol=1)

ZHAZ<-crossprod(PZ,AZ)
ZHAZinv<-solve(ZHAZ)

ZAy<-crossprod(PZ,Ay)
delta<-ZHAZinv%*%ZAy
ZHAZH<-crossprod(PZ,AZH)
VC<-solve(ZHAZH)
model.data <- data.frame(cbind(y,x[,-1]))

spmod <- list(method=method, coefficients=delta, errcomp=NULL, vcov=VC, 
			  vcov.errcomp= NULL, residuals=NULL, fitted.values=NULL,
			  sigma2=NULL, type=type, model= model.data,  N=obs,
			  EQ=eq,K=k, call=cl,terms=mt, Xnames=colnames(x),Ynames=colnames(y), 
			  spec=spec)

class(spmod)<- "splm"
return(spmod)
}

