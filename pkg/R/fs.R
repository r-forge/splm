`fs` <-
function(listw,u,N,T){
	ind<-seq(1,T)
	inde<-rep(ind,each=N)
	NT<-N*T
#ub<-matrix(,NT,1)
#for (i in 1:T) ub[(N*i-N+1):(N*i),]<-lag.listw(listw,u[(N*i-N+1):(N*i)])
#ubb<-matrix(,NT,1)
#for (i in 1:T) ubb[(N*i-N+1):(N*i),]<-lag.listw(listw,ub[(N*i-N+1):(N*i)])
#ubmt<-matrix(,N,1)
#for (i in 1:N) ubmt[i,]<-mean(ub[seq.int(from=i,to=NT,by=N)])
#ubbmt<-matrix(,N,1)		
#for (i in 1:N) ubbmt[i,]<-mean(ubb[seq.int(from=i,to=NT,by=N)])
#umt<-matrix(,N,1)		
#for (i in 1:N) umt[i,]<-mean(u[seq.int(from=i,to=NT,by=N)])
	ub<-unlist(tapply(u,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))###Wu
	ubb<-unlist(tapply(ub,inde, function(TT) lag.listw(listw,TT), simplify=TRUE))###WWu
	ind1<-seq(1,N)
	inde1<-rep(ind1,T)
	umt<-tapply(u, inde1, mean) 
	ubmt<-tapply(ub, inde1, mean) ####mean over time 
	ubbmt<-tapply(ubb, inde1, mean)
	umNT<-rep(umt,T)
	ubmNT<-rep(ubmt,T)
	ubbmNT<-rep(ubbmt,T)
	Q0ub<-ub-ubmNT
	Q0ubb<-ubb-ubbmNT
	Q0u<-u-umNT
	uQ0ub<-crossprod(u,Q0ub)
	ubQ0ubb<-crossprod(ub,Q0ubb)
	ubbQ0ubb<-crossprod(ubb,Q0ubb)
	ubQ0ub<-crossprod(ub,Q0ub)
	uQ0u<-crossprod(u,Q0u)
	ubbQ0ub<-crossprod(ubb,Q0ub)
	uQ0ubb<-crossprod(u,Q0ubb)
	tr <- matrix(0, N, 1)
	for (i in 1:N) {
        tr[i] <- sum(listw$weights[[i]]^2)
	}      
	trwpw <- sum(tr)
	G1c<-(1/(N*(T-1)))*rbind(2*uQ0ub, 2*ubbQ0ub,(uQ0ubb+ubQ0ub) )
	G2c<- (-1/(N*(T-1)))* rbind(ubQ0ub,ubbQ0ubb,ubQ0ubb)
	G3c<- rbind(1,trwpw/N, 0)
	G<-cbind(G1c,G2c,G3c)	
	g<-(1/(N*(T-1)))*rbind(uQ0u,ubQ0ub,uQ0ub)
	output<-list(bigG=G, smallg=g, Q1u=umNT,Q1ub=ubmNT,Q1ubb=ubbmNT, ub=ub,ubb=ubb,TR=trwpw)
}

