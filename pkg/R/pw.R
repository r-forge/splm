`pw` <-
function(bigG, smallg, Q1u,Q1ub,Q1ubb, u, ub,ubb,N, TR){
	uQ1u<-crossprod(u,Q1u)
	uQ1ub<-crossprod(u,Q1ub)
	ubbQ1ub<-crossprod(ubb,Q1ub)
	ubbQ1ubb<-crossprod(ubb,Q1ubb)
	uQ1ubb<-crossprod(u,Q1ubb)
	ubQ1ub<-crossprod(ub,Q1ub)
	ubQ1ubb<-crossprod(ub,Q1ubb)
	G1c1<-rbind(2*uQ1ub, 2*ubbQ1ub,  (uQ1ubb + ubQ1ub))/N
	G1c2<-rbind(ubQ1ub, ubbQ1ubb, ubQ1ubb)/(-N)
	G1c3<-rbind(1,TR/N,0)
	G1c<-cbind(G1c1,G1c2,rep(0,3),G1c3)
	g1<-rbind(uQ1u, ubQ1ub, uQ1ub)/N
	GG<-rbind(cbind(bigG,rep(0,3)),G1c)	
	gg<-rbind(smallg,g1)
	out<-list(GG=GG,gg=gg)
}

