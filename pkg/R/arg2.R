`arg2` <-
function (rhopar, v, ss,SS,T,TW, verbose = verbose) 
{
	Ga<-cbind( (1/(T-1))*ss^2,0)
	Gb<-cbind( 0, SS^2)
	Gc<-rbind(Ga,Gb)
	Gamma<-kronecker(Gc,TW) 
	Gammainv<-solve(Gamma)
    vv <-  v$GG %*% c(rhopar[1], rhopar[1]^2, rhopar[2], rhopar[3]) - v$gg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}



`arg3` <-
function (rhopar, v, ss,T,TW, verbose = verbose) 
{
	Ga<-(1/(T-1))*ss^2
	Gamma<-as.numeric(Ga)*TW 
	Gammainv<-solve(Gamma)
    vv <-  v$bigG %*% c(rhopar[1], rhopar[1]^2, rhopar[2]) - v$smallg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}

