`arg1` <-
function (rhopar, v, ss,SS,T, verbose = verbose) 
{
	Ga<-cbind( (1/(T-1))*ss^2,0)
	Gb<-cbind( 0, SS^2)
	Gc<-rbind(Ga,Gb)
	Gamma<-kronecker(Gc,diag(3)) 
	Gammainv<-solve(Gamma)
    vv <-  v$GG %*% c(rhopar[1], rhopar[1]^2, rhopar[2], rhopar[3]) - v$gg
    value <- t(vv)%*% Gammainv %*% vv
    if (verbose) 
	cat("function:", value, "rho:", rhopar[1], "sig2:", 
		rhopar[2], "\n")
    value
}



