`tw` <-
function(W,N){
## depends on listw2dgCMatrix.R
	Wst<-listw2dgCMatrix(W)
	Ws<-t(Wst)
	WpW<-crossprod(Ws)
	WpWWpW<-WpW%*%WpW
	WppW<-Wst+Ws
	WpWWppW<-WpW%*%WppW
	WW<-Ws%*%Ws
	WWpWpW<-WW+WpW
	tr1<-sum(diag(WpW/N))
	tr2<-sum(diag(WpWWpW/N))
	tr3<-sum(diag(WpWWppW/N))
	tr4<-sum(diag(WWpWpW/N))
	TW1c<-rbind(2,2*tr1, 0)
	TW2c<-rbind(2*tr1, 2*tr2,tr3)
	TW3c<-rbind(0, tr3,tr4)
	TW<-cbind(TW1c,TW2c,TW3c)
	out<-list(TW=matrix(TW,3,3))
}

