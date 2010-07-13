`write.effects.splm` <-
function (x, file = "effects", ncolumns = if (is.character(x)) 1 else 5, 
    append = FALSE, sep = ","){
    	object<-x
	effects<-object$effects
if(effects=="tpfe"){
out<-rbind(object$INTTable,object$TETable)
write(out, file=file, ncolumns=ncolumns, append=append, sep=sep)

}
	
if(effects=="spfe"){
out<-rbind(object$INTTable,object$SETable)
write(out, file=file, ncolumns=ncolumns, append=append, sep=sep)
}
if(effects=="sptpfe"){
out<-rbind(object$INTTable,object$SETable,object$TETable)
write(out, file=file, ncolumns=ncolumns, append=append, sep=sep)
}
	}
