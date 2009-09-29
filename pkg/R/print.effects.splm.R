`print.effects.splm` <-
function(x,digits= max(3, getOption("digits") - 2),
...){
	object<-x
	effects<-object$effects
if(effects=="tpfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits, signif.legend=FALSE)

 cat("\n")  

	  cat("\nTime period fixed effects:\n")
  printCoefmat(object$TETable,digits=digits)

out<-rbind(object$INTTable,object$TETable)
}
	
if(effects=="spfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits,signif.legend=FALSE)

 cat("\n")  
 
	  cat("\nSpatial fixed effects:\n")
  printCoefmat(object$SETable,digits=digits)


out<-rbind(object$INTTable,object$SETable)
}
if(effects=="sptpfe"){
	  cat("\nIntercept:\n")
  printCoefmat(object$INTTable,digits=digits,signif.legend=FALSE)

 cat("\n")  
 
	  cat("\nSpatial fixed effects:\n")
  printCoefmat(object$SETable,digits=digits,signif.legend=FALSE)
 
  cat("\n")  
   
  	  cat("\nTime period fixed effects:\n")
  printCoefmat(object$TETable,digits=digits)

out<-rbind(object$INTTable,object$SETable,object$TETable)
}
	}


`write.effects.splm` <-
function(object,filename='effects.csv',...){
	effects<-object$effects
if(effects=="tpfe"){

out<-rbind(object$INTTable,object$TETable)
write.csv(out, "filename")

}
	
if(effects=="spfe"){
out<-rbind(object$INTTable,object$SETable)
write.csv(out, "filename")
}
if(effects=="sptpfe"){
out<-rbind(object$INTTable,object$SETable,object$TETable)
write.csv(out, "filename")
}
	}
