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


