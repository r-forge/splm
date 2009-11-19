effects.splm<-function(object,...){
	x<-object
	if (class(x) != "splm") stop(paste("methos not implemented for objects of class", class(x)))
	if (class(x) != "splm" && (x$type != "fixed effects lag" || x$type != "fixed effects error")) stop(paste("methos not implemented for objects of type", x$type))
	all.FE<-x$res.eff[[1]]
	effects <- x$effects
if (effects=="pooled") stop("No fixed effects available if effects == pooled")
if(effects=="spfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.se.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,SETable=SETable, effects=effects)
	}
if(effects=="tpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.se.tfe 
	z <- TP.EFF/se.TP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   TETable <- cbind(TP.EFF,se.TP.EFF,z,p)
	colnames(TETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,TETable=TETable,effects=effects)
	}
if(effects=="sptpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.se.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.se.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
		TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.se.tfe 
	z <- TP.EFF/se.TP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   TETable <- cbind(TP.EFF,se.TP.EFF,z,p)
	colnames(TETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,SETable=SETable,TETable=TETable,effects=effects)
	}
res
class(res) <- "effects.splm"
return(res)
	}
	
	
print.effects.splm <-
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


write.effects.splm <-
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
