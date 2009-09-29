`effects.splm`<- function(object,...){
	x<-object
	if (class(x) != "splm") stop(paste("methos not implemented for objects of class", class(x)))
	if (class(x) != "splm" && (x$type != "fixed effects lag" || x$type != "fixed effects error")) stop(paste("methos not implemented for objects of type", x$type))
	all.FE<-x$res.eff[[1]]
	effects <- x$effects
if (effects=="pooled") stop("No fixed effects available if effects == pooled")
if(effects=="spfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.t.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.t.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,SETable=SETable, effects=effects)
	}
if(effects=="tpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.t.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.t.tfe 
	z <- TP.EFF/se.TP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   TETable <- cbind(TP.EFF,se.TP.EFF,z,p)
	colnames(TETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
res<-list(INTTable=INTTable,TETable=TETable,effects=effects)
	}
if(effects=="sptpfe"){
	INT <- all.FE$intercept
	se.INT<- all.FE$res.t.con
	z <- INT/se.INT
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   INTTable <- cbind(INT,se.INT,z,p)
	colnames(INTTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	rownames(INTTable) <- "(Intercept)"
	SP.EFF <- all.FE$res.sfe
	se.SP.EFF <- all.FE$res.t.sfe 
	z <- SP.EFF/se.SP.EFF
   p <- 2*pnorm(abs(z),lower.tail=FALSE)
   SETable <- cbind(SP.EFF,se.SP.EFF,z,p)
	colnames(SETable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
		TP.EFF <- all.FE$res.tfe
	se.TP.EFF <- all.FE$res.t.tfe 
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

