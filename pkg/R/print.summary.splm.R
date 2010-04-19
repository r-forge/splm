`print.summary.splm` <-
function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...) {
    if (x$type=='spsegm') {
        cat("\nSimultaneous Equations Model:\n")
        cat("\nCall:\n")
        print(x$call)
        eq<-x$eq
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))
		  numx<-vector("numeric",eq)
		  
        for (i in 1:eq) {
            cat(" \n" )
            cat(paste('Equation', i, sep = " ", collapse = ""),"\n", sep = "")
            #cat("\n _______________________________________________________ \n")
            tables<- cbind( x$b[[i]], x$se[[i]], x$t[[i]], x$p[[i]] )
 #           print(!is.null(which(x$lags[[i]]==TRUE)))
if(length(which(x$lags[[i]]==TRUE)) !=0 ){
	Wynames<-paste("W",x$ynam[which(x$lags[[i]]==TRUE)])
   rn<-c(Wynames, x$ynam[which(x$endogenous[[i]]==TRUE)], x$xnam[[i]] )
	} 
	else rn<-c(x$ynam[which(x$endogenous[[i]]==TRUE)], x$xnam[[i]] )
        
            dimnames(tables)<-list(rn ,c("Estimate", "Std.Error", "t value", "Pr(>|t|)"))


            if(i==eq) {
                legend=TRUE
            } else {
                legend=FALSE
            }

            printCoefmat(tables,digits=digits, signif.stars=TRUE,signif.legend=legend)
       cat(" \n" )     
            if(x$errors[i] != FALSE) cat(paste('Spatial autoregressive parameter:', round(x$rho[i],4), sep = " ", collapse = ""),"\n", sep = "")

            cat("\n _______________________________________________________ \n")
        }

    } else {

        cat(paste("Spatial panel",x$type,"model\n"))
        cat("\nCall:\n")
        print(x$call)
        cat("\nResiduals:\n")
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))
        print(sumres(x))

        if(!is.null(x$ErrCompTable)) {
            cat("\nError variance parameters:\n")
            printCoefmat(x$ErrCompTable,digits=digits,signif.legend=FALSE)
        }

        if(is.numeric(x$lambda)) {
            cat("\nEstimated spatial coefficient, variance components and theta:\n")
            print(x$lambda)
        }

        if(!is.null(x$ARCoefTable)) {
            cat("\nSpatial autoregressive coefficient:\n")
            printCoefmat(x$ARCoefTable,digits=digits,signif.legend=FALSE)
        }

        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable,digits=digits)
        cat("\n")

    }

    invisible(x)
}

