`print.summary.splm` <-
function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...) {
    if (x$type=='spsegm') {
        cat("\nSymultaneous Equation Model:\n")
        cat("\nCall:\n")
        print(x$call)
        eq<-x$eq
        save.digits <- unlist(options(digits=digits))
        on.exit(options(digits=save.digits))

        for (i in 1:eq) {
            cat(" \n" )
            cat(paste('Equation', i, sep = " ", collapse = ""),"\n", sep = "")
            cat("\n _______________________________________________________ \n")
            tables<- cbind( x$b[[i]], x$se[[i]], x$t[[i]], x$p[[i]] )

            if (is.list(x$sp)) {
                spo<-x$sp[[i]]
                rn<-c(paste("W",x$ynam), x$ynam[-i], x$xnam[spo] )
            } else {
                rn<-c( paste("W",x$ynam), x$ynam[-i], x$xnam )
            }

            dimnames(tables)<-list(rn ,c("Estimate", "Std.Error", "t value", "Pr(>|t|)"))

            if(i==eq) {
                legend=TRUE
            } else {
                legend=FALSE
            }

            printCoefmat(tables,digits=digits, signif.stars=TRUE,signif.legend=legend)
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

