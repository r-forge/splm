
`summary.splm` <-
function(object,...){

  ## summary method for splm objects
  ## adds incrementally to the model object, as summary.plm does
  ## structure remans the same for all type but 'spsegm' (symultaneous equations requires a special printing)

  if (object$type=='spsegm'){
    	coeff<-object$coefficients
		eq<-object$EQ
		var<-as.matrix(object$vcov)
		ser<-sqrt(diag(as.matrix(var)))
		tr<-coeff/ser
		pr<-pnorm(abs(as.matrix(tr)), lower.tail=FALSE)*2
		Xnam<-object$Xnames
		Ynam<-object$Ynames
		reg<-object$K
		sp<-object$spec
		lags<-object$lags
		errors<-object$errors
		endogenous<-object$endogenous
		rho<-object$rho
		nlags<-length(which(unlist(lags)==TRUE))
      nend<-length(which(unlist(endogenous)==TRUE))
		numx<-vector("numeric",eq)
for (i in 1:eq) numx[i]<- sp[i] + length(which(lags[[i]]==TRUE)) + length(which(endogenous[[i]]==TRUE))
		tmp<-seq(1,eq)
      tmp2<-rep(tmp,numx)
		b<-split(coeff,tmp2)
		se<-split(as.matrix(ser), tmp2)
		t<-split(tr,tmp2)
		p<-split(pr,tmp2)
		object$b <- b
		object$se<-se
		object$t<-t
		object$p<-p
		object$eq<-eq
		object$xnam<-Xnam
		object$ynam<-Ynam
		object$sp<-object$spec
		object$lags<-object$lags
		object$errors<-object$errors
		object$endogenous<-object$endogenous

		object$type<-'spsegm'
		class(object)<- c("summary.splm","splm")
		object
  	} else {
            ## to date, only balanced panels are allowed for 'splm'
            balanced <- TRUE #attr(object,"pdim")$balanced
            model.name <- object$type #attr(object,"pmodel")$model
            effect <- "individual" #attr(object,"pmodel")$effect
            ## make coefficients' table if vcov exist
            if (!is.null(object$vcov)) {
                std.err <- sqrt(diag(object$vcov)) #vcov(object) doesn't work
                b <- coefficients(object)
                z <- b/std.err
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                CoefTable <- cbind(b,std.err,z,p)
                colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$CoefTable <- CoefTable
            } else {
                object$CoefTable <- cbind(coefficients(object))
                colnames(object$CoefTable) <- c("Estimate")
            }

            if (object$type == "fixed effects error" && object$method != "eigen") {
                lambda <- object$spat.coef
                object$lambda <- lambda
            }

            if (object$type == "random effects GM" ) {
                lambda <- object$rho
                object$lambda <- lambda
            }

            ## make AR coefficient of y's table
            if(!is.null(object$vcov.arcoef)) {
                std.err1 <- sqrt(diag(object$vcov.arcoef))
                ar <- object$arcoef
                z <- ar/std.err1
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                ARCoefTable <- cbind(ar,std.err1,z,p)
                colnames(ARCoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$ARCoefTable <- ARCoefTable
            }


            ## make error comps' table
            if(!is.null(object$vcov.errcomp)) {
                std.err2 <- sqrt(diag(object$vcov.errcomp))
                ec <- object$errcomp
                z <- ec/std.err2
                p <- 2*pnorm(abs(z),lower.tail=FALSE)
                ErrCompTable <- cbind(ec,std.err2,z,p)
                colnames(ErrCompTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
                object$ErrCompTable <- ErrCompTable
            }

            object$ssr <- sum(residuals(object)^2)
            object$tss <- tss(object$model[[1]])
            object$rsqr <- 1-object$ssr/object$tss
            object$fstatistic <- "nil" #Ftest(object)
            class(object) <- c("summary.splm","splm")
            object
        }
}

