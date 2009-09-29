`spreml` <-
function(formula, data, index = NULL, w, lag=FALSE,
           errors = c("semsrre","semsr","srre","semre","re","sr","sem"),
           pvar = FALSE, hess=FALSE, quiet=TRUE,
           initval = c("zeros", "estimate"),
           x.tol=1.5e-18, rel.tol=1e-15,
           ...) {

  ## set trace parm for optimizer
  trace <- as.numeric(!quiet)

  ## check time variation
  if(pvar) print("<implement pvar>")  ## see what it does!

  ## reorder data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }

  index <- data[,1]
  tindex <- data[,2]

  ## record call
  cl <- match.call()

  require(nlme) # for numerical hessians

  ## check w
  if(!is.matrix(w)) {
      if("listw" %in% class(w)) {
          require(spdep)
          w <- listw2mat(w)
      } else {
          stop("w has to be either a 'matrix' or a 'listw' object")
      }}

  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")

  ## reduce X,y to model matrix values (no NAs)
  X<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  ## reduce index accordingly
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(X))]
  tind<-tindex[which(names(index)%in%row.names(X))]

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  X<-X[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  ## det. number of groups and df
  n<-length(unique(ind))
  k<-dim(X)[[2]]
  ## det. max. group numerosity
  t<-max(tapply(X[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT<-length(ind)

  ## check dim(w)
  if(dim(w)[[1]]!=n) stop("Non conformable spatial weights")

  ## check whether the panel is balanced
  balanced<-n*t==nT
  if(!balanced) stop("Estimation method unavailable for unbalanced panels")

  ## supply vector or estimate model for beta-parms' starting values

  ## consistency check: expected coefficients' vector length
      sv.length <- switch(match.arg(errors), semsrre = 3, semsr = 2, ssrre = 2,
                          semre = 2, re = 1, ssr = 1, sem = 1)

  ## extract 'errors' argument from admissible set
      errors. <- match.arg(errors)

  ## set initial values for parameters
    if(is.numeric(initval)) {

        ## if numeric, just check length
        if(length(initval)!=sv.length) {
            stop("Incorrect number of initial values supplied for error vcov parms")
        }
        coef0 <- initval

        } else {

            ## if char,
            switch(match.arg(initval), zeros = {

                coef0 <- rep(0, sv.length)

                }, estimate = {

                    ## check that the model has >1 parms
                    if(nchar(errors.)<4) {
                        stop("Pre-estimation of unique vcov parm is meaningless: \n please select (default) option 'zeros' or supply a scalar")
                        }

                    coef0 <- NULL

                    if(grepl("re", errors.)) {
                        REmodel<- REmod(X, y, ind, tind, n, k, t, nT, w, coef0=0,
                                        hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...)
                        coef0 <- c(coef0, REmodel$errcomp)}

                    if(grepl("sr", errors.)) {
                        ARmodel<- ssrmod(X, y, ind, tind, n, k, t, nT, w, coef0=0,
                                         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...)
                        coef0 <- c(coef0, ARmodel$errcomp)}

                    if(grepl("sem", errors.)) {
                        SEMmodel<- semmod(X, y, ind, tind, n, k, t, nT, w, coef0=0,
                                          hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...)
                        coef0 <- c(coef0, SEMmodel$errcomp)}

                    })
            }


  ## call the relevant computing engine:
  ## estimator is passed on as a function
  if(lag) stop("Method not yet implemented")
    
    #  est.fun <- switch(match.arg(errors), semsrre = {
     #     saremsrREmod
      #}, semsr = {
       #   saremsrmod
      #}, ssrre = {
       #   sarsrREmod
     # }, semre = {
      #    saremREmod
     # }, re = {
      #    sarREmod
     # }, ssr = {
      #    sarsrmod
     # }, sem = {
      #    saremmod
     # })
 else {
      est.fun <- switch(match.arg(errors), semsrre = {
          semarREmod
      }, semsr = {
          semarmod
      }, ssrre = {
          ssrREmod
      }, semre = {
          semREmod
      }, re = {
          REmod
      }, ssr = {
          ssrmod
      }, sem = {
          semmod
      })
      arcoef <- NULL
  }

  RES <- est.fun(X, y, ind, tind, n, k, t, nT,
                 w=w, coef0=coef0, hess=hess,
                 trace=trace, x.tol=x.tol, rel.tol=rel.tol)

  ## from here on has to be fixed to allow for lag
  ## in calc. of fitted values etc.

  ## calc. fitted and residuals
  y.hat <- as.vector(X%*%RES$betas)
  res <- y - y.hat

  nam.rows <- dimnames(X)[[1]]
  names(y.hat) <- nam.rows
  names(res) <- nam.rows

  ## make model data
  model.data <- data.frame(cbind(y,X[,-1]))
  dimnames(model.data)[[1]] <- nam.rows
  #dimnames(model.data)[[2]] <-  # fix y's name

  type <- "random effects ML"

  sigma2 <- list(one=3,idios=2,id=1) # fix this!

  spmod <- list(coefficients=RES$betas, arcoef=RES$arcoef,
                errcomp=RES$errcomp,
                vcov=RES$covB, vcov.arcoef=RES$covAR,
                vcov.errcomp=RES$covPRL,
                residuals=res, fitted.values=y.hat,
                sigma2=sigma2, model=model.data,
                type=type, call=cl,
                errors=errors, logLik=RES$ll)

  class(spmod) <- "splm"

  return(spmod)
}

