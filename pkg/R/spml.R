spml <- function(formula, data, index=NULL, listw, listw2=listw,
                 model=c("within","random","pooling"),
                 effect=c("individual","time","twoways"),
                 lag=FALSE, spatial.error=c("b","kkp","none"),
                 ...) {

  ## wrapper function for all ML models

  ## check class(listw)
  checklw <- function(x) {
    if(!("listw" %in% class(x))) {
      if("matrix" %in% class(x)) {
        require(spdep)
        x <- listw2mat(x)
      } else {
        stop("'listw' has to be either a 'listw' or a 'matrix' object")
      }}
    return(x)
  }
  checklw(listw)
  checklw(listw2)

  ## dimensions check is moved downstream

  switch(match.arg(model), within={
    if(lag) {
      model <- switch(match.arg(spatial.error), b="sarar",
                      kkp="sarar", none="lag")
    } else {
      model <- "error"
    }
    effects <- switch(match.arg(effect), individual="spfe",
                      time="tpfe", twoways="sptpfe")
    res <- spfeml(formula=formula, data=data, index=index,
                  listw=listw, listw2=listw2,
                  model=model, effects=effects,
                  ...)
  }, random={
    switch(match.arg(effect),
           time={stop("time random effects not implemented")},
           twoways={stop("twoway random effects not implemented")},
           individual={
             errors <- switch(match.arg(spatial.error),
                              b="semre", kkp="sem2re", none="re")})
    res <- spreml(formula=formula, data=data, index=index,
                  w=listw2mat(listw), w2=listw2mat(listw2),
                  lag=lag, errors=errors, ...)
  }, pooling={
           errors <- switch(match.arg(spatial.error),
                              b="sem", kkp="sem", none="ols")
    res <- spreml(formula=formula, data=data, index=index,
                  w=listw2mat(listw), w2=listw2mat(listw2),
                  lag=lag, errors=errors, ...)
         })

  return(res)
}

    