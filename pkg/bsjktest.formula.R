`bsjktest.formula` <-
function(x, data, w, test=c(paste("C",1:3,sep="."),"J"), index=NULL, ...){

  ## transform data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }

  gindex <- dimnames(data)[[2]][1]
  tindex <- dimnames(data)[[2]][2]

  switch(match.arg(test), C.1 = {

    bsjk = pbsjkSDtest(formula=x, data=data, w=w, index=index, ...)

  }, C.2 = {

    bsjk = pbsjkARtest(formula=x, data=data, w=w, index=index, ...)

  }, C.3 = {

    bsjk = pbsjkREtest(formula=x, data=data, w=w, index=index, ...)

  }, J = {

    bsjk = pbsjkJtest(formula=x, data=data, w=w, index=index, ...)

  })

  return(bsjk)

}

