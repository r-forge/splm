ivsplm <-function(formula,data=list(), index=NULL, 
                  endog = NULL, instruments= NULL, 
                  method = c("w2sls", "b2sls", "g2sls", "ec2sls"), 
                  lag = FALSE, listw = NULL, twow = FALSE, 
                  listw2 = NULL, effects = NULL, lag.instruments = FALSE, Durbin = FALSE){

# If the user do not make any choice in terms of method, when effects is Fixed the function calculates the w2sls. On the other hand, when effects is random the function calculates the ec2sls
if(length(method) != 1 && effects == "fixed")  method <- "w2sls" 	
if(length(method) != 1 && effects == "random") method <- "ec2sls" 	
		
 #if(!is.null(index)) {
    #require(plm)
 #   data <- plm.data(data, index)
 #   }

 ## substituted module based on plm.data (deprecated)
 ## with one based on pdata.frame
  if(!("pdata.frame" %in% class(data))) {
    data <- pdata.frame(data, index)
  } 
 
  index <- data[,1]
  tindex <- data[,2]

  names(index) <- row.names(data)
  ind          <- index[which(names(index)%in%row.names(data))]
  tind         <- tindex[which(names(index)%in%row.names(data))]
  spord        <- order(tind, ind)
  data         <-  data[spord,]


  ## record call
  cl <- match.call()

  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")
  
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
 
    N  <- length(unique(ind))
    k  <- dim(x)[[2]]
    t  <- max(tapply(x[,1],ind,length))
    NT <- length(ind)



  balanced <- N*t == NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")

if(isTRUE(lag) | isTRUE(lag.instruments) | isTRUE(Durbin) | inherits(Durbin, "formula")){

#### creating the block diagonal matrix for the lag model and for additional instruments

  I_T <- Diagonal(t)
  Ws  <- kronecker(I_T, listw)
  
  if(twow)  W2  <- kronecker(I_T, listw2)
  else      W2  <- NULL

  
}
else{
  Ws <- NULL
  W2  <- NULL
}
  
	#if not lag, check if there are endogenous 
if(!lag && is.null(endog)) stop("No endogenous variables specified. Please use plm instead of splm")
    
  #check if there are instruments
  if(!is.null(endog) && is.null(instruments)) stop("No instruments specified for the additional variable")

  if(!is.null(endog)){
    
    endog <- as.matrix(lm(endog, data, na.action = na.fail, method = "model.frame"))
    
    
    if(lag.instruments){
      
      instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))  
      winst <- Ws %*% instruments
      wwinst <- Ws %*% winst
      
      if(twow){
      W2 <- kronecker(I_T, listw2)
      w2inst <-   Ws %*% instruments
      w2ws.inst <- W2 %*% winst
      w2ww.inst <- W2 %*% wwinst
      instruments <- cbind(instruments, winst, wwinst, w2inst, w2ws.inst, w2ww.inst)
      }
      
      else instruments <- cbind(instruments, winst, wwinst)
      }
    else instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))

  }
  
  

  
  if(inherits(Durbin, "formula")) xdur <- as.matrix(lm(Durbin, data, na.action = na.fail, method="model.frame"))
  else xdur <- NULL
  
switch(method, 
w2sls = {

  	result <- ivplm.w2sls(Y = y, X = x, H = instruments, endog = endog, 
	                      twow = twow, lag = lag, listw = Ws, listw2 = W2,
	                      lag.instruments = lag.instruments,
	                      t = t, N = N, NT = NT, 
	                      Durbin = Durbin, xdur = xdur)
	},
b2sls = {
	result <- ivplm.b2sls(Y = y, X = x, H = instruments, endog = endog, 
	                      twow = twow, lag = lag, listw = Ws, listw2 = W2,
	                      lag.instruments = lag.instruments,
	                      t = t, N = N, NT = NT, 
	                      Durbin = Durbin, xdur = xdur)
	},
ec2sls = {
	result <- ivplm.ec2sls(Y = y, X = x, H = instruments, endog = endog, 
	                       twow = twow, lag = lag, listw = Ws, listw2 = W2,
	                       lag.instruments = lag.instruments,
	                       t = t, N = N, NT = NT, 
	                       Durbin = Durbin, xdur = xdur)
	},
g2sls = {
	result <-ivplm.g2sls(Y = y, X = x, H = instruments, endog = endog, 
	                     twow = twow, lag = lag, listw = Ws, listw2 = W2,
	                     lag.instruments = lag.instruments,
	                     t = t, N = N, NT = NT, 
	                     Durbin = Durbin, xdur = xdur)
	},
stop("...\nUnknown method\n"))


    result$zero.policy <- FALSE
    result$robust <- FALSE
    result$legacy <- FALSE
    result$listw_style <- NULL
    result$call <- match.call()
    
result
}
