`pbsjkREtest` <-
function(formula, data, w, index=NULL, ...) {

  ## performs Baltagi, Song, Jung and Koh C.3 test
  ## for RE conditional on serial correlation and spatial corr.
  ## Giovanni Millo, Trieste, this version: 19/03/2009

  ## new interface for operation with pbsjktest2.R--> and
  ## -->spreml.R

  ## depends: spreml()>ssrREmod(), fdHess{nlme} for numerical hessians

# PROBLEMS WITH THIS ONE! negative test value on the Bivand\montecarlo.R model
    ## check: negative variance comps
    ##        corner solutions due to Calzolari's Problem

  require(nlme)

  require(nlme)

  ## reorder data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }

  gindex <- data[,1]
  tindex <- data[,2]

  ## for our purpose data has to be (re)ordered
  ## by time, then group
  data <- data[order(tindex, gindex),]

  ## est. MLE SEM-AR model
  mymod <- spreml(formula=formula, data=data, w=w,
                  index=index, lag=FALSE, errors="semsr", ...)

  ## def. 'trace' function
  tr<-function(x) sum(diag(x))

  ## def. 'matrix square' function
  msq<-function(x) x%*%x

  ## make W matrix from listw object, if needed
  if("listw" %in% class(w)) if(require(spdep)) w<-listw2mat(w)

  ## retrieve restricted model's residuals (ordered!)
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula,data))
  beta0 <- mymod$coefficients
  u.hat<-y-X%*%beta0

  ## calc. data numerosities
  nt. <- dim(data)[[1]]
  n. <- length(unique(gindex))
  t. <- length(unique(tindex))

  ## calc. error variance (unique component under H0)
  sigma2e<-as.numeric(crossprod(u.hat)/nt.)

  ## retrieve autoregressive coefficient from model coef
  rho <- mymod$errcomp["rho"]
  ## retrieve SEM coefficient
  lambda <- mymod$errcomp["lambda"]

  ## henceforth notation as in Baltagi, Song, Jung, Koh (JE 2007)
  Jt<-matrix(1,ncol=t.,nrow=t.)

  V1<-matrix(ncol=t.,nrow=t.)
  for(i in 1:t.) V1[i,]<-rho^abs(1:t.-i)

  Vrho <- (1/(1-rho^2)) * V1
  iVrho<-solve(Vrho)
  VrhoJt <- solve(Vrho,Jt)

  ## this is tr(VJt) on original V=sigma2e*Vrho, see below (3.6)
  g. <- (1-rho)/sigma2e^2 * ( 2 + (t.-2)*(1-rho) )

  B<-diag(1,n.)-lambda*w
  BB<-crossprod(B)
  BB.1 <- solve(BB)

  wBBw<-crossprod(w,B)+crossprod(B,w)

  blackspade <- kronecker(VrhoJt %*% iVrho, msq(BB))

  Dhat <- -g./2 * tr(BB) + 1/(2*sigma2e^2) * crossprod(u.hat, blackspade) %*% u.hat

  ## information matrix:
  d3<-tr( wBBw%*%BB.1 )
  d6<-tr( msq( wBBw %*% BB.1 ) )

  j11<-nt./(2*sigma2e^2)
  j12<-g.*tr(BB)/(2*sigma2e)
  j13<-(n.*rho)/(sigma2e*(1-rho^2))
  j14<-t.*d3/(2*sigma2e)
  j22<-g.^2*tr(msq(BB))/2
  j23<-tr(BB)/(sigma2e*(1+rho)) * ( (2-t.)*rho^2 + (t.-1) + rho )
  j24<-g./2*tr(wBBw)
  j33<-n./(1-rho^2)^2 * (3*rho^2 - t.*rho^2 +t.-1)
  j34<-(rho*d3)/(1-rho^2)
  j44<-t.*d6/2

  Jtheta<-matrix(ncol=4,nrow=4)
  Jtheta[1,]<-c(j11,j12,j13,j14)
  Jtheta[2,]<-c(j12,j22,j23,j24)
  Jtheta[3,]<-c(j13,j23,j33,j34)
  Jtheta[4,]<-c(j14,j24,j34,j44)

  J22.1<-solve(Jtheta)[2,2]

  LMm.rl <- (Dhat^2) * J22.1

  df.<-1
  pval <- pchisq(LMm.rl,df=df.,lower.tail=FALSE)

  names(LMm.rl)="LM"
  names(df.)<-"df"

  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = LMm.rl, parameter = df.,
               method = "Baltagi, Song, Jung and Koh C.3 conditional test",
               alternative = "random effects, sub serial corr. and spatial dependence in error terms",
               p.value = pval,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}

