`semarmod` <-
function(X, y, ind, tind, n, k, t, nT, w, coef0=rep(0,3),
         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {
  ## spatial- and time- error autoregressive (pooling) panel model estimation
  ## modified by restriction (phi=0) from Appendix A.2 in Baltagi, Song,
  ## Jung and Koh WP, version May 2004; restricted likelihood in (3.15).

  ## V matrix, as V/sigma.e^2 in BJSK 2007
    Vmat<-function(rho,t) {
      V1<-matrix(ncol=t,nrow=t)
      for(i in 1:t) V1[i,]<-rho^abs(1:t-i)
      V <- (1/(1-rho^2)) * V1
      }

  ## spatial lag operator
    B<-function(lambda) diag(1,n)-lambda*w

  ## determinant of B
    detB<-function(lambda) det(B(lambda)) # use more efficient versions from Elhorst

  ## inverse of Sigma
    invSigma <- function(rho, lambda, n, t) {
                   invVmat<-solve(Vmat(rho,t))
                   BB<-crossprod(B(lambda))
                   invSigma<-kronecker(invVmat,BB)
                   invSigma
                   }


  ## concentrated likelihood
    ll.c<-function(rholambda,y,X,n,t,w) {
            rho<-rholambda[1]
            lambda<-rholambda[2]

            ## beta_hat as GLS estimator and error variance estimator
            sigma.1<-invSigma(rho,lambda,n,t)
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)
            e <- bhat[[2]]
            s2e <- bhat[[3]]

            uno <- n/2*log(1-rho^2)
            tre <- -(n*t)/2*log(s2e)
            quattro <- t*log(detB(lambda))
            cinque <- -1/(2*s2e)*crossprod(e,invSigma(rho, lambda, n, t))%*%e
            const <- -(n*t)/2*log(2*pi)

            ll.c<-const+uno+tre+quattro+cinque
            llc <- -ll.c

            }

  ## iterate (=traballa) until convergence:

  myrholambda0 <- coef0

  optimum<-nlminb(myrholambda0, ll.c,
                  lower=c(-0.999,-0.999), upper=c(0.999,0.999),
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, w=w, ...)


  myrholambda<-optimum$par
  myll <- optimum$objective

  ## optimal values of parms:
  rho<-myrholambda[1]
  lambda<-myrholambda[2]

  ## perform GLS
            sigma.1<-invSigma(rho,lambda,n,t)
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)


  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- c("rho","lambda")

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,invSigma(rho, lambda, n, t))%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  ## calc. cov(phi,rho,lambda) by numerical Hessian
  covPRL <- solve(-fdHess(myrholambda, function(x) -ll.c(x,y,X,n,t,w))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- c(rho, lambda)
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)
  }

