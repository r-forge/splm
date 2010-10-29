semarREmod<-function(X, y, ind, tind, n, k, t, nT, w, coef0=rep(0,3),
                     hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {
  ## spatial- and time- error autoregressive, random effects panel model estimation
  ## following Appendix A.2 in Baltagi, Song, Jung and Koh WP, version May 2004

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

  ## some useful pieces:
    alfa2<-function(rho) (1+rho)/(1-rho)
    d2<-function(rho,t) alfa2(rho)+t-1
    Jt<-matrix(1,ncol=t,nrow=t)
    In<-diag(1,n)

  ## second determinant in (A.17)
    det2<-function(phi,rho,lambda,t) det( d2(rho,t)* (1-rho)^2*phi*In+solve(crossprod(B(lambda))))
  ## Z0 for (A.13)
    Z0 <- function(phi,rho,lambda,t) solve(d2(rho,t)*(1-rho)^2*phi*In + solve(crossprod(B(lambda))) )

  ## inverse of Sigma
    invSigma <- function(phi, rho, lambda, n, t) {
                   invVmat<-solve(Vmat(rho,t))
                   BB<-crossprod(B(lambda))
                   invSi1<-kronecker(invVmat,BB)
                   invSi2<-1/(d2(rho,t)*(1-rho)^2)
                   invSi3<-kronecker( solve(Vmat(rho,t),Jt)%*%invVmat, Z0(phi,rho,lambda,t)-BB )
                   invSigma <- invSi1 + invSi2*invSi3
                   invSigma
                   }


  ## concentrated likelihood
    ll.c<-function(phirholambda,y,X,n,t,w) {
            phi<-phirholambda[1]
            rho<-phirholambda[2]
            lambda<-phirholambda[3]

            ## perform GLS
            sigma.1<-invSigma(phi,rho,lambda,n,t)
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-(crossprod(ehat,sigma.1)%*%ehat)/(n*t)
            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)
            e <- bhat[[2]]
            s2e <- bhat[[3]]

            uno <- n/2*log(1-rho)
            due <- -1/2*log(det2(phi,rho,lambda,t))
            tre <- -(n*t)/2*log(s2e)
            quattro <- (t-1)*log(detB(lambda))
            cinque <- -1/(2*s2e)*crossprod(e,sigma.1)%*%e

            const <- -(n*t)/2*log(2*pi)
            ll.c <- const+uno+due+tre+quattro+cinque
            llc <- - ll.c

            }


  myphirholambda0 <- coef0

  optimum<-nlminb(myphirholambda0, ll.c,
                  lower=c(1e-8,-0.999,-0.999), upper=c(10e8,0.999,0.999),
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, w=w, ...)


  myphirholambda<-optimum$par
  myll <- optimum$objective

  ## optimal values of parms:
  phi<-myphirholambda[1]
  rho<-myphirholambda[2]
  lambda<-myphirholambda[3]

  ## perform GLS
            sigma.1<-invSigma(phi,rho,lambda,n,t)
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)


  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- c("phi","rho","lambda")

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,invSigma(phi, rho, lambda, n, t))%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  ## calc. cov(phi,rho,lambda) by numerical Hessian
  covPRL <- solve(-fdHess(myphirholambda, function(x) -ll.c(x,y,X,n,t,w))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- c(phi, rho, lambda)
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)

  }


