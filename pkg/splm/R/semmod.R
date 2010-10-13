`semmod` <-
function(X, y, ind, tind, n, k, t, nT, w, coef0=0,
         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {

  ## spatial error (pooling) panel model estimation
  ## based on general framework, spatial structure on errors
  ## (see likelihood and Sigmas in Anselin, LeGallo and Jayet, page 25)

  ## some useful pieces:
    It<-diag(1,t)

  ## spatial lag operator
    A<-function(psi) diag(1,n)-psi*w

  ## determinant of A
    detA<-function(psi) det(A(psi)) # use more efficient versions from Elhorst

  ## inverse of Sigma (Spatial error, pooling)
    invSigma <- function(lambda, n, t) {
			 BB <- crossprod(A(lambda))
                   invSigma <- kronecker(It, BB)
                   invSigma
                   }

  ## (log-) determinant of Sigma (Spatial error, pooling)
    detSigma <- function(lambda) {
                  detSigma <- -2*t*log(detA(lambda))
                  detSigma
                  }

  ## concentrated likelihood: pooling SEM
    ll.c<-function(lambda, y, X, n, t, w) {

            ## perform GLS

            ## invert Sigma: given n,t: Sigma.1=Sigma.1(phi)
            sigma.1<-invSigma(lambda,n,t)

            ## GLS step
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)

            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat,sigma.1=sigma.1)

            e <- bhat[[2]]
            s2e <- bhat[[3]]

            due <- t*log(detA(lambda)) #-1/2*detSigma(lambda)
            tre <- -n*t/2*log(s2e)
            quattro <- -1/(2*s2e)*crossprod(e,sigma.1)%*%e

            const <- -(n*t)/2*log(2*pi)
            ll.c <- const+due+tre+quattro
            llc <- - ll.c
            }

  ## iterate (=traballa) until convergence:

  mylambda0 <- coef0

  optimum<-nlminb(mylambda0, ll.c,
                  lower=-0.999, upper=0.999,
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, w=w, ...)


  mylambda<-optimum$par
  myll <- optimum$objective

  ## optimal values of parms:
  lambda<-mylambda

  ## perform GLS
            ## invert Sigma: given n,t: Sigma.1=Sigma.1(phi)
            sigma.1<-invSigma(lambda,n,t)
            ## GLS step
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat,sigma.1=sigma.1)

  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- "lambda"

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,invSigma(lambda, n, t))%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  covPRL <- solve(-fdHess(mylambda, function(x) -ll.c(x,y,X,n,t,w))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- lambda
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)
  }

