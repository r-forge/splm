`REmod` <-
function(X, y, ind, tind, n, k, t, nT, w=NULL, coef0=0,
         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {
  ## random effects panel model 'vanilla' ML estimation
  ## based on general framework, RE structure on errors

  ## MUCH room for optimization! Take framework from ssrRE,
  ## using kinship etc. (notice inversion of order(ind,tind)!)

  ## Giovanni Millo, Trieste; this version: 22/10/2008.

  ## some useful pieces:
    Jt<-matrix(1,ncol=t,nrow=t)
    In<-diag(1,n)
    It<-diag(1,t)
    Jbart<-Jt/t
    Et<-It-Jbart

  ## inverse of Sigma
    invSigma <- function(phi, n, t) {
			 invSigma <- 1/(t*phi+1) * kronecker(Jbart, In) + kronecker(Et, In)
                   invSigma
                   }

  ## outsource determinant of sigma with an algebraically efficient form!

  ## concentrated likelihood
    ll.c<-function(phi, y, X, n, t) {

        ## perform GLS

            sigma.1<-invSigma(phi,n,t)

            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)

            e <- bhat[[2]]
            s2e <- bhat[[3]]

            due <- -n*t/2*log(s2e)
            tre <- -n/2*log(t*phi+1)
            quattro <- -1/(2*s2e)*crossprod(e,sigma.1)%*%e

            const <- -(n*t)/2*log(2*pi)
            ll.c <- const+due+tre+quattro
            llc <- - ll.c
            }

  ## iterate (=traballa) until convergence:

  myphi0 <- coef0

  optimum<-nlminb(myphi0, ll.c,
                  lower=1e-08, upper=1e08,
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, ...)


  myphi<-optimum$par
  myll <- optimum$objective


  ## optimal values of parms:
  phi<-myphi

  ## perform GLS

            sigma.1<-invSigma(phi,n,t)

            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)

  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- "phi"

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,sigma.1)%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  ## calc. cov(phi) by numerical Hessian
  covPRL <- -solve(-fdHess(myphi, function(x) -ll.c(x,y,X,n,t))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)


  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- phi
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)
  }

