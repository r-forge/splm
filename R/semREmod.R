`semREmod` <-
function(X, y, ind, tind, n, k, t, nT, w, coef0=c(0,0),
         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {

  ## spatial error random effects panel model estimation
  ## based on general framework, spatial structure on errors
  ## (see likelihood and Sigmas in Baltagi et al.)

  ## some useful pieces:
    Jt<-matrix(1,ncol=t,nrow=t)
    In<-diag(1,n)
    It<-diag(1,t)
    Jbart<-Jt/t
    Et<-It-Jbart

  ## spatial lag operator
    B<-function(lambda) diag(1,n)-lambda*w

  ## determinant of A
    detB<-function(lambda) det(B(lambda)) # use more efficient versions from Elhorst

  ## inverse of Sigma (Spatial error and Random Effects)
    invSigma <- function(phi, lambda, n, t) { # use more efficient algebra here
			 BB <- crossprod(B(lambda))
                   invSigma <- kronecker(Jbart, solve(t*phi*In + solve(BB))) + kronecker(Et, BB)
                   invSigma
                   }

  ## determinant of Sigma (Spatial error and Random Effects)
    detSigma <- function(phi, lambda, n, t) { # use more efficient algebra here
                  detSigma <- -1/2*log( det( t*phi*In +
                              solve(crossprod(B(lambda))) ) ) +
                              (t-1)*log(detB(lambda))
                  detSigma
                  }

  ## concentrated likelihood: random effects SEM
    ll.c<-function(philambda, y, X, n, t, w) {
            phi<-philambda[1]
            lambda<-philambda[2]

            ## perform GLS

            ## invert Sigma:
            sigma.1<-invSigma(phi, lambda, n, t)

            ## GLS step
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)

            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat,sigma.1=sigma.1)

            e <- bhat[[2]]
            s2e <- bhat[[3]]

            due <- detSigma(phi, lambda, n, t)
            tre <- -n*t/2*log(s2e)
            quattro <- -1/(2*s2e)*crossprod(e,sigma.1)%*%e

            const <- -(n*t)/2*log(2*pi)
            ll.c <- const+due+tre+quattro
            llc <- - ll.c
            }

  ## iterate (=traballa) until convergence:

  myphilambda0 <- coef0

  optimum<-nlminb(myphilambda0, ll.c,
                  lower=c(1e-08, -0.999), upper=c(1e08, 0.999),
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, w=w, ...)


  myphilambda<-optimum$par
  myll <- optimum$objective

  ## optimal values of parms:
  phi<-myphilambda[1]
  lambda<-myphilambda[2]

  ## perform GLS
            ## invert Sigma: given n,t: Sigma.1=Sigma.1(phi)
            sigma.1<-invSigma(phi, lambda, n, t)
            ## GLS step
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%y )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat,sigma.1=sigma.1)

  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- c("phi", "lambda")

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,invSigma(phi, lambda, n, t))%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  covPRL <- solve(-fdHess(myphilambda, function(x) -ll.c(x,y,X,n,t,w))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- c(phi, lambda)
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)
  }

