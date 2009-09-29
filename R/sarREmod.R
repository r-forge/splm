`sarREmod` <-
function(X,y,ind,tind,n,k,t,nT,w,coef0=NULL,
                   quiet=TRUE,hess=FALSE,
                   tol=10e-03,ltol=1000,...) {
  ## spatial autoregressive, random effects panel model estimation
  ## based on general framework, spatial filter on y plus RE structure on errors

  ## based on (spatial/spatialpanel/AR/g/) sarREmodg()
  ## version 3 (efficient lag on y)

  ## Giovanni Millo, Trieste; this version: 25/2/2009.

  require(nlme) # for numerical hessians

  ## preliminary stuff removed

  ## some useful pieces:
    Jt<-matrix(1,ncol=t,nrow=t)
    In<-diag(1,n)
    It<-diag(1,t)
    Jbart<-Jt/t
    Et<-It-Jbart

  ## inverse of Sigma (Random Effects)
    invSigma <- function(phi, n, t) {
			 invSigma <- 1/(t*phi+1) * kronecker(Jt/t, In) + kronecker(Et, In)
                   invSigma
                   }

  ## specific for spatial lag: ##

  ## spatial lag operator
    A<-function(psi) diag(1,n)-psi*w

  ## determinant of A
    detA<-function(psi) det(A(psi)) # use more efficient versions from Elhorst

  ## calc. Wy (spatial lag of y)
    wy<-list()
    for (j in 1:length(unique(tind))) {
      yT<-y[tind==unique(tind)[j]]
      wy[[j]] <- w %*% yT
      }
    wy<-unlist(wy)

  ## beta_hat as GLS estimator and error variance estimator
  ## *spatial-filter-on-y version*
    bhat<-function(y,X,psi,phi,n,t) {

            ## spatial filter on y: Ay=Ay(psi)
            Ay<-y-psi*wy

            ## invert Sigma: given n,t: Sigma.1=Sigma.1(phi)
            sigma.1<-invSigma(phi,n,t)

            ## GLS step
            b.hat<-solve( crossprod(X,sigma.1)%*%X, crossprod(X,sigma.1)%*%Ay )
            ehat<-Ay-X%*%b.hat
            sigma2ehat<-crossprod(ehat,sigma.1)%*%ehat/(n*t)

            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat,sigma.1=sigma.1)
            bhat
            }

  ## end: specific for spatial lag ##

  ## concentrated likelihood: SAR with RE
    ll.c<-function(theta) {
            psi<-theta[1]
            phi<-theta[2]

            ## the concentrated likelihood needs to include
            ## the spatial filtering step, as Ay=Ay(psi)

            beta<-bhat(y,X,psi,phi,n,t)
            b.hat<-beta$b.hat
            e<-beta$e
            s2e<-beta$sigma2
            sigma.1<-beta$sigma.1

            uno <- t*log(detA(psi))
            due <- -n*t/2*log(s2e)
            tre <- -n/2*log(t*phi+1)
            quattro <- -1/(2*s2e)*crossprod(e,sigma.1)%*%e

            ll.c<-uno+due+tre+quattro
            }

  ## iterate (=traballa) until convergence:

  ## init mybhat as beta OLS, ll parms as (0.5,1)
  mybhat<-solve(crossprod(X),crossprod(X,y))
  mytheta<-c(0.5,1)

    ## initialize values
    i<-0
    mytol<-1
    beta<-bhat(y,X,mytheta[1],mytheta[2],n,t)

    ## iterate
    while(abs(mytol)>tol) {
      beta0<-beta
      mytheta0<-mytheta

      beta<-bhat(y,X,mytheta[1],mytheta[2],n,t)

      mytheta<-optim(mytheta,ll.c,method="L-BFGS-B",
                        lower=c(-0.9,1e-08),upper=c(0.9,10e8),
                        control=list(fnscale=-1,factr=ltol))$par

      i<-i+1
      if(!quiet) {
        print(paste("Iteration no.",i))
        print(beta[1])
        print(mytheta)
        }

      ## update tolerance condition
      mytol<-max(max(abs(beta[[1]]-beta0[[1]])),max(abs(mytheta-mytheta0)))
      }

  ## now beta and phirholambda are optimal values
  ## calc. the log-likelihood:
  myll <- ll.c(mytheta)


  ## optimal values of parms:
  psi<-mytheta[1]
  phi<-mytheta[2]

  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.arcoef <- "psi"
  nam.errcomp <- "phi"

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]])*solve(crossprod(X,invSigma(phi, n, t))%*%X)
  dimnames(covB) <- list(nam.beta, nam.beta)

  ## calc. cov(psi, phi) by numerical Hessian
  covTheta <- solve(-fdHess(mytheta,ll.c)$Hessian)
  covAR <- covTheta[1,1,drop=FALSE]
  covPRL <- covTheta[2,2,drop=FALSE]
  dimnames(covAR) <- list(nam.arcoef, nam.arcoef)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  arcoef <- psi
  errcomp <- phi
  names(betas) <- nam.beta
  names(arcoef) <- nam.arcoef
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, arcoef=arcoef, errcomp=errcomp,
              covB=covB, covAR=covAR, covPRL=covPRL, ll=myll)

  return(RES)
  }

