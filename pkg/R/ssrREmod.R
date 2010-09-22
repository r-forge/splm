`ssrREmod` <-
function(X, y, ind, tind, n, k, t, nT, w=NULL, coef0=rep(0,2),
         hess=FALSE, trace=trace, x.tol=1.5e-18, rel.tol=1e-15, ...) {
  ## time- error autoregressive, random effects panel model estimation
  ## following Appendix A.2 in Baltagi, Song, Jung and Koh WP, version May 2004

  ## w included only for parameters compatibility with splm(), but
  ## obviously not used!

  ## this version 5+ works!! cfr. lme() (some diff. in phi only)
  ## corrected typo, missing log() in likelihood (3.7) BSJK's paper

  ## from version 4: elimination of Kron. prods, use bdsmatrices,
  ## reversed data ordering (unaffected, see cfr. ver.3).

  ## timings against ver. 3: 7'' both on Munnell's data 48x5,
  ## 13'' vs. 110'' on 48x10, 20'' vs. 390'' on full 48x17.

  ## Giovanni Millo, Trieste; this 'splm' version: 20/10/2008.

  require(kinship) # for bdsmatrices

  ## reorder data again (were passed as order(tind, ind) for compatibility
  ## with the wrapper splm(), but here it is more convenient to order 'the
  ## usual way')
  oo<-order(ind,tind)
  X<-X[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  ## (re-)determine number of groups and df
  n<-length(unique(ind))
  k<-dim(X)[[2]]
  ## det. max. group numerosity
  t<-max(tapply(X[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT<-length(ind)


  ## V matrix, as V/sigma.e^2 in BJSK 2007
    Vmat<-function(rho,t) {
      V1<-matrix(ncol=t,nrow=t)
      for(i in 1:t) V1[i,]<-rho^abs(1:t-i)
      V <- (1/(1-rho^2)) * V1
      }

  ## some useful pieces:
    alfa2<-function(rho) (1+rho)/(1-rho)
    d2<-function(rho,t) alfa2(rho)+t-1
    Jt<-matrix(1,ncol=t,nrow=t)
    In<-diag(1,n)


  ## typical TxT block, Sigma
    bSigma <- function(phi, rho, n, t) {
			 bSigma <- phi*Jt + Vmat(rho,t)
                   bSigma
                   }

  ## make bdsmatrix of full Sigma
    fullSigma<-function(phi, rho, n, t) {
            sigma.i<-bSigma(phi,rho,n,t)
            fullSigma<-bdsmatrix(rep(t,n),rep(as.numeric(sigma.i),n))
            fullSigma
            }

  ## concentrated likelihood
    ll.c<-function(phirho, y, X, n, t, w) {
            phi<-phirho[1]
            rho<-phirho[2]

            ## perform GLS

            ## use direct matrix here, make it a bdsmatrix obj.
            ## and invert later, just as in REmod{panel10}
            sigma.1<-fullSigma(phi,rho,n,t)
            b.hat<-solve( crossprod(X,solve(sigma.1,X)), crossprod(X,solve(sigma.1,y)) )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat, solve(sigma.1,ehat)) / (n*t)
            bhat<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)
            e <- bhat[[2]]
            s2e <- bhat[[3]]

            uno <- n/2*log(1-rho^2)
            due <- -n/2*log( d2(rho,t) * (1-rho)^2 * phi + 1 )
            tre <- -(n*t)/2*log(s2e)
            ## use only the TxT block (phiJt+Vrho)^(-1) and calculate block by block
            cinque <- -1/(2*s2e)*crossprod(e, solve(sigma.1, e))

            const <- -(n*t)/2*log(2*pi)
            ll.c <- const+uno+due+tre+cinque
            llc <- - ll.c
        }

  ## iterate (=traballa) until convergence:

  myphirho0 <- coef0

  optimum<-nlminb(myphirho0, ll.c,
                  lower=c(1e-8,-0.999), upper=c(10e8,0.999),
                  control=list(x.tol=x.tol, rel.tol=rel.tol, trace=trace),
                  y=y, X=X, n=n, t=t, w=w, ...)


  myphirho <- optimum$par
  myll <- optimum$objective

  ## optimal values of parms:
  phi<-myphirho[1]
  rho<-myphirho[2]

  ## perform GLS
            ## as above: use direct matrix here, make it a bdsmatrix obj.
            ## and invert later, just as in REmod{panel10}
            sigma.1<-fullSigma(phi,rho,n,t)
            b.hat<-solve( crossprod(X,solve(sigma.1,X)), crossprod(X,solve(sigma.1,y)) )
            ehat<-y-X%*%b.hat
            sigma2ehat<-crossprod(ehat, solve(sigma.1,ehat)) / (n*t)
            beta<-list(betahat=b.hat,e=ehat,sigma2=sigma2ehat)

  ## names for coefs and error comp.s
  nam.beta <- dimnames(X)[[2]]
  nam.errcomp <- c("phi","rho")

  ## calc. cov(b) by GLS
  covB<-as.numeric(beta[[3]]) * solve(crossprod(X, solve(fullSigma(phi, rho, n, t), X)))
  dimnames(covB) <- list(nam.beta, nam.beta)

  ## calc. cov(phi,rho) by numerical Hessian
  covPRL <- solve(-fdHess(myphirho, function(x) -ll.c(x,y,X,n,t,w))$Hessian)
  dimnames(covPRL) <- list(nam.errcomp, nam.errcomp)

  ## make (separate) coefficients' vectors
  betas <- as.vector(beta[[1]])
  errcomp <- c(phi, rho)
  names(betas) <- nam.beta
  names(errcomp) <- nam.errcomp

  RES <- list(betas=betas, errcomp=errcomp,
              covB=covB, covPRL=covPRL, ll=myll)

  return(RES)

  }

