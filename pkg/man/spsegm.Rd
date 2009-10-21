\name{spsegm}
\alias{spsegm}
\title{Spatial Simultaneous Equations}

\description{Feasible generalized three stages least square estimator (FGS3SLS)
of symultaneous systems of spatially interrelated cross sectional equations of the form: 

\deqn{Y = Y B + X C +  \bar{Y} L + U}

\deqn{ U = \bar{U} R + E}

with \eqn{Y = (y_1,\dots,y_m)}, 
\eqn{X = (x_1,\dots,x_m)}, \eqn{U = (u_1,\dots,u_m)}, 
\eqn{\bar{Y} = (\bar{y}_1,\dots,\bar{y}_m)} 
and \eqn{\bar{y}_j = W y_j} \eqn{j=1,\dots,m}, 
\eqn{E = (e_1,\dots,e_m)}, \eqn{\bar{U} = (\bar{u}_1,\dots,\bar{u}_m)} 
and \eqn{\bar{u}_j = W u_j} \eqn{j=1,\dots,m}. 
B, C, L and \eqn{R=diag_{j=1}^m (\rho_j)} 
are matrix of parameters


}

\usage{
spsegm(formula,data=list(), spec='default', listw, method='spatialsim',zero.policy = FALSE)
}
\arguments{
  \item{formula}{an object of class \code{formula} (see Details)}
  \item{data}{an object of class \code{\link{data.frame}}. An optional data frame containing the variables
  in the model}
  \item{spec}{a list used to define the specification (see Details) }
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw}}
  \item{method}{\code{spatialsim}}
  \item{zero.policy}{See \code{lagsarlm} for details}
}


\details{
The function can be specified with any number of equations. 
The number of equations is determined through the \code{formula} object.
In particular, all the responses for all equations should be organized in a matrix and
specified on the left hand side of the formula object. Analogously, all the explanatory 
variables from all equations should be specified on the right hand side of the 
formula. The number of equations is then retrieved from the dimension of the matrix 
containing the responses.

The argument \code{spec} should either be \code{default} or a \code{list}.
When \code{spec="default"}, all the equations in the system have the same explanatory variables 
(including the intercept). The intercept can be eliminated (for all the equations at the same time)
as it is done in \code{lm}. 

When \code{spec=list} it should have the same length as the number of equations (i.e.
the number of columns of y in the formula object). The presence of the intercept is managed 
by \code{spec} (see Examples). When trying to eliminate the intercept in the classical way, if 
\code{spec} is a list an error message will occour. 

}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{FG3SLS coefficients estimate of the model parameters (for all equations)}
  \item{vcov}{the variance covariance matrix of the estimated coefficients}
  \item{type}{'spsegm'}
  \item{model}{the matrix of the data used (responses in each equation are reported first, then the 
  explanatory variables)}
  \item{N}{the number of cross-sectional observations}  
  \item{Eq}{the number of equations in the system}  
  \item{k}{the number of columns of the matrix of regressors 
  (i.e. this corresponds to the number of explanatory variables in each equation only when 
  \code{spec="default"})}
  \item{call}{the call used to create the object}
  \item{terms}{the \code{terms} object generated from \code{formula} and some optional arguments}
  \item{Xnames}{the names of the variables in the matrix of explanatory variables}
  \item{Ynames}{the names of the responses}
  \item{spec}{the argument \code{spec}}
}

\references{
  Kelejian, H.H. and Prucha, I.R. (2004) 
Estimation of Simultaneous systems of spatially interrelated cross sectional equations,
    \emph{Journal of Econometrics}, \bold{118}, pages 27--50.

  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.

}

\author{Gianfranco Piras\email{gpiras@mac.com}}

\seealso{
\code{\link{spregm}}
}
\examples{
n<-70
N<-n^2
ny<-3
nx<-3
W<-nb2listw(cell2nb(n,n)) 
y<-matrix(rnorm(N*ny),N,ny)
X<-matrix(runif(N*nx,-10,10),N,nx)
ynseq<-paste("y",seq(1:ny))
colnames(y)<-ynseq
xnseq<-paste("y",seq(1:nx))
colnames(X)<-xnseq	
data<-cbind(y,X)
data<-data.frame(data)
SE<-spsegm( y~ X ,spec='default', listw=W)
summary(SE)
spec=list(c(1,2,3),c(1,3,4),c(2,3,4))
#c(1,2,3): intercept (1) and first and second columns of X (2,3)
#c(1,3,4): intercept (1) and second and third columns of X (3,4)
#c(2,3,4): first, second and third columns of X (2,3,4)
SE<-spsegm(y~ X ,spec=list(c(1,2,3),c(1,3,4),c(2,3,4)), listw=W)
}

\keyword{spatial}