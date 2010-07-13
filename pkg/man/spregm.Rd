\name{spregm}
\alias{spregm}
\title{GM estimator for spatial panel data models}

\description{GM estimator for panel data models with spatially correlated 
errors components of the form: 

\deqn{ y_N(t) = X_N(t) \beta + u_N(t) }

\deqn{ u_N(t) = \rho W_N u_N(t) + \epsilon(t)}

\deqn{ \epsilon_N = (e_T \otimes I_N ) \mu_N + \nu_N }

where \eqn{ \rho}, and the variance components \eqn{\sigma^2_\mu} and \eqn{\sigma^2_\nu} 
are estimated by GM, and the model coefficients by GLS. 
}

\usage{
spregm(formula, data=list(), index=NULL, w, method = c("init", "weigh", "fulweigh"), lag=FALSE)
}

\arguments{
  \item{formula}{a description of the model to be fit. The details of model specification are given
  for \code{lm} }
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{\link{data.frame}}, the first two columns may contain the indexes. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{w}{an object of class \code{listw} created for example by \code{nb2listw}. 
  If a \code{matrix} is given as input the function will transform it using \code{mat2listw}}
  \item{method}{\code{"init"} (default) defines the set of GM estimator to be used. Alternatives are  
  \code{"weigh"} and \code{"fulweigh"} (See Details)}
  \item{lag}{if \code{TRUE} a spatial lag of the dependent variable is added to the regression equation}
}


\details{ When method=init, the initial estimator is calculated. This first set 
of GM estimators is based only on a subset of the moments conditions and assigns 
equal weigths to each of them. 
When method=fulweigh, the second set of GM estimators is calculated. This estimator 
is based on the full set of moments conditions. It also involves the expression for the variance 
covariance matrix of the sample moments calculated under the assumption of 
normally distributed innovations. The calculation of the trace terms in the expression 
of the variance covariance matrix of the sample moments 
uses codes from the \code{\link{Matrix}} package. 
When method=weigh, the third set of GM estimator is used. This is motivated by computational 
issues. The procedure is analogous to the second one but uses a simplified expression for the 
variance covariance matrix of the sample moments

Note that \eqn{\sigma^2_\mu} is not reported. \eqn{\sigma^2_1} is reported instead. 
However, a value for \eqn{\sigma^2_\mu} can easily be obtained from:
\deqn{\sigma^2_1 = \sigma^2_\nu + T \sigma^2_\mu}

The function also produces an estimate for \eqn{\theta} which is a
function of the variance components.
}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{GLS coefficients estimate of the model parameters}
  \item{vcov}{the variance covariance matrix of the estimated coefficients}
  \item{residuals}{the GLS residuals}
  \item{fitted.values}{difference between response variable and residuals}
  \item{sigma2}{GLS residuals variance}
  \item{type}{'random effect GM'}
  \item{rho}{a vector including the spatial parameter and the variance components (see Details)}  
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
}

\references{
  Kapoor, M., Kelejian, H.H. and Prucha, I.R. (2007) 
  Panel data model with spatially correlated error components,
    \emph{Journal of Econometrics}, \bold{140}, pages 97--130. 
    
  Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
    \emph{International Economic Review}, \bold{40}, pages 509--533.
    
      Kelejian, H.H. and Prucha, I.R. (1999) 
A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
Model with Autoregressive Disturbances,
    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.

}

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{
  \code{\link{spreml}}, \code{\link{spsegm}}
}
\examples{
data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")
summary(GM)
}

\keyword{spatial}