\name{spreml}
\alias{spreml}
\title{Spatial Panel Random Effects Model Estimation}
\description{Maximum likelihood (ML) estimation of spatial random effects panel of the form:
\deqn{
 y_it = X'_it \beta + u_it, 
 i=1,\dots,N, t=1,\dots,T
}
with 
\deqn{
u_t=\mu +\epsilon_t
}
and
\deqn{
\epsilon_t=\lambda W \epsilon_t+ \nu_t }
where
\deqn{\nu_t=\rho \nu_t-1 +e_t.
}
Depending on the restrictions on the vector of parameters one can
differently combine error features giving rise to various nested
specifications. 
}


\usage{
spreml(formula, data, index = NULL, w, w2=w, lag=FALSE,
           errors = c("semsrre","semsr","srre","semre",
                      "re", "sr", "sem","ols", "sem2re"),
           pvar = FALSE, hess=FALSE, quiet=TRUE,
           initval = c("zeros", "estimate"),
           x.tol=1.5e-18, rel.tol=1e-15,
           ...)}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated}
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{\link{data.frame}}, the first two columns may contain the indexes. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{w}{an object of class \code{listw} or a \code{matrix}}
  \item{w2}{an object of class \code{listw} or a \code{matrix}}
  \item{lag}{default=\code{FALSE}. If \code{TRUE}, a spatial lag of the dependent variable is added.}
  \item{errors}{one of \code{c("semsrre","semsr","srre","semre","re","sr","sem", "ols", "sem2re")}. See details.}
  \item{pvar}{if \code{TRUE} the \code{pvar} function is called}
  \item{hess}{if \code{TRUE} use numerical Hessian instead of GLS for the standard errors of the estimates}
  \item{quiet}{if \code{FALSE} report function and parameters values during optimization}
  \item{initval}{one of \code{c("zeros", "estimate")}, the initial values for the parameters. If \code{"zeros"} a vector of zeros is used. if \code{"estimate"} the initial 
  values are retreived from the estimation of the nested specifications. Alternatively, a numeric vector can be specified.}
  \item{x.tol}{Tolerance. See \code{\link{nlminb}} for details. }
  \item{rel.tol}{Relative tolerance. See \code{\link{nlminb}} for details. }
  \item{...}{additional argument to pass over to other functions}
}
\details{
  The models are estimated by two-step Maximum Likelihood.

  The covariance structures allowed are:
    \code{"semsrre"} random effects, spatial and serial correlation in the idiosyncratic error
    term; 
  \code{"semsr"} spatial and serial correlation in the idiosyncratic error
    term; 
  \code{"srre"} individual random effects and serial correlation;
  \code{"semre"} individual random effects and spatial correlation;
  \code{"re"} individual random effects;
  \code{"sr"} serial correlation;
  \code{"sem"} spatial correlation.
  \code{"ols"} none (spherical residuals).
  \code{"sem2re"} spatial correlation in both individual random effects
    and idiosncratic errors (Kapoor, Kelejian and Prucha's specification).
}
\value{
  An object of class \code{"splm"}.
  \item{coefficients}{coefficients estimate of the model parameters }
  \item{arcoef}{the coefficient for the spatial lag on \code{y}}
  \item{errcomp}{the estimates of the error variance components}
  \item{vcov}{the asymptotic variance covariance matrix of the estimated coefficients}
  \item{vcov.arcoef}{the asymptotic variance of the
  estimated spatial lag parameter}
  \item{vcov.errcomp}{the asymptotic variance covariance matrix of the
  estimated error covariance parameters}
  \item{type}{'random effects ML'}
  \item{residuals}{the model residuals}
  \item{fitted.values}{the fitted values, calculated as \eqn{\hat{y}=X \hat{\beta}}}
  \item{sigma2}{GLS residuals variance}
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
  \item{logLik}{the value of the log likelihood function at the optimum }
  \item{errors}{the value of the \code{errors} argument}
}
\references{
Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.
}

\author{Giovanni Millo}

\seealso{\code{\link{spregm}}, \code{\link{spfeml}}}

\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1974, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
## the two standard specifications:
## random effects panel with spatial errors
respaterr <- spreml(fm, data = Produc, w = usaww, errors="semre")
summary(respaterr)
## random effects panel with spatial lag
respatlag <- spreml(fm, data = Produc, w = usaww, errors="re", lag=TRUE)
summary(respatlag)
}

\keyword{spatial}