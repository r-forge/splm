\name{spreml}
\alias{spreml}
\title{Spatial Panel Random Effects Model Estimation}
\description{Maximum likelihood (ML) estimation of spatial random effects panel}
\usage{
spreml(formula, data, index = NULL, w, lag=FALSE,
           errors = c("semsrre","semsr","srre","semre","re","sr","sem"),
           pvar = FALSE, hess=FALSE, quiet=TRUE,
           initval = c("zeros", "estimate"),
           x.tol=1.5e-18, rel.tol=1e-15,
           ...)}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated}
  \item{data}{a \code{data.frame}}
  \item{w}{an object of either class \code{listw} or class \code{matrix}}
  \item{index}{}
  \item{lag}{}
  \item{errors}{}
  \item{pvar}{}
  \item{hess}{}
  \item{quiet}{}
  \item{initval}{}
  \item{x.tol}{}
  \item{rel.tol}{}
  \item{...}{}
}
\details{
  The models are estimated by two-step Maximum Likelihood.

  The covariance structures allowed are:
  \emph{semsr:} spatial and serial correlation in the idiosyncratic error
    term; 
  \emph{ssrre:} individual random effects and serial correlation;
  \emph{semre:} individual random effects and spatial correlation;
  \emph{re:} individual random effects;
  \emph{ssr:} serial correlation;
  \emph{sem:} spatial correlation.  
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
  \item{fitted.values}{the fitted values, calculated as \code{\hat{y}=X \hat{\beta}}}
  \item{sigma2}{GLS residuals variance}
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
  \item{logLik}{the value of the log likelihood function at the optimum }
  \item{errors}{the value of the \code{errors} argument}
}
\references{}
\author{Giovanni Millo}
\note{}
\seealso{\code{\link{spregm}}, \code{\link{spfeml}}}

\examples{
data(Produc, package = "Ecdat")
\dontrun{data(usaww)}
Produc <- Produc[Produc$year<1974, ]
\dontrun{fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp}
\dontrun{respaterr <- spreml(fm, data = Produc, w = usaww, errors="semre")}
\dontrun{summary(respaterr)}
}

\keyword{spatial}