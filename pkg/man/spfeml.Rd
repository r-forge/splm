\name{spfeml}
\alias{spfeml}
\title{Spatial Panel Fixed Effects Models Estimation}
\description{Maximum likelihood (ML) estimation of spatial panel
fixed effects models}

\usage{
spfeml(formula, data, index=NULL,listw, model=c("lag","error"),
effects=c('pooled','spfe','tpfe','sptpfe'), method="eigen",
na.action=na.fail,quiet=TRUE,zero.policy = FALSE,
 interval = c(-1, 0.999), tol.solve = 1e-10, tol.opt = .Machine$double.eps^0.5)
}
\arguments{
  \item{formula}{a description of the model to be fit }
  \item{data}{an object of class \code{\link{data.frame}}. An optional data frame containing the variables
  in the model. The first two columns may contain the indexes. see \code{index}}
  \item{index}{if not NULL (default), a vector containing the individual and time indexes}
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw}}
  \item{model}{may be set to "lag" or "error" }
  \item{effects}{one of "pooled" (no spatial effects, i.e. OLS on the pooled model),
  "spfe" (only spatial fixed effects), "tpfe" (only time period fixed effects),
  "sptpfe" (both time period and spatial fixed effects)}
  \item{method}{"eigen" (default) - the Jacobian is computed as
  \eqn{\prod_{i=1}^N (1 - \rho\omega_i)} (with \eqn{\omega_i} the eigenvalues of the spatial
  weigths matrix) using \code{eigenw}, and "spam" or "Matrix" using code from the spam or Matrix
  packages to calculate the determinant (see \code{lagsarlm} for details)}
  \item{na.action}{a function (default na.fail), can also be na.omit or na.exclude with consequences
  for residuals and fitted values - (see \code{lagsarlm} for details)}
  \item{quiet}{default=TRUE; if FALSE, reports function values during optimization}
 \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours,
if FALSE (default) assign NA - causing \code{spfeml} to terminate with an error}
  \item{interval}{search interval for spatial parameter when not using method="eigen";
  default is c(-1,1); method="Matrix" will attempt to search for an appropriate interval}
  \item{tol.solve}{tolerance for detecting linear dependence in the columns of matrices to be inverted
  (default= 1.0e-10)- (see \code{lagsarlm} for details)}
  \item{tol.opt}{desired accuracy for optimization
  (default=square root of double precision machine tolerance)}
}


\details{
The function deals with both fixed effects spatial lag and error models, of the form:

\deqn{ y = \lambda (I_T \otimes W_N)y + (\iota_T \otimes \alpha) + X\beta + \epsilon }

where \eqn{\lambda} is called the spatial autoregressive coefficient, and:
\deqn{ y =  (\iota_T \otimes \alpha) + X\beta + u }
\deqn{ u =  \rho (\iota_T \otimes W_N)u + \epsilon}

respectively, where \eqn{\rho} is the spatial autocorrelation coefficient.
Both \eqn{\lambda} and \eqn{\rho} are found by \code{\link{optimize}},
and the model parameters by a feasible GLS procedure.

The asymptotic standard error of
\eqn{\lambda} is only computed when
\code{method=eigen}, because the full matrix
operations involved would be costly for
dimensions associated with the choice of "spam" or "Matrix".
The same applies to the coefficient covariance matrix.

The model variables are first transformed depending on the
typology of fixed effects. As an example, if time period
fixed effects are considered, the demeaned form is obtained
by subtracting the average for each cross-sectional unit computed over the time
dimension.

The transformed variables are then used in a concentrated likelihood
approach to estimate the spatial parameters.

Once the spatial parameters are obtained, the model parameters are
then estimated by a feasible GLS procedure.

Fixed effects can then be extracted using the function \code{\link{effects}}.
}


\value{
  An object of class \code{"splm"}.
  \item{coefficients}{coefficients estimate of the model parameters }
  \item{vcov}{the asymptotic variance covariance matrix of the estimated coefficients (when calculated)}
  \item{type}{'fixed effects spatial lag (or error) model'}
  \item{spat.coeff}{the spatial coefficients \eqn{\rho} or \eqn{\lambda} }
  \item{residuals}{the GLS residuals}
  \item{fitted.values}{difference between response variable and residuals}
  \item{sigma2}{GLS residuals variance}
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
  \item{logLik}{the value of the log likelihood function at the optimum }
  \item{method}{the method used to calculate the Jacobian}
  \item{effects}{the typology of effects included in the model}
  \item{res.eff}{a list whose elements are the fixed effects (see \code{\link{effects}})}

}

\references{
  Elhorst, J.P. (2003)
Specification and estimation of spatial panel data models,
    \emph{International Regional Science Review}, \bold{26}, pages 244--268.

  Elhorst, J.P. (2009)
Spatial panel data models,
    \emph{In} Fischer, M.M. and Getis, A. (eds),
    \emph{Handbook of Applied Spatial Analysis} Springer, Berlin.

Anselin, L. and Le Gallo, J. and Jayet, H. (2008)
Spatial Panel Econometrics, \emph{In}
Matyas, L. and Sevestre, P. (eds), \emph{The econometrics of Panel
Data, Fundamentals and Recent Developments in Theory and Practice (3rd
Edition)}, pages 624--660. Springer-Verlag, Berlin Heidelberg.
}

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{
  \code{\link{spreml}}, \code{\link{spregm}}, \code{lagsarlm},
  \code{errorsarlm}, \code{\link{effects}}
}
\examples{
data(Produc, package = "Ecdat")
\dontrun{data(usaww)}
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
\dontrun{lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE)}
\dontrun{summary(lag)}
\dontrun{eff <- effects(lag)}
\dontrun{err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)}
\dontrun{summary(err)}
\dontrun{eff <- effects(err)}
\dontrun{summary(eff)}
}

\keyword{spatial}