\name{spfeml}
\alias{spfeml}
\title{Spatial Panel Fixed Effects Models Estimation}
\description{Maximum likelihood (ML) estimation of spatial panel
fixed effects models}

\usage{
spfeml(formula, data = list(), index = NULL,listw,listw2 = NULL, 
model= c("lag","error", "sarar"), effects=c('pooled','spfe','tpfe','sptpfe'), 
method ="eigen", na.action = na.fail, quiet = TRUE, zero.policy = NULL, interval = NULL, 
tol.solve = 1e-10, control = list(), legacy = FALSE, llprof = NULL)
}


\arguments{
  \item{formula}{a description of the model to be fit }
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{\link{data.frame}}, the first two columns may contain the indexes. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{listw}{an object of class \code{listw} created for example by \code{nb2listw}}
  \item{listw2}{only if \code{model} is \code{"sarar"}. an object of class \code{listw} created for example by \code{nb2listw}. if not given, set to the same spatial weights as the listw argument}
  \item{model}{one of \code{c("lag", "error", "sarar")} }
  \item{effects}{one of \code{"pooled"} (no spatial effects, i.e. OLS on the pooled model),
  \code{"spfe"} (only spatial fixed effects), \code{"tpfe"} (only time period fixed effects),
  \code{"sptpfe"} (both time period and spatial fixed effects)}
  
  \item{method}{\code{"eigen"} (default) - the Jacobian is computed as
  \eqn{\prod_{i=1}^N (1 - \rho \omega_i) } (with \eqn{\omega_i} the eigenvalues of the spatial
  weigths matrix) using \code{eigenw} - (see \code{lagsarlm} for details). 
  Other available methods are \code{"spam"} or \code{"Matrix"} for strictly symmetric weights lists of styles "B" and "C", 
  or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", 
  using code from the \pkg{spam} or \pkg{Matrix} packages to calculate the determinant; \code{"LU"}
  provides an alternative sparse matrix decomposition approach. 
  Finally, \code{"Chebyshev"} and Monte Carlo \code{"MC"} provides approximations to the computation of the log-determinant.}
\item{na.action}{see \code{lagsarlm} for details } 
  \item{quiet}{default=\code{TRUE}; if \code{FALSE}, reports function values during optimization}
  \item{zero.policy}{see \code{lagsarlm} for details}
  \item{interval}{search interval for spatial parameters}
  \item{tol.solve}{tolerance for detecting linear dependence in the columns of matrices to be inverted (default= 1.0e-10)- (see \code{lagsarlm} for details)}
  \item{control}{a list of control parameters for the optimization - see \code{lagsarlm} for details}
  \item{legacy}{FALSE When TRUE and model is \code{"lag"} calculate goodness of fit measures}
  \item{llprof}{default NULL, can either be an integer, to divide the feasible ranges into a grid of points, or a two-column matrix of spatial coefficient values, at which to evaluate the likelihood function}
}


\details{
The function deals with fixed effects spatial sarar, lag and error models, of the general form:

\deqn{ y = \lambda (I_T \otimes W_N) y + (\iota_T \otimes \alpha) + X\beta + \epsilon }
\deqn{ u =  \rho (\iota_T \otimes W_N)u + \epsilon}

where \eqn{\rho} is called the spatial autoregressive coefficient, and
\eqn{\lambda} is the spatial autocorrelation coefficient.

Both \eqn{\rho} and \eqn{\lambda} are found by \code{\link{optimize}},
and the model parameters by a feasible GLS procedure.

The asymptotic standard error of 
\eqn{\rho} in the error model is only computed when
\code{method=eigen}, because the full matrix
operations involved would be costly for
dimensions associated with the choice of \code{"spam"} or \code{"Matrix"}.
The same applies to the coefficient covariance matrix.

The asymptotic standard error of 
\eqn{\rho} and \eqn{lambda} in the full model 
are calculated by a numerical Hessian.

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
  \code{\link{spreml}}, \code{\link{spgm}}, \code{\link{effects}}
}
\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE)
summary(lag)
eff <- effects(lag)
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err)
print(eff)
write.effects.splm(eff)
}

\keyword{spatial}