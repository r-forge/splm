\name{effects.splm}
\alias{effects.splm}
\title{method for extracting fixed effects}

\description{
Methods used for extracting fixed effects from objects of class 'splm' where type is one of 
"fixed effects lag" or "fixed effects error"
}

\usage{
\method{effects}{splm}(object,...)
}

\arguments{
\item{object}{an object of class 'splm'}
\item{...}{additional arguments to be passed over}
}

\details{
If the argument \code{object} is not of class 'splm' 
the function will terminate with an error message.

If the argument \code{object} is of class 'splm' 
type is not one of 
"fixed effects lag" or "fixed effects error",
the function will terminate with an error message.

}


\value{
  An object of class \code{"effects.splm"}
  \item{res}{a list whose elements are various type of fixed effects and the intercept (when present)}
}

\references{
  Elhorst, J.P. (2003) 
Specification and estimation of spatial panel data models,
    \emph{International Regional Science Review}, \bold{26}, pages 244--268.

  Elhorst, J.P. (2009) 
Spatial panel data models,
    \emph{In} Fischer, M.M. and Getis, A. (eds),
    \emph{Handbook of Applied Spatial Analysis} Springer, Berlin.

}

\author{ Gianfranco Piras \email{gpiras@mac.com}}

\seealso{
\code{\link{spfeml}} \code{print.effects.splm} \code{summary.effects.splm}
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
summary(eff)
}
\keyword{spatial}