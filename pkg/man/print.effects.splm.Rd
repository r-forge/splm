\name{print.effects.splm}
\alias{print.effects.splm}
\title{method for printing fixed effects from objects of class effects.splm}

\description{
Methods used for printing fixed effects from objects of class \code{effects.splm} generated using 
\link{effects.splm}
}

\usage{
\method{print}{effects.splm}(x,digits= max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{x}{an object of class \code{effects.splm}}
\item{digits}{specifies the minimum number of significant digits to be printed in values. See \code{\link{print.default}} for details}
\item{...}{additional arguments to be passed}
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
print(eff)
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err) 
print(eff)
}

\keyword{spatial}