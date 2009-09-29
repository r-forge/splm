\name{write.effects.splm}
\alias{write.effects.splm}
\title{method for writing a table with fixed effects}

\description{
Methods used for extracting fixed effects from objects of class 'splm' where type is one of 
"fixed effects lag" or "fixed effects error" and writing the corresponding table to disk
}

\usage{
\method{write}{effects.splm}(object,filename='effects.csv',...)
}

\arguments{
\item{object}{an object of class 'splm'}
\item{filename}{name of the file, default set to effects.csv}
\item{...}{additional arguments}
}

\details{
If the argument \code{x} is not of class 'splm' 
the function will terminate with an error message.

If the argument \code{x} is of class 'splm' 
type is not one of 
"fixed effects lag" or "fixed effects error",
the function will terminate with an error message.

}


\value{
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
\dontrun{data(usaww)}
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
\dontrun{lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE) }
\dontrun{summary(lag)}
\dontrun{eff <- effects(lag) }
\dontrun{err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)}
\dontrun{summary(err)}
\dontrun{eff <- effects(err) }
\dontrun{summary(eff)}
}

\keyword{spatial}