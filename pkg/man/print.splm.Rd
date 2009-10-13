\name{print.splm}
\alias{print.splm}
\title{print method for class splm}

\description{
Method used to print objects of class \code{'summary.splm'} and \code{'splm'}
}

\usage{
\method{print}{splm}(x, digits = max(3,getOption("digits") -3), ...)
}

\arguments{
\item{x}{an object of class \code{'summary.splm'} and \code{splm}}
\item{digits}{minimal number of significant digits, see \code{\link{print.default}}}
\item{...}{additional arguments to be passed}
}


\details{
The summary function summary.splm 
returns an objects of class 'splm' 
organized in a coefficient matrix.

Also a matrix for the error components, 
or the spatial coefficient will be generated 
depending on the estimated model. 

When the 'splm' objects comes from 
the function 'spsegm', the summary 
will be generated looping over the number
of equations in the system.
}

\value{

}

\references{
}

\seealso{
  \code{\link{spreml}}, \code{\link{spregm}}, \code{lagsarlm},
  \code{errorsarlm}, \code{\link{spreml}}, \code{\link{effects}}
}
\author{ Giovanni Millo \email{Giovanni\_Millo@Generali.com}, Gianfranco Piras\email{gpiras@mac.com}}

\examples{
data(Produc, package = "Ecdat") 
\dontrun{data(usaww)}
Produc <- Produc[Produc$year<1975, ] 
\dontrun{GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")}
\dontrun{summary(GM)}
}
\keyword{spatial}