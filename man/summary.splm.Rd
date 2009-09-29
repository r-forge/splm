\name{summary.splm}
\alias{summary.splm}
\title{summary method for class splm}

\description{
Methods used for presenting the results of objects of class 'splm'
}

\usage{

\method{summary}{splm}(object,...)
}

\arguments{
\item{object}{an object of class 'splm'}
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
An object of class 'summary.splm'

}

\references{

}

\seealso{
  \code{\link{spreml}}, \code{\link{spregm}}, 
 \code{\link{spreml}}, \code{\link{effects}}
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