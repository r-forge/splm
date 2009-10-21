\name{bsjktest}
\alias{bsjktest.splm}
\alias{bsjktest.formula}
\alias{bsjktest}
\title{Baltagi, Song, Jung and Koh LM test for spatial panels}
\description{
   Baltagi, Song, Jung and Koh marginal or conditional LM test for
   spatial error correlation, serial error correlation or random effects in panel models
}
\usage{
bsjktest(x,...)

\method{bsjktest}{formula}(x, data, w, test=c(paste("C",1:3,sep="."),"J"), index=NULL, ...)
\method{bsjktest}{splm}(x, x2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{an object of class  \code{formula} or \code{splm}}
  \item{data}{A data frame}
  \item{w}{An object of either class \code{listw} or class \code{matrix}}
  \item{test}{One of c(paste("C",1:3,sep="."),"J") the test to be performed}
  \item{index}{if not NULL (default), a vector containing the individual and time indexes}
  \item{x2}{}
  \item{...}{additional arguments to be passed}
}
\details{
If \code{test=J} the joint test for spatial or serial error correlation or random individual 
effects is returned; if \code{test} is one of \code{C.1}, \code{C.2} or \code{C.3} 
the conditional test for, respectively, spatial errror correlation; serial error correlation; 
random individual effects; allowing for the other two effects to be present is returned.
}
\value{
  An object of class \code{htest}
}
\references{Baltagi, B.H., Song, S.H., Jung xx and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{xxx}, xxx--xxx.}
\author{Giovanni Millo}
\seealso{\code{\link{bsktest}}}

\examples{
data(Produc, package="Ecdat")
Produc <- Produc[Produc$year<1975, ]
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
\dontrun{test1<-bsjktest(fm,data=Produc, w=usaww,
  test="C.1")}
\dontrun{test1}
}
\keyword{htest}
