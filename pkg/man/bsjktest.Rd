\name{bsjktest}
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

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{an object of class  \code{formula} or \code{splm}}
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. An optional data frame containing the variables
  in the model. When the obect is a \code{\link{data.frame}}, the first two columns may contain the indexes. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{w}{An object of class \code{listw} or a \code{matrix}}
  \item{test}{One of \code{c(paste("C",1:3,sep="."),"J")} the test to be performed}
  \item{...}{additional arguments to be passed}
}
\details{
If \code{test="J"} the joint test for spatial or serial error correlation or random individual 
effects is returned; if \code{test} is one of \code{C.1}, \code{C.2} or \code{C.3} 
the conditional test for, respectively, spatial errror correlation; serial error correlation; 
random individual effects (with the other two effects possibly present) is returned.
}
\value{
  An object of class \code{htest}
}
\references{
Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.
}
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
