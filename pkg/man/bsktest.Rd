\name{bsktest}
\alias{bsktest}
\alias{bsktest.formula}
\alias{bsktest.lm}
\alias{bsktest.splm}

\title{Baltagi, Song and Koh LM test for spatial panels}
\description{
   Baltagi, Song and Koh marginal or conditional LM test for
   spatial error correlation or random effects in panel models
}
\usage{
bsktest(x,...)
\method{bsktest}{formula}(x, data, w, test=c("SLM1","SLM2","LMJOINT","CLMlambda","CLMmu"), index=NULL, ...)
\method{bsktest}{lm}(x, w, index=NULL, test=c("SLM1","SLM2","LMJOINT"), ...)
\method{bsktest}{splm}(x, w, index=NULL, test=c("CLMlambda","CLMmu"), ...)
}
\arguments{
\item{x}{an object of class  \code{formula} or \code{lm} or \code{splm} }
\item{data}{an object of class \code{\link{data.frame}}. An optional
  data frame containing the variables in the model.The first two columns may contain the indexes. See \code{index}}
\item{w}{an object of class \code{listw} created for example by \code{nb2listw} }
\item{index}{if not NULL (default), a vector containing the individual and time indexes}
\item{test}{one of "SLM1","SLM2","LMJOINT","CLMlambda","CLMmu", the test to be performed}
\item{...}{additional arguments to be passed}
}

\value{
an object of class \code{htest}
}
\references{Baltagi, B.H., Song, S.H. and Koh, W. (2003)
Testing panel data regression models with spatial error correlation.
\emph{Journal of Econometrics}, \bold{117}, 123--150.}
\author{Gianfranco Piras}
\seealso{\code{\link{bsjktest}}}

\examples{
data(Produc, package="Ecdat")
Produc <- Produc[Produc$year<1975, ]
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
test1<-bsktest(fm,data=Produc, w=mat2listw(usaww),
  test="SLM1")
class(test1)
test1
ml2 <- spfeml(fm, data = Produc, , mat2listw(usaww), model = "error", effects = "pooled")
class(ml2)
test5bis<-bsktest(ml2, w=mat2listw(usaww),index=Produc[,c(1,2)] ,test="CLMmu")
summary(test5bis)
}
\keyword{htest}