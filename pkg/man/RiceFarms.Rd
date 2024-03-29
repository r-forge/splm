\name{RiceFarms}
\docType{data}
\alias{RiceFarms}
\title{Production of Rice in India }
\description{
  yearly observations of 171 farms

 \emph{number of observations} :  1026

   \emph{country} :  Indonesia

\emph{economic topic} : producer behavior

\emph{econometrics topic} : error component

}
\usage{data(RiceFarms)}
\format{A dataframe containing :
  \describe{
    \item{id}{the farm identifier}
    \item{time}{the growing season}
    \item{size}{the total area cultivated with rice, measured in hectares}
    \item{status}{land status, on of \code{'owner'} (non sharecroppers, owner operators or leasholders or both), \code{'share'} (sharecroppers), \code{'mixed'} (mixed of the two previous status)}
    \item{varieties}{one of \code{'trad'} (traditional varieties), \code{'high'} (high yielding varieties) and \code{'mixed'} (mixed varieties)}
    \item{bimas}{bIMAS is an intensification program ; one of \code{'no'} (non-bimas famer), \code{'yes'} (bimas farmer) or \code{'mixed'} (part but not all of farmer's land was registered to be in the bimas program)}
    \item{seed}{seed in kilogram}
    \item{urea}{urea in kilogram}
    \item{phosphate}{phosphate in kilogram}
    \item{pesticide}{pesticide cost in Rupiah}
    \item{pseed}{price of seed in Rupiah per kg}
    \item{purea}{price of urea in Rupiah per kg}
    \item{pphosph}{price of phosphate in Rupiah per kg}
    \item{hiredlabor}{hired labor in hours}
    \item{famlabor}{family labor in hours}
    \item{totlabor}{total labor (excluding harvest labor)}
    \item{wage}{labor wage in Rupiah per hour}
    \item{goutput}{gross output of rice in kg}
    \item{noutput}{net output, gross output minus harvesting cost (paid in terms of rice)}
    \item{price}{price of rough rice in Rupiah per Kg}
    \item{region}{one of \code{'wargabinangun'}, \code{'langan'}, \code{'gunungwangi'}, \code{'malausma'}, \code{'sukaambit'}, \code{'ciwangi'}}
  }
}
\source{
    Journal of Applied Econometrics Data Archive.
}
\references{
     Qu Feng and William C. Horrace, (2012) \dQuote{Alternative Measures of Technical Efficiency: Skew, Bias and Scale}, \emph{Journal of Applied Econometrics}, \bold{forthcoming}.

    Horrace, W.C.  and  P.  Schmidt (1996) \dQuote{Confidence statements for efficiency estimates from stochastic frontier models}, \emph{Journal of Productivity Analysis}, \bold{ 7}, 257--282.

}

\keyword{datasets}
