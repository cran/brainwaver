\name{wave.trans}
\alias{wave.trans}
\title{Computation of the wavelet transform}
\description{
Uses the wavelet decomposition implemented by Whitcher in the library \code{waveslim}. See all the details there.
}
\usage{
wave.trans(x, method = "modwt", wf = "la8", n.levels = 4, boundary = "periodic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{original vector to be decomposed}
 \item{method}{wavelet decomposition to be used, algorithm implemented in the \code{waveslim} package (Whitcher, 2000). By default, the Maximal Overlap Discrete Wavelet Transform is used \code{"modwt"}. It is also possible to use the classical Discrete Wavelet Transform \code{"dwt"}.}
  \item{wf}{name of the wavelet filter to use in the decomposition.  By default
    this is set to \code{"la8"}, the Daubechies orthonormal compactly
    supported wavelet of length \eqn{L=8} (Daubechies, 1992), least
    asymmetric family.}
  \item{n.levels}{specifies the depth of the decomposition.  This must be a number
    less than or equal to
    \eqn{\log_2(\mbox{length}(x))}{log(length(x),2)}.}
  \item{boundary}{Character string specifying the boundary condition.   If
    \code{boundary=="periodic"} the default, then the vector you
    decompose is assumed to be periodic on its defined interval,\cr
    if \code{boundary=="reflection"}, the vector beyond its boundaries
    is assumed to be a symmetric reflection of itself.}
}
\details{
See the library package \code{waveslim} (Whitcher, 2000). 
}
\value{
 Object of class \code{"modwt"}, basically, a list with the following
  components 
  \item{d?}{Wavelet coefficient vectors.}
  \item{s?}{Scaling coefficient vector.}
  \item{wavelet}{Name of the wavelet filter used.}
  \item{boundary}{How the boundaries were handled.}
}
\references{ 
R. Gencay, F. Selcuk and B. Whitcher (2001)
  \emph{An Introduction to Wavelets and Other Filtering Methods in
    Finance and Economics},
  Academic Press.

 D. B. Percival and A. T. Walden (2000)
  \emph{Wavelet Methods for Time Series Analysis},
  Cambridge University Press.


}
\author{ S. Achard }
\examples{

data(brain) # the result brain is a matrix
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]

PreCG.R<-brain[,1]
# LA(8)
PreCG.R.la8 <- wave.trans(PreCG.R, wf="la8")
names(PreCG.R.la8) <- c("w1", "w2", "w3", "w4", "v4")
## plot partial MODWT for PreCG.R data
par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
plot.ts(PreCG.R, axes=FALSE, ylab="", main="(a)")
for(i in 1:5)
  plot.ts(PreCG.R.la8[[i]], axes=FALSE, ylab=names(PreCG.R.la8)[i])
axis(side=1, at=seq(0,518,by=50),
  labels=c(0,"",100,"",200,"",300,"",400,"",500))

}
\keyword{ts}
