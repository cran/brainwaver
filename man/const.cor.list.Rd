\name{const.cor.list}
\alias{const.cor.list}
\title{Computation of the list of correlation matrices}
\description{
 Computes the list of the correlation matrices in terms of the scale of the wavelet decomposition.
}
\usage{
const.cor.list(data.mat, names.data = 0, method = "modwt", wf = "la8", 
                n.levels = 4, boundary = "periodic", p.corr = 0.975, 
                save.wave = FALSE, export.data = FALSE)
}
\arguments{
  \item{data.mat}{ matrix containing the data time series. Each column of the matrix represents one time series.}
  \item{names.data}{optional character vector containing the name associated to the column of the matrix \code{data.mat}.}
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
  \item{p.corr}{(one minus the) two-sided p-value for the confidence interval}
  \item{save.wave}{logical. If TRUE all the wavelet coefficient are saved.}
  \item{export.data}{logical. If TRUE the correlation matrices with the upper and lower bound are exported to text file.}
}
\details{
  This function uses the wavelet decomposition implemented in the R package \code{waveslim}, (whitcher, 2000).
}
\value{
 Object of class \code{"Wave Correlation"}, basically, a list with the following
  components
 \item{d?}{Correlation matrix for each scale of the wavelet decomposition.}
 \item{lowerd?}{matrix containing the lower bound of the correlation for each scale of the wavelet decomposition.}
 \item{upperd?}{matrix containing the upper bound of the correlation for each scale of the wavelet decomposition.}
}
\references{
R. Gencay, F. Selcuk and B. Whitcher (2001)
  \emph{An Introduction to Wavelets and Other Filtering Methods in
    Finance and Economics},
  Academic Press.

 D. B. Percival and A. T. Walden (2000)
  \emph{Wavelet Methods for Time Series Analysis},
  Cambridge University Press.

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 

}
\author{S. Achard}
\note{change between version 1.1 and 1.2, now the length of the time series is saved with 
the values of the correlation.}
\seealso{\code{\link{const.var.list}}, \code{\link{read.cor.txt}}, \code{\link{save.cor.txt}}}
\examples{
data(brain)
brain<-as.matrix(brain)
# WARNING : To process only the first five regions
brain<-brain[,1:5]


n.levels<-4
wave.cor.list<-const.cor.list(brain,n.levels=n.levels)

tot.regions <- dim(brain)[2]
n.series <- dim(brain)[1]
col.regions<-1

nb.comp.regions <- 8
    comp.regions <- round(runif(nb.comp.regions,2,tot.regions))
    sym.region <- col.regions+1

comp.regions <- c(sym.region,comp.regions)

    name.ps <- "example-1.ps"
    postscript(name.ps)
    par(mfrow=c(3,3))
    for(k in 1:(nb.comp.regions+1)){
	
        reg <- comp.regions[k]


       cor.vector<-matrix(0,4,3)
	for(i in 1:n.levels){

       cor.vector[i,1]<-(wave.cor.list[[i]])[1,reg]
       cor.vector[i,2]<-(wave.cor.list[[i+n.levels]])[1,reg]
       cor.vector[i,3]<-(wave.cor.list[[i+2*n.levels]])[1,reg]

}

         title <- paste("1 -- ",comp.regions[k],sep="")
        matplot(2^(0:(n.levels-1)),cor.vector,main=title,type="b", 
                 log="x", pch="*LU", xaxt="n", lty=1, col=c(1,4,4), 
                 xlab="Wavelet Scale",ylab="Wavelet Covariance",ylim=c(-0.5,1))
        axis(side=1, at=2^(0:7))
        abline(h=0)
      }
    dev.off()






}
\keyword{multivariate}
\keyword{ts}
