\name{const.var.list}
\alias{const.var.list}
\title{Computation of the list of variance vectors}
\description{
  Computes the list of the variance vectors in terms of the scale of the wavelet decomposition.
}
\usage{
const.var.list(data.mat, names.data = 0, method = "modwt", wf = "la8", 
                n.levels = 4, boundary = "periodic", save.wave = FALSE, 
                export.data = FALSE)
}

\arguments{

  \item{data.mat}{ matrix containing the data time series. Each column of the matrix represents one time series. }
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
   \item{save.wave}{logical. If TRUE all the wavelet coefficient are saved.}
  \item{export.data}{logical. If TRUE the variance vectors with the upper and lower bound are exported to text file.}
}
\details{
  This function uses the wavelet decomposition implemented in the R package \code{waveslim}, (whitcher, 2000).
}
\value{
 Object of class \code{"Wave Variance"}, basically, a list with the following
  components
 \item{d?}{Variance vectors for each scale of the wavelet decomposition.}
 \item{lowerd?}{vector containing the lower bound of the variance for each scale of the wavelet decomposition.}
 \item{upperd?}{vector containing the upper bound of the variance for each scale of the wavelet decomposition.}
}

\references{
Gencay, R., F. Selcuk and B. Whitcher (2001)
  \emph{An Introduction to Wavelets and Other Filtering Methods in
    Finance and Economics},
  Academic Press.

 Percival, D. B. and A. T. Walden (2000)
  \emph{Wavelet Methods for Time Series Analysis},
  Cambridge University Press.

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 

}
\author{S. Achard}
\seealso{\code{\link{const.cor.list}}, \code{\link{read.var.txt}}, \code{\link{save.var.txt}}}
\examples{
data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]



n.levels<-4
wave.var.list<-const.var.list(brain,n.levels=n.levels)

tot.regions <- dim(brain)[2]
n.series <- dim(brain)[1]

nb.num.regions <- 9
   num.regions <- round(runif(nb.num.regions,2,tot.regions))
    par(mfrow=c(3,3))
    for(k in 1:(nb.num.regions)){
	
        reg <- num.regions[k]

       var.vector<-matrix(0,4,3)
for(i in 1:n.levels){

	var.vector[i,1]<-(wave.var.list[[i]])[reg]
       	var.vector[i,2]<-(wave.var.list[[i+n.levels]])[reg]
       	var.vector[i,3]<-(wave.var.list[[i+2*n.levels]])[reg]
      }
        title <- num.regions[k]
        matplot(2^(0:(n.levels-1)),var.vector,main=title,type="b", 
                log="x", pch="*LU", xaxt="n", lty=1, col=c(1,4,4), 
                xlab="Wavelet Scale",ylab="Wavelet Variance")
        axis(side=1, at=2^(0:7))
        abline(h=0)
      }
  

}
\keyword{multivariate}
\keyword{ts}

